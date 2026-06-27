// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"encoding/json"
	"fmt"
	"math"
	"strconv"

	"oblikovati.org/traceon/core/constants"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/solver"
	"oblikovati.org/traceon/core/tracing"
)

// StudyResult summarizes a completed study for the status bar / CLI.
type StudyResult struct {
	ElectrodeCount   int
	CoilCount        int
	ElementCount     int
	RayCount         int
	FocusZ           float64 // axial focus position (cm), NaN if the beam does not cross the axis
	GraphicsClientID string
}

// boundsMargin pads the tracing/field box around the geometry extent (metres).
const boundsMargin = 0.02

// driftFactor extends the downstream trace region to this multiple of the geometry's axial
// span, so a focus that forms past the lens is captured. driftRadius widens the radial bound
// (× rMax) so a converging ray is not clipped before it reaches the axis.
const (
	driftFactor = 8.0
	driftRadius = 2.0
)

// beamAperture is the paraxial launch radius as a fraction of the geometry's max radius:
// rays start near the axis (well inside any bore) so they sample the paraxial focusing
// field rather than grazing the electrode wall.
const beamAperture = 0.25

// cmToMetres converts the host DB unit (cm, ADR-0042) to SI metres, the unit the BEM solve
// and the tracer (which use SI constants) require for physically-correct beam dynamics.
// metresToCm converts results back for rendering in the host's coordinate system.
const (
	cmToMetres = 0.01
	metresToCm = 100.0
)

// attrSet / attrVoltages name the document attribute holding per-electrode voltages: a JSON
// object mapping body index (as a string) to volts, e.g. {"0":0,"1":5000,"2":0}.
const (
	attrSet      = "traceon"
	attrVoltages = "voltages"
)

// electrode is one sectioned body: its (r, z) meridian (cm, for rendering) and the voltage
// applied to it.
type electrode struct {
	prof    *profile
	voltage float64
}

// electronChargeOverMass is q/m for an electron (C/kg): -e / m_e.
var electronChargeOverMass = -constants.ElementaryCharge / constants.ElectronMass

// RunStudy is the whole add-in study: section every solid body in the active part into an
// axisymmetric (r, z) electrode, solve the combined electrostatic BEM for all of them
// together (so the electrodes interact), trace a beam through the resulting field, and push
// the electrodes + trajectories + potential map back into the viewport.
//
// Per-electrode voltages come from the document attribute traceon/voltages (a JSON
// {bodyIndex: volts} map); bodies absent from it use the panel's default voltage. Geometry is
// converted from the host DB unit (cm) to SI metres for the physics, and trajectories are
// converted back to cm for rendering. The bodyIndex parameter is ignored (kept for the
// command path); the study always sections the whole part.
func (e *Engine) RunStudy(int) (*StudyResult, error) {
	e.mu.Lock()
	params := e.params
	e.mu.Unlock()

	electrodes, coils, err := e.collectBodies(params)
	if err != nil {
		return nil, err
	}
	if len(electrodes) == 0 && len(coils) == 0 {
		return nil, fmt.Errorf("no solid bodies could be sectioned into electrodes or coils")
	}

	// Electrostatic charges from the electrode boundaries, magnetic current rings from the
	// coils — combined into one field the beam is traced through (Lorentz force E + μ₀ v×H).
	var elec solver.EffectivePointCharges
	lines, types, values := assembleElements(electrodes)
	if len(lines) > 0 {
		elec, err = solver.SolveElectrostatic(lines, types, values)
		if err != nil {
			return nil, fmt.Errorf("solve electrostatic: %w", err)
		}
	}
	current := buildCoilCharges(coils)
	bem := field.NewFieldRadialBEMFull(elec, solver.EffectivePointCharges{}, current)

	rays := e.traceBeam(bem, electrodes, coils, params)

	nodes := renderNodes(electrodes, coils, bem, rays)
	if err := e.pushGraphics(nodes); err != nil {
		return nil, err
	}
	return &StudyResult{
		ElectrodeCount:   len(electrodes),
		CoilCount:        len(coils),
		ElementCount:     len(lines),
		RayCount:         len(rays),
		FocusZ:           focusZcm(rays),
		GraphicsClientID: graphicsClientID,
	}, nil
}

// collectBodies sections every solid body in the active part and sorts each into an electrode
// (voltage boundary) or a current coil. A body is a coil if it carries a current in the
// traceon/currents attribute or its name contains "coil"; otherwise it is an electrode.
//
// Electrode voltages come from traceon/voltages (a JSON {bodyIndex: volts} map); when that is
// unset, the einzel convention applies — the panel voltage biases the axially central
// electrode, the others are grounded. Coil currents come from traceon/currents, else the
// panel coil current.
func (e *Engine) collectBodies(params studyParams) ([]electrode, []coil, error) {
	list, err := e.api.Body().List()
	if err != nil {
		return nil, nil, fmt.Errorf("list bodies: %w", err)
	}
	currents := e.coilCurrents()
	voltages := e.electrodeVoltages()

	var elecProfs []*profile
	var elecIdx []int
	var coils []coil
	for _, b := range list.Bodies {
		if !b.Solid {
			continue
		}
		prof, err := e.extractProfile(b.Index)
		if err != nil {
			continue // a body that cannot be sectioned (e.g. non-axisymmetric) is skipped
		}
		if amps, ok := isCoil(b.Index, b.Name, currents, params.coilCurrent); ok {
			coils = append(coils, coil{prof: prof, current: amps})
			continue
		}
		elecProfs = append(elecProfs, prof)
		elecIdx = append(elecIdx, b.Index)
	}

	central := -1
	if len(voltages) == 0 {
		central = centralElectrode(elecProfs)
	}
	electrodes := make([]electrode, len(elecProfs))
	for i, prof := range elecProfs {
		var v float64
		switch {
		case len(voltages) > 0:
			v = voltages[elecIdx[i]]
		case i == central:
			v = params.voltage
		}
		electrodes[i] = electrode{prof: prof, voltage: v}
	}
	return electrodes, coils, nil
}

// centralElectrode returns the index of the electrode whose axial mid-point is closest to
// the centroid of all electrode mid-points — the one the einzel default biases.
func centralElectrode(profs []*profile) int {
	if len(profs) == 0 {
		return -1
	}
	mids := make([]float64, len(profs))
	centroid := 0.0
	for i, p := range profs {
		_, lo, hi := p.extent()
		mids[i] = (lo + hi) / 2
		centroid += mids[i]
	}
	centroid /= float64(len(profs))
	best, bestDist := 0, math.Inf(1)
	for i, m := range mids {
		if d := math.Abs(m - centroid); d < bestDist {
			best, bestDist = i, d
		}
	}
	return best
}

// assembleElements flattens every electrode's profile into one combined BEM element set
// (lines in metres) so all electrodes are solved together in a single influence matrix.
func assembleElements(electrodes []electrode) ([]radial.Line, []radial.ExcitationType, []float64) {
	var lines []radial.Line
	var types []radial.ExcitationType
	var values []float64
	for _, el := range electrodes {
		l, t, v := el.prof.lineElements(el.voltage, cmToMetres)
		lines = append(lines, l...)
		types = append(types, t...)
		values = append(values, v...)
	}
	return lines, types, values
}

// electrodeVoltages reads the per-electrode voltage map (body index → volts) from the active
// document's traceon/voltages attribute. Returns an empty map when unset or unreadable.
func (e *Engine) electrodeVoltages() map[int]float64 {
	out := map[int]float64{}
	docID, ok := e.activeDocID()
	if !ok {
		return out
	}
	res, err := e.api.Attributes().Get(docID, attrSet, attrVoltages)
	if err != nil || !res.Found {
		return out
	}
	s, ok := res.Attribute.Value.Str()
	if !ok {
		return out
	}
	var byKey map[string]float64
	if json.Unmarshal([]byte(s), &byKey) != nil {
		return out
	}
	for k, v := range byKey {
		if i, err := strconv.Atoi(k); err == nil {
			out[i] = v
		}
	}
	return out
}

// activeDocID returns the active document's id, if any.
func (e *Engine) activeDocID() (uint64, bool) {
	docs, err := e.api.Documents().List()
	if err != nil {
		return 0, false
	}
	for _, d := range docs.Documents {
		if d.Active {
			return d.ID, true
		}
	}
	return 0, false
}

// traceBeam launches params.numRays parallel rays from below the geometry (in metres), spread
// across the radial aperture, along +z, and traces each through the BEM field. Returns the
// trajectories (in metres).
func (e *Engine) traceBeam(bem field.FieldRadialBEM, electrodes []electrode, coils []coil, params studyParams) [][]tracing.State {
	rMaxCm, zMinCm, zMaxCm := studyExtent(electrodes, coils)
	rMax, zMin, zMax := rMaxCm*cmToMetres, zMinCm*cmToMetres, zMaxCm*cmToMetres
	// Trace through a generous downstream drift region so the focus (which forms past the
	// lens) is captured. The radial bound is loose so a converging ray is not clipped early.
	drift := driftFactor * (zMax - zMin)
	bounds := tracing.Bounds{
		{-driftRadius * rMax, driftRadius * rMax},
		{-boundsMargin, boundsMargin},
		{zMin - boundsMargin, zMax + drift},
	}
	// The field the particle feels: E from the electrodes and H from the coils. The tracer
	// applies the Lorentz force a = q/m·(E + μ₀·v×H).
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		p := geom2d.Vertex{pos[0], pos[1], pos[2]}
		ef := bem.FieldAtPoint(p)
		hf := bem.CurrentFieldAtPoint(p)
		return geom3d.Vec3{ef[0], ef[1], ef[2]}, geom3d.Vec3{hf[0], hf[1], hf[2]}
	}
	v0 := tracing.VelocityVec(params.energyEV, geom3d.Vec3{0, 0, 1}, constants.ElectronMass)
	startZ := zMin - 0.5*boundsMargin

	var rays [][]tracing.State
	// Launch a PARAXIAL beam well inside the bore (the focus is a paraxial property, and rays
	// near the electrode wall would graze the conductor). beamAperture is a fraction of rMax.
	aperture := beamAperture * rMax
	for i := 0; i < params.numRays; i++ {
		frac := 0.0
		if params.numRays > 1 {
			frac = float64(i) / float64(params.numRays-1)
		}
		r0 := 0.04*rMax + frac*aperture
		_, states := tracing.TraceParticle(geom3d.Vec3{r0, 0, startZ}, v0, electronChargeOverMass, fieldFn, bounds, 1e-8)
		if len(states) > 1 {
			rays = append(rays, states)
		}
	}
	return rays
}

// combinedExtent returns the (r, z) bounding box (cm) spanning every electrode.
func combinedExtent(electrodes []electrode) (rMax, zMin, zMax float64) {
	rMax, zMin, zMax = 0, math.Inf(1), math.Inf(-1)
	for _, el := range electrodes {
		r, lo, hi := el.prof.extent()
		rMax = math.Max(rMax, r)
		zMin = math.Min(zMin, lo)
		zMax = math.Max(zMax, hi)
	}
	return rMax, zMin, zMax
}

// studyExtent returns the (r, z) bounding box (cm) spanning every electrode and coil.
func studyExtent(electrodes []electrode, coils []coil) (rMax, zMin, zMax float64) {
	er, ez0, ez1 := combinedExtent(electrodes)
	cr, cz0, cz1 := coilExtent(coils)
	return math.Max(er, cr), math.Min(ez0, cz0), math.Max(ez1, cz1)
}

// focusZcm returns the axial (z) focus of the ray bundle in cm, or NaN if it cannot be
// determined (fewer than two rays, or no axis crossing).
func focusZcm(rays [][]tracing.State) float64 {
	if len(rays) < 2 {
		return math.NaN()
	}
	// Average the per-ray axis crossing; simple and robust for a paraxial bundle.
	sum, n := 0.0, 0
	for _, ray := range rays {
		if z, ok := tracing.AxisIntersection(ray); ok {
			sum += z
			n++
		}
	}
	if n == 0 {
		return math.NaN()
	}
	return sum / float64(n) * metresToCm
}
