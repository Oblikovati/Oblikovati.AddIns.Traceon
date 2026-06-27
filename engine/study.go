// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"fmt"

	"oblikovati.org/traceon/core/constants"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/solver"
	"oblikovati.org/traceon/core/tracing"
)

// StudyResult summarizes a completed study for the status bar / CLI.
type StudyResult struct {
	ElementCount     int
	RayCount         int
	GraphicsClientID string
}

// boundsMargin pads the tracing/field box around the geometry extent (cm).
const boundsMargin = 1.0

// electronChargeOverMass is q/m for an electron (C/kg): -e / m_e.
var electronChargeOverMass = -constants.ElementaryCharge / constants.ElectronMass

// RunStudy is the whole add-in study on a host body: section the body into an axisymmetric
// (r, z) electrode profile, solve the electrostatic BEM for its surface charges, trace a beam
// of rays through the resulting field, and push the electrode + trajectories + potential map
// back into the viewport as client graphics.
//
// Units note: the geometry arrives in the host DB unit (cm) and the study traces in those
// units; absolute beam dynamics depend on a cm→m calibration that a later milestone pins, but
// the trajectory shape (focusing/deflection) is already faithful.
func (e *Engine) RunStudy(bodyIndex int) (*StudyResult, error) {
	e.mu.Lock()
	params := e.params
	e.mu.Unlock()

	prof, err := e.extractProfile(bodyIndex)
	if err != nil {
		return nil, err
	}
	lines, types, values := prof.lineElements(params.voltage)
	if len(lines) == 0 {
		return nil, fmt.Errorf("section produced no BEM elements")
	}

	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		return nil, fmt.Errorf("solve electrostatic: %w", err)
	}
	bem := field.NewFieldRadialBEM(charges)

	rays := e.traceBeam(bem, prof, params)

	nodes := renderNodes(prof, bem, rays)
	if err := e.pushGraphics(nodes); err != nil {
		return nil, err
	}
	return &StudyResult{ElementCount: len(lines), RayCount: len(rays), GraphicsClientID: graphicsClientID}, nil
}

// traceBeam launches params.numRays parallel rays from below the geometry, spread across the
// radial aperture, along +z, and traces each through the BEM field. Returns the trajectories.
func (e *Engine) traceBeam(bem field.FieldRadialBEM, prof *profile, params studyParams) [][]tracing.State {
	rMax, zMin, zMax := prof.extent()
	bounds := tracing.Bounds{
		{-rMax - boundsMargin, rMax + boundsMargin},
		{-boundsMargin, boundsMargin},
		{zMin - boundsMargin, zMax + boundsMargin},
	}
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		ef := bem.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{ef[0], ef[1], ef[2]}, geom3d.Vec3{}
	}
	v0 := tracing.VelocityVec(params.energyEV, geom3d.Vec3{0, 0, 1}, constants.ElectronMass)
	startZ := zMin - 0.5*boundsMargin

	var rays [][]tracing.State
	aperture := 0.8 * rMax
	for i := 0; i < params.numRays; i++ {
		// Spread launch radii across the aperture; the innermost ray is near the axis.
		frac := 0.0
		if params.numRays > 1 {
			frac = float64(i) / float64(params.numRays-1)
		}
		r0 := 0.05*rMax + frac*aperture
		_, states := tracing.TraceParticle(geom3d.Vec3{r0, 0, startZ}, v0, electronChargeOverMass, fieldFn, bounds, 1e-8)
		if len(states) > 1 {
			rays = append(rays, states)
		}
	}
	return rays
}
