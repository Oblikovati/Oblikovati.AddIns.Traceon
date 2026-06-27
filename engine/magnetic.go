// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"encoding/json"
	"math"
	"strconv"
	"strings"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/solver"
)

// attrCurrents / attrMagnets name the document attributes holding per-coil currents and
// per-magnet axial magnetisations: JSON objects mapping body index (as a string) to amperes
// and to A/m respectively, e.g. {"3": 2.5}.
const (
	attrCurrents     = "currents"
	attrMagnets      = "magnets"
	attrPermeability = "permeability"
)

// coilTriGrid is the (r, z) triangulation resolution for a coil cross-section: the bounding
// region is split into this many bands each way, then two triangles per cell. Finer grids
// integrate the current ring more accurately at more cost.
const coilTriGrid = 6

// coil is a sectioned current-carrying body: its (r, z) cross-section (cm, for rendering) and
// the total azimuthal current it carries (amperes).
type coil struct {
	prof    *profile
	current float64
}

// coilCurrents reads the per-coil current map (body index → amperes) from the active
// document's traceon/currents attribute. Returns an empty map when unset or unreadable.
func (e *Engine) coilCurrents() map[int]float64 {
	return e.floatAttributeMap(attrCurrents)
}

// magnet is a sectioned permanent-magnet body: its (r, z) cross-section (cm, for rendering)
// and its axial magnetisation (A/m).
type magnet struct {
	prof          *profile
	magnetisation float64
}

// buildMagnetCharges turns permanent magnets into magnetostatic surface charges: each
// boundary element carries a magnetic charge equal to the axial magnetisation projected onto
// the element normal (n_z · M). Mirrors MagnetostaticSolverRadial.get_permanent_magnet_field.
func buildMagnetCharges(magnets []magnet) solver.EffectivePointCharges {
	var lines []radial.Line
	var charges []float64
	for _, m := range magnets {
		ml, _, _ := m.prof.lineElements(0, cmToMetres)
		for _, l := range ml {
			lines = append(lines, l)
			n := radial.ElementNormal(l) // (n_r, n_z) unit normal
			charges = append(charges, n[1]*m.magnetisation)
		}
	}
	if len(lines) == 0 {
		return solver.EffectivePointCharges{}
	}
	jac, pos := radial.FillJacobianBufferRadial(lines)
	return solver.EffectivePointCharges{Charges: charges, Jac: jac, Pos: pos}
}

// isCoil reports whether a body is a current coil rather than an electrode: it carries a
// current in the traceon/currents attribute, or its name contains "coil" (the drivable
// convention when the attribute cannot be set). The resolved current (amperes) is returned.
func isCoil(index int, name string, currents map[int]float64, defaultCurrent float64) (float64, bool) {
	if c, ok := currents[index]; ok {
		return c, true
	}
	if strings.Contains(strings.ToLower(name), "coil") {
		return defaultCurrent, true
	}
	return 0, false
}

// isMagnet reports whether a body is an axially-magnetised permanent magnet: it carries a
// magnetisation in the traceon/magnets attribute, or its name contains "magnet". The resolved
// axial magnetisation (A/m, along the optical axis) is returned.
func isMagnet(index int, name string, magnets map[int]float64, defaultMag float64) (float64, bool) {
	if m, ok := magnets[index]; ok {
		return m, true
	}
	if strings.Contains(strings.ToLower(name), "magnet") {
		return defaultMag, true
	}
	return 0, false
}

// magnetMagnetisations reads the per-magnet axial magnetisation map (body index → A/m) from
// the active document's traceon/magnets attribute.
func (e *Engine) magnetMagnetisations() map[int]float64 {
	return e.floatAttributeMap(attrMagnets)
}

// ironPermeabilities reads the per-iron relative-permeability map (body index → μr) from the
// active document's traceon/permeability attribute.
func (e *Engine) ironPermeabilities() map[int]float64 {
	return e.floatAttributeMap(attrPermeability)
}

// isIron reports whether a body is a magnetizable (soft-magnetic) body that RESPONDS to the
// field — it has a permeability in the traceon/permeability attribute, or its name contains
// "iron" or "magnetizable". The resolved relative permeability is returned.
func isIron(index int, name string, perms map[int]float64, defaultMu float64) (float64, bool) {
	if mu, ok := perms[index]; ok {
		return mu, true
	}
	n := strings.ToLower(name)
	if strings.Contains(n, "iron") || strings.Contains(n, "magnetizable") {
		return defaultMu, true
	}
	return 0, false
}

// iron is a sectioned magnetizable body: its (r, z) cross-section (cm) and relative permeability.
type iron struct {
	prof         *profile
	permeability float64
}

// floatAttributeMap reads a JSON {bodyIndex: value} attribute in the traceon set into a map.
func (e *Engine) floatAttributeMap(name string) map[int]float64 {
	out := map[int]float64{}
	docID, ok := e.activeDocID()
	if !ok {
		return out
	}
	res, err := e.api.Attributes().Get(docID, attrSet, name)
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

// buildCoilCharges turns coil cross-sections into the effective current rings the magnetic
// field is evaluated from: each cross-section is triangulated in the (r, z) meridian plane
// (scaled to metres), carrying a uniform current density = total current / cross-section
// area. Mirrors MagnetostaticSolverRadial.get_current_field.
func buildCoilCharges(coils []coil) solver.CurrentCharges {
	var tris []geom3d.Triangle
	var currents []float64
	for _, c := range coils {
		ct := triangulateCrossSection(c.prof, cmToMetres)
		if len(ct) == 0 {
			continue
		}
		area := 0.0
		for _, t := range ct {
			area += geom3d.TriangleArea(t[0], t[1], t[2])
		}
		if area == 0 {
			continue
		}
		density := c.current / area
		for _, t := range ct {
			tris = append(tris, t)
			currents = append(currents, density)
		}
	}
	if len(tris) == 0 {
		return solver.CurrentCharges{}
	}
	jac, pos := geom3d.FillJacobianBuffer3D(tris)
	return solver.CurrentCharges{Currents: currents, Jac: jac, Pos: pos}
}

// buildIronCharges solves for the magnetisation induced in the magnetizable bodies by the
// pre-existing field (the coils' and permanent magnets' field, supplied as preField), and
// returns the resulting magnetostatic surface charges. Mirrors MagnetostaticSolverRadial's
// magnetizable solve: each iron boundary element is a Magnetizable row whose right-hand side
// is the negated pre-field flux through its normal.
func buildIronCharges(irons []iron, preField solver.PreField) (solver.EffectivePointCharges, error) {
	var lines []radial.Line
	var perms []float64
	for _, ir := range irons {
		il, _, _ := ir.prof.lineElements(0, cmToMetres)
		for _, l := range il {
			lines = append(lines, l)
			perms = append(perms, ir.permeability)
		}
	}
	if len(lines) == 0 {
		return solver.EffectivePointCharges{}, nil
	}
	types := make([]radial.ExcitationType, len(lines))
	for i := range types {
		types[i] = radial.Magnetizable
	}
	return solver.SolveMagnetostatic(lines, types, perms, preField)
}

// combineCharges concatenates two effective-charge sets (their charges, Jacobian and position
// buffers) into one — used to merge the permanent-magnet and induced-iron magnetostatic charges.
func combineCharges(a, b solver.EffectivePointCharges) solver.EffectivePointCharges {
	if len(a.Charges) == 0 {
		return b
	}
	if len(b.Charges) == 0 {
		return a
	}
	return solver.EffectivePointCharges{
		Charges: append(append([]float64{}, a.Charges...), b.Charges...),
		Jac:     append(append(radial.JacobianBuffer{}, a.Jac...), b.Jac...),
		Pos:     append(append(radial.PositionBuffer{}, a.Pos...), b.Pos...),
	}
}

// ironExtent returns the (r, z) bounding box (cm) spanning every iron body.
func ironExtent(irons []iron) (rMax, zMin, zMax float64) {
	rMax, zMin, zMax = 0, math.Inf(1), math.Inf(-1)
	for _, ir := range irons {
		r, lo, hi := ir.prof.extent()
		rMax = math.Max(rMax, r)
		zMin = math.Min(zMin, lo)
		zMax = math.Max(zMax, hi)
	}
	return rMax, zMin, zMax
}

// ironNode draws each iron cross-section outline (cm).
func ironNode(irons []iron) (graphicsLines, bool) {
	return profilesNode(func(yield func(*profile)) {
		for _, ir := range irons {
			yield(ir.prof)
		}
	})
}

// triangulateCrossSection fills the coil's (r, z) cross-section with triangles in the meridian
// plane (vertices at (r, 0, z), scaled). It grids the cross-section's bounding region — the
// common revolved-rectangle coil fills its bounds, so a uniform grid captures it exactly.
func triangulateCrossSection(p *profile, scale float64) []geom3d.Triangle {
	rMin, rMax, zMin, zMax := math.Inf(1), math.Inf(-1), math.Inf(1), math.Inf(-1)
	for _, loop := range p.loops {
		for _, pt := range loop {
			rMin, rMax = math.Min(rMin, pt[0]), math.Max(rMax, pt[0])
			zMin, zMax = math.Min(zMin, pt[1]), math.Max(zMax, pt[1])
		}
	}
	if rMax <= rMin || zMax <= zMin {
		return nil
	}
	vert := func(r, z float64) geom3d.Vec3 { return geom3d.Vec3{r * scale, 0, z * scale} }
	var tris []geom3d.Triangle
	for i := 0; i < coilTriGrid; i++ {
		for j := 0; j < coilTriGrid; j++ {
			r0 := rMin + (rMax-rMin)*float64(i)/coilTriGrid
			r1 := rMin + (rMax-rMin)*float64(i+1)/coilTriGrid
			z0 := zMin + (zMax-zMin)*float64(j)/coilTriGrid
			z1 := zMin + (zMax-zMin)*float64(j+1)/coilTriGrid
			tris = append(tris,
				geom3d.Triangle{vert(r0, z0), vert(r1, z0), vert(r1, z1)},
				geom3d.Triangle{vert(r0, z0), vert(r1, z1), vert(r0, z1)})
		}
	}
	return tris
}

// coilExtent returns the (r, z) bounding box (cm) spanning every coil — for placing the
// magnetic-field sampling grid.
func coilExtent(coils []coil) (rMax, zMin, zMax float64) {
	rMax, zMin, zMax = 0, math.Inf(1), math.Inf(-1)
	for _, c := range coils {
		r, lo, hi := c.prof.extent()
		rMax = math.Max(rMax, r)
		zMin = math.Min(zMin, lo)
		zMax = math.Max(zMax, hi)
	}
	return rMax, zMin, zMax
}

// magnetExtent returns the (r, z) bounding box (cm) spanning every permanent magnet.
func magnetExtent(magnets []magnet) (rMax, zMin, zMax float64) {
	rMax, zMin, zMax = 0, math.Inf(1), math.Inf(-1)
	for _, m := range magnets {
		r, lo, hi := m.prof.extent()
		rMax = math.Max(rMax, r)
		zMin = math.Min(zMin, lo)
		zMax = math.Max(zMax, hi)
	}
	return rMax, zMin, zMax
}

// profilesNode draws each profile outline (cm) as line segments — used for the coil and
// magnet overlays (distinguished by colour at the call site).
func profilesNode(loops func(yield func(*profile))) (node graphicsLines, has bool) {
	loops(func(p *profile) {
		for _, loop := range p.loops {
			for i := 0; i+1 < len(loop); i++ {
				node.add(loop[i], loop[i+1])
				has = true
			}
		}
	})
	return node, has
}

// coilNode draws each coil cross-section outline (cm).
func coilNode(coils []coil) (graphicsLines, bool) {
	return profilesNode(func(yield func(*profile)) {
		for _, c := range coils {
			yield(c.prof)
		}
	})
}

// magnetNode draws each permanent-magnet cross-section outline (cm).
func magnetNode(magnets []magnet) (graphicsLines, bool) {
	return profilesNode(func(yield func(*profile)) {
		for _, m := range magnets {
			yield(m.prof)
		}
	})
}

// graphicsLines accumulates (r, z) line segments for a graphics node (drawn in the xz-plane).
type graphicsLines struct {
	coords  []float64
	indices []int
}

func (g *graphicsLines) add(a, b geom2d.Point2) {
	base := len(g.coords) / 3
	g.coords = append(g.coords, a[0], 0, a[1], b[0], 0, b[1])
	g.indices = append(g.indices, base, base+1)
}
