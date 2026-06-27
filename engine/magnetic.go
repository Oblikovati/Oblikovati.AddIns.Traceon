// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"encoding/json"
	"math"
	"strconv"
	"strings"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/solver"
)

// attrCurrents names the document attribute holding per-coil currents: a JSON object mapping
// body index (as a string) to amperes, e.g. {"3": 2.5}.
const attrCurrents = "currents"

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
	out := map[int]float64{}
	docID, ok := e.activeDocID()
	if !ok {
		return out
	}
	res, err := e.api.Attributes().Get(docID, attrSet, attrCurrents)
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

// coilNode draws each coil cross-section outline (cm) in a copper colour so coils are
// visually distinct from electrodes.
func coilNode(coils []coil) (node graphicsLines, has bool) {
	for _, c := range coils {
		for _, loop := range c.prof.loops {
			for i := 0; i+1 < len(loop); i++ {
				node.add(loop[i], loop[i+1])
				has = true
			}
		}
	}
	return node, has
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
