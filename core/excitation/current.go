// SPDX-License-Identifier: MPL-2.0

package excitation

import (
	"fmt"

	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/mesh"
	"oblikovati.org/traceon/core/solver"
)

// CoilExcitation assigns total currents to the named coils of a triangle mesh and builds the
// resulting current-ring charges — the magnetic pre-field for a magnetostatic solve. Currents
// live on triangle groups (a coil is a solid cross-section), separate from the line-element
// excitations, mirroring add_current's requirement that radial coils be triangle meshes.
type CoilExcitation struct {
	tm       *mesh.TriangleMesh
	currents map[string]float64
}

// NewCoils creates a coil excitation over the triangle mesh tm.
func NewCoils(tm *mesh.TriangleMesh) *CoilExcitation {
	return &CoilExcitation{tm: tm, currents: map[string]float64{}}
}

// AddCurrent assigns a total current (amperes) to the named coil. The current is spread as a
// uniform density over the coil's cross-section. Port of add_current for a radial coil. Panics
// if name is not a coil group in the mesh.
func (c *CoilExcitation) AddCurrent(name string, amps float64) {
	if _, ok := c.tm.PhysicalToTriangles[name]; !ok {
		names := make([]string, 0, len(c.tm.PhysicalToTriangles))
		for n := range c.tm.PhysicalToTriangles {
			names = append(names, n)
		}
		panic(fmt.Sprintf("excitation: %q is not a coil in the mesh; have %v", name, names))
	}
	c.currents[name] = amps
}

// Charges builds the combined current-ring charges from all assigned coils. Each coil's
// triangles carry a uniform current density (its total current divided by its total
// cross-sectional area), so the integrated current equals the assigned value. Coils with zero
// area are skipped. Feed the result to solver.NewFieldRadialBEMFull / as the magnetostatic
// pre-field. Mirrors the CURRENT branch of the radial solver assembly.
func (c *CoilExcitation) Charges() solver.CurrentCharges {
	var tris []geom3d.Triangle
	var currents []float64
	for name, amps := range c.currents {
		group := c.tm.Group(name)
		area := 0.0
		for _, t := range group {
			area += geom3d.TriangleArea(t[0], t[1], t[2])
		}
		if area == 0 {
			continue
		}
		density := amps / area
		for _, t := range group {
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
