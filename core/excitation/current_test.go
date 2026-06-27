// SPDX-License-Identifier: MPL-2.0

package excitation

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/mesh"
)

// twoTriangleCoil is a unit-square coil cross-section (two triangles, total area 1) in the
// xz-plane, named "coil".
func twoTriangleCoil() *mesh.TriangleMesh {
	pts := []geom3d.Vec3{{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}}
	tris := [][3]int{{0, 1, 2}, {0, 2, 3}}
	return mesh.NewTriangleMesh(pts, tris, map[string][]int{"coil": {0, 1}})
}

// TestAddCurrentUniformDensity checks the assigned current is spread as a uniform density
// (total current / total area) over the coil's triangles.
func TestAddCurrentUniformDensity(t *testing.T) {
	tm := twoTriangleCoil() // total area 1
	c := NewCoils(tm)
	c.AddCurrent("coil", 10)

	ch := c.Charges()
	if len(ch.Currents) != 2 {
		t.Fatalf("currents = %d, want 2 (one per triangle)", len(ch.Currents))
	}
	for i, d := range ch.Currents {
		if math.Abs(d-10.0) > 1e-12 { // density = 10 A / area 1 = 10
			t.Errorf("density[%d] = %g, want 10 (current/area)", i, d)
		}
	}
}

// TestAddCurrentScalesWithArea checks two coils with equal current but different areas get
// inversely-scaled densities (a smaller coil is denser).
func TestAddCurrentScalesWithArea(t *testing.T) {
	// "big" has area 1, "small" has area 0.25.
	pts := []geom3d.Vec3{
		{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}, // big unit square
		{5, 0, 0}, {5.5, 0, 0}, {5.5, 0, 0.5}, {5, 0, 0.5}, // small 0.5×0.5 square
	}
	tris := [][3]int{{0, 1, 2}, {0, 2, 3}, {4, 5, 6}, {4, 6, 7}}
	tm := mesh.NewTriangleMesh(pts, tris, map[string][]int{"big": {0, 1}, "small": {2, 3}})

	c := NewCoils(tm)
	c.AddCurrent("big", 1)
	c.AddCurrent("small", 1)
	ch := c.Charges()

	// 4 triangles total; densities are 1/1=1 (big) and 1/0.25=4 (small).
	var ones, fours int
	for _, d := range ch.Currents {
		switch {
		case math.Abs(d-1) < 1e-9:
			ones++
		case math.Abs(d-4) < 1e-9:
			fours++
		default:
			t.Errorf("unexpected density %g", d)
		}
	}
	if ones != 2 || fours != 2 {
		t.Errorf("densities: got %d×1 and %d×4, want 2 and 2", ones, fours)
	}
}

// TestAddCurrentUnknownPanics checks assigning current to a non-coil is rejected.
func TestAddCurrentUnknownPanics(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("AddCurrent on an unknown coil did not panic")
		}
	}()
	NewCoils(twoTriangleCoil()).AddCurrent("nope", 1)
}

// TestChargesEmpty checks a coil excitation with no assigned currents yields empty charges.
func TestChargesEmpty(t *testing.T) {
	if ch := NewCoils(twoTriangleCoil()).Charges(); len(ch.Currents) != 0 {
		t.Errorf("currents = %d, want 0", len(ch.Currents))
	}
}
