// SPDX-License-Identifier: MPL-2.0

package excitation

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/mesh"
	"oblikovati.org/traceon/core/radial"
)

// twoElectrodeMesh builds a tiny radial mesh with two named electrodes: "a" (two straight
// segments up the r=1 wall) and "b" (one segment across the top). Points are distinct so
// dedup leaves the line order intact.
func twoElectrodeMesh() *mesh.Mesh {
	pts := []geom2d.Vertex{{1, 0, 0}, {1, 0, 1}, {1, 0, 2}, {2, 0, 2}}
	lines := [][]int{{0, 1}, {1, 2}, {2, 3}}
	phys := map[string][]int{"a": {0, 1}, "b": {2}}
	return mesh.New(pts, lines, phys, false)
}

// TestAddVoltageFixed checks a fixed voltage is applied to every element of its electrode and
// only the assigned electrode's elements are active.
func TestAddVoltageFixed(t *testing.T) {
	exc := New(twoElectrodeMesh())
	exc.AddVoltage("a", 5)

	lines, types, values := exc.Electrostatic()
	if len(lines) != 2 {
		t.Fatalf("active lines = %d, want 2 (electrode b is unassigned)", len(lines))
	}
	for i := range lines {
		if types[i] != radial.VoltageFixed {
			t.Errorf("type[%d] = %d, want VoltageFixed", i, types[i])
		}
		if values[i] != 5 {
			t.Errorf("value[%d] = %g, want 5", i, values[i])
		}
	}
}

// TestAddVoltageFunc checks a position-dependent voltage is sampled at each element centre:
// electrode "b" spans (1,2)→(2,2) in (r,z), so its centre is r=1.5, z=2.
func TestAddVoltageFunc(t *testing.T) {
	exc := New(twoElectrodeMesh())
	exc.AddVoltageFunc("b", func(x, _, z float64) float64 { return x*10 + z })

	lines, types, values := exc.Electrostatic()
	if len(lines) != 1 {
		t.Fatalf("active lines = %d, want 1", len(lines))
	}
	if types[0] != radial.VoltageFun {
		t.Errorf("type = %d, want VoltageFun", types[0])
	}
	if want := 1.5*10 + 2; math.Abs(values[0]-want) > 1e-9 {
		t.Errorf("sampled voltage = %g, want %g (centre r=1.5, z=2)", values[0], want)
	}
}

// TestElectrostaticBoundary checks a boundary becomes a zero-permittivity dielectric.
func TestElectrostaticBoundary(t *testing.T) {
	exc := New(twoElectrodeMesh())
	exc.AddVoltage("a", 1)
	exc.AddElectrostaticBoundary("b")

	_, types, values := exc.Electrostatic()
	found := false
	for i, ty := range types {
		if ty == radial.Dielectric {
			found = true
			if values[i] != 0 {
				t.Errorf("boundary dielectric value = %g, want 0", values[i])
			}
		}
	}
	if !found {
		t.Error("no Dielectric element produced for the electrostatic boundary")
	}
}

// TestOrderPreserved checks active elements come out in mesh-line order.
func TestOrderPreserved(t *testing.T) {
	exc := New(twoElectrodeMesh())
	exc.AddVoltage("a", 1)
	exc.AddVoltage("b", 2)
	_, _, values := exc.Electrostatic()
	want := []float64{1, 1, 2}
	if len(values) != len(want) {
		t.Fatalf("values = %v, want %v", values, want)
	}
	for i := range want {
		if values[i] != want[i] {
			t.Errorf("value[%d] = %g, want %g", i, values[i], want[i])
		}
	}
}

// TestMagnetostaticActiveElements checks the magnetostatic builder selects scalar-potential
// and magnetizable elements (and nothing else), with the assigned values.
func TestMagnetostaticActiveElements(t *testing.T) {
	exc := New(twoElectrodeMesh())
	exc.AddMagnetostaticPotential("a", 50)
	exc.AddMagnetizable("b", 1000)

	_, types, values := exc.Magnetostatic()
	if len(types) != 3 { // electrode a = 2 lines, b = 1 line
		t.Fatalf("magnetostatic elements = %d, want 3", len(types))
	}
	for i, ty := range types {
		switch ty {
		case radial.MagnetostaticPot:
			if values[i] != 50 {
				t.Errorf("potential value[%d] = %g, want 50", i, values[i])
			}
		case radial.Magnetizable:
			if values[i] != 1000 {
				t.Errorf("permeability value[%d] = %g, want 1000", i, values[i])
			}
		default:
			t.Errorf("unexpected magnetostatic type %d at %d", ty, i)
		}
	}
	// The electrostatic builder must ignore magnetostatic excitations.
	if lines, _, _ := exc.Electrostatic(); len(lines) != 0 {
		t.Errorf("electrostatic builder returned %d magnetostatic elements, want 0", len(lines))
	}
}

// TestMagnetostaticBoundary checks a magnetostatic boundary becomes a zero-permeability
// magnetizable element.
func TestMagnetostaticBoundary(t *testing.T) {
	exc := New(twoElectrodeMesh())
	exc.AddMagnetostaticPotential("a", 1)
	exc.AddMagnetostaticBoundary("b")

	_, types, values := exc.Magnetostatic()
	found := false
	for i, ty := range types {
		if ty == radial.Magnetizable {
			found = true
			if values[i] != 0 {
				t.Errorf("boundary magnetizable value = %g, want 0", values[i])
			}
		}
	}
	if !found {
		t.Error("no Magnetizable element produced for the magnetostatic boundary")
	}
}

// TestElectrostaticGroupIndices checks the name→active-index map covers every active element
// exactly once, in active order, matching the slice Electrostatic returns.
func TestElectrostaticGroupIndices(t *testing.T) {
	exc := New(twoElectrodeMesh())
	exc.AddVoltage("a", 1)
	exc.AddVoltage("b", 2)

	lines, _, _ := exc.Electrostatic()
	groups := exc.ElectrostaticGroupIndices()

	// Every active element index appears in exactly one group, covering [0, len(lines)).
	seen := make([]bool, len(lines))
	for _, idxs := range groups {
		for _, i := range idxs {
			if i < 0 || i >= len(lines) {
				t.Fatalf("group index %d out of range [0,%d)", i, len(lines))
			}
			if seen[i] {
				t.Errorf("active index %d appears in more than one group", i)
			}
			seen[i] = true
		}
	}
	for i, ok := range seen {
		if !ok {
			t.Errorf("active index %d missing from all groups", i)
		}
	}
	// Electrode "a" owns the first two elements, "b" the third.
	if got := groups["a"]; len(got) != 2 {
		t.Errorf("group a = %v, want 2 indices", got)
	}
	if got := groups["b"]; len(got) != 1 {
		t.Errorf("group b = %v, want 1 index", got)
	}
}

// TestAssignUnknownElectrodePanics checks assigning to a non-existent electrode is rejected.
func TestAssignUnknownElectrodePanics(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Error("assigning to an unknown electrode did not panic")
		}
	}()
	New(twoElectrodeMesh()).AddVoltage("nope", 1)
}
