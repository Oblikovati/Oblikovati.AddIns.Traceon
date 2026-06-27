// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/solver"
)

const sqrt2 = math.Sqrt2

func close(a, b float64) bool { return math.Abs(a-b) <= 1e-7+1e-7*math.Abs(b) }

func closeVec(a, b Point) bool {
	return close(a[0], b[0]) && close(a[1], b[1]) && close(a[2], b[2])
}

// TestArc ports geometric_tests.test_arc: arc length is a quarter circumference and the
// 1/8-circumference sample lands on the 45° point, for the plain and reversed sweeps.
func TestArc(t *testing.T) {
	r := 2.0
	p := Arc(Point{0, 0, 0}, Point{r, 0, 0}, Point{0, r, 0}, false)
	if !close(p.Length, 0.25*2*math.Pi*r) {
		t.Errorf("arc length = %g, want %g", p.Length, 0.25*2*math.Pi*r)
	}
	if got := p.At(0.125 * 2 * math.Pi * r); !closeVec(got, Point{r / sqrt2, r / sqrt2, 0}) {
		t.Errorf("arc(1/8) = %v, want [%g %g 0]", got, r/sqrt2, r/sqrt2)
	}

	// Reversed quarter arc takes the long way: 3/4 of the circle.
	r = 3
	center := Point{0, r, 0}
	rev := Arc(center, Point{0, 0, 0}, Point{0, r, r}, true)
	if !close(rev.Length, 0.75*2*math.Pi*r) {
		t.Errorf("reversed arc length = %g, want %g", rev.Length, 0.75*2*math.Pi*r)
	}
	want := add(center, Point{0, r / sqrt2, -r / sqrt2})
	if got := rev.At(0.5 * 0.75 * 2 * math.Pi * r); !closeVec(got, want) {
		t.Errorf("reversed arc(mid) = %v, want %v", got, want)
	}
}

// TestLineEndpoints checks the straight-line builder's endpoints and arc-length midpoint.
func TestLineEndpoints(t *testing.T) {
	p := Line(Point{1, 0, 0}, Point{1, 0, 2})
	if !close(p.Length, 2) {
		t.Errorf("length = %g, want 2", p.Length)
	}
	if got := p.MiddlePoint(); !closeVec(got, Point{1, 0, 1}) {
		t.Errorf("middle = %v, want [1 0 1]", got)
	}
}

// TestRectangleClosesAndBreakpoints checks the rectangle is a closed loop (start==end) with
// the three interior corners recorded as breakpoints (the fourth corner is the start).
func TestRectangleClosesAndBreakpoints(t *testing.T) {
	p := RectangleXZ(0.5, 1.0, -0.5, 0.5)
	if !closeVec(p.StartingPoint(), p.Endpoint()) {
		t.Errorf("rectangle not closed: start %v end %v", p.StartingPoint(), p.Endpoint())
	}
	if len(p.Breakpoints) != 3 {
		t.Errorf("breakpoints = %v, want 3 corners", p.Breakpoints)
	}
}

// TestCircleXZ checks the full-circle builder: circumference length and the quarter-turn
// sample, in the meridian (x, z) plane (y stays 0).
func TestCircleXZ(t *testing.T) {
	c := CircleXZ(0, 0, 2, 2*math.Pi)
	if !close(c.Length, 2*math.Pi*2) {
		t.Errorf("circle length = %g, want %g", c.Length, 2*math.Pi*2)
	}
	if got := c.At(0.25 * c.Length); !closeVec(got, Point{0, 0, 2}) {
		t.Errorf("circle(1/4) = %v, want [0 0 2]", got)
	}
}

// TestReverse checks Reverse swaps the endpoints and mirrors breakpoints.
func TestReverse(t *testing.T) {
	p := Line(Point{1, 0, 0}, Point{1, 0, 2}).
		ExtendWithLine(Point{2, 0, 2}).Reverse()
	if !closeVec(p.StartingPoint(), Point{2, 0, 2}) {
		t.Errorf("reversed start = %v, want [2 0 2]", p.StartingPoint())
	}
	if !closeVec(p.Endpoint(), Point{1, 0, 0}) {
		t.Errorf("reversed end = %v, want [1 0 0]", p.Endpoint())
	}
	if len(p.Breakpoints) != 1 {
		t.Errorf("breakpoints = %v, want 1", p.Breakpoints)
	}
}

// TestExtendWithArcAndName checks an arc extension joins smoothly and WithName tags the path.
func TestExtendWithArcAndName(t *testing.T) {
	p := Line(Point{0, 0, 0}, Point{2, 0, 0}).
		ExtendWithArc(Point{2, 0, 2}, Point{4, 0, 2}, false).
		WithName("lens")
	if p.Name != "lens" {
		t.Errorf("name = %q, want lens", p.Name)
	}
	if !closeVec(p.Endpoint(), Point{4, 0, 2}) {
		t.Errorf("arc-extended end = %v, want [4 0 2]", p.Endpoint())
	}
}

// TestDiscretizeMembership ports test_discretize_path: 0, the length, and every breakpoint
// appear in the sample parameters.
func TestDiscretizeMembership(t *testing.T) {
	u := discretizePath(10, []float64{3.33, 5, 9}, 1.0, 0, 1)
	for _, want := range []float64{0, 10, 3.33, 5, 9} {
		if !contains(u, want) {
			t.Errorf("sample %g missing from %v", want, u)
		}
	}
}

func contains(xs []float64, v float64) bool {
	for _, x := range xs {
		if close(x, v) {
			return true
		}
	}
	return false
}

// TestMeshGroupNamedGroups checks MeshGroup meshes several named paths into one mesh with a
// physical group per path, and that the shared junction node between abutting paths is merged.
func TestMeshGroupNamedGroups(t *testing.T) {
	lower := Line(Point{1, 0, 0}, Point{1, 0, 1}).WithName("lower")
	upper := Line(Point{1, 0, 1}, Point{1, 0, 2}).WithName("upper") // starts where lower ends
	m := MeshGroup([]Path{lower, upper}, MeshOptions{MeshSize: 0.5})

	if _, ok := m.PhysicalToLines["lower"]; !ok {
		t.Errorf("missing physical group lower; have %v", m.PhysicalToLines)
	}
	if _, ok := m.PhysicalToLines["upper"]; !ok {
		t.Errorf("missing physical group upper; have %v", m.PhysicalToLines)
	}
	// Every line is assigned to exactly one group (no leaks, no double-counting).
	total := len(m.PhysicalToLines["lower"]) + len(m.PhysicalToLines["upper"])
	if total != len(m.Lines) {
		t.Errorf("group line count %d != total lines %d", total, len(m.Lines))
	}
	// Each length-1 path meshes to the 3-element floor → 4 nodes; the shared junction node
	// (1,0,1) is merged across the two paths, so 4+4−1 = 7 points, not 8.
	if len(m.Points) != 7 {
		t.Errorf("points = %d, want 7 (junction node merged across paths)", len(m.Points))
	}
}

// TestRadialLinesSolves checks the meshed parametric electrode feeds the radial solver: a
// charged aperture electrode produces a non-trivial surface-charge solution.
func TestRadialLinesSolves(t *testing.T) {
	m := Aperture(0.5, 0.3, 1.5, 0.0).Mesh(MeshOptions{MeshSize: 0.2, HigherOrder: true, Name: "ap", EnsureOutwardNormals: true})
	lines := RadialLines(m)
	if len(lines) == 0 {
		t.Fatal("aperture meshed to zero lines")
	}

	// Bias every element to 1 V and solve the electrostatic BEM.
	types := make([]radial.ExcitationType, len(lines))
	values := make([]float64, len(lines))
	for i := range lines {
		types[i] = radial.VoltageFixed
		values[i] = 1.0
	}
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve: %v", err)
	}
	if len(charges.Charges) != len(lines) {
		t.Errorf("charges = %d, want %d", len(charges.Charges), len(lines))
	}
	nonzero := false
	for _, c := range charges.Charges {
		if c != 0 {
			nonzero = true
		}
	}
	if !nonzero {
		t.Error("all surface charges zero — electrode did not respond to its bias")
	}
}
