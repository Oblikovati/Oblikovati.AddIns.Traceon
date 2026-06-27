// SPDX-License-Identifier: MPL-2.0

package validation

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/excitation"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geometry"
	"oblikovati.org/traceon/core/solver"
)

// Edwards 2007 geometry g5, from validation/edwards2007.py (D. Edwards, "High precision
// electrostatic potential calculations for cylindrically symmetric lenses", 2007). Two coaxial
// stepped cylinders at 0 V and 10 V; the value of interest is the potential at (r=12, z=4).
const edw2007MSF = 32

// edw2007PaperValue is the high-precision potential from the paper.
const edw2007PaperValue = 6.69099430708

// edw2007UpstreamMSF32 is the potential Traceon computes at MSF=32 (higher-order) — the
// port-equivalence oracle (rel err to the paper ~6.5e-4 at this resolution).
const edw2007UpstreamMSF32 = 6.686658832380274

// TestEdwards2007 reproduces the stepped-cylinder potential: solve the two electrodes and read
// the potential at (12, 4). A pure electrostatic-potential check over a multi-segment geometry.
func TestEdwards2007(t *testing.T) {
	p := func(x, z float64) geometry.Point { return geometry.Point{x, 0, z} }
	inner := geometry.Line(p(0, 5), p(12, 5)).ExtendWithLine(p(12, 15)).ExtendWithLine(p(0, 15)).WithName("inner")
	boundary := geometry.Line(p(0, 0), p(20, 0)).ExtendWithLine(p(20, 20)).ExtendWithLine(p(0, 20)).WithName("boundary")

	m := geometry.MeshGroup([]geometry.Path{inner, boundary},
		geometry.MeshOptions{MeshSizeFactor: edw2007MSF, HigherOrder: true, EnsureOutwardNormals: true})

	exc := excitation.New(m)
	exc.AddVoltage("boundary", 0)
	exc.AddVoltage("inner", 10)

	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve electrostatic: %v", err)
	}

	pot := field.NewFieldRadialBEM(charges).PotentialAtPoint(geom2d.Vertex{12, 0, 4})

	t.Logf("potential = %.12g (upstream %.12g, paper %.12g)", pot, edw2007UpstreamMSF32, edw2007PaperValue)
	if rel := math.Abs(pot-edw2007UpstreamMSF32) / edw2007UpstreamMSF32; rel > 1e-5 {
		t.Errorf("potential = %.12g, want %.12g (upstream MSF=32); rel err %.2e > 1e-5", pot, edw2007UpstreamMSF32, rel)
	}
	if rel := math.Abs(pot-edw2007PaperValue) / edw2007PaperValue; rel > 2e-3 {
		t.Errorf("potential = %.12g vs paper %.12g; rel err %.2e > 2e-3", pot, edw2007PaperValue, rel)
	}
}
