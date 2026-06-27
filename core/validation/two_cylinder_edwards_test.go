// SPDX-License-Identifier: MPL-2.0

// Package validation reproduces Traceon's /validation cases end-to-end in pure Go and checks
// them against their published reference values. Unlike the per-module oracle tests, these
// exercise the whole stack as one unit — parametric geometry → mesh → named excitation →
// BEM solve → field evaluation — so a regression anywhere shows up against a physics oracle.
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

// Two-cylinder electrostatic lens geometry, from validation/two_cylinder_edwards.py.
const (
	tcGapSize        = 0.2
	tcRadius         = 1.0
	tcBoundaryLength = 20.0
	tcMSF            = 16 // upstream mesh_size_factor; gives 240 higher-order line elements
)

// tcCylinderLength is the length of each tube: (boundary − gap) / 2 = 9.9.
const tcCylinderLength = (tcBoundaryLength - tcGapSize) / 2

// tcPaperValue is the reference potential at (r=0, z=9.6) from Edwards Jr. 2007, "Accurate
// Potential Calculations For The Two Tube Electrostatic Lens Using A Multiregion FDM Method".
const tcPaperValue = 2.5966375108359858

// tcUpstreamMSF16 is the potential Traceon itself computes at MSF=16, higher_order — the
// port-equivalence oracle (regenerate with validation/two_cylinder_edwards.py if the mesher or
// solver changes). Traceon's own relative error to the paper value at this resolution is 1.4e-3.
const tcUpstreamMSF16 = 2.5930493947741513

// tcGapVoltage ramps linearly from 0 V at the bottom of the gap (z=9.9) to 10 V at the top
// (z=10.1): (z − 9.9)/0.2 · 10. Port of two_cylinder_edwards.gap_voltage.
func tcGapVoltage(_, _, z float64) float64 { return (z - 9.9) / tcGapSize * 10 }

// TestTwoCylinderEdwards reproduces the two-cylinder lens potential and checks it against both
// Traceon's own computed value (tight, port equivalence) and the paper value (loose, physics).
func TestTwoCylinderEdwards(t *testing.T) {
	bottom := geometry.Line(geometry.Point{0, 0, 0}, geometry.Point{tcRadius, 0, 0}).
		ExtendWithLine(geometry.Point{tcRadius, 0, tcCylinderLength}).WithName("v1")
	gap := geometry.Line(geometry.Point{tcRadius, 0, tcCylinderLength},
		geometry.Point{tcRadius, 0, tcCylinderLength + tcGapSize}).WithName("gap")
	top := geometry.Line(geometry.Point{tcRadius, 0, tcCylinderLength + tcGapSize},
		geometry.Point{tcRadius, 0, tcBoundaryLength}).
		ExtendWithLine(geometry.Point{0, 0, tcBoundaryLength}).WithName("v2")

	m := geometry.MeshGroup([]geometry.Path{bottom, gap, top},
		geometry.MeshOptions{MeshSizeFactor: tcMSF, HigherOrder: true})

	exc := excitation.New(m)
	exc.AddVoltage("v1", 0)
	exc.AddVoltage("v2", 10)
	exc.AddVoltageFunc("gap", tcGapVoltage)

	lines, types, values := exc.Electrostatic()
	if len(lines) != 240 {
		t.Errorf("active line elements = %d, want 240 (matches upstream at MSF=16)", len(lines))
	}

	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve electrostatic: %v", err)
	}

	pot := field.NewFieldRadialBEM(charges).PotentialAtPoint(geom2d.Vertex{0, 0, tcBoundaryLength/2 - 0.4})

	// Port equivalence: Go must reproduce Traceon's own number to np.isclose tolerance.
	if rel := math.Abs(pot-tcUpstreamMSF16) / math.Abs(tcUpstreamMSF16); rel > 1e-5 {
		t.Errorf("potential = %.16g, want %.16g (upstream MSF=16); rel err %.2e > 1e-5",
			pot, tcUpstreamMSF16, rel)
	}
	// Physics validation: within the BEM's discretization error of the paper value.
	if rel := math.Abs(pot-tcPaperValue) / math.Abs(tcPaperValue); rel > 2e-3 {
		t.Errorf("potential = %.16g vs paper %.16g; rel err %.2e > 2e-3", pot, tcPaperValue, rel)
	}
}
