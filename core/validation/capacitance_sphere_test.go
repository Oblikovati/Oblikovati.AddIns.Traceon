// SPDX-License-Identifier: MPL-2.0

package validation

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/excitation"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geometry"
	"oblikovati.org/traceon/core/solver"
)

// Concentric spherical capacitor with a dielectric layer, from validation/capacitance_sphere.py.
// Inner sphere at 1 V, outer at 0 V, with a dielectric shell (εr=2) between r3 and r4; the value
// of interest is the capacitance, integrated from the surface charge on the two conductors.
const (
	csR1  = 0.5 // inner conductor
	csR2  = 1.0 // outer conductor
	csR3  = 0.6 // dielectric inner surface
	csR4  = 0.9 // dielectric outer surface
	csK   = 2.0 // dielectric relative permittivity
	csMSF = 10
)

// csUpstreamMSF10 is the capacitance Traceon computes at MSF=10 (higher-order) — the
// port-equivalence oracle (rel err to the analytic value ~3e-8 at this resolution).
const csUpstreamMSF10 = 17.3995894992897

// csAnalytic is the exact capacitance 4π / ((1/r1−1/r3) + (1/r3−1/r4)/K + (1/r4−1/r2)).
func csAnalytic() float64 {
	return 4 * math.Pi / ((1/csR1 - 1/csR3) + (1/csR3-1/csR4)/csK + (1/csR4 - 1/csR2))
}

// csShell builds one spherical shell as two quarter arcs (−z pole → +r equator → +z pole).
func csShell(radius float64, name string) geometry.Path {
	return geometry.Arc(geometry.Point{0, 0, 0}, geometry.Point{0, 0, -radius}, geometry.Point{radius, 0, 0}, false).
		ExtendWithArc(geometry.Point{0, 0, 0}, geometry.Point{0, 0, radius}, false).
		WithName(name)
}

// TestCapacitanceSphere reproduces the dielectric spherical capacitor: solve the two conductors
// plus the dielectric layer, integrate the surface charge on each conductor, and form the
// capacitance. This exercises a real dielectric (εr=2, not a boundary) and charge integration.
func TestCapacitanceSphere(t *testing.T) {
	m := geometry.MeshGroup(
		[]geometry.Path{csShell(csR1, "inner"), csShell(csR2, "outer"),
			csShell(csR3, "dielectric_inner"), csShell(csR4, "dielectric_outer")},
		geometry.MeshOptions{MeshSizeFactor: csMSF, HigherOrder: true, EnsureOutwardNormals: true})
	// The dielectric's inner surface faces inward (the outer surface keeps the outward default).
	m.EnsureInwardNormals("dielectric_inner")

	exc := excitation.New(m)
	exc.AddVoltage("inner", 1)
	exc.AddVoltage("outer", 0)
	exc.AddDielectric("dielectric_inner", csK)
	exc.AddDielectric("dielectric_outer", csK)

	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve electrostatic: %v", err)
	}

	bem := field.NewFieldRadialBEM(charges)
	groups := exc.ElectrostaticGroupIndices()
	qInner := bem.ChargeOnElements(groups["inner"])
	qOuter := bem.ChargeOnElements(groups["outer"])
	capacitance := (math.Abs(qOuter) + math.Abs(qInner)) / 2

	t.Logf("capacitance = %.12g (upstream %.12g, analytic %.12g)", capacitance, csUpstreamMSF10, csAnalytic())
	if rel := math.Abs(capacitance-csUpstreamMSF10) / csUpstreamMSF10; rel > 1e-5 {
		t.Errorf("capacitance = %.12g, want %.12g (upstream MSF=10); rel err %.2e > 1e-5", capacitance, csUpstreamMSF10, rel)
	}
	if rel := math.Abs(capacitance-csAnalytic()) / csAnalytic(); rel > 1e-3 {
		t.Errorf("capacitance = %.12g vs analytic %.12g; rel err %.2e > 1e-3", capacitance, csAnalytic(), rel)
	}
}
