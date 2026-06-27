// SPDX-License-Identifier: MPL-2.0

package ring

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
)

// ringGolden mirrors core/ring/testdata/ring.golden.json (tools/gen_fixtures.py _ring).
type ringGolden struct {
	Samples       [][]oracle.F `json:"samples"` // (r0, z0, delta_r, delta_z)
	Potential     []oracle.F   `json:"potential"`
	Dr1           []oracle.F   `json:"dr1"`
	Dz1           []oracle.F   `json:"dz1"`
	AxialDerivs   [][]oracle.F `json:"axial_derivs"` // [n][9]
	CurSamples    [][]oracle.F `json:"cur_samples"`  // (x0, y0, x, y)
	CurPotential  []oracle.F   `json:"cur_potential"`
	CurField      [][]oracle.F `json:"cur_field"`        // [n][2]
	CurAxialDeriv [][]oracle.F `json:"cur_axial_derivs"` // [n][9]
}

// TestRingAgainstGolden verifies every ring kernel reproduces the upstream C output over
// the sampled grid — direct port-equivalence for the BEM Green's functions.
func TestRingAgainstGolden(t *testing.T) {
	var fx ringGolden
	oracle.LoadGolden(t, "ring", &fx)
	if len(fx.Samples) == 0 {
		t.Fatal("no ring samples loaded")
	}

	for i, s := range fx.Samples {
		r0, z0, dr, dz := s[0].Float(), s[1].Float(), s[2].Float(), s[3].Float()
		oracle.CheckClose(t, "potential", PotentialRadialRing(r0, z0, dr, dz), fx.Potential[i].Float())
		oracle.CheckClose(t, "dr1", Dr1PotentialRadialRing(r0, z0, dr, dz), fx.Dr1[i].Float())
		oracle.CheckClose(t, "dz1", Dz1PotentialRadialRing(r0, z0, dr, dz), fx.Dz1[i].Float())

		got := AxialDerivativesRadialRing(z0, math.Max(dr, 1e-3), z0+dz)
		for k := 0; k < Deriv2DMax; k++ {
			oracle.CheckClose(t, "axial_deriv", got[k], fx.AxialDerivs[i][k].Float())
		}
	}

	for i, s := range fx.CurSamples {
		x0, y0, x, y := s[0].Float(), s[1].Float(), s[2].Float(), s[3].Float()
		oracle.CheckClose(t, "cur_potential", CurrentPotentialAxialRadialRing(y0, x, y), fx.CurPotential[i].Float())
		f := CurrentFieldRadialRing(x0, y0, x, y)
		oracle.CheckClose(t, "cur_field.r", f[0], fx.CurField[i][0].Float())
		oracle.CheckClose(t, "cur_field.z", f[1], fx.CurField[i][1].Float())

		got := CurrentAxialDerivativesRadialRing(y0, x, y)
		for k := 0; k < Deriv2DMax; k++ {
			oracle.CheckClose(t, "cur_axial_deriv", got[k], fx.CurAxialDeriv[i][k].Float())
		}
	}
}

// TestCurrentFieldAtCenter ports test_current_field_at_center: the axial field at the
// centre of a unit current ring of radius R is 1/(2R). Analytic oracle, fixture-independent.
func TestCurrentFieldAtCenter(t *testing.T) {
	for r := 0.5; r <= 5.0; r += 0.25 {
		f := CurrentFieldRadialRing(0, 0, r, 0)
		oracle.CheckClose(t, "B_z(center)", f[1], 1/(2*r))
	}
}

// TestCurrentPotentialAxial ports the closed form -dz/(2*sqrt(dz^2+r^2)) and verifies the
// axial field is its negative derivative (the radial-ring identity the upstream relies on).
func TestCurrentPotentialAxial(t *testing.T) {
	const rRing = 2.0
	for _, z := range []float64{-5, -2, -0.5, 0, 0.5, 2, 5} {
		dz := z
		wantPot := -dz / (2 * math.Sqrt(dz*dz+rRing*rRing))
		oracle.CheckClose(t, "cur_pot_axial", CurrentPotentialAxialRadialRing(z, rRing, 0), wantPot)
		wantField := rRing * rRing / (2 * math.Pow(z*z+rRing*rRing, 1.5))
		f := CurrentFieldRadialRing(0, z, rRing, 0)
		oracle.CheckClose(t, "cur_field_axial", f[1], wantField)
	}
}

// TestFluxDensityFactor checks the dielectric factor 2(K-1)/(K+1): 0 for vacuum (K=1),
// →2 as K→∞.
func TestFluxDensityFactor(t *testing.T) {
	oracle.CheckClose(t, "K=1", FluxDensityToChargeFactor(1), 0.0)
	oracle.CheckClose(t, "K=3", FluxDensityToChargeFactor(3), 1.0)
}
