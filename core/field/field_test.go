// SPDX-License-Identifier: MPL-2.0

package field

import (
	"testing"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/solver"
)

type fieldGolden struct {
	Lines       [][][]oracle.F `json:"lines"`
	Charges     []oracle.F     `json:"charges"`
	Z           []oracle.F     `json:"z"`
	Derivs      [][]oracle.F   `json:"derivs"`      // (N, 9)
	Coeffs      [][][]oracle.F `json:"coeffs"`      // (N-1, 9, 6)
	EvalPoints  [][]oracle.F   `json:"eval_points"` // (m, 3)
	PotDirect   []oracle.F     `json:"pot_direct"`
	FieldDirect [][]oracle.F   `json:"field_direct"`
	PotInterp   []oracle.F     `json:"pot_interp"`
	FieldInterp [][]oracle.F   `json:"field_interp"`
}

func (g *fieldGolden) lines() []radial.Line {
	out := make([]radial.Line, len(g.Lines))
	for i := range g.Lines {
		for v := 0; v < 4; v++ {
			out[i][v] = geom2d.Vertex{g.Lines[i][v][0].Float(), g.Lines[i][v][1].Float(), g.Lines[i][v][2].Float()}
		}
	}
	return out
}

func (g *fieldGolden) floats(fs []oracle.F) []float64 {
	out := make([]float64, len(fs))
	for i, f := range fs {
		out[i] = f.Float()
	}
	return out
}

func (g *fieldGolden) point(p []oracle.F) geom2d.Vertex {
	return geom2d.Vertex{p[0].Float(), p[1].Float(), p[2].Float()}
}

func load(t *testing.T) (fieldGolden, solver.EffectivePointCharges) {
	t.Helper()
	var fx fieldGolden
	oracle.LoadGolden(t, "field", &fx)
	lines := fx.lines()
	jac, pos := radial.FillJacobianBufferRadial(lines)
	return fx, solver.EffectivePointCharges{Charges: fx.floats(fx.Charges), Jac: jac, Pos: pos}
}

// TestAxialDerivatives verifies the accumulated on-axis derivatives match upstream.
func TestAxialDerivatives(t *testing.T) {
	fx, epc := load(t)
	z := fx.floats(fx.Z)
	got := radial.AxialDerivativesRadial(z, epc.Charges, epc.Jac, epc.Pos)
	for i := range z {
		for l := 0; l < radial.Deriv2DMax; l++ {
			oracle.CheckClose(t, "deriv", got[i][l], fx.Derivs[i][l].Float())
		}
	}
}

// TestQuinticSplineCoefficients verifies the assembled interpolation coefficients match
// upstream _quintic_spline_coefficients (white-box: combines core/interp quintic + cubic).
func TestQuinticSplineCoefficients(t *testing.T) {
	fx, epc := load(t)
	z := fx.floats(fx.Z)
	derivs := radial.AxialDerivativesRadial(z, epc.Charges, epc.Jac, epc.Pos)
	coeffs, err := quinticSplineCoefficients(z, derivs)
	if err != nil {
		t.Fatalf("quinticSplineCoefficients: %v", err)
	}
	for iv := range coeffs {
		for order := 0; order < radial.Deriv2DMax; order++ {
			for c := 0; c < 6; c++ {
				oracle.CheckClose(t, "coeff", coeffs[iv][order][c], fx.Coeffs[iv][order][c].Float())
			}
		}
	}
}

// TestFieldRadialBEMDirect verifies the direct (boundary-integral) potential and field.
func TestFieldRadialBEMDirect(t *testing.T) {
	fx, epc := load(t)
	f := NewFieldRadialBEM(epc)
	for m, ep := range fx.EvalPoints {
		p := fx.point(ep)
		oracle.CheckClose(t, "pot_direct", f.PotentialAtPoint(p), fx.PotDirect[m].Float())
		fv := f.FieldAtPoint(p)
		for c := 0; c < 3; c++ {
			oracle.CheckClose(t, "field_direct", fv[c], fx.FieldDirect[m][c].Float())
		}
	}
}

// TestFieldRadialAxialInterp verifies the axial-series interpolated potential and field,
// including the zero returned outside [zmin, zmax].
func TestFieldRadialAxialInterp(t *testing.T) {
	fx, epc := load(t)
	zmin, zmax, n := fx.Z[0].Float(), fx.Z[len(fx.Z)-1].Float(), len(fx.Z)
	fa, err := NewFieldRadialAxial(epc, zmin, zmax, n)
	if err != nil {
		t.Fatalf("NewFieldRadialAxial: %v", err)
	}
	for m, ep := range fx.EvalPoints {
		p := fx.point(ep)
		oracle.CheckClose(t, "pot_interp", fa.PotentialAtPoint(p), fx.PotInterp[m].Float())
		fv := fa.FieldAtPoint(p)
		for c := 0; c < 3; c++ {
			oracle.CheckClose(t, "field_interp", fv[c], fx.FieldInterp[m][c].Float())
		}
	}
}
