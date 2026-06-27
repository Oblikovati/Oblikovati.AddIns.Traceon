// SPDX-License-Identifier: MPL-2.0

package radial

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/internal/oracle"
)

// radialGolden mirrors core/radial/testdata/radial.golden.json (gen_fixtures _radial).
type radialGolden struct {
	Lines         [][][]oracle.F `json:"lines"` // [n][4][3]
	Jac           [][]oracle.F   `json:"jac"`   // [n][16]
	Pos           [][][]oracle.F `json:"pos"`   // [n][16][2]
	Charges       []oracle.F     `json:"charges"`
	ChargeRadial  []oracle.F     `json:"charge_radial"`
	EvalPoints    [][]oracle.F   `json:"eval_points"` // [m][3]
	Potential     []oracle.F     `json:"potential"`
	Field         [][]oracle.F   `json:"field"` // [m][3]
	Alphas        []oracle.F     `json:"alphas"`
	K             oracle.F       `json:"K"`
	SelfPotential [][]oracle.F   `json:"self_potential"` // [n][len(alphas)]
	SelfField     [][]oracle.F   `json:"self_field"`     // [n][len(alphas)]
	ExcTypes      []int          `json:"exc_types"`
	ExcValues     []oracle.F     `json:"exc_values"`
	Matrix        [][]oracle.F   `json:"matrix"` // [n][n]
}

func (g *radialGolden) line(i int) Line {
	var l Line
	for v := 0; v < 4; v++ {
		l[v] = geom2d.Vertex{g.Lines[i][v][0].Float(), g.Lines[i][v][1].Float(), g.Lines[i][v][2].Float()}
	}
	return l
}

func (g *radialGolden) lines() []Line {
	out := make([]Line, len(g.Lines))
	for i := range g.Lines {
		out[i] = g.line(i)
	}
	return out
}

func (g *radialGolden) chargesF() []float64 {
	out := make([]float64, len(g.Charges))
	for i, c := range g.Charges {
		out[i] = c.Float()
	}
	return out
}

func (g *radialGolden) point(p []oracle.F) geom2d.Vertex {
	return geom2d.Vertex{p[0].Float(), p[1].Float(), p[2].Float()}
}

func TestRadialAgainstGolden(t *testing.T) {
	var fx radialGolden
	oracle.LoadGolden(t, "radial", &fx)
	lines := fx.lines()
	if len(lines) == 0 {
		t.Fatal("no lines loaded")
	}

	// Jacobian/position buffers.
	jac, pos := FillJacobianBufferRadial(lines)
	for i := range lines {
		for k := 0; k < NQuad2D; k++ {
			oracle.CheckClose(t, "jac", jac[i][k], fx.Jac[i][k].Float())
			oracle.CheckClose(t, "pos.r", pos[i][k][0], fx.Pos[i][k][0].Float())
			oracle.CheckClose(t, "pos.z", pos[i][k][1], fx.Pos[i][k][1].Float())
		}
	}

	// Charge per element.
	for i := range lines {
		oracle.CheckClose(t, "charge_radial", ChargeRadial(lines[i], 1.0), fx.ChargeRadial[i].Float())
	}

	// Potential & field evaluation from effective point charges.
	charges := fx.chargesF()
	for m, ep := range fx.EvalPoints {
		p := fx.point(ep)
		oracle.CheckClose(t, "potential", PotentialRadial(p, charges, jac, pos), fx.Potential[m].Float())
		f := FieldRadial(p, charges, jac, pos)
		for c := 0; c < 3; c++ {
			oracle.CheckClose(t, "field", f[c], fx.Field[m][c].Float())
		}
	}

	// Singular self-term integrands at sampled α.
	for i := range lines {
		for a, alpha := range fx.Alphas {
			oracle.CheckClose(t, "self_potential", SelfPotentialRadialIntegrand(alpha.Float(), lines[i]), fx.SelfPotential[i][a].Float())
			oracle.CheckClose(t, "self_field", SelfFieldDotNormalRadialIntegrand(alpha.Float(), lines[i], fx.K.Float()), fx.SelfField[i][a].Float())
		}
	}

	// Dense matrix assembly.
	n := len(lines)
	types := make([]ExcitationType, n)
	values := make([]float64, n)
	for i := range types {
		types[i] = ExcitationType(fx.ExcTypes[i])
		values[i] = fx.ExcValues[i].Float()
	}
	matrix := make([]float64, n*n)
	FillMatrixRadial(matrix, lines, types, values, jac, pos, 0, n-1)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			oracle.CheckClose(t, "matrix", matrix[i*n+j], fx.Matrix[i][j].Float())
		}
	}
}

// TestChargeRadialVertical ports test_charge_radial_vertical: a unit-density vertical line
// element of length 1 at r=1 carries total charge 2π (= surface area of the revolved band).
func TestChargeRadialVertical(t *testing.T) {
	// GMSH line4 ordering: [start, end, 1/3, 2/3].
	line := Line{{1, 0, 0}, {1, 0, 1}, {1, 0, 1.0 / 3}, {1, 0, 2.0 / 3}}
	oracle.CheckClose(t, "charge", ChargeRadial(line, 1.0), 2*math.Pi)
}
