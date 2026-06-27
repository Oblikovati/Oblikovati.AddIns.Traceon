// SPDX-License-Identifier: MPL-2.0

package solver

import (
	"testing"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/radial"
)

// solverGolden mirrors core/solver/testdata/solver.golden.json (gen_fixtures _solver).
type solverGolden struct {
	Lines         [][][]oracle.F `json:"lines"` // [n][4][3]
	Types         []int          `json:"types"`
	Values        []oracle.F     `json:"values"`
	SelfPotential []oracle.F     `json:"self_potential"`
	SelfField     []oracle.F     `json:"self_field"`
	Matrix        [][]oracle.F   `json:"matrix"`
	RHS           []oracle.F     `json:"rhs"`
	Charges       []oracle.F     `json:"charges"`
}

func (g *solverGolden) lines() []radial.Line {
	out := make([]radial.Line, len(g.Lines))
	for i := range g.Lines {
		for v := 0; v < 4; v++ {
			out[i][v] = geom2d.Vertex{g.Lines[i][v][0].Float(), g.Lines[i][v][1].Float(), g.Lines[i][v][2].Float()}
		}
	}
	return out
}

func loadSolverGolden(t *testing.T) (solverGolden, []radial.Line, []radial.ExcitationType, []float64) {
	t.Helper()
	var fx solverGolden
	oracle.LoadGolden(t, "solver", &fx)
	lines := fx.lines()
	types := make([]radial.ExcitationType, len(fx.Types))
	values := make([]float64, len(fx.Values))
	for i := range types {
		types[i] = radial.ExcitationType(fx.Types[i])
		values[i] = fx.Values[i].Float()
	}
	return fx, lines, types, values
}

// TestSelfTermsMatchOracle is the critical check: the integrated singular self-term diagonal
// entries reproduce the upstream scipy-quad values. These gate the whole solve.
func TestSelfTermsMatchOracle(t *testing.T) {
	fx, lines, _, values := loadSolverGolden(t)
	for i := range lines {
		oracle.CheckClose(t, "self_potential", SelfPotentialRadial(lines[i]), fx.SelfPotential[i].Float())
		oracle.CheckClose(t, "self_field", SelfFieldDotNormalRadial(lines[i], values[i]), fx.SelfField[i].Float())
	}
}

// TestAssembleMatrix verifies the full assembled matrix (off-diagonal + singular diagonal)
// reproduces the upstream matrix.
func TestAssembleMatrix(t *testing.T) {
	fx, lines, types, values := loadSolverGolden(t)
	jac, pos := radial.FillJacobianBufferRadial(lines)
	m := AssembleMatrix(lines, types, values, jac, pos)
	n := len(lines)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			oracle.CheckClose(t, "matrix", m.At(i, j), fx.Matrix[i][j].Float())
		}
	}
}

// TestSolveElectrostatic is the end-to-end check: the solved effective charges reproduce
// the upstream np.linalg.solve result for the assembled system.
func TestSolveElectrostatic(t *testing.T) {
	fx, lines, types, values := loadSolverGolden(t)
	epc, err := SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("SolveElectrostatic: %v", err)
	}
	if len(epc.Charges) != len(fx.Charges) {
		t.Fatalf("got %d charges, want %d", len(epc.Charges), len(fx.Charges))
	}
	for i := range epc.Charges {
		oracle.CheckClose(t, "charge", epc.Charges[i], fx.Charges[i].Float())
	}
}

// TestSolveEmpty checks the degenerate empty-geometry case returns no charges, no panic.
func TestSolveEmpty(t *testing.T) {
	epc, err := SolveElectrostatic(nil, nil, nil)
	if err != nil {
		t.Fatalf("empty solve: %v", err)
	}
	if len(epc.Charges) != 0 {
		t.Errorf("expected 0 charges, got %d", len(epc.Charges))
	}
}
