// SPDX-License-Identifier: MPL-2.0

package solver

import (
	"testing"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/radial"
)

type magnetostaticGolden struct {
	Lines   [][][]oracle.F `json:"lines"`
	Types   []int          `json:"types"`
	Values  []oracle.F     `json:"values"`
	Current oracle.F       `json:"current"`
	RRing   oracle.F       `json:"r_ring"`
	Matrix  [][]oracle.F   `json:"matrix"`
	RHS     []oracle.F     `json:"rhs"`
	Charges []oracle.F     `json:"charges"`
}

// ringCurrentCharges builds the ideal single current ring (jac=[1,0,…], pos=12×[r,0,0]).
func ringCurrentCharges(current, r float64) CurrentCharges {
	jac := make(geom3d.CurrentJacobianBuffer, 1)
	pos := make(geom3d.CurrentPositionBuffer, 1)
	jac[0][0] = 1.0
	for k := 0; k < geom3d.N_TRIANGLE_QUAD; k++ {
		pos[0][k] = geom3d.Vec3{r, 0, 0}
	}
	return CurrentCharges{Currents: []float64{current}, Jac: jac, Pos: pos}
}

// TestSolveMagnetostatic verifies the magnetostatic solve end to end: a magnetizable
// element and a magnetic-scalar-potential element responding to a current ring's pre-field.
// The assembled matrix, the right-hand side (incl. the pre-field flux), and the solved
// charges must all match the upstream MagnetostaticSolverRadial computation.
func TestSolveMagnetostatic(t *testing.T) {
	var fx magnetostaticGolden
	oracle.LoadGolden(t, "magnetostatic", &fx)

	lines := make([]radial.Line, len(fx.Lines))
	for i := range fx.Lines {
		for v := 0; v < 4; v++ {
			lines[i][v] = geom2d.Vertex{fx.Lines[i][v][0].Float(), fx.Lines[i][v][1].Float(), fx.Lines[i][v][2].Float()}
		}
	}
	types := make([]radial.ExcitationType, len(fx.Types))
	values := make([]float64, len(fx.Values))
	for i := range types {
		types[i] = radial.ExcitationType(fx.Types[i])
		values[i] = fx.Values[i].Float()
	}

	current := ringCurrentCharges(fx.Current.Float(), fx.RRing.Float())
	preField := current.CurrentFieldAt

	// Matrix and RHS.
	jac, pos := radial.FillJacobianBufferRadial(lines)
	m := AssembleMatrix(lines, types, values, jac, pos)
	n := len(lines)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			oracle.CheckClose(t, "matrix", m.At(i, j), fx.Matrix[i][j].Float())
		}
	}
	rhs := rightHandSideMagnetostatic(lines, types, values, preField)
	for i := range rhs {
		oracle.CheckClose(t, "rhs", rhs[i], fx.RHS[i].Float())
	}

	// Full solve.
	epc, err := SolveMagnetostatic(lines, types, values, preField)
	if err != nil {
		t.Fatalf("SolveMagnetostatic: %v", err)
	}
	for i := range epc.Charges {
		oracle.CheckClose(t, "charge", epc.Charges[i], fx.Charges[i].Float())
	}
}

// TestSolveMagnetostaticNoPreField checks a current-free magnetostatic solve (nil preField):
// a lone magnetizable element gets a zero right-hand side, so its charge is zero.
func TestSolveMagnetostaticNoPreField(t *testing.T) {
	line := radial.Line{{0.5, 0, -0.2}, {0.5, 0, 0.2}, {0.5, 0, -0.2 + 0.4/3}, {0.5, 0, -0.2 + 0.8/3}}
	epc, err := SolveMagnetostatic([]radial.Line{line}, []radial.ExcitationType{radial.Magnetizable}, []float64{500}, nil)
	if err != nil {
		t.Fatalf("SolveMagnetostatic: %v", err)
	}
	oracle.CheckClose(t, "charge", epc.Charges[0], 0.0)
}
