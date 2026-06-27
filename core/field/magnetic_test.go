// SPDX-License-Identifier: MPL-2.0

package field

import (
	"testing"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/solver"
)

type magneticGolden struct {
	MagLines   [][][]oracle.F `json:"mag_lines"`
	MagCharges []oracle.F     `json:"mag_charges"`
	Current    oracle.F       `json:"current"`
	RRing      oracle.F       `json:"r_ring"`
	Points     [][]oracle.F   `json:"points"`
	CurField   [][]oracle.F   `json:"cur_field"`
	MagField   [][]oracle.F   `json:"mag_field"`
}

// TestMagnetostaticFieldEvaluation verifies FieldRadialBEM's current and total magnetostatic
// field evaluators reproduce the upstream FieldRadialBEM.{current,magnetostatic}_field_at_point.
func TestMagnetostaticFieldEvaluation(t *testing.T) {
	var fx magneticGolden
	oracle.LoadGolden(t, "magnetic_field", &fx)

	magLines := make([]radial.Line, len(fx.MagLines))
	for i := range fx.MagLines {
		for v := 0; v < 4; v++ {
			magLines[i][v] = geom2d.Vertex{fx.MagLines[i][v][0].Float(), fx.MagLines[i][v][1].Float(), fx.MagLines[i][v][2].Float()}
		}
	}
	magJac, magPos := radial.FillJacobianBufferRadial(magLines)
	magCharges := make([]float64, len(fx.MagCharges))
	for i := range magCharges {
		magCharges[i] = fx.MagCharges[i].Float()
	}
	mag := solver.EffectivePointCharges{Charges: magCharges, Jac: magJac, Pos: magPos}

	// Ideal single current ring.
	r := fx.RRing.Float()
	curJac := make(geom3d.CurrentJacobianBuffer, 1)
	curPos := make(geom3d.CurrentPositionBuffer, 1)
	curJac[0][0] = 1.0
	for k := 0; k < geom3d.N_TRIANGLE_QUAD; k++ {
		curPos[0][k] = geom3d.Vec3{r, 0, 0}
	}
	current := solver.CurrentCharges{Currents: []float64{fx.Current.Float()}, Jac: curJac, Pos: curPos}

	f := NewFieldRadialBEMFull(solver.EffectivePointCharges{}, mag, current)
	for m, pf := range fx.Points {
		p := geom2d.Vertex{pf[0].Float(), pf[1].Float(), pf[2].Float()}
		cur := f.CurrentFieldAtPoint(p)
		tot := f.MagnetostaticFieldAtPoint(p)
		for c := 0; c < 3; c++ {
			oracle.CheckClose(t, "cur_field", cur[c], fx.CurField[m][c].Float())
			oracle.CheckClose(t, "mag_field", tot[c], fx.MagField[m][c].Float())
		}
	}
}
