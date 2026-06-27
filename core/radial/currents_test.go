// SPDX-License-Identifier: MPL-2.0

package radial

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/constants"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/quad"
)

type currentsGolden struct {
	Tris          [][][]oracle.F `json:"tris"`
	Jac3D         [][]oracle.F   `json:"jac3d"`
	Pos3D         [][][]oracle.F `json:"pos3d"`
	Current       oracle.F       `json:"current"`
	R             oracle.F       `json:"r"`
	FieldPoints   [][]oracle.F   `json:"field_points"`
	CurField      [][]oracle.F   `json:"cur_field"`
	ZAxis         []oracle.F     `json:"z_axis"`
	CurPotential  []oracle.F     `json:"cur_potential"`
	CurAxialDeriv [][]oracle.F   `json:"cur_axial_derivs"`
}

// ringEffectiveCharges builds the ideal single current ring effective point charge
// (get_ring_effective_point_charges): all quadrature weight in the first point, every
// position the ring point (r, 0, 0).
func ringEffectiveCharges(current, r float64) ([]float64, geom3d.CurrentJacobianBuffer, geom3d.CurrentPositionBuffer) {
	jac := make(geom3d.CurrentJacobianBuffer, 1)
	pos := make(geom3d.CurrentPositionBuffer, 1)
	jac[0][0] = 1.0
	for k := 0; k < geom3d.N_TRIANGLE_QUAD; k++ {
		pos[0][k] = geom3d.Vec3{r, 0, 0}
	}
	return []float64{current}, jac, pos
}

func TestCurrentsAgainstGolden(t *testing.T) {
	var fx currentsGolden
	oracle.LoadGolden(t, "currents", &fx)

	// Triangle jacobian buffer.
	tris := make([]geom3d.Triangle, len(fx.Tris))
	for i := range fx.Tris {
		for v := 0; v < 3; v++ {
			tris[i][v] = geom3d.Vec3{fx.Tris[i][v][0].Float(), fx.Tris[i][v][1].Float(), fx.Tris[i][v][2].Float()}
		}
	}
	jac3d, pos3d := geom3d.FillJacobianBuffer3D(tris)
	for i := range tris {
		for k := 0; k < geom3d.N_TRIANGLE_QUAD; k++ {
			oracle.CheckClose(t, "jac3d", jac3d[i][k], fx.Jac3D[i][k].Float())
			for c := 0; c < 3; c++ {
				oracle.CheckClose(t, "pos3d", pos3d[i][k][c], fx.Pos3D[i][k][c].Float())
			}
		}
	}

	// Current ring field / potential / axial derivatives.
	currents, jac, pos := ringEffectiveCharges(fx.Current.Float(), fx.R.Float())
	for m, fp := range fx.FieldPoints {
		p := geom3d.Vec3{fp[0].Float(), fp[1].Float(), fp[2].Float()}
		f := CurrentFieldRadial(p, currents, jac, pos)
		for c := 0; c < 3; c++ {
			oracle.CheckClose(t, "cur_field", f[c], fx.CurField[m][c].Float())
		}
	}
	z := make([]float64, len(fx.ZAxis))
	for i := range z {
		z[i] = fx.ZAxis[i].Float()
		oracle.CheckClose(t, "cur_potential", CurrentPotentialAxial(z[i], currents, jac, pos), fx.CurPotential[i].Float())
	}
	derivs := CurrentAxialDerivativesRadial(z, currents, jac, pos)
	for i := range z {
		for l := 0; l < Deriv2DMax; l++ {
			oracle.CheckClose(t, "cur_axial_deriv", derivs[i][l], fx.CurAxialDeriv[i][l].Float())
		}
	}
}

// TestCurrentFieldVsBiotSavart cross-checks the current ring field against the analytic
// Biot-Savart law for a circular loop (ported from test_radial_ring.biot_savart_loop): for
// a unit-radius loop, μ₀·CurrentFieldRadial must equal the Biot-Savart field.
func TestCurrentFieldVsBiotSavart(t *testing.T) {
	const current = 2.5
	currents, jac, pos := ringEffectiveCharges(current, 1.0)
	points := []geom3d.Vec3{{0.5, 0, 0.5}, {0.5, 0, -0.25}, {0.2, 0, 0.1}, {1.2, 0, 0.3}}
	for _, p := range points {
		got := CurrentFieldRadial(p, currents, jac, pos)
		// μ₀·H is the B field; compare to Biot-Savart (which already includes μ₀).
		gotB := geom3d.Vec3{constants.VacuumPermeability * got[0], constants.VacuumPermeability * got[1], constants.VacuumPermeability * got[2]}
		want := biotSavartLoop(current, p)
		// Biot-Savart returns [Bx, By, Bz]; the radial field's x-component is Bx, z is Bz.
		if !oracle.IsClose(gotB[0], want[0], 1e-6, 1e-9) || !oracle.IsClose(gotB[2], want[2], 1e-6, 1e-9) {
			t.Errorf("point %v: got B=(%.10g,%.10g), want Biot-Savart (%.10g,%.10g)", p, gotB[0], gotB[2], want[0], want[2])
		}
	}
}

// biotSavartLoop computes the magnetic field of a unit-radius circular current loop in the
// xy-plane at r_point, by quadrature over the loop. Port of test_radial_ring.biot_savart_loop.
func biotSavartLoop(current float64, rPoint geom3d.Vec3) geom3d.Vec3 {
	integrand := func(axis int) quad.Func {
		return func(tt float64) float64 {
			rLoop := geom3d.Vec3{math.Cos(tt), math.Sin(tt), 0}
			dl := geom3d.Vec3{-math.Sin(tt), math.Cos(tt), 0}
			r := geom3d.Vec3{rPoint[0] - rLoop[0], rPoint[1] - rLoop[1], rPoint[2] - rLoop[2]}
			db := geom3d.CrossProduct3D(dl, r)
			norm := geom3d.Norm3D(r[0], r[1], r[2])
			return db[axis] / (norm * norm * norm)
		}
	}
	scale := current * constants.VacuumPermeability / (4 * math.Pi)
	var out geom3d.Vec3
	for axis := 0; axis < 3; axis++ {
		out[axis] = scale * quad.AdaptiveRecursive(integrand(axis), 0, 2*math.Pi, 1e-11, 1e-11)
	}
	return out
}
