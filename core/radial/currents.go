// SPDX-License-Identifier: MPL-2.0

package radial

import (
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/ring"
)

// CurrentFieldRadial returns the magnetic field (Hx, Hy, Hz) at point produced by a set of
// axisymmetric current rings (coil-cross-section triangle quadrature points carrying a
// current density). On the axis (r < MinDistanceAxis) the in-plane components are zeroed.
// Port of current_field_radial.
func CurrentFieldRadial(point geom3d.Vec3, currents []float64, jac geom3d.CurrentJacobianBuffer, pos geom3d.CurrentPositionBuffer) geom3d.Vec3 {
	r := geom2d.Norm2D(point[0], point[1])
	br, bz := 0.0, 0.0
	for i := range currents {
		for k := 0; k < geom3d.N_TRIANGLE_QUAD; k++ {
			p := pos[i][k]
			f := ring.CurrentFieldRadialRing(r, point[2], p[0], p[2])
			br += currents[i] * jac[i][k] * f[0]
			bz += currents[i] * jac[i][k] * f[1]
		}
	}
	var result geom3d.Vec3
	if r >= ring.MinDistanceAxis {
		result[0] = point[0] / r * br
		result[1] = point[1] / r * br
	}
	result[2] = bz
	return result
}

// CurrentPotentialAxial returns the on-axis magnetic scalar potential at z0 produced by the
// current rings. Port of current_potential_axial.
func CurrentPotentialAxial(z0 float64, currents []float64, jac geom3d.CurrentJacobianBuffer, pos geom3d.CurrentPositionBuffer) float64 {
	sum := 0.0
	for i := range currents {
		for k := 0; k < geom3d.N_TRIANGLE_QUAD; k++ {
			p := pos[i][k]
			sum += currents[i] * jac[i][k] * ring.CurrentPotentialAxialRadialRing(z0, p[0], p[2])
		}
	}
	return sum
}

// CurrentAxialDerivativesRadial accumulates, for each axis sample z[i], the first Deriv2DMax
// on-axis z-derivatives of the current rings' axial potential. Port of
// current_axial_derivatives_radial.
func CurrentAxialDerivativesRadial(z []float64, currents []float64, jac geom3d.CurrentJacobianBuffer, pos geom3d.CurrentPositionBuffer) [][Deriv2DMax]float64 {
	derivs := make([][Deriv2DMax]float64, len(z))
	for i := range z {
		for j := range currents {
			for k := 0; k < geom3d.N_TRIANGLE_QUAD; k++ {
				p := pos[j][k]
				d := ring.CurrentAxialDerivativesRadialRing(z[i], p[0], p[2])
				w := jac[j][k] * currents[j]
				for l := 0; l < Deriv2DMax; l++ {
					derivs[i][l] += w * d[l]
				}
			}
		}
	}
	return derivs
}
