// SPDX-License-Identifier: MPL-2.0

package radial

import (
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/ring"
)

// Deriv2DMax is the number of on-axis derivative terms (== ring.Deriv2DMax).
const Deriv2DMax = ring.Deriv2DMax

// nCoeff is the number of polynomial coefficients per derivative order in the axial
// interpolation: a quintic (degree 5) has 6 coefficients (matches the C coeff[..][..][6]).
const nCoeff = 6

// AxialDerivativesRadial accumulates, for each axis sample z[i], the first Deriv2DMax
// on-axis z-derivatives of the potential produced by the effective point charges. Port of
// axial_derivatives_radial: derivs[i][l] = Σ_j Σ_k jac[j][k]·charges[j]·D_l(z[i], pos[j][k]).
func AxialDerivativesRadial(z []float64, charges []float64, jac JacobianBuffer, pos PositionBuffer) [][Deriv2DMax]float64 {
	derivs := make([][Deriv2DMax]float64, len(z))
	for i := range z {
		for j := range charges {
			for k := 0; k < NQuad2D; k++ {
				p := pos[j][k]
				d := ring.AxialDerivativesRadialRing(z[i], p[0], p[1])
				w := jac[j][k] * charges[j]
				for l := 0; l < Deriv2DMax; l++ {
					derivs[i][l] += w * d[l]
				}
			}
		}
	}
	return derivs
}

// AxialCoeffs is the piecewise-quintic interpolation of the Deriv2DMax axial derivatives:
// one [Deriv2DMax][nCoeff] block of polynomial coefficients per z-interval (there are
// len(z)-1 intervals). Coefficients are in descending powers of the local variable
// diffz = z - z[index], matching the C coeff[interval][order][6] layout.
type AxialCoeffs [][Deriv2DMax][nCoeff]float64

// evalDerivs evaluates the Deriv2DMax interpolated derivatives at point z, given the axis
// grid zInter and the quintic coefficients. Mirrors the shared head of potential/field_radial_derivs.
// Returns the derivatives and whether z is in range (out-of-range → caller returns zero field).
func evalDerivs(r, z float64, zInter []float64, coeff AxialCoeffs) ([Deriv2DMax]float64, bool) {
	var derivs [Deriv2DMax]float64
	nz := len(zInter)
	z0, zlast := zInter[0], zInter[nz-1]
	if !(z0 <= z && z <= zlast) {
		return derivs, false
	}
	dz := zInter[1] - zInter[0]
	index := int((z - z0) / dz)
	// Guard against roundoff: there are nz grid points but nz-1 coefficient intervals.
	if index < 0 {
		index = 0
	}
	if index > nz-2 {
		index = nz - 2
	}
	diffz := z - zInter[index]
	c := &coeff[index]
	for i := 0; i < Deriv2DMax; i++ {
		ci := c[i]
		derivs[i] = ci[0]*pow5(diffz) + ci[1]*pow4(diffz) + ci[2]*diffz*diffz*diffz + ci[3]*diffz*diffz + ci[4]*diffz + ci[5]
	}
	return derivs, true
}

func pow4(x float64) float64 { return x * x * x * x }
func pow5(x float64) float64 { return x * x * x * x * x }

// PotentialRadialDerivs evaluates the interpolated potential at point using the radial
// series of axial derivatives. Port of potential_radial_derivs. Returns 0 outside [z0, zlast].
func PotentialRadialDerivs(point geom2d.Vertex, zInter []float64, coeff AxialCoeffs) float64 {
	r := geom2d.Norm2D(point[0], point[1])
	d, ok := evalDerivs(r, point[2], zInter, coeff)
	if !ok {
		return 0.0
	}
	r2 := r * r
	r4 := r2 * r2
	r6 := r4 * r2
	r8 := r4 * r4
	return d[0] - r2/4*d[2] + r4/64*d[4] - r6/2304*d[6] + r8/147456*d[8]
}

// FieldRadialDerivs evaluates the interpolated electric field at point using the radial
// series of axial derivatives. Port of field_radial_derivs. Returns the zero field outside
// [z0, zlast]. The radial component is already divided by r so the axis is safe.
func FieldRadialDerivs(point geom2d.Vertex, zInter []float64, coeff AxialCoeffs) geom2d.Vertex {
	r := geom2d.Norm2D(point[0], point[1])
	d, ok := evalDerivs(r, point[2], zInter, coeff)
	if !ok {
		return geom2d.Vertex{}
	}
	r2 := r * r
	r4 := r2 * r2
	r6 := r4 * r2
	fieldRadial := 0.5 * (d[2] - r2/8*d[4] + r4/192*d[6] - r6/9216*d[8])
	fieldZ := -d[1] + r2/4*d[3] - r4/64*d[5] + r6/2304*d[7]
	return geom2d.Vertex{point[0] * fieldRadial, point[1] * fieldRadial, fieldZ}
}
