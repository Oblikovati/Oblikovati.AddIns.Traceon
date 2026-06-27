// SPDX-License-Identifier: MPL-2.0

// Package field evaluates the radial field from a solved charge distribution, porting
// traceon/field.py's FieldRadialBEM (direct boundary integration) and FieldRadialAxial
// (the fast axial-series interpolation). The direct evaluators are thin wrappers over
// core/radial; the axial field builds a piecewise-quintic interpolation of the on-axis
// derivatives (core/interp) and reconstructs the off-axis field from the radial series.
//
// Scope: electrostatic radial. Magnetostatic/current contributions are added in their PBI.
package field

import (
	"fmt"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/interp"
	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/solver"
)

// FieldRadialBEM evaluates the field directly from the effective point charges by
// integrating over every element — accurate everywhere, slower per evaluation. Port of
// FieldRadialBEM. Carries the electrostatic, magnetostatic, and current charge sets; any
// may be empty.
type FieldRadialBEM struct {
	elec    solver.EffectivePointCharges
	mag     solver.EffectivePointCharges
	current solver.CurrentCharges
}

// NewFieldRadialBEM wraps a solved electrostatic charge distribution for direct evaluation.
func NewFieldRadialBEM(epc solver.EffectivePointCharges) FieldRadialBEM {
	return FieldRadialBEM{elec: epc}
}

// NewFieldRadialBEMFull wraps electrostatic, magnetostatic, and current charge sets.
func NewFieldRadialBEMFull(elec, mag solver.EffectivePointCharges, current solver.CurrentCharges) FieldRadialBEM {
	return FieldRadialBEM{elec: elec, mag: mag, current: current}
}

// PotentialAtPoint returns the electrostatic potential at point. Port of
// electrostatic_potential_at_local_point (= backend.potential_radial).
func (f FieldRadialBEM) PotentialAtPoint(point geom2d.Vertex) float64 {
	return radial.PotentialRadial(point, f.elec.Charges, f.elec.Jac, f.elec.Pos)
}

// FieldAtPoint returns the electric field (Ex, Ey, Ez) at point. Port of
// electrostatic_field_at_local_point (= backend.field_radial).
func (f FieldRadialBEM) FieldAtPoint(point geom2d.Vertex) geom2d.Vertex {
	return radial.FieldRadial(point, f.elec.Charges, f.elec.Jac, f.elec.Pos)
}

// CurrentFieldAtPoint returns the magnetic field (Hx, Hy, Hz) at point from the current
// rings alone. Port of current_field_at_local_point.
func (f FieldRadialBEM) CurrentFieldAtPoint(point geom2d.Vertex) geom2d.Vertex {
	h := radial.CurrentFieldRadial(geom3d.Vec3{point[0], point[1], point[2]}, f.current.Currents, f.current.Jac, f.current.Pos)
	return geom2d.Vertex{h[0], h[1], h[2]}
}

// MagnetostaticFieldAtPoint returns the total magnetic field (Hx, Hy, Hz): the current
// field plus the field of the magnetostatic surface charges. Port of
// magnetostatic_field_at_local_point (= current_field + field_radial(mag charges)).
func (f FieldRadialBEM) MagnetostaticFieldAtPoint(point geom2d.Vertex) geom2d.Vertex {
	cur := f.CurrentFieldAtPoint(point)
	mag := radial.FieldRadial(point, f.mag.Charges, f.mag.Jac, f.mag.Pos)
	return geom2d.Vertex{cur[0] + mag[0], cur[1] + mag[1], cur[2] + mag[2]}
}

// FieldRadialAxial evaluates the field via the axial-series interpolation: it samples the
// on-axis derivatives on a uniform grid, fits a piecewise-quintic interpolation, and
// reconstructs the off-axis field from the radial series. Fast near the axis; accuracy
// degrades with radius and the field is zero outside [zmin, zmax]. Port of FieldRadialAxial.
type FieldRadialAxial struct {
	z      []float64
	coeffs radial.AxialCoeffs
}

// NewFieldRadialAxial builds the axial interpolation of an electrostatic charge
// distribution over N uniformly-spaced samples in [zmin, zmax]. Port of
// FieldRadialAxial._get_interpolation_coefficients (electrostatic).
func NewFieldRadialAxial(epc solver.EffectivePointCharges, zmin, zmax float64, n int) (FieldRadialAxial, error) {
	if zmax <= zmin {
		return FieldRadialAxial{}, fmt.Errorf("field.NewFieldRadialAxial: zmax %g must exceed zmin %g", zmax, zmin)
	}
	if n < 3 {
		return FieldRadialAxial{}, fmt.Errorf("field.NewFieldRadialAxial: need >=3 samples, got %d", n)
	}
	z := linspace(zmin, zmax, n)
	derivs := radial.AxialDerivativesRadial(z, epc.Charges, epc.Jac, epc.Pos)
	coeffs, err := quinticSplineCoefficients(z, derivs)
	if err != nil {
		return FieldRadialAxial{}, err
	}
	return FieldRadialAxial{z: z, coeffs: coeffs}, nil
}

// PotentialAtPoint returns the interpolated potential at point (zero outside [zmin, zmax]).
func (f FieldRadialAxial) PotentialAtPoint(point geom2d.Vertex) float64 {
	return radial.PotentialRadialDerivs(point, f.z, f.coeffs)
}

// FieldAtPoint returns the interpolated electric field at point (zero outside [zmin, zmax]).
func (f FieldRadialAxial) FieldAtPoint(point geom2d.Vertex) geom2d.Vertex {
	return radial.FieldRadialDerivs(point, f.z, f.coeffs)
}

// quinticSplineCoefficients assembles the per-interval [Deriv2DMax][6] coefficient blocks
// from the sampled axial derivatives, porting field.py _quintic_spline_coefficients: the
// lower derivative orders (i+2 < Deriv2DMax) get a quintic Hermite using (d_i, d_{i+1},
// d_{i+2}); the top two orders get a not-a-knot cubic placed in the cubic slots [2:6].
func quinticSplineCoefficients(z []float64, derivsByZ [][radial.Deriv2DMax]float64) (radial.AxialCoeffs, error) {
	n := len(z)
	coeffs := make(radial.AxialCoeffs, n-1)
	column := func(order int) []float64 {
		out := make([]float64, n)
		for k := range z {
			out[k] = derivsByZ[k][order]
		}
		return out
	}
	for i := 0; i < radial.Deriv2DMax; i++ {
		if i+2 < radial.Deriv2DMax { // high order → quintic Hermite (6 coeffs, slots [0:6])
			q, err := interp.QuinticHermite(z, column(i), column(i+1), column(i+2))
			if err != nil {
				return nil, err
			}
			for iv := 0; iv < n-1; iv++ {
				coeffs[iv][i] = q[iv]
			}
		} else { // top orders → not-a-knot cubic (4 coeffs, slots [2:6]; [0:2] stay zero)
			cu, err := interp.NotAKnotCubic(z, column(i))
			if err != nil {
				return nil, err
			}
			for iv := 0; iv < n-1; iv++ {
				for c := 0; c < 4; c++ {
					coeffs[iv][i][2+c] = cu[iv][c]
				}
			}
		}
	}
	return coeffs, nil
}

// linspace returns n evenly-spaced samples over [start, stop] inclusive (numpy.linspace).
func linspace(start, stop float64, n int) []float64 {
	out := make([]float64, n)
	if n == 1 {
		out[0] = start
		return out
	}
	step := (stop - start) / float64(n-1)
	for i := 0; i < n; i++ {
		out[i] = start + float64(i)*step
	}
	out[n-1] = stop // exact endpoint (numpy pins it)
	return out
}
