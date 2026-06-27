// SPDX-License-Identifier: MPL-2.0

// Package ring ports traceon/backend/radial_ring.c: the closed-form kernels for a single
// charged or current-carrying ring in the (r, z) half-plane. These are the Green's
// functions the radial BEM integrates over each line element — the potential of a ring
// of charge, its r/z field derivatives, the axial current potential/field, and the
// on-axis derivative recurrences used by the axial field interpolation.
//
// The expressions are written in terms of the complete elliptic integrals via the
// km1/em1 forms (accurate near the m→1 self-interaction singularity), matching the C.
package ring

import (
	"math"

	"oblikovati.org/traceon/core/elliptic"
	"oblikovati.org/traceon/core/geom2d"
)

// MinDistanceAxis guards the on-axis singularity (MIN_DISTANCE_AXIS in three_d.c).
const MinDistanceAxis = 1e-10

// Deriv2DMax is the number of on-axis derivative terms produced by the recurrences
// (DERIV_2D_MAX in radial_ring.c). The axial field interpolation consumes these.
const Deriv2DMax = 9

// invPi is 1/π; the C uses 1./M_PI with M_PI = 3.14159265358979323846 (== math.Pi).
const invPi = 1.0 / math.Pi

// PotentialRadialRing returns the potential at (r0, z0) of a unit ring offset by
// (deltaR, deltaZ), i.e. the ring passing through (r0+deltaR, z0+deltaZ). Port of
// potential_radial_ring. (Result is in units where it must be divided by epsilon_0 to
// get SI potential — see the tests.)
func PotentialRadialRing(r0, z0, deltaR, deltaZ float64) float64 {
	t := (deltaZ*deltaZ + deltaR*deltaR) / (deltaZ*deltaZ + deltaR*deltaR + 4*r0*deltaR + 4*r0*r0)
	return invPi * elliptic.Ellipkm1(t) * (deltaR + r0) / math.Sqrt(deltaZ*deltaZ+(deltaR+2*r0)*(deltaR+2*r0))
}

// Dr1PotentialRadialRing returns ∂/∂r of the ring potential. Port of dr1_potential_radial_ring.
// Returns 0 on the axis (r0 < MinDistanceAxis), matching the C guard.
func Dr1PotentialRadialRing(r0, z0, deltaR, deltaZ float64) float64 {
	if r0 < MinDistanceAxis {
		return 0.0
	}
	r := r0 + deltaR
	commonArg := (deltaZ*deltaZ + deltaR*deltaR) / (4*r*r - 4*deltaR*r + deltaZ*deltaZ + deltaR*deltaR)
	denominator := ((-2 * deltaR * deltaR * r) + deltaZ*deltaZ*(2*deltaR-2*r) + 2*deltaR*deltaR*deltaR) * math.Sqrt(4*r*r-4*deltaR*r+deltaZ*deltaZ+deltaR*deltaR)
	ellipkm1Term := (deltaZ*deltaZ*r + deltaR*deltaR*r) * elliptic.Ellipkm1(commonArg)
	ellipem1Term := ((-2 * deltaR * r * r) - deltaZ*deltaZ*r + deltaR*deltaR*r) * elliptic.Ellipem1(commonArg)
	return invPi * (ellipkm1Term + ellipem1Term) / denominator
}

// Dz1PotentialRadialRing returns ∂/∂z of the ring potential. Port of dz1_potential_radial_ring.
func Dz1PotentialRadialRing(r0, z0, deltaR, deltaZ float64) float64 {
	commonArg := (deltaZ*deltaZ + deltaR*deltaR) / (4*r0*r0 + 4*deltaR*r0 + deltaZ*deltaZ + deltaR*deltaR)
	denominator := (deltaZ*deltaZ + deltaR*deltaR) * math.Sqrt(4*r0*r0+4*deltaR*r0+deltaZ*deltaZ+deltaR*deltaR)
	ellipem1Term := -deltaZ * (r0 + deltaR) * elliptic.Ellipem1(commonArg)
	return invPi * -ellipem1Term / denominator
}

// FluxDensityToChargeFactor maps relative permittivity K to the effective dielectric
// charge factor 2(K-1)/(K+1). Port of flux_density_to_charge_factor.
func FluxDensityToChargeFactor(k float64) float64 { return 2.0 * (k - 1) / (1 + k) }

// FieldDotNormalRadial returns the dielectric-weighted normal field at (r0,z0) for a ring
// at the offset, with surface normal (nr, nz) and relative permittivity K. Port of
// field_dot_normal_radial.
func FieldDotNormalRadial(r0, z0, deltaR, deltaZ float64, normal geom2d.Point2, k float64) float64 {
	factor := FluxDensityToChargeFactor(k)
	er := -Dr1PotentialRadialRing(r0, z0, deltaR, deltaZ)
	ez := -Dz1PotentialRadialRing(r0, z0, deltaR, deltaZ)
	return factor * (normal[0]*er + normal[1]*ez)
}

// CurrentPotentialAxialRadialRing returns the on-axis magnetic scalar potential of a
// current ring of radius r at height z, evaluated at z0. Port of
// current_potential_axial_radial_ring.
func CurrentPotentialAxialRadialRing(z0, r, z float64) float64 {
	dz := z0 - z
	return -dz / (2 * math.Sqrt(dz*dz+r*r))
}

// CurrentFieldRadialRing returns the (r, z) magnetic field components at (x0, y0) of a
// unit current ring of radius x at height y. Port of current_field_radial_ring (result[0]
// = radial, result[1] = axial). Uses the full elliptic integrals (not the km1 forms).
func CurrentFieldRadialRing(x0, y0, x, y float64) geom2d.Point2 {
	a := x
	r := x0
	z := y0 - y

	A := z*z + (a+r)*(a+r)
	B := z*z + (r-a)*(r-a)
	k := 4 * r * a / A

	var result geom2d.Point2
	if x < MinDistanceAxis {
		// Unphysical: infinitely small ring.
		return result
	}
	if x0 < MinDistanceAxis {
		result[0] = 0.0
	} else {
		result[0] = invPi * z / (2 * r * math.Sqrt(A)) * ((z*z+r*r+a*a)/B*elliptic.Ellipe(k) - elliptic.Ellipk(k))
	}
	result[1] = invPi * 1 / (2 * math.Sqrt(A)) * ((a*a-z*z-r*r)/B*elliptic.Ellipe(k) + elliptic.Ellipk(k))
	return result
}

// AxialDerivativesRadialRing fills the first Deriv2DMax on-axis z-derivatives of the ring
// potential via the upward recurrence in radial_ring.c (electrostatic case). Port of
// axial_derivatives_radial_ring.
func AxialDerivativesRadialRing(z0, r, z float64) [Deriv2DMax]float64 {
	var d [Deriv2DMax]float64
	R := geom2d.Norm2D(z0-z, r)
	d[0] = 1 / R
	d[1] = -(z0 - z) / (R * R * R)
	for n := 1; n+1 < Deriv2DMax; n++ {
		d[n+1] = -1.0 / (R * R) * ((2*float64(n)+1)*(z0-z)*d[n] + float64(n)*float64(n)*d[n-1])
	}
	for n := 0; n < Deriv2DMax; n++ {
		d[n] *= r / 2
	}
	return d
}

// CurrentAxialDerivativesRadialRing fills the first Deriv2DMax on-axis z-derivatives of the
// current ring's axial potential. Port of current_axial_derivatives_radial_ring.
func CurrentAxialDerivativesRadialRing(z0, r, z float64) [Deriv2DMax]float64 {
	var d [Deriv2DMax]float64
	dz := z0 - z
	R := geom2d.Norm2D(dz, r)
	mu := dz / R
	d[0] = -dz / (2 * math.Sqrt(dz*dz+r*r))
	d[1] = -r * r / (2 * math.Pow(dz*dz+r*r, 1.5))
	for n := 2; n < Deriv2DMax; n++ {
		d[n] = -(2*float64(n)-1)*mu/R*d[n-1] - (float64(n)*float64(n)-2*float64(n))/(R*R)*d[n-2]
	}
	return d
}
