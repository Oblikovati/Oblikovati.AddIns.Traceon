// SPDX-License-Identifier: MPL-2.0

package tracing

import (
	"math"

	"oblikovati.org/traceon/core/constants"
	"oblikovati.org/traceon/core/geom3d"
)

// ConvertVelocityToSI converts a velocity vector whose magnitude is expressed in electron-
// volts of kinetic energy into one in m/s for a particle of the given mass (kg). Port of
// _convert_velocity_to_SI.
func ConvertVelocityToSI(velocity geom3d.Vec3, mass float64) geom3d.Vec3 {
	speedEV := geom3d.Norm3D(velocity[0], velocity[1], velocity[2])
	speed := math.Sqrt(2 * speedEV * constants.ElementaryCharge / mass)
	scale := speed / speedEV
	return geom3d.Vec3{velocity[0] * scale, velocity[1] * scale, velocity[2] * scale}
}

// VelocityVec builds an initial velocity (m/s) of the given kinetic energy (eV) along
// direction (normalized internally) for a particle of mass kg. Port of velocity_vec.
func VelocityVec(eV float64, direction geom3d.Vec3, mass float64) geom3d.Vec3 {
	norm := geom3d.Norm3D(direction[0], direction[1], direction[2])
	scaled := geom3d.Vec3{eV * direction[0] / norm, eV * direction[1] / norm, eV * direction[2] / norm}
	return ConvertVelocityToSI(scaled, mass)
}

// VelocityVecSpherical builds the initial velocity from energy (eV) and spherical angles
// (theta from the z-axis, phi from the x-axis). Port of velocity_vec_spherical.
func VelocityVecSpherical(eV, theta, phi, mass float64) geom3d.Vec3 {
	dir := geom3d.Vec3{math.Sin(theta) * math.Cos(phi), math.Sin(theta) * math.Sin(phi), math.Cos(theta)}
	return VelocityVec(eV, dir, mass)
}

// VelocityVecXZPlane builds the initial velocity in the xz-plane at the given angle to the
// z-axis (downward → negative z). Port of velocity_vec_xz_plane.
func VelocityVecXZPlane(eV, angle float64, downward bool, mass float64) geom3d.Vec3 {
	sign := 1.0
	if downward {
		sign = -1.0
	}
	return VelocityVec(eV, geom3d.Vec3{math.Sin(angle), 0.0, sign * math.Cos(angle)}, mass)
}

// ZToBounds returns a z-extent enclosing [z1, z2] with a unit margin, matching the
// upstream _z_to_bounds convention (used to build tracing bounds around an optical span).
func ZToBounds(z1, z2 float64) (float64, float64) {
	switch {
	case z1 < 0 && z2 < 0:
		return math.Min(z1, z2) - 1, 1.0
	case z1 > 0 && z2 > 0:
		return -1.0, math.Max(z1, z2) + 1
	default:
		return math.Min(z1, z2) - 1, math.Max(z1, z2) + 1
	}
}

// PlaneIntersection finds where a trajectory crosses the plane through p0 with the given
// normal, scanning from the end backward for the first sign change of the signed distance
// and linearly interpolating the 6-state there. Returns the state and true, or false if the
// trajectory never crosses. Port of plane_intersection (backend) — note it scans from the
// LAST point backward, matching the C.
func PlaneIntersection(positions []State, p0, normal geom3d.Vec3) (State, bool) {
	var out State
	n := len(positions)
	if n < 2 {
		return out, false
	}
	nn := geom3d.Norm3D(normal[0], normal[1], normal[2])
	kappa := func(s State) float64 {
		return (normal[2]*p0[2] - s[2]*normal[2] + normal[1]*p0[1] - s[1]*normal[1] + normal[0]*p0[0] - s[0]*normal[0]) / nn
	}
	prev := kappa(positions[n-1])
	for i := n - 2; i >= 0; i-- {
		k := kappa(positions[i])
		if sign(k) != sign(prev) {
			diff := k - prev
			factor := -prev / diff
			prevFactor := k / diff
			for c := 0; c < 6; c++ {
				out[c] = prevFactor*positions[i+1][c] + factor*positions[i][c]
			}
			return out, true
		}
		prev = k
	}
	return out, false
}

func sign(x float64) int {
	if x > 0 {
		return 1
	}
	return -1
}

// XYPlaneIntersection finds the crossing of the trajectory with the plane z = const.
func XYPlaneIntersection(positions []State, z float64) (State, bool) {
	return PlaneIntersection(positions, geom3d.Vec3{0, 0, z}, geom3d.Vec3{0, 0, 1})
}

// XZPlaneIntersection finds the crossing with the plane y = const.
func XZPlaneIntersection(positions []State, y float64) (State, bool) {
	return PlaneIntersection(positions, geom3d.Vec3{0, y, 0}, geom3d.Vec3{0, 1, 0})
}

// YZPlaneIntersection finds the crossing with the plane x = const.
func YZPlaneIntersection(positions []State, x float64) (State, bool) {
	return PlaneIntersection(positions, geom3d.Vec3{x, 0, 0}, geom3d.Vec3{1, 0, 0})
}

// AxisIntersection returns the z-value where the trajectory crosses the x=0 plane (the
// optical axis crossing). Port of axis_intersection. Returns false if there is no crossing.
func AxisIntersection(positions []State) (float64, bool) {
	s, ok := YZPlaneIntersection(positions, 0)
	return s[2], ok
}
