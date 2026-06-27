// SPDX-License-Identifier: MPL-2.0

// Package geom2d ports traceon/backend/utilities_2d.c: the 2D/axisymmetric geometry
// helpers for the radial BEM. A radial line element is a cubic (4-vertex) curve in the
// (r, z) half-plane; each vertex is stored as a 3D point and the element uses only its
// r = v[0] and z = v[2] components (the upstream reaches v[2] through a numpy slice that
// is a view onto the original 3-element buffer — see the port reference).
package geom2d

import "math"

// Vertex is a 3D point; radial elements use components 0 (r) and 2 (z).
type Vertex = [3]float64

// Point2 is a point in the (r, z) half-plane.
type Point2 = [2]float64

// N_QUAD_2D is the order of the Gauss-Legendre rule used to integrate radial elements
// (matches N_QUAD_2D in utilities_2d.c).
const N_QUAD_2D = 16

// GaussQuadWeights / GaussQuadPoints are the 16-point Gauss-Legendre weights and nodes
// on [-1, 1], verbatim from utilities_2d.c. Consumed by the radial Jacobian buffer fill.
var (
	GaussQuadWeights = [N_QUAD_2D]float64{
		0.1894506104550685, 0.1894506104550685, 0.1826034150449236, 0.1826034150449236,
		0.1691565193950025, 0.1691565193950025, 0.1495959888165767, 0.1495959888165767,
		0.1246289712555339, 0.1246289712555339, 0.0951585116824928, 0.0951585116824928,
		0.0622535239386479, 0.0622535239386479, 0.0271524594117541, 0.0271524594117541,
	}
	GaussQuadPoints = [N_QUAD_2D]float64{
		-0.0950125098376374, 0.0950125098376374, -0.2816035507792589, 0.2816035507792589,
		-0.4580167776572274, 0.4580167776572274, -0.6178762444026438, 0.6178762444026438,
		-0.7554044083550030, 0.7554044083550030, -0.8656312023878318, 0.8656312023878318,
		-0.9445750230732326, 0.9445750230732326, -0.9894009349916499, 0.9894009349916499,
	}
)

// Norm2D returns sqrt(x^2 + y^2).
func Norm2D(x, y float64) float64 { return math.Sqrt(x*x + y*y) }

// Length2D returns the distance between two (r, z) points.
func Length2D(v1, v2 Point2) float64 { return Norm2D(v2[0]-v1[0], v2[1]-v1[1]) }

// Normal2D returns the unit normal to the segment p1→p2, rotated (tx,ty)→(ty,-tx).
// Port of normal_2d.
func Normal2D(p1, p2 Point2) Point2 {
	tangentX, tangentY := p2[0]-p1[0], p2[1]-p1[1]
	nx, ny := tangentY, -tangentX
	length := Norm2D(nx, ny)
	return Point2{nx / length, ny / length}
}

// derivRZ evaluates the analytic derivative (d/dα) of the cubic element parametrization
// in r and z at α, used by both the normal and the Jacobian. Returns (dr, dz). This is
// the bracketed expression shared by higher_order_normal_radial and the jacobian term.
func derivRZ(alpha float64, v1, v2, v3, v4 Vertex) (dr, dz float64) {
	a2 := alpha * alpha
	r1, z1 := v1[0], v1[2]
	r2, z2 := v2[0], v2[2]
	r3, z3 := v3[0], v3[2]
	r4, z4 := v4[0], v4[2]
	dr = (2*alpha*(9*r4-9*r3-9*r2+9*r1) + 3*a2*(9*r4-27*r3+27*r2-9*r1) - r4 + 27*r3 - 27*r2 + r1) / 16
	dz = (2*alpha*(9*z4-9*z3-9*z2+9*z1) + 3*a2*(9*z4-27*z3+27*z2-9*z1) - z4 + 27*z3 - 27*z2 + z1) / 16
	return dr, dz
}

// HigherOrderNormalRadial returns the unit normal of the cubic radial element at α.
// Port of higher_order_normal_radial.
func HigherOrderNormalRadial(alpha float64, v1, v2, v3, v4 Vertex) Point2 {
	dr, dz := derivRZ(alpha, v1, v2, v3, v4)
	return Normal2D(Point2{0, 0}, Point2{dr, dz})
}

// PositionAndJacobianRadial returns the Jacobian (|dr/dα, dz/dα| times 1) and the (r, z)
// position of the cubic radial element at α ∈ [-1, 1]. Port of position_and_jacobian_radial;
// returns (jac, pos) to match the Python wrapper's tuple order.
func PositionAndJacobianRadial(alpha float64, v1, v2, v3, v4 Vertex) (jac float64, pos Point2) {
	a2 := alpha * alpha
	a3 := a2 * alpha
	r1, z1 := v1[0], v1[2]
	r2, z2 := v2[0], v2[2]
	r3, z3 := v3[0], v3[2]
	r4, z4 := v4[0], v4[2]
	pos[0] = (a2*(9*r4-9*r3-9*r2+9*r1) + a3*(9*r4-27*r3+27*r2-9*r1) - r4 + alpha*(-r4+27*r3-27*r2+r1) + 9*r3 + 9*r2 - r1) / 16
	pos[1] = (a2*(9*z4-9*z3-9*z2+9*z1) + a3*(9*z4-27*z3+27*z2-9*z1) - z4 + alpha*(-z4+27*z3-27*z2+z1) + 9*z3 + 9*z2 - z1) / 16
	dr, dz := derivRZ(alpha, v1, v2, v3, v4)
	jac = Norm2D(dr, dz)
	return jac, pos
}

// DeltaPositionAndJacobianRadial is PositionAndJacobianRadial with the position taken
// relative to α=0 (drops the constant terms), used for the singular self-potential
// integrals. Port of delta_position_and_jacobian_radial.
func DeltaPositionAndJacobianRadial(alpha float64, v1, v2, v3, v4 Vertex) (jac float64, pos Point2) {
	a2 := alpha * alpha
	a3 := a2 * alpha
	r1, z1 := v1[0], v1[2]
	r2, z2 := v2[0], v2[2]
	r3, z3 := v3[0], v3[2]
	r4, z4 := v4[0], v4[2]
	pos[0] = (a2*(9*r4-9*r3-9*r2+9*r1) + a3*(9*r4-27*r3+27*r2-9*r1) + alpha*(-r4+27*r3-27*r2+r1)) / 16
	pos[1] = (a2*(9*z4-9*z3-9*z2+9*z1) + a3*(9*z4-27*z3+27*z2-9*z1) + alpha*(-z4+27*z3-27*z2+z1)) / 16
	dr, dz := derivRZ(alpha, v1, v2, v3, v4)
	jac = Norm2D(dr, dz)
	return jac, pos
}
