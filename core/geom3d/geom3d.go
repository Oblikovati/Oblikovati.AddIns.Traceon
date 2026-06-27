// SPDX-License-Identifier: MPL-2.0

// Package geom3d ports traceon/backend/utilities_3d.c and triangle.c: 3D vector
// primitives, triangle area/normal, the flat-triangle parametrization, and the
// triangle-orientation test the mesher uses to make surface normals consistent. The
// radial path needs only a few of these (triangle areas/normals for current rings and
// effective point charges); the full 3D BEM is Traceon Pro and out of scope.
package geom3d

import "math"

// Vec3 is a 3D vector/point.
type Vec3 = [3]float64

// Triangle is three 3D vertices.
type Triangle = [3]Vec3

// N_TRIANGLE_QUAD is the order of the symmetric triangle quadrature rule (matches
// N_TRIANGLE_QUAD in utilities_3d.c). Used when reducing a triangle panel to effective
// point charges.
const N_TRIANGLE_QUAD = 12

// QuadWeights / QuadB1 / QuadB2 are the 12-point symmetric triangle-quadrature weights
// and barycentric nodes, verbatim from utilities_3d.c.
var (
	QuadWeights = [N_TRIANGLE_QUAD]float64{
		0.0254224531851035, 0.0254224531851035, 0.0254224531851035,
		0.0583931378631895, 0.0583931378631895, 0.0583931378631895,
		0.041425537809187, 0.041425537809187, 0.041425537809187,
		0.041425537809187, 0.041425537809187, 0.041425537809187,
	}
	QuadB1 = [N_TRIANGLE_QUAD]float64{
		0.873821971016996, 0.063089014491502, 0.063089014491502,
		0.501426509658179, 0.249286745170910, 0.249286745170910,
		0.636502499121399, 0.636502499121399, 0.310352451033785,
		0.310352451033785, 0.053145049844816, 0.053145049844816,
	}
	QuadB2 = [N_TRIANGLE_QUAD]float64{
		0.063089014491502, 0.873821971016996, 0.063089014491502,
		0.249286745170910, 0.501426509658179, 0.249286745170910,
		0.310352451033785, 0.053145049844816, 0.636502499121399,
		0.053145049844816, 0.636502499121399, 0.310352451033785,
	}
)

// Dot3D returns the dot product v1·v2.
func Dot3D(v1, v2 Vec3) float64 { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] }

// Norm3D returns sqrt(x^2 + y^2 + z^2).
func Norm3D(x, y, z float64) float64 { return math.Sqrt(x*x + y*y + z*z) }

// Distance3D returns the Euclidean distance between v0 and v1.
func Distance3D(v0, v1 Vec3) float64 { return Norm3D(v0[0]-v1[0], v0[1]-v1[1], v0[2]-v1[2]) }

// Normalize3D returns v scaled to unit length.
func Normalize3D(v Vec3) Vec3 {
	l := Norm3D(v[0], v[1], v[2])
	return Vec3{v[0] / l, v[1] / l, v[2] / l}
}

// CrossProduct3D returns v1 × v2.
func CrossProduct3D(v1, v2 Vec3) Vec3 {
	return Vec3{
		v1[1]*v2[2] - v1[2]*v2[1],
		v1[2]*v2[0] - v1[0]*v2[2],
		v1[0]*v2[1] - v1[1]*v2[0],
	}
}

// TriangleArea returns the area of the triangle (v0, v1, v2).
func TriangleArea(v0, v1, v2 Vec3) float64 {
	vec1 := Vec3{v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]}
	vec2 := Vec3{v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]}
	out := CrossProduct3D(vec1, vec2)
	return Norm3D(out[0], out[1], out[2]) / 2.0
}

// TriangleAreas returns the area of each triangle (batch form of three_d.c triangle_areas).
func TriangleAreas(tris []Triangle) []float64 {
	out := make([]float64, len(tris))
	for i, t := range tris {
		out[i] = TriangleArea(t[0], t[1], t[2])
	}
	return out
}

// Normal3D returns the unit normal of a triangle. Port of normal_3d (note: this uses the
// raw edge-cross formula and a different component layout than CrossProduct3D, kept verbatim
// from the C so the sign convention matches the upstream exactly).
func Normal3D(t Triangle) Vec3 {
	x1, y1, z1 := t[0][0], t[0][1], t[0][2]
	x2, y2, z2 := t[1][0], t[1][1], t[1][2]
	x3, y3, z3 := t[2][0], t[2][1], t[2][2]
	nx := (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1)
	ny := (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1)
	nz := (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
	l := Norm3D(nx, ny, nz)
	return Vec3{nx / l, ny / l, nz / l}
}

// CurrentJacobianBuffer holds the triangle quadrature Jacobians (2·w·area) at the
// N_TRIANGLE_QUAD points; CurrentPositionBuffer the corresponding 3D positions. Indexed
// [triangleIndex][quadIndex]. Used to treat axisymmetric coil-cross-section triangles as
// current rings.
type (
	CurrentJacobianBuffer [][N_TRIANGLE_QUAD]float64
	CurrentPositionBuffer [][N_TRIANGLE_QUAD]Vec3
)

// FillJacobianBuffer3D precomputes, for every triangle, the quadrature Jacobian (2·w·area)
// and barycentric position at the N_TRIANGLE_QUAD points. Port of fill_jacobian_buffer_3d.
func FillJacobianBuffer3D(triangles []Triangle) (CurrentJacobianBuffer, CurrentPositionBuffer) {
	jac := make(CurrentJacobianBuffer, len(triangles))
	pos := make(CurrentPositionBuffer, len(triangles))
	for i, t := range triangles {
		area := TriangleArea(t[0], t[1], t[2])
		for k := 0; k < N_TRIANGLE_QUAD; k++ {
			b1, b2 := QuadB1[k], QuadB2[k]
			jac[i][k] = 2 * QuadWeights[k] * area
			pos[i][k] = Vec3{
				t[0][0] + b1*(t[1][0]-t[0][0]) + b2*(t[2][0]-t[0][0]),
				t[0][1] + b1*(t[1][1]-t[0][1]) + b2*(t[2][1]-t[0][1]),
				t[0][2] + b1*(t[1][2]-t[0][2]) + b2*(t[2][2]-t[0][2]),
			}
		}
	}
	return jac, pos
}

// PositionAndJacobian3D returns (jac, pos) for the flat-triangle parametrization at
// barycentric (alpha, beta): pos = t0 + alpha*(t1-t0) + beta*(t2-t0); jac = 2*area.
// Port of position_and_jacobian_3d (returns (jac, pos) to match the Python tuple order).
func PositionAndJacobian3D(alpha, beta float64, t Triangle) (jac float64, pos Vec3) {
	v1 := Vec3{t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]}
	v2 := Vec3{t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]}
	pos = Vec3{
		t[0][0] + alpha*v1[0] + beta*v2[0],
		t[0][1] + alpha*v1[1] + beta*v2[1],
		t[0][2] + alpha*v1[2] + beta*v2[2],
	}
	jac = 2 * TriangleArea(t[0], t[1], t[2])
	return jac, pos
}
