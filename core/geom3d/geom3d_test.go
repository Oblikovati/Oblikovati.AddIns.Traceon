// SPDX-License-Identifier: MPL-2.0

package geom3d

import (
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
)

// TestNormal3D ports test_normal_3d: the normal is unit-length and orthogonal to both edges.
func TestNormal3D(t *testing.T) {
	tri := Triangle{{0, 0, 0}, {2, 0, 0}, {0, 3, 0}}
	n := Normal3D(tri)
	vec1 := sub(tri[1], tri[0])
	vec2 := sub(tri[2], tri[0])
	oracle.CheckClose(t, "|normal|", Norm3D(n[0], n[1], n[2]), 1.0)
	oracle.CheckClose(t, "normal·vec1", Dot3D(n, vec1), 0.0)
	oracle.CheckClose(t, "normal·vec2", Dot3D(n, vec2), 0.0)
}

// TestPositionAndJacobian3D ports test_position_and_jacobian: pos = v0 + alpha*vec1 + beta*vec2.
func TestPositionAndJacobian3D(t *testing.T) {
	tri := Triangle{{1, 1, 1}, {2, 0, 0}, {0, 3, 0}}
	vec1 := sub(tri[1], tri[0])
	vec2 := sub(tri[2], tri[0])
	alpha, beta := 1.0/3, 1.0/4
	_, pos := PositionAndJacobian3D(alpha, beta, tri)
	for k := 0; k < 3; k++ {
		want := tri[0][k] + vec1[k]*alpha + beta*vec2[k]
		oracle.CheckClose(t, "pos", pos[k], want)
	}
}

// TestTriangleArea checks a right triangle of legs 2 and 3 has area 3.
func TestTriangleArea(t *testing.T) {
	oracle.CheckClose(t, "area", TriangleArea(Vec3{0, 0, 0}, Vec3{2, 0, 0}, Vec3{0, 3, 0}), 3.0)
}

// TestFillJacobianBuffer3D checks the triangle quadrature buffer: the Jacobians (2·w·area)
// sum to the triangle area (since the reference-triangle weights sum to 1/2), and each
// quadrature position is the correct barycentric combination of the vertices.
func TestFillJacobianBuffer3D(t *testing.T) {
	tri := Triangle{{1, 0, 0}, {1.5, 0, 0}, {1.25, 0, 0.5}}
	jac, pos := FillJacobianBuffer3D([]Triangle{tri})
	area := TriangleArea(tri[0], tri[1], tri[2])
	sum := 0.0
	for k := 0; k < N_TRIANGLE_QUAD; k++ {
		sum += jac[0][k]
		// Position must equal v0 + b1·(v1-v0) + b2·(v2-v0).
		for c := 0; c < 3; c++ {
			want := tri[0][c] + QuadB1[k]*(tri[1][c]-tri[0][c]) + QuadB2[k]*(tri[2][c]-tri[0][c])
			oracle.CheckClose(t, "quad pos", pos[0][k][c], want)
		}
	}
	oracle.CheckClose(t, "sum(jac)=area", sum, area)
}

// TestOrientation checks the triangle-orientation test on a shared-edge pair: two
// triangles wound the same way agree; flipping one disagrees; a disjoint pair has no edge.
func TestOrientation(t *testing.T) {
	// Two triangles sharing edge (0,1) in the z=0 plane, both wound counter-clockwise.
	points := []Vec3{{0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, -1, 0}}
	// tri A = (0,1,2) CCW; tri B shares edge (0,1) but on the other side: (1,0,3) keeps a
	// consistent outward winding across the shared edge.
	tris := [][3]uint64{{0, 1, 2}, {1, 0, 3}}
	if got := TriangleOrientationIsEqual(0, 1, tris, points); got != OrientationEqual {
		t.Errorf("consistent pair: got %d, want OrientationEqual", got)
	}
	// Flip triangle B's winding → orientations disagree.
	trisFlipped := [][3]uint64{{0, 1, 2}, {0, 1, 3}}
	if got := TriangleOrientationIsEqual(0, 1, trisFlipped, points); got != OrientationOpposite {
		t.Errorf("flipped pair: got %d, want OrientationOpposite", got)
	}
	// Disjoint triangles share no edge.
	disjoint := []Vec3{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {9, 9, 9}, {9, 9, 8}, {9, 8, 9}}
	trisDisjoint := [][3]uint64{{0, 1, 2}, {3, 4, 5}}
	if got := TriangleOrientationIsEqual(0, 1, trisDisjoint, disjoint); got != OrientationNoCommonEdge {
		t.Errorf("disjoint pair: got %d, want OrientationNoCommonEdge", got)
	}
}
