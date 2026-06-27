// SPDX-License-Identifier: MPL-2.0

package geom3d

import "math"

// OrientationResult is the tri-state outcome of TriangleOrientationIsEqual.
type OrientationResult int

const (
	// OrientationNoCommonEdge means the two triangles share no common edge (C return -1).
	OrientationNoCommonEdge OrientationResult = -1
	// OrientationOpposite means their winding orders disagree across the shared edge (C return 0).
	OrientationOpposite OrientationResult = 0
	// OrientationEqual means their winding orders agree (C return 1).
	OrientationEqual OrientationResult = 1
)

// TriangleOrientationIsEqual reports whether two triangles (given as index triples into
// points) wind consistently across their shared edge. Port of triangle.c
// triangle_orientation_is_equal, used by the mesher to make surface normals consistent.
//
// The method: at the common vertex, build the in-plane "upward" normal from the vector to
// each triangle's centroid crossed with the shared edge, and compare its sign against each
// triangle's face normal. Equal signs ⇒ equal orientation.
func TriangleOrientationIsEqual(idx1, idx2 int, triangles [][3]uint64, points []Vec3) OrientationResult {
	t1 := triangles[idx1]
	t2 := triangles[idx2]
	tri1 := Triangle{points[t1[0]], points[t1[1]], points[t1[2]]}
	tri2 := Triangle{points[t2[0]], points[t2[1]], points[t2[2]]}

	mid1 := centroid(tri1)
	mid2 := centroid(tri2)
	normal1 := Normal3D(tri1)
	normal2 := Normal3D(tri2)

	// Find a shared vertex (the C keeps the LAST matching pair, so iterate in the same order).
	commonI, commonJ, found := -1, -1, false
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if t1[i] == t2[j] {
				commonI, commonJ, found = i, j, true
			}
		}
	}
	if !found {
		return OrientationNoCommonEdge
	}
	commonVertex := points[t1[commonI]]

	// The four edges from the common vertex: two from each triangle (excluding the common vertex).
	var edges [4]Vec3
	n := 0
	for i := 0; i < 3; i++ {
		if i == commonI {
			continue
		}
		edges[n] = sub(points[t1[i]], commonVertex)
		n++
	}
	for j := 0; j < 3; j++ {
		if j == commonJ {
			continue
		}
		edges[n] = sub(points[t2[j]], commonVertex)
		n++
	}
	for e := range edges {
		edges[e] = Normalize3D(edges[e])
	}

	// The common edge: a tri-1 edge parallel to a tri-2 edge (dot product ≈ 1).
	var commonEdge Vec3
	haveCommon := false
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			if math.Abs(Dot3D(edges[i], edges[2+j])-1.0) < 1e-14 {
				commonEdge = edges[i]
				haveCommon = true
			}
		}
	}
	if !haveCommon {
		return OrientationNoCommonEdge
	}

	toMid1 := sub(mid1, commonVertex)
	toMid2 := sub(mid2, commonVertex)
	up1 := CrossProduct3D(toMid1, commonEdge)
	up2 := CrossProduct3D(commonEdge, toMid2)

	dot1 := Dot3D(up1, normal1)
	dot2 := Dot3D(up2, normal2)
	if (dot1 > 0.0) == (dot2 > 0.0) {
		return OrientationEqual
	}
	return OrientationOpposite
}

func centroid(t Triangle) Vec3 {
	return Vec3{
		(t[0][0] + t[1][0] + t[2][0]) / 3.0,
		(t[0][1] + t[1][1] + t[2][1]) / 3.0,
		(t[0][2] + t[1][2] + t[2][2]) / 3.0,
	}
}

func sub(a, b Vec3) Vec3 { return Vec3{a[0] - b[0], a[1] - b[1], a[2] - b[2]} }
