// SPDX-License-Identifier: MPL-2.0

package mesh

import (
	"testing"

	"oblikovati.org/traceon/core/geom3d"
)

// TestNewTriangleMeshDropsDegenerate checks a repeated-vertex (zero-area) triangle is removed
// and the physical group is remapped to the surviving triangle indices.
func TestNewTriangleMeshDropsDegenerate(t *testing.T) {
	points := []geom3d.Vec3{{0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {1, 0, 1}}
	triangles := [][3]int{
		{0, 1, 2}, // good
		{0, 1, 1}, // degenerate (repeated vertex → zero area)
		{1, 3, 2}, // good
	}
	m := NewTriangleMesh(points, triangles, map[string][]int{"coil": {0, 1, 2}})

	if len(m.Triangles) != 2 {
		t.Fatalf("triangles = %d, want 2 (degenerate dropped)", len(m.Triangles))
	}
	// Group must point only at the two surviving triangles, remapped to 0 and 1.
	if got := m.PhysicalToTriangles["coil"]; len(got) != 2 || got[0] != 0 || got[1] != 1 {
		t.Errorf("coil group = %v, want [0 1]", got)
	}
}

// TestNewTriangleMeshDeduplicates checks coincident points are merged and triangle indices
// remapped, so a duplicated point does not survive.
func TestNewTriangleMeshDeduplicates(t *testing.T) {
	// Point 3 duplicates point 0; both triangles should end up sharing the merged point.
	points := []geom3d.Vec3{{0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {0, 0, 0}}
	triangles := [][3]int{{0, 1, 2}, {3, 1, 2}}
	m := NewTriangleMesh(points, triangles, map[string][]int{"c": {0, 1}})

	if len(m.Points) != 3 {
		t.Fatalf("points = %d, want 3 (duplicate merged)", len(m.Points))
	}
	if len(m.Triangles) != 2 {
		t.Fatalf("triangles = %d, want 2", len(m.Triangles))
	}
	// Both triangles now reference the same three distinct points.
	for _, tri := range m.Triangles {
		if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
			t.Errorf("triangle %v has a repeated index after dedup", tri)
		}
	}
}

// TestTriangleMeshGroup checks Group materializes the named triangles as coordinates, and an
// absent group returns nil.
func TestTriangleMeshGroup(t *testing.T) {
	points := []geom3d.Vec3{{0, 0, 0}, {2, 0, 0}, {0, 0, 2}}
	m := NewTriangleMesh(points, [][3]int{{0, 1, 2}}, map[string][]int{"coil": {0}})

	tris := m.Group("coil")
	if len(tris) != 1 {
		t.Fatalf("coil triangles = %d, want 1", len(tris))
	}
	if area := geom3d.TriangleArea(tris[0][0], tris[0][1], tris[0][2]); area != 2 {
		t.Errorf("triangle area = %g, want 2", area)
	}
	if m.Group("missing") != nil {
		t.Error("absent group should return nil")
	}
}
