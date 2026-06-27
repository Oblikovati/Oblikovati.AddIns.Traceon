// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"math"
	"testing"
)

// TestRevolveYDiskMeridian checks RevolveY's circumferential length: revolving the radius line
// gives length2 = 2π·r_avg with r_avg = radius/2 (the arc-length mean of r along the line).
func TestRevolveYDiskMeridian(t *testing.T) {
	radius := 2.0
	s := Line(Point{0, 0, 0}, Point{radius, 0, 0}).RevolveY(2 * math.Pi)
	if !close(s.PathLength1, radius) {
		t.Errorf("PathLength1 = %g, want %g (the radius)", s.PathLength1, radius)
	}
	if want := 2 * math.Pi * (radius / 2); !close(s.PathLength2, want) {
		t.Errorf("PathLength2 = %g, want %g (2π·r_avg)", s.PathLength2, want)
	}
	// The far edge (u = radius) revolved a quarter turn lands on the +z axis at r = radius.
	p := s.At(radius, s.PathLength2/4)
	if !closeVec(p, Point{0, 0, radius}) {
		t.Errorf("disk rim at quarter turn = %v, want [0 0 %g]", p, radius)
	}
}

// flatPlane is a flat rectangular surface in the xz-plane with a breakpoint in each direction,
// so meshing it exercises the multi-section path: per-section grids plus edge welding.
func flatPlane(w, h float64) Surface {
	return Surface{
		Fun:          func(u, v float64) Point { return Point{u, 0, v} },
		PathLength1:  w,
		PathLength2:  h,
		Breakpoints1: []float64{w / 2},
		Breakpoints2: []float64{h / 2},
		Name:         "plate",
	}
}

// TestMeshSurfaceSectionsWeld checks a surface with breakpoints (four sections) meshes into a
// valid, non-empty triangle mesh, exercising the section-welding path. Exact equivalence to
// upstream is asserted by the plate_sections case in TestSurfaceMesherGolden.
func TestMeshSurfaceSectionsWeld(t *testing.T) {
	m := flatPlane(2.0, 1.0).Mesh(0.5)

	if len(m.Triangles) == 0 {
		t.Fatal("plate meshed to zero triangles")
	}
	for _, tri := range m.Triangles {
		for _, idx := range tri {
			if idx < 0 || idx >= len(m.Points) {
				t.Fatalf("triangle index %d out of range [0,%d)", idx, len(m.Points))
			}
		}
	}
}

// TestMeshSurfaceGroupNamed checks two surfaces mesh into one triangle mesh with a group each.
func TestMeshSurfaceGroupNamed(t *testing.T) {
	a := DiskXZ(0, 5e-3, 1e-3).WithName("coil1")
	b := DiskXZ(0, -5e-3, 1e-3).WithName("coil2")
	m := MeshSurfaceGroup([]Surface{a, b}, 0.25e-3)

	for _, name := range []string{"coil1", "coil2"} {
		if len(m.PhysicalToTriangles[name]) == 0 {
			t.Errorf("group %q has no triangles; have %v", name, keysOf(m.PhysicalToTriangles))
		}
	}
	total := len(m.PhysicalToTriangles["coil1"]) + len(m.PhysicalToTriangles["coil2"])
	if total != len(m.Triangles) {
		t.Errorf("grouped triangles %d != total %d", total, len(m.Triangles))
	}
}

func keysOf(m map[string][]int) []string {
	out := make([]string, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	return out
}
