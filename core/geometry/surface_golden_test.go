// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"sort"
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/mesh"
)

// surfaceMeshGolden mirrors the JSON emitted by tools/gen_fixtures.py::_surface_mesh.
type surfaceMeshGolden struct {
	Cases []struct {
		Name      string       `json:"name"`
		Points    [][]oracle.F `json:"points"`
		Triangles [][]int      `json:"triangles"`
	} `json:"cases"`
}

// buildSurfaceCase reconstructs the Go disk mesh for a fixture case, mirroring the upstream
// construction in _surface_mesh.
func buildSurfaceCase(name string) *mesh.TriangleMesh {
	switch name {
	case "disk_coil":
		return DiskXZ(10e-3, 5e-3, 1e-3).Mesh(0.25e-3)
	case "disk_origin":
		return DiskXZ(0, 0, 1.0).Mesh(0.5)
	case "plate_sections":
		return flatPlane(2.0, 1.0).Mesh(0.5)
	default:
		return nil
	}
}

// TestSurfaceMesherGolden verifies the surface mesher reproduces upstream Surface.mesh: the same
// post-deduplication point coordinates (ordering included) and the same triangles. Triangles are
// compared as sorted vertex-index sets, since the outward-normal reorientation may flip winding.
func TestSurfaceMesherGolden(t *testing.T) {
	var fx surfaceMeshGolden
	oracle.LoadGolden(t, "surface_mesh", &fx)

	for _, c := range fx.Cases {
		t.Run(c.Name, func(t *testing.T) {
			m := buildSurfaceCase(c.Name)
			if m == nil {
				t.Fatalf("no Go builder for case %q", c.Name)
			}
			if len(m.Points) != len(c.Points) {
				t.Fatalf("point count = %d, want %d", len(m.Points), len(c.Points))
			}
			for i, p := range c.Points {
				for k := 0; k < 3; k++ {
					oracle.CheckClose(t, c.Name+" point", m.Points[i][k], p[k].Float())
				}
			}
			if len(m.Triangles) != len(c.Triangles) {
				t.Fatalf("triangle count = %d, want %d", len(m.Triangles), len(c.Triangles))
			}
			if g, w := triangleKeySet(m.Triangles), triangleKeySetInts(c.Triangles); !equalKeySet(g, w) {
				t.Errorf("triangle vertex-index sets differ from upstream")
			}
		})
	}
}

// triangleKeySet builds a multiset of triangles keyed by their sorted vertex indices, so two
// meshes match regardless of per-triangle winding.
func triangleKeySet(tris [][3]int) map[[3]int]int {
	out := map[[3]int]int{}
	for _, t := range tris {
		out[sortedTriple(t[0], t[1], t[2])]++
	}
	return out
}

func triangleKeySetInts(tris [][]int) map[[3]int]int {
	out := map[[3]int]int{}
	for _, t := range tris {
		out[sortedTriple(t[0], t[1], t[2])]++
	}
	return out
}

func sortedTriple(a, b, c int) [3]int {
	s := []int{a, b, c}
	sort.Ints(s)
	return [3]int{s[0], s[1], s[2]}
}

func equalKeySet(a, b map[[3]int]int) bool {
	if len(a) != len(b) {
		return false
	}
	for k, v := range a {
		if b[k] != v {
			return false
		}
	}
	return true
}
