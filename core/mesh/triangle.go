// SPDX-License-Identifier: MPL-2.0

package mesh

import "oblikovati.org/traceon/core/geom3d"

// degenerateAreaTol drops a triangle whose area falls below it (in the mesh length unit²) as
// degenerate — the disk centre and revolution seams collapse to slivers there. Mirrors
// mesher._remove_degenerate_triangles (areas < 1e-12).
const degenerateAreaTol = 1e-12

// TriangleMesh is a planar surface mesh: 3D points plus triangle index triples, grouped by
// name. Radial current coils are meshed this way (a solid cross-section in the xz-plane). It
// is the triangle counterpart of the line Mesh; the two are kept separate because coil
// triangles feed the current pre-field while boundary/material lines feed the BEM matrix.
// Port of the triangle half of mesher.Mesh.
type TriangleMesh struct {
	Points              []geom3d.Vec3
	Triangles           [][3]int
	PhysicalToTriangles map[string][]int
}

// NewTriangleMesh builds a triangle mesh: it drops degenerate (near-zero-area) triangles —
// remapping the physical groups — then deduplicates coincident points (remapping triangles).
// Mirrors mesher.Mesh.__init__ for the triangle case (_remove_degenerate_triangles then
// _deduplicate_points).
func NewTriangleMesh(points []geom3d.Vec3, triangles [][3]int, physical map[string][]int) *TriangleMesh {
	tris, phys := removeDegenerateTriangles(points, triangles, physical)
	dedupPoints, dedupTris := deduplicateTrianglePoints(points, tris)
	return &TriangleMesh{Points: dedupPoints, Triangles: dedupTris, PhysicalToTriangles: phys}
}

// Group returns the named physical group's triangles as geom3d.Triangle (actual coordinates),
// or nil if the group is absent — the form the current-ring builder consumes.
func (m *TriangleMesh) Group(name string) []geom3d.Triangle {
	idxs, ok := m.PhysicalToTriangles[name]
	if !ok {
		return nil
	}
	out := make([]geom3d.Triangle, len(idxs))
	for k, i := range idxs {
		t := m.Triangles[i]
		out[k] = geom3d.Triangle{m.Points[t[0]], m.Points[t[1]], m.Points[t[2]]}
	}
	return out
}

// removeDegenerateTriangles drops triangles with area below degenerateAreaTol and remaps each
// physical group to the surviving triangles' new indices. Mirrors _remove_degenerate_triangles.
func removeDegenerateTriangles(points []geom3d.Vec3, triangles [][3]int, physical map[string][]int) ([][3]int, map[string][]int) {
	mapIndex := make([]int, len(triangles)) // old index → new index, or -1 if dropped
	newTris := make([][3]int, 0, len(triangles))
	for i, t := range triangles {
		if geom3d.TriangleArea(points[t[0]], points[t[1]], points[t[2]]) < degenerateAreaTol {
			mapIndex[i] = -1
			continue
		}
		mapIndex[i] = len(newTris)
		newTris = append(newTris, t)
	}

	newPhys := make(map[string][]int, len(physical))
	for k, idxs := range physical {
		var kept []int
		for _, i := range idxs {
			if mapIndex[i] != -1 {
				kept = append(kept, mapIndex[i])
			}
		}
		newPhys[k] = kept
	}
	return newTris, newPhys
}

// deduplicateTrianglePoints merges coincident points (reusing the line-mesh dedup, since a
// triangle is just a 3-index element) and remaps the triangle indices.
func deduplicateTrianglePoints(points []geom3d.Vec3, triangles [][3]int) ([]geom3d.Vec3, [][3]int) {
	elems := make([][]int, len(triangles))
	for i, t := range triangles {
		elems[i] = []int{t[0], t[1], t[2]}
	}
	newPoints, newElems := DeduplicatePoints(points, elems)
	newTris := make([][3]int, len(newElems))
	for i, e := range newElems {
		newTris[i] = [3]int{e[0], e[1], e[2]}
	}
	return newPoints, newTris
}
