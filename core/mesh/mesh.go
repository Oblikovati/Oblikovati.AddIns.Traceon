// SPDX-License-Identifier: MPL-2.0

// Package mesh is the pure-Go port of Traceon's radial (2D line-element) mesh container
// and the topology operations the parametric mesher relies on: point deduplication
// (so adjacent path segments share nodes), connected-component grouping, and line
// orientation/outward-normal consistency. The 3D triangle machinery in upstream
// mesher.py is closed-source "Traceon Pro" and is intentionally not ported.
//
// A Mesh holds a flat point array and a line-element index array; each element is two
// indices (a straight "line") or four ("line4": [start, end, node@1/3, node@2/3]) — the
// same node ordering the radial solver consumes (see core/radial.reorder).
package mesh

import (
	"math"
	"sort"

	"oblikovati.org/traceon/core/geom2d"
)

// Mesh is a radial boundary mesh: points in the (x=r, y=0, z) plane plus line elements
// indexing into them. PhysicalToLines maps an electrode name to the indices of the lines
// that belong to it (mirrors Traceon's physical_to_lines).
type Mesh struct {
	Points          []geom2d.Vertex
	Lines           [][]int
	PhysicalToLines map[string][]int
}

// mantissaMask zeroes the low 16 bits of a float64's 52-bit mantissa, snapping points
// that agree to ~2^-37 relative precision to identical bit patterns so they deduplicate.
// This is Traceon's _deduplicate_points quantization (mesher.py): reinterpret each double
// as uint64 and AND with 0xFFFFFFFFFFFF0000.
const mantissaMask = uint64(0xFFFFFFFFFFFF0000)

// snap quantizes one coordinate by masking off the low mantissa bits.
func snap(x float64) float64 {
	return math.Float64frombits(math.Float64bits(x) & mantissaMask)
}

// DeduplicatePoints merges coincident points (after mantissa-snapping) and remaps the
// element index array, returning the surviving points and the remapped elements. Points
// are returned in lexicographic (z, y, x) ascending order — the upstream lexsort key —
// and carry the snapped coordinates of the survivor. Elements may be any arity (2 for
// lines, 4 for line4); their integer indices are simply remapped.
func DeduplicatePoints(points []geom2d.Vertex, elements [][]int) ([]geom2d.Vertex, [][]int) {
	if len(points) == 0 {
		return points, elements
	}
	snapped := make([]geom2d.Vertex, len(points))
	for i, p := range points {
		snapped[i] = geom2d.Vertex{snap(p[0]), snap(p[1]), snap(p[2])}
	}

	order := make([]int, len(points))
	for i := range order {
		order[i] = i
	}
	// Stable lexicographic sort with z as the primary key, then y, then x — numpy
	// lexsort(points.T) uses the LAST row (z) as primary. Stable keeps duplicates in
	// their original relative order, matching numpy's stable lexsort.
	sort.SliceStable(order, func(a, b int) bool {
		return lessZYX(snapped[order[a]], snapped[order[b]])
	})

	newPoints := make([]geom2d.Vertex, 0, len(points))
	oldToNew := make([]int, len(points))
	for rank, oldIdx := range order {
		if rank == 0 || !equalVertex(snapped[oldIdx], snapped[order[rank-1]]) {
			newPoints = append(newPoints, snapped[oldIdx])
		}
		oldToNew[oldIdx] = len(newPoints) - 1
	}

	newElems := make([][]int, len(elements))
	for i, el := range elements {
		ne := make([]int, len(el))
		for j, idx := range el {
			ne[j] = oldToNew[idx]
		}
		newElems[i] = ne
	}
	return newPoints, newElems
}

// lessZYX orders vertices by z, then y, then x (all ascending).
func lessZYX(a, b geom2d.Vertex) bool {
	if a[2] != b[2] {
		return a[2] < b[2]
	}
	if a[1] != b[1] {
		return a[1] < b[1]
	}
	return a[0] < b[0]
}

func equalVertex(a, b geom2d.Vertex) bool {
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]
}

// ConnectedElements groups elements into connected components: two elements are connected
// when they share at least one vertex index. Returns one ascending index slice per
// component, in order of first appearance (mirrors mesher._get_connected_elements). Works
// for any element arity (lines or triangles).
func ConnectedElements(elements [][]int) [][]int {
	vertexToElems := map[int][]int{}
	for i, el := range elements {
		for _, v := range el {
			vertexToElems[v] = append(vertexToElems[v], i)
		}
	}

	labels := make([]int, len(elements))
	for i := range labels {
		labels[i] = -1
	}
	var components [][]int
	for start := range elements {
		if labels[start] != -1 {
			continue
		}
		comp := len(components)
		labels[start] = comp
		queue := []int{start}
		for len(queue) > 0 {
			cur := queue[0]
			queue = queue[1:]
			for _, v := range elements[cur] {
				for _, nb := range vertexToElems[v] {
					if labels[nb] == -1 {
						labels[nb] = comp
						queue = append(queue, nb)
					}
				}
			}
		}
		components = append(components, nil)
	}

	for i, lbl := range labels {
		components[lbl] = append(components[lbl], i)
	}
	return components
}

// LineOrientationEqual reports whether lines[i] and lines[j] are head-to-tail consistent
// (so their 2D normals agree): true when line i's end equals line j's start, or line j's
// end equals line i's start. Only the first two (endpoint) indices are used, so it applies
// to both 2-node and 4-node lines. Mirrors mesher._line_orientation_equal.
func LineOrientationEqual(i, j int, lines [][]int) bool {
	p1, p2 := lines[i][0], lines[i][1]
	n1, n2 := lines[j][0], lines[j][1]
	return p2 == n1 || n2 == p1
}
