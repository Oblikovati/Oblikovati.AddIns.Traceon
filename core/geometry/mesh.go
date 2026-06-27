// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"fmt"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/mesh"
	"oblikovati.org/traceon/core/radial"
)

// MeshOptions controls how a Path is discretized into line elements. Exactly one of
// MeshSize / MeshSizeFactor must be > 0 (MeshSize takes precedence). HigherOrder produces
// curved 4-node "line4" elements (the radial solver's higher-order elements). Name tags
// the produced electrode (falling back to the Path's own name).
type MeshOptions struct {
	MeshSize             float64
	MeshSizeFactor       float64
	HigherOrder          bool
	Name                 string
	EnsureOutwardNormals bool
}

// Mesh discretizes the path into line elements and returns a deduplicated radial mesh.
// Flat meshing yields straight 2-node lines; HigherOrder yields curved 4-node line4
// elements ordered [start, end, node@1/3, node@2/3]. Mirrors geometry.Path.mesh.
func (p Path) Mesh(opts MeshOptions) *mesh.Mesh {
	nFactor := 1
	if opts.HigherOrder {
		nFactor = 3
	}
	u := discretizePath(p.Length, p.Breakpoints, opts.MeshSize, opts.MeshSizeFactor, nFactor)

	points := make([]geom2d.Vertex, len(u))
	for i, t := range u {
		points[i] = p.Fun(t)
	}
	lines := elementIndices(len(u), opts.HigherOrder)

	name := opts.Name
	if name == "" {
		name = p.Name
	}
	physical := map[string][]int{}
	if name != "" {
		idxs := make([]int, len(lines))
		for i := range idxs {
			idxs[i] = i
		}
		physical[name] = idxs
	}
	return mesh.New(points, lines, physical, opts.EnsureOutwardNormals)
}

// elementIndices builds the line index array for n sampled points: straight [k, k+1] pairs
// when flat, or line4 [start, end, node@1/3, node@2/3] quads when higher order (requiring
// n ≡ 1 mod 3, the invariant discretizePath guarantees for nFactor=3).
func elementIndices(n int, higherOrder bool) [][]int {
	if !higherOrder {
		lines := make([][]int, 0, n-1)
		for k := 0; k+1 < n; k++ {
			lines = append(lines, []int{k, k + 1})
		}
		return lines
	}
	if n%3 != 1 {
		panic(fmt.Sprintf("higher-order mesh needs (#points) ≡ 1 mod 3, got %d", n))
	}
	lines := make([][]int, 0, (n-1)/3)
	for s := 0; s+3 < n; s += 3 {
		lines = append(lines, []int{s, s + 3, s + 1, s + 2})
	}
	return lines
}

// RadialLines converts a meshed path into the 4-node line elements the radial BEM solver
// consumes ([start, end, node@1/3, node@2/3]). A flat (2-node) mesh is promoted by placing
// the two interior nodes at 1/3 and 2/3 along each straight segment, matching upstream
// mesher._lines_to_higher_order.
func RadialLines(m *mesh.Mesh) []radial.Line {
	out := make([]radial.Line, len(m.Lines))
	for i, l := range m.Lines {
		start, end := m.Points[l[0]], m.Points[l[1]]
		var n1, n2 geom2d.Vertex
		if len(l) == 4 {
			n1, n2 = m.Points[l[2]], m.Points[l[3]]
		} else {
			n1 = lerp(start, end, 1.0/3.0)
			n2 = lerp(start, end, 2.0/3.0)
		}
		out[i] = radial.Line{start, end, n1, n2}
	}
	return out
}

// lerp linearly interpolates between two vertices at parameter f.
func lerp(a, b geom2d.Vertex, f float64) geom2d.Vertex {
	return geom2d.Vertex{
		a[0] + (b[0]-a[0])*f,
		a[1] + (b[1]-a[1])*f,
		a[2] + (b[2]-a[2])*f,
	}
}
