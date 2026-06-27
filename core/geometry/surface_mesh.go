// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"math"

	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/mesh"
)

// startDepth is the initial uniform subdivision level (2^startDepth cells per side) the
// adaptive mesher refines from. Mirrors the _mesh default start_depth=2.
const startDepth = 2

// Mesh discretizes the surface into triangles at the given element size and returns a
// deduplicated triangle mesh (one physical group named after the surface). Mirrors
// geometry.Surface.mesh / mesher._mesh for an explicit mesh size.
func (s Surface) Mesh(meshSize float64) *mesh.TriangleMesh {
	return meshSurfaces([]Surface{s}, meshSize)
}

// MeshByFactor meshes the surface at an element size derived from the mesh-size factor:
// min(PathLength1, PathLength2)/4 / √factor — the dimension-relative sizing upstream uses when
// only mesh_size_factor is given. Mirrors Surface.mesh's mesh_size_factor branch.
func (s Surface) MeshByFactor(factor float64) *mesh.TriangleMesh {
	meshSize := math.Min(s.PathLength1, s.PathLength2) / 4
	if factor > 0 {
		meshSize /= math.Sqrt(factor)
	}
	return s.Mesh(meshSize)
}

// MeshSurfaceGroup meshes several named surfaces into one triangle mesh, each surface its own
// physical group. Coincident points across surfaces are merged by the dedup. Mirrors meshing a
// surface collection: (a + b).mesh(...).
func MeshSurfaceGroup(surfaces []Surface, meshSize float64) *mesh.TriangleMesh {
	return meshSurfaces(surfaces, meshSize)
}

// meshSurfaces meshes each surface independently, concatenates the triangles into one
// point-shared mesh (group per surface name), and deduplicates. Each surface goes through the
// adaptive quad-subdivision → triangle pipeline (mesher._mesh).
func meshSurfaces(surfaces []Surface, meshSize float64) *mesh.TriangleMesh {
	var allPoints []geom3d.Vec3
	var allTris [][3]int
	physical := map[string][]int{}
	for _, s := range surfaces {
		points, tris := meshOneSurface(s, meshSize, len(allPoints))
		base := len(allTris)
		allPoints = append(allPoints, points...)
		allTris = append(allTris, tris...)
		if s.Name != "" {
			idxs := make([]int, len(tris))
			for i := range idxs {
				idxs[i] = base + i
			}
			physical[s.Name] = append(physical[s.Name], idxs...)
		}
	}
	return mesh.NewTriangleMesh(allPoints, allTris, physical)
}

// meshOneSurface runs the adaptive quad-subdivision mesher on one surface and returns its
// freshly-created points (to be appended at offset pointOffset in the combined mesh) and its
// triangles (already offset). Mirrors mesher._mesh for a single surface.
func meshOneSurface(s Surface, meshSize float64, pointOffset int) ([]geom3d.Vec3, [][3]int) {
	points, stacks, quadsList := meshSubsectionsToQuads(s, meshSize)

	maxDepth := 0
	for _, ps := range stacks {
		if d := ps.depth(); d > maxDepth {
			maxDepth = d
		}
	}
	pqs := make([]*pointsWithQuads, len(stacks))
	for i, ps := range stacks {
		pqs[i] = ps.normalizeToDepth(maxDepth, quadsList[i])
	}

	// Weld shared section edges so abutting sub-surface grids reference one point per node.
	nx := len(s.Breakpoints1) + 1
	ny := len(s.Breakpoints2) + 1
	for i := 0; i < nx-1; i++ {
		for j := 0; j < ny; j++ {
			copyOverEdge(rowPtrs(pqs[j*nx+i], -1), rowPtrs(pqs[j*nx+i+1], 0))
		}
	}
	for i := 0; i < nx; i++ {
		for j := 0; j < ny-1; j++ {
			copyOverEdge(colPtrs(pqs[j*nx+i], -1), colPtrs(pqs[(j+1)*nx+i], 0))
		}
	}

	var tris [][3]int
	for _, pq := range pqs {
		tris = append(tris, pq.toTriangles()...)
	}
	// Offset the triangle indices to the combined point array.
	for i := range tris {
		tris[i] = [3]int{tris[i][0] + pointOffset, tris[i][1] + pointOffset, tris[i][2] + pointOffset}
	}
	return points, tris
}

// quad is one un-split cell pending triangulation: (depth, i0, i1, j0, j1) on the depth grid.
type cell = [5]int

// pointStack lazily materialises a surface's sample points on a hierarchy of 2^depth+1 grids,
// creating a point the first time a grid node is touched. Points are shared (by pointer) across
// a surface's sections so the combined array grows once. Port of mesher.PointStack.
type pointStack struct {
	points  *[]geom3d.Vec3
	surf    Surface
	indices [][][]int // indices[depth][i][j] = point index, or -1 if not yet created
}

func (ps *pointStack) numIndices(depth int) int { return (1 << depth) + 1 }

func (ps *pointStack) indexToPoint(depth, i, j int) geom3d.Vec3 {
	n := float64(ps.numIndices(depth) - 1)
	u := ps.surf.PathLength1 / n * float64(i)
	v := ps.surf.PathLength2 / n * float64(j)
	return ps.surf.Fun(u, v)
}

func (ps *pointStack) addLevel() {
	newDepth := len(ps.indices)
	n := ps.numIndices(newDepth)
	grid := make([][]int, n)
	for i := range grid {
		grid[i] = make([]int, n)
		for j := range grid[i] {
			grid[i][j] = -1
		}
	}
	if newDepth != 0 {
		prev := ps.indices[newDepth-1]
		for i := range prev {
			for j := range prev[i] {
				grid[2*i][2*j] = prev[i][j]
			}
		}
	}
	ps.indices = append(ps.indices, grid)
}

func (ps *pointStack) toPointIndex(depth, i, j int) int {
	for depth >= len(ps.indices) {
		ps.addLevel()
	}
	if ps.indices[depth][i][j] == -1 {
		*ps.points = append(*ps.points, ps.indexToPoint(depth, i, j))
		ps.indices[depth][i][j] = len(*ps.points) - 1
	}
	return ps.indices[depth][i][j]
}

func (ps *pointStack) point(depth, i, j int) geom3d.Vec3 {
	return (*ps.points)[ps.toPointIndex(depth, i, j)]
}

func (ps *pointStack) depth() int { return len(ps.indices) - 1 }

// normalizeToDepth brings the stack to the common depth, propagates every coarse-level point
// down into the finest grid, and rescales the quads onto that grid — so all sub-surfaces share
// one index resolution before triangulation. Port of PointStack.normalize_to_depth.
func (ps *pointStack) normalizeToDepth(depth int, quads []cell) *pointsWithQuads {
	for ps.depth() < depth {
		ps.addLevel()
	}
	for d := startDepth; d < len(ps.indices)-1; d++ {
		prev := ps.indices[d]
		for x := range prev {
			for y := range prev[x] {
				if prev[x][y] != -1 {
					ps.indices[d+1][2*x][2*y] = prev[x][y]
				}
			}
		}
	}
	scaled := make([]cell, len(quads))
	for i, q := range quads {
		qd, i0, i1, j0, j1 := q[0], q[1], q[2], q[3], q[4]
		for qd < depth {
			i0, i1, j0, j1, qd = 2*i0, 2*i1, 2*j0, 2*j1, qd+1
		}
		scaled[i] = cell{qd, i0, i1, j0, j1}
	}
	return &pointsWithQuads{indices: ps.indices[len(ps.indices)-1], quads: scaled}
}

// subdivideQuads adaptively splits a seed cell until each leaf's edges are within the mesh
// size, accumulating leaves into quads. A cell splits along a direction whose edge length
// exceeds the mesh size (or is far more elongated than the cross direction). Port of
// mesher._subdivide_quads with a constant mesh size.
func subdivideQuads(ps *pointStack, meshSize float64, seed cell, quads *[]cell) {
	toSubdivide := []cell{seed}
	for len(toSubdivide) > 0 {
		q := toSubdivide[len(toSubdivide)-1]
		toSubdivide = toSubdivide[:len(toSubdivide)-1]
		depth, i0, i1, j0, j1 := q[0], q[1], q[2], q[3], q[4]

		p1 := ps.point(depth, i0, j0)
		p2 := ps.point(depth, i0, j1)
		p3 := ps.point(depth, i1, j0)
		p4 := ps.point(depth, i1, j1)
		horizontal := math.Max(dist3(p1, p2), dist3(p3, p4))
		vertical := math.Max(dist3(p1, p3), dist3(p2, p4))

		h := horizontal > meshSize || (horizontal > 2.5*vertical && horizontal > meshSize/8)
		v := vertical > meshSize || (vertical > 2.5*horizontal && vertical > meshSize/8)

		switch {
		case h && v:
			toSubdivide = append(toSubdivide,
				cell{depth + 1, 2 * i0, 2*i0 + 1, 2 * j0, 2*j0 + 1},
				cell{depth + 1, 2 * i0, 2*i0 + 1, 2*j0 + 1, 2*j0 + 2},
				cell{depth + 1, 2*i0 + 1, 2*i0 + 2, 2 * j0, 2*j0 + 1},
				cell{depth + 1, 2*i0 + 1, 2*i0 + 2, 2*j0 + 1, 2*j0 + 2})
		case h:
			toSubdivide = append(toSubdivide,
				cell{depth + 1, 2 * i0, 2 * i1, 2 * j0, 2*j0 + 1},
				cell{depth + 1, 2 * i0, 2 * i1, 2*j0 + 1, 2*j0 + 2})
		case v:
			toSubdivide = append(toSubdivide,
				cell{depth + 1, 2 * i0, 2*i0 + 1, 2 * j0, 2 * j1},
				cell{depth + 1, 2*i0 + 1, 2*i0 + 2, 2 * j0, 2 * j1})
		default:
			*quads = append(*quads, q)
		}
	}
}

// meshSubsectionsToQuads runs the adaptive subdivision over every section of the surface,
// sharing one growing point array, and returns the points plus per-section stacks and quads.
// Port of mesher._mesh_subsections_to_quads.
func meshSubsectionsToQuads(surface Surface, meshSize float64) ([]geom3d.Vec3, []*pointStack, [][]cell) {
	points := []geom3d.Vec3{}
	var stacks []*pointStack
	var allQuads [][]cell
	for _, s := range surface.sections() {
		ps := &pointStack{points: &points, surf: s}
		var quads []cell
		n := ps.numIndices(startDepth)
		for i := 0; i < n-1; i++ {
			for j := 0; j < n-1; j++ {
				subdivideQuads(ps, meshSize, cell{startDepth, i, i + 1, j, j + 1}, &quads)
			}
		}
		stacks = append(stacks, ps)
		allQuads = append(allQuads, quads)
	}
	return points, stacks, allQuads
}

// pointsWithQuads is one section's finest-resolution index grid plus its leaf quads, the input
// to triangulation. Port of mesher.Points3DWithQuads.
type pointsWithQuads struct {
	indices [][]int
	quads   []cell
}

// toTriangles converts each leaf quad into two triangles, or three when a finer neighbour put a
// point on one edge (a T-junction), so the mesh stays watertight. Port of
// Points3DWithQuads.to_triangles.
func (pq *pointsWithQuads) toTriangles() [][3]int {
	var tris [][3]int
	add := func(a, b, c [2]int) {
		tris = append(tris, [3]int{pq.indices[a[0]][a[1]], pq.indices[b[0]][b[1]], pq.indices[c[0]][c[1]]})
	}
	for _, q := range pq.quads {
		i0, i1, j0, j1 := q[1], q[2], q[3], q[4]
		p0, p1, p2, p3 := [2]int{i0, j0}, [2]int{i0, j1}, [2]int{i1, j1}, [2]int{i1, j0}

		split := false
		for edge := 0; edge < 4; edge++ {
			mid := [2]int{(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2}
			if (absInt(p0[0]-p1[0]) > 1 || absInt(p0[1]-p1[1]) > 1) && pq.indices[mid[0]][mid[1]] != -1 {
				add(p0, mid, p3)
				add(mid, p2, p3)
				add(mid, p1, p2)
				split = true
				break
			}
			p0, p1, p2, p3 = p1, p2, p3, p0
		}
		if !split {
			add(p0, p1, p2)
			add(p0, p2, p3)
		}
	}
	return tris
}

// rowPtrs returns pointers to row i of the index grid (i = -1 means the last row), and colPtrs
// the same for a column — the borders welded by copyOverEdge.
func rowPtrs(pq *pointsWithQuads, i int) []*int {
	if i < 0 {
		i += len(pq.indices)
	}
	ptrs := make([]*int, len(pq.indices))
	for j := range pq.indices[i] {
		ptrs[j] = &pq.indices[i][j]
	}
	return ptrs
}

func colPtrs(pq *pointsWithQuads, j int) []*int {
	if j < 0 {
		j += len(pq.indices)
	}
	ptrs := make([]*int, len(pq.indices))
	for i := range pq.indices {
		ptrs[i] = &pq.indices[i][j]
	}
	return ptrs
}

// copyOverEdge welds two abutting section borders: a node set on either side is shared by both.
// Port of mesher._copy_over_edge.
func copyOverEdge(e1, e2 []*int) {
	for k := range e1 {
		if *e2[k] != -1 {
			*e1[k] = *e2[k]
		}
	}
	for k := range e1 {
		if *e1[k] != -1 {
			*e2[k] = *e1[k]
		}
	}
}

func dist3(a, b geom3d.Vec3) float64 {
	return math.Sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]))
}

func absInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
