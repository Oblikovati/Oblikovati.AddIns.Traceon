// SPDX-License-Identifier: MPL-2.0

package mesh

import "oblikovati.org/traceon/core/geom2d"

// New builds a radial mesh from raw points and line elements: it deduplicates coincident
// points (so adjacent segments share nodes) and, when ensureOutward is set, makes each
// named electrode's line normals consistent and outward-pointing — the sign convention the
// BEM solver assumes. Mirrors mesher.Mesh.__init__ for the radial (line-only) case.
func New(points []geom2d.Vertex, lines [][]int, physicalToLines map[string][]int, ensureOutward bool) *Mesh {
	p, l := DeduplicatePoints(points, lines)
	m := &Mesh{Points: p, Lines: l, PhysicalToLines: physicalToLines}
	if ensureOutward {
		for _, idxs := range physicalToLines {
			m.ensureLineOrientation(idxs, true)
		}
	}
	return m
}

// EnsureInwardNormals makes the named electrode's line normals consistent and inward-pointing
// — the convention an electrostatic/magnetostatic boundary (a zero-constant dielectric/
// magnetizable element) uses. A no-op if the electrode is absent. Mirrors
// mesher.Mesh.ensure_inward_normals for the radial (line-only) case.
func (m *Mesh) EnsureInwardNormals(electrode string) {
	if idxs, ok := m.PhysicalToLines[electrode]; ok {
		m.ensureLineOrientation(idxs, false)
	}
}

// ensureLineOrientation makes the given electrode's lines consistently oriented and then
// flips the whole group when its winding does not match the wanted direction, so the normals
// end up outward (outward=true) or inward (outward=false). Mirrors
// mesher._ensure_line_orientation with should_be_outwards = outward.
func (m *Mesh) ensureLineOrientation(group []int, outward bool) {
	if len(group) == 0 {
		return
	}
	m.reorientGroup(group)
	if m.normalsOutward(group) != outward {
		for _, i := range group {
			m.Lines[i] = flipLine(m.Lines[i])
		}
	}
}

// reorientGroup flood-fills the group from its first line, flipping each shared-vertex
// neighbour that is head-to-head rather than head-to-tail, so every line in a connected
// electrode winds the same way. Mirrors mesher._reorient_lines restricted to one group.
func (m *Mesh) reorientGroup(group []int) {
	oriented := make(map[int]bool, len(group))
	vertexToLines := map[int][]int{}
	for _, i := range group {
		vertexToLines[m.Lines[i][0]] = append(vertexToLines[m.Lines[i][0]], i)
		vertexToLines[m.Lines[i][1]] = append(vertexToLines[m.Lines[i][1]], i)
	}
	start := group[0]
	oriented[start] = true
	queue := []int{start}
	for len(queue) > 0 {
		cur := queue[0]
		queue = queue[1:]
		for _, v := range m.Lines[cur][:2] {
			for _, nb := range vertexToLines[v] {
				if oriented[nb] {
					continue
				}
				if !LineOrientationEqual(cur, nb, m.Lines) {
					m.Lines[nb] = flipLine(m.Lines[nb])
				}
				oriented[nb] = true
				queue = append(queue, nb)
			}
		}
	}
}

// normalsOutward reports whether the group's line normals point outward, using the
// divergence-theorem flux of the field (x, 0): ∮ x·n_x dl equals the enclosed area and is
// positive exactly when the normals point outward. Mirrors
// mesher._are_line_normals_pointing_outwards.
func (m *Mesh) normalsOutward(group []int) bool {
	sum := 0.0
	for _, i := range group {
		a, b := m.Points[m.Lines[i][0]], m.Points[m.Lines[i][1]]
		midX := (a[0] + b[0]) / 2
		n := geom2d.Normal2D(geom2d.Point2{a[0], a[2]}, geom2d.Point2{b[0], b[2]})
		length := geom2d.Length2D(geom2d.Point2{a[0], a[2]}, geom2d.Point2{b[0], b[2]})
		sum += midX * n[0] * length
	}
	return sum > 0
}

// flipLine reverses a line's traversal direction: [start, end] for a 2-node line, or
// [start, end, node@1/3, node@2/3] → [end, start, node@2/3, node@1/3] for a line4 (so the
// interior nodes stay in along-curve order after the flip). Mirrors the flips in
// mesher._reorient_lines.
func flipLine(line []int) []int {
	if len(line) == 4 {
		return []int{line[1], line[0], line[3], line[2]}
	}
	return []int{line[1], line[0]}
}
