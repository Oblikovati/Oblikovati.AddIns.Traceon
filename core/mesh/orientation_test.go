// SPDX-License-Identifier: MPL-2.0

package mesh

import (
	"reflect"
	"testing"

	"oblikovati.org/traceon/core/geom2d"
)

// squareLoop returns the four corner points (in the x,z meridian plane) and the line index
// loop of a unit square electrode, wound in the given direction.
func squareLoop(ccw bool) ([]geom2d.Vertex, [][]int) {
	pts := []geom2d.Vertex{vtx(1, 0, -1), vtx(2, 0, -1), vtx(2, 0, 1), vtx(1, 0, 1)}
	if ccw {
		return pts, [][]int{{0, 1}, {1, 2}, {2, 3}, {3, 0}}
	}
	return pts, [][]int{{0, 3}, {3, 2}, {2, 1}, {1, 0}}
}

// TestNewOrientsOutward checks that New makes a named electrode's normals point outward
// regardless of the input winding direction (the BEM sign convention).
func TestNewOrientsOutward(t *testing.T) {
	for _, ccw := range []bool{true, false} {
		pts, lines := squareLoop(ccw)
		m := New(pts, lines, map[string][]int{"sq": {0, 1, 2, 3}}, true)
		if !m.normalsOutward([]int{0, 1, 2, 3}) {
			t.Errorf("ccw=%v: normals not outward after New", ccw)
		}
	}
}

// TestNewWithoutEnsureKeepsWinding checks that with ensureOutward off, New only deduplicates
// and leaves the line winding untouched.
func TestNewWithoutEnsureKeepsWinding(t *testing.T) {
	pts, lines := squareLoop(false) // clockwise → normals point inward
	m := New(pts, lines, map[string][]int{"sq": {0, 1, 2, 3}}, false)
	if m.normalsOutward([]int{0, 1, 2, 3}) {
		t.Error("ensureOutward=false should not have reoriented the inward winding")
	}
}

// TestFlipLine checks the 2-node and 4-node line reversals (interior nodes swap so they
// stay in along-curve order).
func TestFlipLine(t *testing.T) {
	if got := flipLine([]int{2, 5}); !reflect.DeepEqual(got, []int{5, 2}) {
		t.Errorf("flipLine([2 5]) = %v, want [5 2]", got)
	}
	if got := flipLine([]int{0, 3, 1, 2}); !reflect.DeepEqual(got, []int{3, 0, 2, 1}) {
		t.Errorf("flipLine([0 3 1 2]) = %v, want [3 0 2 1]", got)
	}
}

// TestEnsureLineOrientationEmpty checks an empty group is a no-op (no panic).
func TestEnsureLineOrientationEmpty(t *testing.T) {
	m := &Mesh{}
	m.ensureLineOrientation(nil, true)
}

// TestEnsureInwardNormals checks that EnsureInwardNormals makes a named electrode's normals
// point inward regardless of the input winding, and is a no-op for an absent electrode.
func TestEnsureInwardNormals(t *testing.T) {
	for _, ccw := range []bool{true, false} {
		pts, lines := squareLoop(ccw)
		m := New(pts, lines, map[string][]int{"sq": {0, 1, 2, 3}}, false)
		m.EnsureInwardNormals("sq")
		if m.normalsOutward([]int{0, 1, 2, 3}) {
			t.Errorf("ccw=%v: normals still outward after EnsureInwardNormals", ccw)
		}
	}
	// Absent electrode is a no-op (must not panic).
	pts, lines := squareLoop(true)
	New(pts, lines, nil, false).EnsureInwardNormals("missing")
}

// TestNewEmpty checks New tolerates an empty mesh.
func TestNewEmpty(t *testing.T) {
	m := New(nil, nil, nil, true)
	if len(m.Points) != 0 || len(m.Lines) != 0 {
		t.Errorf("empty mesh = %d points / %d lines, want 0/0", len(m.Points), len(m.Lines))
	}
}
