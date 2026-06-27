// SPDX-License-Identifier: MPL-2.0

package mesh

import (
	"reflect"
	"testing"

	"oblikovati.org/traceon/core/geom2d"
)

// These tests port the assertions in upstream tests/test_mesher.py so the radial mesh
// topology operations stand alone against the same oracle inputs/outputs.

func vtx(x, y, z float64) geom2d.Vertex { return geom2d.Vertex{x, y, z} }

// TestDeduplicateNoDuplicates ports test_no_duplicates: three distinct points survive in
// (z,y,x) order, which for these inputs is the original order.
func TestDeduplicateNoDuplicates(t *testing.T) {
	points := []geom2d.Vertex{vtx(0, 0, 0), vtx(1, 0, 0), vtx(0, 1, 0)}
	elems := [][]int{{0, 1, 2}}
	got, gotElems := DeduplicatePoints(points, elems)
	want := []geom2d.Vertex{vtx(0, 0, 0), vtx(1, 0, 0), vtx(0, 1, 0)}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("points = %v, want %v", got, want)
	}
	if !reflect.DeepEqual(gotElems, [][]int{{0, 1, 2}}) {
		t.Errorf("elements = %v, want [[0 1 2]]", gotElems)
	}
}

// TestDeduplicateWithDuplicates ports test_with_duplicates: point 2 duplicates point 0.
func TestDeduplicateWithDuplicates(t *testing.T) {
	points := []geom2d.Vertex{vtx(0, 0, 0), vtx(1, 0, 0), vtx(0, 0, 0), vtx(0, 1, 0)}
	elems := [][]int{{0, 1, 3}, {2, 1, 3}}
	got, gotElems := DeduplicatePoints(points, elems)
	want := []geom2d.Vertex{vtx(0, 0, 0), vtx(1, 0, 0), vtx(0, 1, 0)}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("points = %v, want %v", got, want)
	}
	if !reflect.DeepEqual(gotElems, [][]int{{0, 1, 2}, {0, 1, 2}}) {
		t.Errorf("elements = %v, want [[0 1 2] [0 1 2]]", gotElems)
	}
}

// TestDeduplicateClosePoints ports test_close_points: two points 1e-12 apart merge; the
// survivors come out in (z,y,x) ascending order with the snapped coordinate.
func TestDeduplicateClosePoints(t *testing.T) {
	points := []geom2d.Vertex{vtx(0, 2, 3), vtx(0, 0, 1), vtx(0, 0, 1+1e-12), vtx(1, 0, 0)}
	elems := [][]int{{0, 1, 2}}
	got, gotElems := DeduplicatePoints(points, elems)
	if len(got) != 3 {
		t.Fatalf("len(points) = %d, want 3", len(got))
	}
	want := []geom2d.Vertex{vtx(1, 0, 0), vtx(0, 0, 1), vtx(0, 2, 3)}
	for i := range want {
		for k := 0; k < 3; k++ {
			if d := got[i][k] - want[i][k]; d > 1e-5 || d < -1e-5 {
				t.Errorf("points[%d] = %v, want ~%v", i, got[i], want[i])
				break
			}
		}
	}
	if !reflect.DeepEqual(gotElems, [][]int{{2, 1, 1}}) {
		t.Errorf("elements = %v, want [[2 1 1]]", gotElems)
	}
}

// TestConnectedElements ports TestConnectedElements: shared-vertex grouping into components.
func TestConnectedElements(t *testing.T) {
	cases := []struct {
		name     string
		elements [][]int
		want     [][]int
	}{
		{"single", [][]int{{1, 2}}, [][]int{{0}}},
		{"two_connected", [][]int{{1, 2}, {2, 3}}, [][]int{{0, 1}}},
		{"two_disconnected", [][]int{{1, 2}, {3, 4}}, [][]int{{0}, {1}}},
		{"multiple", [][]int{{1, 2}, {2, 3}, {3, 4}, {5, 6}}, [][]int{{0, 1, 2}, {3}}},
		{"triangles", [][]int{{1, 2, 3}, {3, 4, 5}, {6, 7, 8}}, [][]int{{0, 1}, {2}}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := ConnectedElements(c.elements)
			if !reflect.DeepEqual(got, c.want) {
				t.Errorf("ConnectedElements(%v) = %v, want %v", c.elements, got, c.want)
			}
		})
	}
}

// TestLineOrientationEqual ports TestLineOrientation: head-to-tail consistency.
func TestLineOrientationEqual(t *testing.T) {
	same := [][]int{{0, 1}, {1, 2}}
	if !LineOrientationEqual(0, 1, same) {
		t.Error("[[0 1] [1 2]] should be orientation-equal (head-to-tail)")
	}
	opposite := [][]int{{0, 1}, {2, 1}}
	if LineOrientationEqual(0, 1, opposite) {
		t.Error("[[0 1] [2 1]] should NOT be orientation-equal")
	}
	// Large-angle V shape: still head-to-tail, still equal.
	vshape := [][]int{{0, 1}, {1, 2}}
	if !LineOrientationEqual(0, 1, vshape) {
		t.Error("V-shape [[0 1] [1 2]] should be orientation-equal")
	}
}
