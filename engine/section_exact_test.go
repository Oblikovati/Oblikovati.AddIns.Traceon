// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"math"
	"testing"

	"oblikovati.org/api/wire"
	"oblikovati.org/traceon/core/geom2d"
)

// closedWire builds a closed section WirePolyline from (x, y) corners (z = 0), as the kernel
// returns it — flattened xyz triplets with Closed set.
func closedWire(corners ...geom2d.Point2) wire.WirePolyline {
	pts := make([]float64, 0, len(corners)*3)
	for _, c := range corners {
		pts = append(pts, c[0], c[1], 0)
	}
	return wire.WirePolyline{Points: pts, Closed: true}
}

// TestFoldWireSolidStraddle folds the section of a solid cylinder (a rectangle spanning −R…R) into
// one open meridian arc: top cap, outer wall, bottom cap, with both ends meeting the axis at r = 0.
func TestFoldWireSolidStraddle(t *testing.T) {
	w := closedWire(
		geom2d.Point2{1, 2}, geom2d.Point2{1, -2}, // right wall (x = R)
		geom2d.Point2{-1, -2}, geom2d.Point2{-1, 2}, // left wall (x = −R)
	)
	loops := foldWire(w)
	if len(loops) != 1 {
		t.Fatalf("solid section: want 1 arc, got %d", len(loops))
	}
	arc := loops[0]
	for _, p := range arc {
		if p[0] < -1e-9 {
			t.Errorf("folded point has negative r: %v", p)
		}
	}
	// The arc must reach the axis exactly twice (the two end-caps) and the outer wall (r = 1).
	if got := countNearR(arc, 0); got != 2 {
		t.Errorf("want 2 on-axis points (the caps), got %d in %v", got, arc)
	}
	if got := countNearR(arc, 1); got < 2 {
		t.Errorf("want ≥2 outer-wall points at r=1, got %d in %v", got, arc)
	}
}

// TestFoldWireTubeKeepsNearWire keeps the tube wall (a rectangle wholly at x > 0) as a closed loop
// and discards its mirror at x < 0, so the meridian carries the wall once, not twice.
func TestFoldWireTubeKeepsNearWire(t *testing.T) {
	near := closedWire(
		geom2d.Point2{0.6, 0.5}, geom2d.Point2{0.6, -0.5},
		geom2d.Point2{0.3, -0.5}, geom2d.Point2{0.3, 0.5},
	)
	far := closedWire(
		geom2d.Point2{-0.6, 0.5}, geom2d.Point2{-0.6, -0.5},
		geom2d.Point2{-0.3, -0.5}, geom2d.Point2{-0.3, 0.5},
	)
	loops := append(foldWire(near), foldWire(far)...)
	if len(loops) != 1 {
		t.Fatalf("tube section: want 1 loop (near wall, mirror discarded), got %d", len(loops))
	}
	loop := loops[0]
	if !samePoint(loop[0], loop[len(loop)-1]) {
		t.Errorf("tube wall loop must be closed (first==last), got %v … %v", loop[0], loop[len(loop)-1])
	}
	for _, p := range loop {
		if p[0] < 0.3-1e-9 || p[0] > 0.6+1e-9 {
			t.Errorf("tube wall r out of [0.3,0.6]: %v", p)
		}
	}
}

// TestFoldWireMirrorOnly discards a wire lying wholly on the far side of the axis.
func TestFoldWireMirrorOnly(t *testing.T) {
	w := closedWire(geom2d.Point2{-1, 1}, geom2d.Point2{-1, -1}, geom2d.Point2{-2, -1}, geom2d.Point2{-2, 1})
	if loops := foldWire(w); loops != nil {
		t.Errorf("far-side wire should fold to nothing, got %v", loops)
	}
}

// TestAxisCrossing checks the r = 0 crossing of a straddling segment is interpolated in z.
func TestAxisCrossing(t *testing.T) {
	p := axisCrossing(geom2d.Point2{2, 0}, geom2d.Point2{-2, 4})
	if math.Abs(p[0]) > 1e-12 || math.Abs(p[1]-2) > 1e-12 {
		t.Errorf("axisCrossing = %v, want (0, 2)", p)
	}
}

// TestExtractProfilePrefersSection proves extractProfile uses the EXACT brep section when one is
// available (and the body passes the axisymmetry guard), not the faceted envelope.
func TestExtractProfilePrefersSection(t *testing.T) {
	h := cylinderHost() // facets describe an axisymmetric r=1 cylinder → passes the guard
	h.sectionWires = []wire.WirePolyline{closedWire(
		geom2d.Point2{1, 1}, geom2d.Point2{1, -1},
		geom2d.Point2{-1, -1}, geom2d.Point2{-1, 1},
	)}
	prof, err := NewEngine(h).extractProfile(0)
	if err != nil {
		t.Fatalf("extractProfile: %v", err)
	}
	if !h.sawCall(wire.MethodBrepSectionWithPlane) {
		t.Error("expected extractProfile to attempt an exact section")
	}
	if !h.sawCall(wire.MethodBrepDelete) {
		t.Error("expected the transient section body to be freed")
	}
	// The exact meridian is a single arc with on-axis caps; the facet envelope (60-band) would
	// instead yield a long polyline with no exact r=0 cap points.
	if len(prof.loops) != 1 {
		t.Fatalf("want 1 exact meridian loop, got %d", len(prof.loops))
	}
	if got := countNearR(prof.loops[0], 0); got != 2 {
		t.Errorf("exact meridian should have 2 on-axis cap points, got %d", got)
	}
}

// TestExtractProfileFallsBackToFacets confirms that with no section wires the extractor still
// produces a meridian from the facet envelope (so a body the kernel cannot section still solves).
func TestExtractProfileFallsBackToFacets(t *testing.T) {
	h := cylinderHost() // sectionWires nil → section returns no wires
	prof, err := NewEngine(h).extractProfile(0)
	if err != nil {
		t.Fatalf("extractProfile: %v", err)
	}
	if len(prof.loops) == 0 || len(prof.loops[0]) < 2 {
		t.Fatalf("facet fallback produced no meridian: %v", prof.loops)
	}
}

// countNearR counts points whose r is within coincidentTol of target.
func countNearR(loop []geom2d.Point2, target float64) int {
	n := 0
	for _, p := range loop {
		if math.Abs(p[0]-target) < coincidentTol {
			n++
		}
	}
	return n
}
