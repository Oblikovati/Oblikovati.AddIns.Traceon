// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"fmt"

	"oblikovati.org/api/wire"
	"oblikovati.org/traceon/core/geom2d"
)

// axisFoldEps is the half-width (cm) of the on-axis band: a section point with |x| below this is
// treated as lying on the optical axis (r = 0). It absorbs the kernel's sampling noise so an
// end-cap that reaches the axis closes cleanly instead of leaving a sub-micron sliver.
const axisFoldEps = 1e-7

// sectionMeridian extracts a body's EXACT axisymmetric (r, z) meridian by intersecting it with the
// world XY plane — the plane that contains the optical (Y) axis — and folding the resulting section
// wires onto the r ≥ 0 half-plane (r = |x|, z = y). Unlike the faceted axial-band envelope this
// captures the true profile: curved walls (sampled to the kernel's chord tolerance), overhangs, and
// fine bore detail a 60-bin envelope smears. Returns the meridian as ordered (r, z) loops, or an
// empty slice when the section yields nothing usable so the caller falls back to the facet envelope.
func (e *Engine) sectionMeridian(bodyIndex int) ([][]geom2d.Point2, error) {
	idx := bodyIndex
	res, err := e.api.TransientBRep().CreateIntersectionWithPlane(
		wire.BrepBodyRef{BodyIndex: &idx}, []float64{0, 0, 0}, []float64{0, 0, 1})
	if err != nil {
		return nil, fmt.Errorf("section body %d with plane: %w", bodyIndex, err)
	}
	// The section lands on a transient body too; we only need the sampled wires, so free it.
	if res.Handle != 0 {
		_ = e.api.TransientBRep().Delete(res.Handle)
	}
	var loops [][]geom2d.Point2
	for _, w := range res.Wires {
		loops = append(loops, foldWire(w)...)
	}
	return loops, nil
}

// foldWire folds one z = 0 section wire onto the r ≥ 0 half-plane. A wire lying wholly on the near
// side (x ≥ 0) is a closed meridian loop (a tube wall band) and is kept closed; one wholly on the
// far side (x < 0) is the mirror image and is discarded; one straddling the axis (a solid whose
// section spans −R…R) is split at its axis crossings into open arcs whose end-caps meet r = 0.
func foldWire(w wire.WirePolyline) [][]geom2d.Point2 {
	seq := wireXY(w)
	if len(seq) < 2 {
		return nil
	}
	lo, hi := xRange(seq)
	if lo >= -axisFoldEps {
		return [][]geom2d.Point2{closeLoop(clampAxis(seq))} // wholly near side: a closed wall loop
	}
	if hi < axisFoldEps {
		return nil // wholly far side: the mirror half, already represented by the near side
	}
	return insideArcs(rotateToFarSide(seq))
}

// wireXY flattens a sampled wire to its (x, y) points in the section (z = 0) plane, dropping a
// trailing point coincident with the first so a closed wire is a clean cycle of distinct vertices.
func wireXY(w wire.WirePolyline) []geom2d.Point2 {
	n := len(w.Points) / 3
	pts := make([]geom2d.Point2, 0, n)
	for i := 0; i < n; i++ {
		pts = append(pts, geom2d.Point2{w.Points[i*3], w.Points[i*3+1]})
	}
	if len(pts) >= 2 && samePoint(pts[0], pts[len(pts)-1]) {
		pts = pts[:len(pts)-1]
	}
	return pts
}

// xRange returns the minimum and maximum x of a point set.
func xRange(pts []geom2d.Point2) (lo, hi float64) {
	lo, hi = pts[0][0], pts[0][0]
	for _, p := range pts {
		if p[0] < lo {
			lo = p[0]
		}
		if p[0] > hi {
			hi = p[0]
		}
	}
	return lo, hi
}

// rotateToFarSide cyclically shifts a closed wire's vertices so it starts at a far-side (x < 0)
// point. Then a near-side run can never straddle the seam between index 0 and the wrap edge, so
// insideArcs extracts each run as one contiguous arc instead of a split head and tail.
func rotateToFarSide(seq []geom2d.Point2) []geom2d.Point2 {
	start := 0
	for i, p := range seq {
		if p[0] < -axisFoldEps {
			start = i
			break
		}
	}
	out := make([]geom2d.Point2, len(seq))
	for i := range seq {
		out[i] = seq[(start+i)%len(seq)]
	}
	return out
}

// insideArcs walks a closed wire (already rotated to start on the far side) and returns each
// maximal run of near-side (x ≥ 0) vertices as one open (r, z) arc, inserting the r = 0 crossing
// point where the wire enters or leaves the axis so an end-cap that reaches the axis is closed.
func insideArcs(seq []geom2d.Point2) [][]geom2d.Point2 {
	n := len(seq)
	var arcs [][]geom2d.Point2
	var cur []geom2d.Point2
	for i := 0; i < n; i++ {
		a, b := seq[i], seq[(i+1)%n]
		inA, inB := a[0] >= -axisFoldEps, b[0] >= -axisFoldEps
		switch {
		case inA && inB:
			if len(cur) == 0 {
				cur = append(cur, clampAxisPoint(a))
			}
			cur = append(cur, clampAxisPoint(b))
		case inA && !inB: // leaving the near side: close the arc on the axis
			if len(cur) == 0 {
				cur = append(cur, clampAxisPoint(a))
			}
			cur = append(cur, axisCrossing(a, b))
			arcs = appendArc(arcs, cur)
			cur = nil
		case !inA && inB: // entering the near side: open a fresh arc on the axis
			cur = append(cur[:0], axisCrossing(a, b))
			cur = append(cur, clampAxisPoint(b))
		}
	}
	return appendArc(arcs, cur)
}

// axisCrossing returns the (r = 0, z) point where segment a→b crosses x = 0, by linear
// interpolation. Callers guarantee a and b straddle the axis, so a[0] − b[0] is non-zero.
func axisCrossing(a, b geom2d.Point2) geom2d.Point2 {
	t := a[0] / (a[0] - b[0])
	return geom2d.Point2{0, a[1] + (b[1]-a[1])*t}
}

// appendArc keeps an arc only if it has at least two points (a single point makes no element).
func appendArc(arcs [][]geom2d.Point2, arc []geom2d.Point2) [][]geom2d.Point2 {
	if len(arc) >= 2 {
		arcs = append(arcs, arc)
	}
	return arcs
}

// clampAxis maps a whole loop onto the r ≥ 0 half-plane (r = |x|), snapping near-axis noise to 0.
func clampAxis(seq []geom2d.Point2) []geom2d.Point2 {
	out := make([]geom2d.Point2, len(seq))
	for i, p := range seq {
		out[i] = clampAxisPoint(p)
	}
	return out
}

// clampAxisPoint folds one point onto r = |x| (z = y), snapping |x| < axisFoldEps to exactly 0.
func clampAxisPoint(p geom2d.Point2) geom2d.Point2 {
	r := p[0]
	if r < 0 {
		r = -r
	}
	if r < axisFoldEps {
		r = 0
	}
	return geom2d.Point2{r, p[1]}
}

// closeLoop appends the first point to a loop so its closing segment (last → first) is included
// when lineElements walks consecutive pairs; a no-op if already closed or too short.
func closeLoop(loop []geom2d.Point2) []geom2d.Point2 {
	if len(loop) < 3 || samePoint(loop[0], loop[len(loop)-1]) {
		return loop
	}
	return append(loop, loop[0])
}

// samePoint reports whether two (r, z) points coincide within coincidentTol.
func samePoint(a, b geom2d.Point2) bool {
	return abs(a[0]-b[0]) < coincidentTol && abs(a[1]-b[1]) < coincidentTol
}

// abs is the float absolute value (kept local so this file has no math import beyond what it needs).
func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
