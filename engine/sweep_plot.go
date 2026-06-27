// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"fmt"
	"math"

	"oblikovati.org/api/types"
	"oblikovati.org/api/wire"
)

// The focus-vs-parameter plot is drawn in the xz-plane at x < 0 — the half-space the axisymmetric
// model never occupies (geometry lives at x = r ≥ 0) — so it sits clearly beside the study without
// overlapping it. Dimensions are in cm (the host DB unit).
const (
	sweepPlotWidth  = 4.0 // box width (parameter axis)
	sweepPlotHeight = 4.0 // box height (focus axis)
	sweepPlotMargin = 1.0 // gap from the optical axis (x = 0)
)

// plotBox is the world-space rectangle (xz-plane, y = 0) the sweep curve is drawn into.
type plotBox struct {
	x0, x1, z0, z1 float64
}

// sweepBox places the plot box just left of the optical axis.
func sweepBox() plotBox {
	return plotBox{x0: -(sweepPlotMargin + sweepPlotWidth), x1: -sweepPlotMargin, z0: 0, z1: sweepPlotHeight}
}

// mapX maps a parameter value in [lo, hi] to the box's horizontal span (degenerate range → centre).
func (b plotBox) mapX(value, lo, hi float64) float64 {
	if hi == lo {
		return (b.x0 + b.x1) / 2
	}
	return b.x0 + (b.x1-b.x0)*(value-lo)/(hi-lo)
}

// mapZ maps a focus value in [lo, hi] to the box's vertical span (degenerate range → centre).
func (b plotBox) mapZ(value, lo, hi float64) float64 {
	if hi == lo {
		return (b.z0 + b.z1) / 2
	}
	return b.z0 + (b.z1-b.z0)*(value-lo)/(hi-lo)
}

// pushSweepPlot renders the focus-vs-parameter plot into its own client-graphics group (leaving the
// last study's field + trajectories untouched).
func (e *Engine) pushSweepPlot(cfg sweepConfig, points []sweepPoint) error {
	_, err := e.api.Graphics().Set(wire.SetClientGraphicsArgs{
		ClientId: sweepClientID,
		Lane:     string(types.GraphicsLanePersistent),
		Nodes:    sweepNodes(cfg, points),
	})
	return err
}

// sweepNodes assembles the plot overlay: an axis frame with labels, plus the focus curve and its
// sample markers (omitted when no focus formed anywhere — only the labelled frame is shown).
func sweepNodes(cfg sweepConfig, points []sweepPoint) []wire.GraphicsNode {
	box := sweepBox()
	fLo, fHi, found := focusSpan(points)
	nodes := []wire.GraphicsNode{frameNode(box), labelsNode(cfg, box, fLo, fHi, found)}
	if !found {
		return nodes
	}
	nodes = append(nodes, curveNodes(cfg, points, box, fLo, fHi)...)
	nodes = append(nodes, markersNode(cfg, points, box, fLo, fHi))
	return nodes
}

// frameNode draws the plot box border (the two axes plus the closing top/right edges). Like the
// study overlay, the plot is drawn in engine frame (parameter on x, focus on the optical axis) and
// remapped to world axes so it stands upright in the lens plane rather than lying flat.
func frameNode(b plotBox) wire.GraphicsNode {
	var coords []float64
	for _, e := range [][4]float64{
		{b.x0, b.z0, b.x1, b.z0}, // bottom (parameter axis)
		{b.x0, b.z0, b.x0, b.z1}, // left (focus axis)
		{b.x0, b.z1, b.x1, b.z1}, // top
		{b.x1, b.z0, b.x1, b.z1}, // right
	} {
		coords = appendWorld(coords, e[0], 0, e[1])
		coords = appendWorld(coords, e[2], 0, e[3])
	}
	return wire.GraphicsNode{Id: "traceon.sweep.frame", Primitives: []wire.GraphicsPrimitive{{
		Kind: string(types.GraphicsLines), Coordinates: coords,
		Indices: []int{0, 1, 2, 3, 4, 5, 6, 7},
		Color:   []float32{0.6, 0.6, 0.65, 1}, OnTop: true,
	}}}
}

// curveNodes draws the focus(parameter) curve as line strips, breaking the strip wherever the beam
// did not focus (NaN) so a gap reads as "no focus here" rather than a false straight segment.
func curveNodes(cfg sweepConfig, points []sweepPoint, b plotBox, fLo, fHi float64) []wire.GraphicsNode {
	var nodes []wire.GraphicsNode
	var coords []float64
	flush := func() {
		if len(coords) >= 6 {
			nodes = append(nodes, wire.GraphicsNode{
				Id: "traceon.sweep.curve." + itoa(len(nodes)),
				Primitives: []wire.GraphicsPrimitive{{
					Kind: string(types.GraphicsLineStrip), Coordinates: coords,
					Color: []float32{0.2, 0.9, 1.0, 1}, LineWeight: 2, OnTop: true,
				}},
			})
		}
		coords = nil
	}
	for _, p := range points {
		if isGap(p.focusZ) {
			flush()
			continue
		}
		coords = appendWorld(coords, b.mapX(p.value, cfg.start, cfg.stop), 0, b.mapZ(p.focusZ, fLo, fHi))
	}
	flush()
	return nodes
}

// markersNode draws a point glyph at every focused sample.
func markersNode(cfg sweepConfig, points []sweepPoint, b plotBox, fLo, fHi float64) wire.GraphicsNode {
	var coords []float64
	for _, p := range points {
		if isGap(p.focusZ) {
			continue
		}
		coords = appendWorld(coords, b.mapX(p.value, cfg.start, cfg.stop), 0, b.mapZ(p.focusZ, fLo, fHi))
	}
	return wire.GraphicsNode{Id: "traceon.sweep.markers", Primitives: []wire.GraphicsPrimitive{{
		Kind: string(types.GraphicsPoints), Coordinates: coords,
		Color: []float32{1, 1, 0.3, 1}, PointSize: 7, OnTop: true,
	}}}
}

// labelsNode anchors text at the plot corners: the parameter name + range along the bottom and the
// focus range up the left, so the curve's axes are readable in the viewport.
func labelsNode(cfg sweepConfig, b plotBox, fLo, fHi float64, found bool) wire.GraphicsNode {
	prims := []wire.GraphicsPrimitive{
		textPrim(cfg.param+" ("+cfg.unit+") →", b.x0, b.z0-0.4),
		textPrim(fmt.Sprintf("%g", cfg.start), b.x0, b.z0-0.2),
		textPrim(fmt.Sprintf("%g", cfg.stop), b.x1, b.z0-0.2),
		textPrim("focus z (cm) ↑", b.x0-0.6, b.z1+0.2),
	}
	if found {
		prims = append(prims,
			textPrim(fmt.Sprintf("%.2f", fHi), b.x0-0.6, b.z1),
			textPrim(fmt.Sprintf("%.2f", fLo), b.x0-0.6, b.z0),
		)
	}
	return wire.GraphicsNode{Id: "traceon.sweep.labels", Primitives: prims}
}

// textPrim builds a text label anchored at an engine-frame (x, focus) point, remapped to world axes.
func textPrim(s string, x, z float64) wire.GraphicsPrimitive {
	wx, wy, wz := worldFromEngine(x, 0, z)
	return wire.GraphicsPrimitive{
		Kind: string(types.GraphicsText), Text: s, Anchor: []float64{wx, wy, wz},
		Color: []float32{0.85, 0.85, 0.9, 1}, FontSize: 14, OnTop: true,
	}
}

// isGap reports whether a focus value is a sweep gap (no axis crossing → NaN, or non-finite).
func isGap(focusZ float64) bool {
	return math.IsNaN(focusZ) || math.IsInf(focusZ, 0)
}
