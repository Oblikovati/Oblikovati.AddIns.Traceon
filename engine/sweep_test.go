// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"math"
	"testing"

	"oblikovati.org/api/wire"
)

func TestSweepConfigValueAt(t *testing.T) {
	cfg := sweepConfig{start: 0.3, stop: 0.6, steps: 4}
	want := []float64{0.3, 0.4, 0.5, 0.6}
	for i, w := range want {
		if got := cfg.valueAt(i); math.Abs(got-w) > 1e-12 {
			t.Errorf("valueAt(%d) = %g, want %g", i, got, w)
		}
	}
}

func TestSweepConfigExpression(t *testing.T) {
	withUnit := sweepConfig{unit: "cm"}
	if got := withUnit.expression(0.45); got != "0.45 cm" {
		t.Errorf("expression with unit = %q, want %q", got, "0.45 cm")
	}
	bare := sweepConfig{}
	if got := bare.expression(0.45); got != "0.45" {
		t.Errorf("bare expression = %q, want %q", got, "0.45")
	}
}

func TestFocusSpanIgnoresGaps(t *testing.T) {
	points := []sweepPoint{
		{value: 0, focusZ: math.NaN()},
		{value: 1, focusZ: 8.0},
		{value: 2, focusZ: 4.0},
		{value: 3, focusZ: math.Inf(1)},
	}
	lo, hi, found := focusSpan(points)
	if !found || lo != 4.0 || hi != 8.0 {
		t.Errorf("focusSpan = (%g, %g, %v), want (4, 8, true)", lo, hi, found)
	}
}

// TestSweepNodesBreaksAtGaps proves the plot frame + labels are always drawn, the focus curve is
// split into separate strips at NaN samples, and curve points stay inside the plot box.
func TestSweepNodesBreaksAtGaps(t *testing.T) {
	cfg := sweepConfig{param: "bore", unit: "cm", start: 0.3, stop: 0.7, steps: 5}
	points := []sweepPoint{
		{value: 0.3, focusZ: 8},
		{value: 0.4, focusZ: 6},
		{value: 0.5, focusZ: math.NaN()}, // gap splits the curve
		{value: 0.6, focusZ: 5},
		{value: 0.7, focusZ: 4},
	}
	nodes := sweepNodes(cfg, points)
	curves, markers := 0, 0
	box := sweepBox()
	for _, n := range nodes {
		switch {
		case n.Id == "traceon.sweep.markers":
			markers++
		case len(n.Id) > 19 && n.Id[:20] == "traceon.sweep.curve.":
			curves++
			for _, p := range n.Primitives {
				// Coordinates are remapped to world axes: parameter → world x, focus → world y.
				for i := 0; i+2 < len(p.Coordinates); i += 3 {
					x, y := p.Coordinates[i], p.Coordinates[i+1]
					if x < box.x0-1e-9 || x > box.x1+1e-9 || y < box.z0-1e-9 || y > box.z1+1e-9 {
						t.Errorf("curve point (%g,%g) outside plot box", x, y)
					}
				}
			}
		}
	}
	if curves != 2 {
		t.Errorf("want 2 curve strips (split at the NaN), got %d", curves)
	}
	if markers != 1 {
		t.Errorf("want 1 markers node, got %d", markers)
	}
	if !hasNode(nodes, "traceon.sweep.frame") || !hasNode(nodes, "traceon.sweep.labels") {
		t.Error("plot must always include a frame and labels")
	}
}

// sweepFakeHost is a cylinder fake wired with a host parameter "bore" = 0.3 cm (displayed 3 mm).
func sweepFakeHost() *fakeHost {
	h := cylinderHost()
	h.paramExpr = "0.3 cm"
	h.paramUnit = "mm"
	h.paramValue = "3 mm"
	h.paramList = []wire.ParameterInfo{{Name: "bore", Kind: "user", Expression: "0.3 cm", Value: "3 mm"}}
	return h
}

// TestRunSweepRestoresParameter runs a full sweep and asserts it sampled every step, pushed a
// focus-vs-parameter plot, and restored the swept parameter to its original expression.
func TestRunSweepRestoresParameter(t *testing.T) {
	h := sweepFakeHost()
	e := NewEngine(h)
	e.params.sweepParam = "bore"
	e.params.sweepStart = 0.3
	e.params.sweepStop = 0.6
	e.params.sweepSteps = 4

	msg, err := e.runSweep()
	if err != nil {
		t.Fatalf("runSweep: %v", err)
	}
	// One set per step plus the final restore.
	if len(h.paramSets) != 5 {
		t.Fatalf("want 5 parameters.set calls (4 steps + restore), got %d", len(h.paramSets))
	}
	if last := h.paramSets[len(h.paramSets)-1]; last.Expression != "0.3 cm" {
		t.Errorf("parameter not restored: last set = %q, want original %q", last.Expression, "0.3 cm")
	}
	if h.sweepGraph.ClientId != sweepClientID {
		t.Errorf("expected a plot push to %q, got %q", sweepClientID, h.sweepGraph.ClientId)
	}
	if msg == "" {
		t.Error("expected a sweep summary message")
	}
}

// TestRunSweepNeedsHostGeometry refuses a sweep in parametric-lens mode (no host parameter to vary).
func TestRunSweepNeedsHostGeometry(t *testing.T) {
	h := sweepFakeHost()
	e := NewEngine(h)
	e.params.lens = lensEinzel
	e.params.sweepParam = "bore"
	e.params.sweepStop = 0.6
	msg, err := e.runSweep()
	if err != nil {
		t.Fatalf("runSweep: %v", err)
	}
	if msg == "" || h.sawCall(wire.MethodParametersSet) {
		t.Errorf("parametric-lens sweep should be refused without touching parameters; msg=%q", msg)
	}
}

// TestRunSweepAutoPicksParameter runs with no panel name or range: it should auto-select the
// model's user parameter and default the range to ±50 % of its current value (3 mm → [1.5, 4.5]).
func TestRunSweepAutoPicksParameter(t *testing.T) {
	h := sweepFakeHost()
	e := NewEngine(h) // sweepParam "" , range unset, steps default 9
	msg, err := e.runSweep()
	if err != nil {
		t.Fatalf("runSweep: %v", err)
	}
	if !h.sawCall(wire.MethodParametersList) {
		t.Error("expected the sweep to auto-pick a parameter via parameters.list")
	}
	first := h.paramSets[0].Expression
	if first != "1.5 mm" {
		t.Errorf("default range should start at 50%% of 3 mm = 1.5 mm; first set = %q", first)
	}
	if last := h.paramSets[len(h.paramSets)-1].Expression; last != "0.3 cm" {
		t.Errorf("parameter not restored: last set = %q", last)
	}
	if msg == "" {
		t.Error("expected a summary")
	}
}

// TestRunSweepRejectsNoParameter reports a friendly message when the model has no user parameter
// to sweep, without erroring or mutating the model.
func TestRunSweepRejectsNoParameter(t *testing.T) {
	h := cylinderHost() // no paramList
	e := NewEngine(h)
	msg, err := e.runSweep()
	if err != nil {
		t.Fatalf("runSweep: %v", err)
	}
	if msg == "" || h.sawCall(wire.MethodParametersSet) {
		t.Errorf("a model with no parameter should be rejected without setting parameters; msg=%q", msg)
	}
}

func hasNode(nodes []wire.GraphicsNode, id string) bool {
	for _, n := range nodes {
		if n.Id == id {
			return true
		}
	}
	return false
}
