// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"encoding/json"
	"testing"
	"time"

	"oblikovati.org/api/wire"
)

// TestParseLens checks the panel string maps to the lens mode, defaulting to host.
func TestParseLens(t *testing.T) {
	cases := map[string]paramLens{
		"einzel": lensEinzel, "cylinder": lensCylinder, "host": lensHost, "": lensHost, "bogus": lensHost,
	}
	for in, want := range cases {
		if got := parseLens(in); got != want {
			t.Errorf("parseLens(%q) = %q, want %q", in, got, want)
		}
	}
}

// TestEinzelElectrodes checks the einzel template builds three electrodes (grounded, biased,
// grounded), each meshed into BEM elements with a render profile.
func TestEinzelElectrodes(t *testing.T) {
	p := defaultParams()
	p.lens = lensEinzel
	p.voltage = 1500

	els, err := buildParametricLens(p)
	if err != nil {
		t.Fatalf("buildParametricLens: %v", err)
	}
	if len(els) != 3 {
		t.Fatalf("electrodes = %d, want 3", len(els))
	}
	wantV := []float64{0, 1500, 0}
	for i, el := range els {
		if el.voltage != wantV[i] {
			t.Errorf("electrode %d voltage = %g, want %g", i, el.voltage, wantV[i])
		}
		if len(el.lines) == 0 {
			t.Errorf("electrode %d has no meshed BEM elements", i)
		}
		if el.prof == nil || len(el.prof.loops) == 0 || len(el.prof.loops[0]) == 0 {
			t.Errorf("electrode %d has no render profile", i)
		}
	}
}

// TestCylinderElectrodes checks the two-cylinder template builds two electrodes at 0 V and the
// bias voltage.
func TestCylinderElectrodes(t *testing.T) {
	p := defaultParams()
	p.lens = lensCylinder
	p.voltage = 800

	els, err := buildParametricLens(p)
	if err != nil {
		t.Fatalf("buildParametricLens: %v", err)
	}
	if len(els) != 2 {
		t.Fatalf("electrodes = %d, want 2", len(els))
	}
	if els[0].voltage != 0 || els[1].voltage != 800 {
		t.Errorf("voltages = [%g %g], want [0 800]", els[0].voltage, els[1].voltage)
	}
	for i, el := range els {
		if len(el.lines) == 0 {
			t.Errorf("cylinder %d has no meshed BEM elements", i)
		}
	}
}

// TestBuildParametricLensUnknown checks an unrecognised lens is an error.
func TestBuildParametricLensUnknown(t *testing.T) {
	p := defaultParams()
	p.lens = paramLens("nope")
	if _, err := buildParametricLens(p); err == nil {
		t.Error("expected error for unknown lens, got nil")
	}
}

// TestParametricEinzelStudy runs a whole study from a parametric einzel lens with no host
// geometry: it must solve, trace, and push the overlay without ever listing host bodies.
func TestParametricEinzelStudy(t *testing.T) {
	h := &fakeHost{}
	e := NewEngine(h)
	e.params.lens = lensEinzel

	res, err := e.RunStudy(0)
	if err != nil {
		t.Fatalf("parametric study: %v", err)
	}
	if res.ElectrodeCount != 3 {
		t.Errorf("electrodes = %d, want 3", res.ElectrodeCount)
	}
	if res.ElementCount == 0 {
		t.Error("no BEM elements assembled from the parametric lens")
	}
	if res.RayCount == 0 {
		t.Error("no rays traced through the parametric lens")
	}
	// The study must not touch host geometry in parametric mode.
	if h.sawCall(wire.MethodBodyList) || h.sawCall(wire.MethodBodyCalculateFacets) {
		t.Errorf("parametric study listed/sectioned host bodies; calls = %v", h.calls)
	}
	// The overlay (electrodes + potential + rays) must be pushed.
	if len(h.lastGraph.Nodes) == 0 {
		t.Error("no client-graphics overlay pushed")
	}
}

// TestParametricStudyNoHostCalls confirms the cylinder template likewise needs no host model.
func TestParametricStudyNoHostCalls(t *testing.T) {
	h := &fakeHost{}
	e := NewEngine(h)
	e.params.lens = lensCylinder

	if _, err := e.RunStudy(0); err != nil {
		t.Fatalf("cylinder study: %v", err)
	}
	if h.sawCall(wire.MethodBodyList) {
		t.Errorf("cylinder study listed host bodies; calls = %v", h.calls)
	}
}

// TestLensForCommand checks the parametric-lens commands map to their lens, and other commands
// are not treated as parametric.
func TestLensForCommand(t *testing.T) {
	if l, ok := lensForCommand(EinzelLensCommandID); !ok || l != lensEinzel {
		t.Errorf("einzel command → (%q, %v), want (einzel, true)", l, ok)
	}
	if l, ok := lensForCommand(CylinderLensCommandID); !ok || l != lensCylinder {
		t.Errorf("cylinder command → (%q, %v), want (cylinder, true)", l, ok)
	}
	if _, ok := lensForCommand(RunStudyCommandID); ok {
		t.Error("RunStudy must not be a parametric-lens command")
	}
}

// TestNotifyEinzelCommandRunsParametric checks firing the einzel command selects the einzel lens
// and runs a host-free study (the overlay is pushed without listing host bodies).
func TestNotifyEinzelCommandRunsParametric(t *testing.T) {
	h := &fakeHost{}
	e := NewEngine(h)

	ev, _ := json.Marshal(map[string]string{"type": wire.EventCommandStarted, "command": EinzelLensCommandID})
	e.Notify(ev)
	// launchStudy runs on a goroutine; wait for it to settle. The deadline is generous because the
	// race detector slows the study past a 2s cap (flaked in CI); the loop exits as soon as the
	// overlay lands, so passing runs never wait it out.
	deadline := time.Now().Add(30 * time.Second)
	for time.Now().Before(deadline) {
		e.mu.Lock()
		done := !e.running && e.params.lens == lensEinzel
		e.mu.Unlock()
		if done && len(h.lastGraph.Nodes) > 0 {
			break
		}
		time.Sleep(10 * time.Millisecond)
	}
	if e.params.lens != lensEinzel {
		t.Errorf("lens = %q, want einzel after the command", e.params.lens)
	}
	if len(h.lastGraph.Nodes) == 0 {
		t.Error("einzel command did not push a study overlay")
	}
	if h.sawCall(wire.MethodBodyList) {
		t.Errorf("parametric command listed host bodies; calls = %v", h.calls)
	}
}

// TestApplyPanelEditLens checks the lens-definition controls write through to study params.
func TestApplyPanelEditLens(t *testing.T) {
	e := NewEngine(&fakeHost{})
	e.applyPanelEdit("lens", "einzel")
	e.applyPanelEdit("lens_radius", "0.4")
	e.applyPanelEdit("lens_thickness", "0.6")
	e.applyPanelEdit("lens_spacing", "0.7")

	if e.params.lens != lensEinzel {
		t.Errorf("lens = %q, want einzel", e.params.lens)
	}
	if e.params.lensRadius != 0.4 || e.params.lensThickness != 0.6 || e.params.lensSpacing != 0.7 {
		t.Errorf("dims = (%g, %g, %g), want (0.4, 0.6, 0.7)", e.params.lensRadius, e.params.lensThickness, e.params.lensSpacing)
	}
	// A non-positive dimension is rejected (keeps the prior value).
	e.applyPanelEdit("lens_radius", "0")
	if e.params.lensRadius != 0.4 {
		t.Errorf("lensRadius after invalid edit = %g, want 0.4", e.params.lensRadius)
	}
}
