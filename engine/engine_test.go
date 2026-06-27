// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"encoding/json"
	"errors"
	"math"
	"strings"
	"testing"

	"oblikovati.org/api/types"
	"oblikovati.org/api/wire"

	"oblikovati.org/traceon/core/geom2d"
)

// fakeHost is a named fake HostCaller (no live host): it answers the wire methods a study
// issues with canned JSON, records the methods it saw, and captures the client-graphics
// payload so a test can assert the section→solve→trace→render pipeline ran end to end.
type fakeHost struct {
	calls        []string
	failOn       string // method to fail, for error-path tests ("" = none)
	facets       wire.FacetSetResult
	voltagesJSON string // traceon/voltages attribute payload ("" = unset)
	lastGraph    wire.SetClientGraphicsArgs
}

func (h *fakeHost) Call(method string, req []byte) ([]byte, error) {
	h.calls = append(h.calls, method)
	if method == h.failOn {
		return nil, errors.New("forced failure")
	}
	switch method {
	case wire.MethodBodyList:
		return json.Marshal(wire.BodyListResult{Bodies: []wire.BodyInfo{{Index: 0, Name: "Solid1", Solid: true, Key: "k0"}}})
	case wire.MethodBodyCalculateFacets:
		return json.Marshal(h.facets)
	case wire.MethodDocumentsList:
		return json.Marshal(wire.ListDocumentsResult{Documents: []wire.DocumentInfo{{ID: 1, Active: true}}})
	case wire.MethodAttributesGet:
		if h.voltagesJSON == "" {
			return json.Marshal(wire.AttributeResult{Found: false})
		}
		return json.Marshal(wire.AttributeResult{Found: true, Attribute: wire.AttributeInfo{
			Set: attrSet, Name: attrVoltages, Value: types.StringVariant(h.voltagesJSON)}})
	case wire.MethodClientGraphicsSet:
		_ = json.Unmarshal(req, &h.lastGraph)
		return []byte("{}"), nil
	default:
		return []byte("{}"), nil
	}
}

func (h *fakeHost) sawCall(method string) bool {
	for _, c := range h.calls {
		if c == method {
			return true
		}
	}
	return false
}

// cylinderHost is a fake whose body facets describe a cylinder of radius 1 about the Y axis
// spanning y∈[-1, 1]; the meridian extractor turns it into a vertical r=1 electrode profile.
// Vertices are given as (x, y, z) at a few axial levels and angles (r = √(x²+z²) = 1).
func cylinderHost() *fakeHost {
	var coords []float64
	for _, y := range []float64{-1, -0.5, 0, 0.5, 1} {
		for _, ang := range []float64{0, 1.57, 3.14, 4.71} {
			coords = append(coords, math.Cos(ang), y, math.Sin(ang)) // radius 1
		}
	}
	return &fakeHost{facets: wire.FacetSetResult{VertexCount: len(coords) / 3, VertexCoordinates: coords}}
}

func TestRunStudyDrivesPipeline(t *testing.T) {
	h := cylinderHost()
	res, err := NewEngine(h).RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	if !h.sawCall(wire.MethodBodyCalculateFacets) {
		t.Error("expected the study to section the body via calculateFacets")
	}
	if !h.sawCall(wire.MethodClientGraphicsSet) {
		t.Error("expected the study to push client graphics")
	}
	if res.ElementCount < 1 {
		t.Errorf("ElementCount = %d, want >=1", res.ElementCount)
	}
	if res.RayCount == 0 {
		t.Error("expected at least one traced ray")
	}
	if res.GraphicsClientID != graphicsClientID {
		t.Errorf("GraphicsClientID = %q, want %q", res.GraphicsClientID, graphicsClientID)
	}
}

// TestRunStudyRendersExpectedNodes checks the pushed overlay contains the potential map, the
// electrode profile, and one node per traced ray.
func TestRunStudyRendersExpectedNodes(t *testing.T) {
	h := cylinderHost()
	res, err := NewEngine(h).RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	ids := map[string]bool{}
	for _, n := range h.lastGraph.Nodes {
		ids[n.Id] = true
	}
	if !ids["traceon.potential"] {
		t.Error("missing potential heatmap node")
	}
	if !ids["traceon.electrode"] {
		t.Error("missing electrode profile node")
	}
	rayNodes := 0
	for id := range ids {
		if strings.HasPrefix(id, "traceon.ray.") {
			rayNodes++
		}
	}
	if rayNodes != res.RayCount {
		t.Errorf("ray nodes = %d, want %d (RayCount)", rayNodes, res.RayCount)
	}
}

// TestRunStudySectionError surfaces a host section failure rather than rendering an empty study.
func TestRunStudySectionError(t *testing.T) {
	h := cylinderHost()
	h.failOn = wire.MethodBodyCalculateFacets
	if _, err := NewEngine(h).RunStudy(0); err == nil {
		t.Error("expected RunStudy to fail when sectioning fails")
	}
}

// TestPerElectrodeVoltages checks the traceon/voltages document attribute is parsed into a
// per-body voltage map (body index → volts).
func TestPerElectrodeVoltages(t *testing.T) {
	h := &fakeHost{voltagesJSON: `{"0":0,"1":5000,"2":-250.5}`}
	v := NewEngine(h).electrodeVoltages()
	want := map[int]float64{0: 0, 1: 5000, 2: -250.5}
	if len(v) != len(want) {
		t.Fatalf("got %d voltages, want %d", len(v), len(want))
	}
	for k, w := range want {
		if v[k] != w {
			t.Errorf("voltage[%d] = %v, want %v", k, v[k], w)
		}
	}
	// Unset attribute → empty map (falls back to the panel default).
	if got := NewEngine(&fakeHost{}).electrodeVoltages(); len(got) != 0 {
		t.Errorf("unset voltages → %v, want empty", got)
	}
}

// TestCentralElectrode checks the einzel default picks the axially-central electrode.
func TestCentralElectrode(t *testing.T) {
	band := func(lo, hi float64) *profile {
		return &profile{loops: [][]geom2d.Point2{{{1, lo}, {1, hi}}}}
	}
	profs := []*profile{band(-3, -1), band(-0.5, 0.5), band(1, 3)}
	if got := centralElectrode(profs); got != 1 {
		t.Errorf("central electrode = %d, want 1 (the middle one)", got)
	}
}

// TestStudyReportsElectrodeCount checks the result reports one electrode for the single body.
func TestStudyReportsElectrodeCount(t *testing.T) {
	res, err := NewEngine(cylinderHost()).RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	if res.ElectrodeCount != 1 {
		t.Errorf("ElectrodeCount = %d, want 1", res.ElectrodeCount)
	}
}

// TestSetupRegistersUI checks Setup registers the command and shows the panel.
func TestSetupRegistersUI(t *testing.T) {
	h := &fakeHost{}
	if err := NewEngine(h).Setup(); err != nil {
		t.Fatalf("Setup: %v", err)
	}
	if !h.sawCall(wire.MethodCommandsCreate) {
		t.Error("Setup did not register the study command")
	}
	if !h.sawCall(wire.MethodDockableWindowsSet) {
		t.Error("Setup did not show the dockable panel")
	}
}

// TestNotifyPanelEdit checks a panel.valueChanged event for our panel updates params inline.
func TestNotifyPanelEdit(t *testing.T) {
	e := NewEngine(&fakeHost{})
	ev, _ := json.Marshal(map[string]string{
		"type": wire.EventPanelValueChanged, "windowId": TraceonPanelID, "controlId": "voltage", "value": "4200",
	})
	e.Notify(ev)
	e.mu.Lock()
	defer e.mu.Unlock()
	if e.params.voltage != 4200 {
		t.Errorf("voltage = %v, want 4200 after panel edit event", e.params.voltage)
	}
}

// TestNotifyIgnoresOtherPanels checks edits to a different window are ignored.
func TestNotifyIgnoresOtherPanels(t *testing.T) {
	e := NewEngine(&fakeHost{})
	before := e.params.voltage
	ev, _ := json.Marshal(map[string]string{
		"type": wire.EventPanelValueChanged, "windowId": "some.other.panel", "controlId": "voltage", "value": "4200",
	})
	e.Notify(ev)
	if e.params.voltage != before {
		t.Error("edit to a foreign panel should not change our params")
	}
}

// TestSimNumFallback checks the lenient form-value parsing.
func TestSimNumFallback(t *testing.T) {
	if got := simNum("", 7); got != 7 {
		t.Errorf("empty → %v, want fallback 7", got)
	}
	if got := simNum("not a number", 7); got != 7 {
		t.Errorf("garbage → %v, want fallback 7", got)
	}
	if got := simNum("3.5 kV", 7); got != 3.5 {
		t.Errorf("'3.5 kV' → %v, want 3.5", got)
	}
}

// TestPanelEditUpdatesParams checks a panel.valueChanged event updates the study parameters.
func TestPanelEditUpdatesParams(t *testing.T) {
	e := NewEngine(&fakeHost{})
	e.applyPanelEdit("voltage", "2500 V")
	e.applyPanelEdit("energy", "300")
	e.applyPanelEdit("rays", "11")
	e.mu.Lock()
	defer e.mu.Unlock()
	if e.params.voltage != 2500 {
		t.Errorf("voltage = %v, want 2500", e.params.voltage)
	}
	if e.params.energyEV != 300 {
		t.Errorf("energy = %v, want 300", e.params.energyEV)
	}
	if e.params.numRays != 11 {
		t.Errorf("numRays = %d, want 11", e.params.numRays)
	}
}
