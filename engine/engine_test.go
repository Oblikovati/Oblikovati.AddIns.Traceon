// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"encoding/base64"
	"encoding/json"
	"errors"
	"math"
	"strings"
	"testing"

	"oblikovati.org/api/types"
	"oblikovati.org/api/wire"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
)

// fakeHost is a named fake HostCaller (no live host): it answers the wire methods a study
// issues with canned JSON, records the methods it saw, and captures the client-graphics
// payload so a test can assert the section→solve→trace→render pipeline ran end to end.
type fakeHost struct {
	calls        []string
	failOn       string // method to fail, for error-path tests ("" = none)
	facets       wire.FacetSetResult
	bodies       []wire.BodyInfo               // nil → one solid electrode body
	voltagesJSON string                        // traceon/voltages attribute payload ("" = unset)
	currentsJSON string                        // traceon/currents attribute payload ("" = unset)
	magnetsJSON  string                        // traceon/magnets attribute payload ("" = unset)
	permJSON     string                        // traceon/permeability attribute payload ("" = unset)
	materials    map[string]wire.MaterialInfo  // id → material, for materials.get
	selRefs      []string                      // model.selection result
	refBodies    []wire.BodyTopology           // model.referenceKeys per-body topology
	attrStore    map[string]wire.AttributeInfo // per-target attributes ("target\x00name" → info)
	lastGraph    wire.SetClientGraphicsArgs
}

func (h *fakeHost) Call(method string, req []byte) ([]byte, error) {
	h.calls = append(h.calls, method)
	if method == h.failOn {
		return nil, errors.New("forced failure")
	}
	switch method {
	case wire.MethodBodyList:
		bodies := h.bodies
		if bodies == nil {
			bodies = []wire.BodyInfo{{Index: 0, Name: "Solid1", Solid: true, Key: "k0"}}
		}
		return json.Marshal(wire.BodyListResult{Bodies: bodies})
	case wire.MethodBodyCalculateFacets:
		return json.Marshal(h.facets)
	case wire.MethodDocumentsList:
		return json.Marshal(wire.ListDocumentsResult{Documents: []wire.DocumentInfo{{ID: 1, Active: true}}})
	case wire.MethodModelSelection:
		return json.Marshal(wire.SelectionResult{Refs: h.selRefs})
	case wire.MethodModelReferenceKeys:
		return json.Marshal(wire.ReferenceKeysResult{Bodies: h.refBodies})
	case wire.MethodAttributesSet:
		var a wire.SetAttributeArgs
		_ = json.Unmarshal(req, &a)
		if h.attrStore == nil {
			h.attrStore = map[string]wire.AttributeInfo{}
		}
		info := wire.AttributeInfo{Set: a.Set, Name: a.Name, Value: a.Value, Target: a.Target}
		h.attrStore[a.Target+"\x00"+a.Name] = info
		return json.Marshal(wire.AttributeResult{Attribute: info, Found: true})
	case wire.MethodAttributesGet:
		var args wire.GetAttributeArgs
		_ = json.Unmarshal(req, &args)
		if args.Target != "" {
			if info, ok := h.attrStore[args.Target+"\x00"+args.Name]; ok {
				return json.Marshal(wire.AttributeResult{Attribute: info, Found: true})
			}
			return json.Marshal(wire.AttributeResult{Found: false})
		}
		payload := h.voltagesJSON
		switch args.Name {
		case attrCurrents:
			payload = h.currentsJSON
		case attrMagnets:
			payload = h.magnetsJSON
		case attrPermeability:
			payload = h.permJSON
		}
		if payload == "" {
			return json.Marshal(wire.AttributeResult{Found: false})
		}
		return json.Marshal(wire.AttributeResult{Found: true, Attribute: wire.AttributeInfo{
			Set: attrSet, Name: args.Name, Value: types.StringVariant(payload)}})
	case wire.MethodMaterialsGet:
		var a wire.AssetRefArgs
		_ = json.Unmarshal(req, &a)
		return json.Marshal(h.materials[a.ID])
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

// ringHost is a fake whose single body is a current coil: a rectangular cross-section ring at
// r∈[2,2.5] cm, y∈[-0.25,0.25] cm, named so isCoil treats it as a coil.
func ringHost() *fakeHost {
	var coords []float64
	for _, r := range []float64{2.0, 2.25, 2.5} {
		for _, y := range []float64{-0.25, 0, 0.25} {
			// place around the ring at a few azimuths so √(x²+z²)=r
			for _, ang := range []float64{0, 1.57, 3.14, 4.71} {
				coords = append(coords, r*math.Cos(ang), y, r*math.Sin(ang))
			}
		}
	}
	return &fakeHost{
		facets: wire.FacetSetResult{VertexCount: len(coords) / 3, VertexCoordinates: coords},
		bodies: []wire.BodyInfo{{Index: 0, Name: "Coil1", Solid: true, Key: "c0"}},
	}
}

// TestCoilStudy checks a coil body is recognised, produces a current field, and the study
// reports it. The on-axis magnetic field of the energised coil must be non-zero (a current
// produces an axial field), confirming the magnetostatic path is wired through the engine.
func TestCoilStudy(t *testing.T) {
	res, err := NewEngine(ringHost()).RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	if res.CoilCount != 1 {
		t.Errorf("CoilCount = %d, want 1", res.CoilCount)
	}
	if res.ElectrodeCount != 0 {
		t.Errorf("ElectrodeCount = %d, want 0 (the body is a coil)", res.ElectrodeCount)
	}
}

// TestBuildCoilCharges checks the coil current field is non-zero on the axis at the coil's
// midplane (an energised loop produces an axial field there).
func TestBuildCoilCharges(t *testing.T) {
	h := ringHost()
	prof, err := NewEngine(h).extractProfile(0)
	if err != nil {
		t.Fatalf("extractProfile: %v", err)
	}
	cc := buildCoilCharges([]coil{{prof: prof, current: 1000}})
	if len(cc.Currents) == 0 {
		t.Fatal("no current rings built")
	}
	hz := cc.CurrentFieldAt(geom3d.Vec3{0, 0, 0}) // on-axis, mid-plane (metres)
	if hz[2] == 0 {
		t.Error("expected a non-zero on-axis axial field from the coil")
	}
}

// TestCoilByCurrentAttribute checks the traceon/currents attribute marks a body as a coil.
func TestCoilByCurrentAttribute(t *testing.T) {
	amps, ok := isCoil(3, "Solid7", map[int]float64{3: 2.5}, 1000)
	if !ok || amps != 2.5 {
		t.Errorf("isCoil(attr) = (%v, %v), want (2.5, true)", amps, ok)
	}
	if _, ok := isCoil(0, "Plain", nil, 1000); ok {
		t.Error("a plain body should not be a coil")
	}
}

// magnetHost is a fake whose single body is an axially-magnetised permanent magnet: a ring
// cross-section named so isMagnet treats it as a magnet.
func magnetHost() *fakeHost {
	h := ringHost()
	h.bodies = []wire.BodyInfo{{Index: 0, Name: "Magnet1", Solid: true, Key: "m0"}}
	return h
}

// TestMagnetStudy checks a permanent-magnet body is recognised, builds magnetic surface
// charges, and the study reports it with a non-zero on-axis field.
func TestMagnetStudy(t *testing.T) {
	res, err := NewEngine(magnetHost()).RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	if res.MagnetCount != 1 || res.CoilCount != 0 || res.ElectrodeCount != 0 {
		t.Errorf("counts = (e=%d,c=%d,m=%d), want (0,0,1)", res.ElectrodeCount, res.CoilCount, res.MagnetCount)
	}
}

// TestBuildMagnetCharges checks an axially-magnetised magnet produces non-zero magnetostatic
// surface charges (n_z·M is non-zero on the end caps).
func TestBuildMagnetCharges(t *testing.T) {
	prof, err := NewEngine(magnetHost()).extractProfile(0)
	if err != nil {
		t.Fatalf("extractProfile: %v", err)
	}
	mc := buildMagnetCharges([]magnet{{prof: prof, magnetisation: 1e6}})
	nonzero := false
	for _, c := range mc.Charges {
		if c != 0 {
			nonzero = true
		}
	}
	if !nonzero {
		t.Error("expected non-zero magnetic surface charges from an axial magnet")
	}
}

// TestMagnetByAttribute checks the traceon/magnets attribute marks a body as a magnet.
func TestMagnetByAttribute(t *testing.T) {
	m, ok := isMagnet(2, "Block", map[int]float64{2: 5e5}, 1e6)
	if !ok || m != 5e5 {
		t.Errorf("isMagnet(attr) = (%v, %v), want (5e5, true)", m, ok)
	}
}

// ironHost is a fake whose single body is magnetizable iron (named so isIron matches), shaped
// like the ring so it sections cleanly.
func ironHost() *fakeHost {
	h := ringHost()
	h.bodies = []wire.BodyInfo{{Index: 0, Name: "Iron1", Solid: true, Key: "i0"}}
	return h
}

// TestIronStudy checks a magnetizable iron body is recognised and the study reports it. With
// no coils or magnets the pre-field is zero, so the iron's induced charges are zero, but the
// classification and the magnetostatic solve path must run without error.
func TestIronStudy(t *testing.T) {
	res, err := NewEngine(ironHost()).RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	if res.IronCount != 1 || res.ElectrodeCount != 0 {
		t.Errorf("counts = (e=%d,iron=%d), want (0,1)", res.ElectrodeCount, res.IronCount)
	}
}

// materialHost is a ring body carrying the given magnetic material (no role-naming): the study
// must classify it from the material alone.
func materialHost(mat wire.MaterialInfo) *fakeHost {
	h := ringHost()
	h.bodies = []wire.BodyInfo{{Index: 0, Name: "Pole", Solid: true, Key: "m0", MaterialID: "mat"}}
	h.materials = map[string]wire.MaterialInfo{"mat": mat}
	return h
}

// TestMaterialDrivenIron checks a body whose assigned material is soft-magnetic is classified as
// magnetizable iron (with the material's μr), without any "iron" in its name.
func TestMaterialDrivenIron(t *testing.T) {
	mat := wire.MaterialInfo{ID: "mat", Magnetic: types.Magnetic{Class: types.SoftMagnetic, RelativePermeability: 2500}}
	res, err := NewEngine(materialHost(mat)).RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	if res.IronCount != 1 || res.ElectrodeCount != 0 {
		t.Errorf("counts = (e=%d, iron=%d), want (0, 1)", res.ElectrodeCount, res.IronCount)
	}
}

// TestMaterialRole maps the magnetic class to a role + value: soft → iron(μr), hard → magnet
// (M = Br/μ0), non-magnetic → no role.
func TestMaterialRole(t *testing.T) {
	e := NewEngine(&fakeHost{materials: map[string]wire.MaterialInfo{
		"soft": {Magnetic: types.Magnetic{Class: types.SoftMagnetic, RelativePermeability: 1000}},
		"hard": {Magnetic: types.Magnetic{Class: types.HardMagnetic, Remanence: 1.3}},
		"alum": {Magnetic: types.Magnetic{Class: types.NonMagnetic}},
	}})
	if role, v, ok := e.materialRole("soft"); !ok || role != "iron" || v != 1000 {
		t.Errorf("soft → (%q, %g, %v), want (iron, 1000, true)", role, v, ok)
	}
	if role, v, ok := e.materialRole("hard"); !ok || role != "magnet" || math.Abs(v-1.3/1.25663706127e-06) > 1 {
		t.Errorf("hard → (%q, %g, %v), want (magnet, ~1.03e6, true)", role, v, ok)
	}
	if _, _, ok := e.materialRole("alum"); ok {
		t.Error("non-magnetic material should yield no role")
	}
	if _, _, ok := e.materialRole(""); ok {
		t.Error("empty material id should yield no role")
	}
}

// TestIronByAttribute checks the traceon/permeability attribute marks a body as iron.
func TestIronByAttribute(t *testing.T) {
	mu, ok := isIron(4, "Plate", map[int]float64{4: 4000}, 1000)
	if !ok || mu != 4000 {
		t.Errorf("isIron(attr) = (%v, %v), want (4000, true)", mu, ok)
	}
	if _, ok := isIron(0, "Plain", nil, 1000); ok {
		t.Error("a plain body should not be iron")
	}
}

// TestCombineCharges checks two charge sets concatenate their charges and buffers.
func TestCombineCharges(t *testing.T) {
	prof, err := NewEngine(ringHost()).extractProfile(0)
	if err != nil {
		t.Fatalf("extractProfile: %v", err)
	}
	a := buildMagnetCharges([]magnet{{prof: prof, magnetisation: 1e6}})
	combined := combineCharges(a, a)
	if len(combined.Charges) != 2*len(a.Charges) || len(combined.Jac) != 2*len(a.Jac) {
		t.Errorf("combined lengths = (%d charges, %d jac), want double of %d/%d",
			len(combined.Charges), len(combined.Jac), len(a.Charges), len(a.Jac))
	}
}

// TestFastTrace checks the fast axial-interpolation tracing path runs end to end and still
// produces a rendered study (the field evaluator swaps to the axial series).
func TestFastTrace(t *testing.T) {
	e := NewEngine(cylinderHost())
	e.applyPanelEdit("fast_trace", "1")
	res, err := e.RunStudy(0)
	if err != nil {
		t.Fatalf("fast-trace RunStudy: %v", err)
	}
	if res.RayCount == 0 {
		t.Error("fast trace produced no rays")
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

// boxFacets builds the facet set of a unit cube spanning x,y,z ∈ [-1,1] (8 corner vertices, 6
// quad faces) — a body that is NOT a surface of revolution about Y.
func boxFacets() wire.FacetSetResult {
	v := []float64{
		-1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, 1, // bottom (y=-1)
		-1, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1, // top (y=1)
	}
	faces := [][]int{
		{0, 1, 2, 3}, {4, 5, 6, 7}, // y caps
		{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, // sides (parallel to Y)
	}
	var idx, counts []int
	for _, f := range faces {
		idx = append(idx, f...)
		counts = append(counts, len(f))
	}
	return wire.FacetSetResult{VertexCount: 8, VertexCoordinates: v, VertexIndices: idx, IndexCountPerFace: counts}
}

// TestNonAxisymmetricBox checks a box is detected as non-axisymmetric (its faces dip to the
// inscribed radius at the edge midpoints while the corners stay at √2).
func TestNonAxisymmetricBox(t *testing.T) {
	if !nonAxisymmetric(boxFacets()) {
		t.Error("a box should be detected as non-axisymmetric")
	}
}

// TestAxisymmetricCylinderAccepted checks the revolved cylinder host's facets pass the check.
func TestAxisymmetricCylinderAccepted(t *testing.T) {
	h := cylinderHost()
	if nonAxisymmetric(h.facets) {
		t.Error("an axisymmetric cylinder was wrongly rejected")
	}
}

// TestExtractProfileRejectsBox checks sectioning a box returns an error instead of garbage.
func TestExtractProfileRejectsBox(t *testing.T) {
	h := ringHost()
	h.facets = boxFacets()
	if _, err := NewEngine(h).extractProfile(0); err == nil {
		t.Error("extractProfile accepted a non-axisymmetric box")
	}
}

// selectionHost is a ring body (axisymmetric) named without any role convention, whose face "f0"
// is owned by body "b0" (per the reference-keys topology), with the given face selected.
func selectionHost() *fakeHost {
	h := ringHost()
	h.bodies = []wire.BodyInfo{{Index: 0, Name: "Electrode", Solid: true, Key: "b0"}}
	// referenceKeys returns RAW keys; the viewport selection canonicalises to "<kind>/<base64>".
	rawFace := "\x03Revolution1:wall#1"
	h.refBodies = []wire.BodyTopology{{Faces: []wire.TopologyRef{{Key: rawFace}, {Key: "f1"}}}}
	h.selRefs = []string{"face/" + base64.RawURLEncoding.EncodeToString([]byte(rawFace))}
	return h
}

// TestNormalizeRef checks a canonical "<kind>/<base64>" selection ref decodes to the raw
// reference-key form (and a plain ref passes through).
func TestNormalizeRef(t *testing.T) {
	raw := "\x03Revolution1:wall#1"
	canonical := "face/" + base64.RawURLEncoding.EncodeToString([]byte(raw))
	if got := normalizeRef(canonical); got != raw {
		t.Errorf("normalizeRef(%q) = %q, want the raw key", canonical, got)
	}
	if got := normalizeRef("already-raw"); got != "already-raw" {
		t.Errorf("normalizeRef passthrough = %q, want unchanged", got)
	}
}

// TestFaceToBody checks a selected face resolves to its owning body's reference key.
func TestFaceToBody(t *testing.T) {
	owner, err := NewEngine(selectionHost()).faceToBody()
	if err != nil {
		t.Fatalf("faceToBody: %v", err)
	}
	// Every face/edge/vertex and the body's own key resolve to the body "b0".
	if len(owner) != 3 || owner["f1"] != "b0" || owner["b0"] != "b0" {
		t.Errorf("owner map = %v, want every entry → b0", owner)
	}
	for ref, body := range owner {
		if body != "b0" {
			t.Errorf("owner[%q] = %q, want b0", ref, body)
		}
	}
}

// TestAssignToSelectionThenClassify is the end-to-end face-pick flow: assign a coil current to the
// body owning the selected face, then a study reads the per-body attribute and classifies it.
func TestAssignToSelectionThenClassify(t *testing.T) {
	h := selectionHost()
	e := NewEngine(h)
	e.params.assignRole = "coil"
	e.params.assignValue = 7

	msg, err := e.assignToSelection()
	if err != nil {
		t.Fatalf("assignToSelection: %v", err)
	}
	if !strings.Contains(msg, "1 body") {
		t.Errorf("assign status = %q, want 1 body assigned", msg)
	}
	// The body now reads back as a 7 A coil.
	role, value, ok := e.bodyAssignment(1, "b0")
	if !ok || role != "coil" || value != 7 {
		t.Errorf("bodyAssignment = (%q, %g, %v), want (coil, 7, true)", role, value, ok)
	}
	// A study classifies it as a coil (overriding the "Electrode" name).
	res, err := e.RunStudy(0)
	if err != nil {
		t.Fatalf("RunStudy: %v", err)
	}
	if res.CoilCount != 1 || res.ElectrodeCount != 0 {
		t.Errorf("counts = (e=%d, coil=%d), want (0, 1)", res.ElectrodeCount, res.CoilCount)
	}
}

// TestAssignNothingSelected checks assigning with an empty selection is a friendly no-op.
func TestAssignNothingSelected(t *testing.T) {
	h := selectionHost()
	h.selRefs = nil
	msg, err := NewEngine(h).assignToSelection()
	if err != nil || !strings.Contains(msg, "select") {
		t.Errorf("assign with no selection = (%q, %v), want a 'select…' hint", msg, err)
	}
}
