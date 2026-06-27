// SPDX-License-Identifier: MPL-2.0

// Package engine is the host-facing core of the Traceon electron-optics add-in: it turns a
// host body into a radially-symmetric BEM study (section → solve → trace → render) using
// only the Apache-2.0 oblikovati.org/api client and the pure-Go numerics in ../core. The
// cgo c-shared shell (../export.go) owns the C ABI; this package owns the study pipeline and
// stays cgo-free so it unit-tests on every platform.
package engine

import (
	"encoding/json"
	"fmt"
	"math"
	"sync"

	"oblikovati.org/api/client"
	"oblikovati.org/api/wire"
)

// HostCaller is the transport the engine talks to the host through — exactly the
// api/client Caller contract, supplied by the cgo shell at Activate (or a fake in tests).
type HostCaller interface {
	Call(method string, req []byte) ([]byte, error)
}

// studyParams are the user-editable simulation parameters (set from the dockable panel).
type studyParams struct {
	voltage       float64 // potential biasing the central electrode (volts)
	energyEV      float64 // initial kinetic energy of the traced beam (eV)
	numRays       int     // number of parallel rays launched
	coilCurrent   float64 // default current for coil bodies (amperes)
	magnetisation float64 // default axial magnetisation for magnet bodies (A/m)
	permeability  float64 // default relative permeability for iron bodies
	fastTrace     bool    // trace through the fast axial-series interpolation of the field

	// Parametric-lens definition: when lens != lensHost the study ignores host geometry and
	// builds the electrodes from core/geometry (dimensions in cm). lensRadius is the aperture/
	// bore radius, lensThickness each electrode's axial thickness, lensSpacing the inter-
	// electrode gap.
	lens          paramLens
	lensRadius    float64
	lensThickness float64
	lensSpacing   float64
}

func defaultParams() studyParams {
	return studyParams{
		voltage: 1000, energyEV: 1000, numRays: 7, coilCurrent: 1000, magnetisation: 1e6, permeability: 1000,
		lens: lensHost, lensRadius: 0.3, lensThickness: 0.5, lensSpacing: 0.5,
	}
}

// Engine runs electron-optics studies against a live host.
type Engine struct {
	host HostCaller
	api  *client.Client

	mu      sync.Mutex // guards params + running
	params  studyParams
	running bool // a study is in flight (coalesces overlapping command triggers)
}

// NewEngine binds the engine to the host transport with default simulation parameters.
func NewEngine(host HostCaller) *Engine {
	return &Engine{host: host, api: client.New(host), params: defaultParams()}
}

// RunStudyCommandID is the host command the add-in registers; firing it (a ribbon click or
// the MCP bridge's execute_command) runs the electron-optics study on the active part.
// EinzelLensCommandID / CylinderLensCommandID run a built-in parametric lens with no host
// geometry — they select the lens template (the panel dimensions still apply) and run, so the
// parametric study is invokable from the command palette and over the MCP bridge.
const (
	RunStudyCommandID     = "Traceon.RunStudy"
	EinzelLensCommandID   = "Traceon.RunEinzelLens"
	CylinderLensCommandID = "Traceon.RunCylinderLens"
)

// lensForCommand maps a parametric-lens command to the lens it selects, and reports whether the
// command is one of them.
func lensForCommand(id string) (paramLens, bool) {
	switch id {
	case EinzelLensCommandID:
		return lensEinzel, true
	case CylinderLensCommandID:
		return lensCylinder, true
	default:
		return lensHost, false
	}
}

// studySummary formats the one-line status reported after a completed study, including the
// axial focus when the beam crosses the optical axis.
func studySummary(res *StudyResult) string {
	s := fmt.Sprintf("Traceon: %d electrode(s), %d coil(s), %d magnet(s), %d iron, %d elements, %d rays",
		res.ElectrodeCount, res.CoilCount, res.MagnetCount, res.IronCount, res.ElementCount, res.RayCount)
	if !math.IsNaN(res.FocusZ) {
		s += fmt.Sprintf(" — focus z = %.3f cm", res.FocusZ)
	}
	return s
}

// RegisterCommands registers the study command with the host so it is invokable the same way a
// ribbon click is — including over the MCP bridge's execute_command. The host action is a no-op;
// executing the command fires command.started, which Notify turns into a study run.
func (e *Engine) RegisterCommands() error {
	cmds := []wire.CreateCommandArgs{
		{ID: RunStudyCommandID, DisplayName: "Run Electron-Optics Study", Category: "Traceon",
			Tooltip: "Solve the radial BEM field for the active geometry and trace particle trajectories."},
		{ID: EinzelLensCommandID, DisplayName: "Run Einzel Lens (parametric)", Category: "Traceon",
			Tooltip: "Build a three-aperture einzel lens parametrically (no host geometry) and run the study."},
		{ID: CylinderLensCommandID, DisplayName: "Run Two-Cylinder Lens (parametric)", Category: "Traceon",
			Tooltip: "Build a two-cylinder immersion lens parametrically (no host geometry) and run the study."},
	}
	for _, c := range cmds {
		if _, err := e.api.Commands().Create(c); err != nil {
			return err
		}
	}
	return nil
}

// Setup performs the one-time host-facing initialization: register the study command and show
// the simulation-parameters panel. It MUST NOT run on the host's session goroutine (host calls
// there block until the frame loop drains the dispatcher, deadlocking the head) — the cgo shell
// runs it on its own goroutine.
func (e *Engine) Setup() error {
	if err := e.RegisterCommands(); err != nil {
		return err
	}
	_, err := e.ShowPanel()
	return err
}

// Notify receives host event bytes. A command.started carrying RunStudyCommandID runs the study
// on a SEPARATE goroutine — never inline, because Notify is invoked on the host's session
// goroutine and a host call from there blocks until the frame loop drains the dispatcher (which
// cannot happen while we're inside it), deadlocking every host call. A panel.valueChanged only
// mutates engine state (no host call) so it is handled inline.
func (e *Engine) Notify(ev []byte) {
	var hdr struct {
		Type string `json:"type"`
	}
	if json.Unmarshal(ev, &hdr) != nil {
		return
	}
	switch hdr.Type {
	case wire.EventCommandStarted:
		var c struct {
			Command string `json:"command"`
		}
		if json.Unmarshal(ev, &c) != nil {
			return
		}
		if lens, ok := lensForCommand(c.Command); ok {
			e.mu.Lock()
			e.params.lens = lens
			e.mu.Unlock()
			e.launchStudy()
		} else if c.Command == RunStudyCommandID {
			e.launchStudy()
		}
	case wire.EventPanelValueChanged:
		var p struct {
			WindowId  string `json:"windowId"`
			ControlId string `json:"controlId"`
			Value     string `json:"value"`
		}
		if json.Unmarshal(ev, &p) == nil && p.WindowId == TraceonPanelID {
			e.applyPanelEdit(p.ControlId, p.Value)
		}
	}
}

// launchStudy starts one study goroutine, coalescing overlapping triggers, and reports the
// outcome to the host status bar so a failed study is visible rather than silently empty.
func (e *Engine) launchStudy() {
	e.mu.Lock()
	if e.running {
		e.mu.Unlock()
		return
	}
	e.running = true
	e.mu.Unlock()

	go func() {
		defer func() {
			e.mu.Lock()
			e.running = false
			e.mu.Unlock()
		}()
		res, err := e.RunStudy(0)
		if err != nil {
			_, _ = e.api.Status().SetText("Traceon study failed: " + err.Error())
			return
		}
		_, _ = e.api.Status().SetText(studySummary(res))
	}()
}
