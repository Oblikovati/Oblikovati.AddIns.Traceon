// SPDX-License-Identifier: MPL-2.0

// Package engine is the host-facing core of the Traceon electron-optics add-in: it
// turns a host body into a radially-symmetric BEM study (section → solve → trace →
// render) using only the Apache-2.0 oblikovati.org/api client and the pure-Go
// numerics in ../core. The cgo c-shared shell (../export.go) owns the C ABI; this
// package owns the study pipeline and stays cgo-free so it unit-tests on every
// platform.
package engine

import (
	"encoding/json"
	"sync"

	"oblikovati.org/api/client"
	"oblikovati.org/api/wire"
)

// HostCaller is the transport the engine talks to the host through — exactly the
// api/client Caller contract, supplied by the cgo shell at Activate (or a fake in
// tests). Keeping it an interface here keeps this package cgo-free and testable.
type HostCaller interface {
	Call(method string, req []byte) ([]byte, error)
}

// Engine runs electron-optics studies against a live host.
type Engine struct {
	host HostCaller
	api  *client.Client

	mu      sync.Mutex // guards running
	running bool       // a study is in flight (coalesces overlapping command triggers)
}

// NewEngine binds the engine to the host transport.
func NewEngine(host HostCaller) *Engine {
	return &Engine{host: host, api: client.New(host)}
}

// RunStudyCommandID is the host command the add-in registers; firing it (a ribbon click or
// the MCP bridge's execute_command) runs the electron-optics study on the active part.
const RunStudyCommandID = "Traceon.RunStudy"

// RegisterCommands registers the study command with the host so it is invokable the same way a
// ribbon click is — including over the MCP bridge's execute_command. The host action is a no-op;
// executing the command fires command.started, which Notify turns into a study run.
func (e *Engine) RegisterCommands() error {
	_, err := e.api.Commands().Create(wire.CreateCommandArgs{
		ID:          RunStudyCommandID,
		DisplayName: "Run Electron-Optics Study",
		Category:    "Traceon",
		Tooltip:     "Solve the radial BEM field for the active geometry and trace particle trajectories.",
	})
	return err
}

// Setup performs the one-time host-facing initialization. It MUST NOT run on the host's session
// goroutine (host calls there block until the frame loop drains the dispatcher, deadlocking the
// head) — the cgo shell runs it on its own goroutine.
func (e *Engine) Setup() error {
	return e.RegisterCommands()
}

// Notify receives host event bytes. A command.started carrying RunStudyCommandID runs the study
// on a SEPARATE goroutine — never inline, because Notify is invoked on the host's session
// goroutine and a host call from there blocks until the frame loop drains the dispatcher (which
// cannot happen while we're inside it), deadlocking every host call. A guard coalesces
// overlapping triggers so one study is in flight at a time.
func (e *Engine) Notify(ev []byte) {
	var hdr struct {
		Type string `json:"type"`
	}
	if json.Unmarshal(ev, &hdr) != nil {
		return
	}
	if hdr.Type != wire.EventCommandStarted {
		return
	}
	var c struct {
		Command string `json:"command"`
	}
	if json.Unmarshal(ev, &c) == nil && c.Command == RunStudyCommandID {
		e.launchStudy()
	}
}

// launchStudy starts one study goroutine, coalescing overlapping triggers. The study pipeline
// itself (section → solve → trace → render) is filled in over M2; M0 wires the command path.
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
		_ = e.runStudy()
	}()
}

// runStudy is the study pipeline entry point. M0 stub: reports readiness to the status bar so a
// command trigger is observably handled. Replaced by the real section→solve→trace→render flow in M2.
func (e *Engine) runStudy() error {
	_, err := e.api.Status().SetText("Traceon: study pipeline not yet implemented")
	return err
}
