//go:build js && wasm

// SPDX-License-Identifier: MPL-2.0

package main

import (
	"sync"

	"oblikovati.org/api/client"
	"oblikovati.org/api/wasmaddin"
	"oblikovati.org/traceon/engine"
)

const addInID = "com.oblikovati.traceon"

type wasmTraceon struct {
	mu     sync.Mutex
	active *engine.Engine
}

func (a *wasmTraceon) ID() string       { return addInID }
func (a *wasmTraceon) Manifest() string { return manifestJSON }

func (a *wasmTraceon) Activate(caller client.Caller) error {
	a.mu.Lock()
	defer a.mu.Unlock()
	eng := engine.NewEngine(caller)
	a.active = eng
	go func() { _ = eng.Setup() }()
	return nil
}

func (a *wasmTraceon) Deactivate() error {
	a.mu.Lock()
	a.active = nil
	a.mu.Unlock()
	return nil
}

func (a *wasmTraceon) Notify(event []byte) {
	a.mu.Lock()
	eng := a.active
	a.mu.Unlock()
	if eng != nil {
		eng.Notify(event)
	}
}

func main() {
	wasmaddin.Serve(&wasmTraceon{})
}
