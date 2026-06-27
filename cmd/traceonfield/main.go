// SPDX-License-Identifier: MPL-2.0

// Command traceonfield runs the Traceon study engine against a built-in capture host and
// writes the exact client-graphics payload the add-in would push to the viewport. It needs
// no live app or render backend, so it is a replayable integration harness: drive the real
// section→solve→trace→render pipeline on a canned electrode profile and inspect (or replay)
// the resulting electrode / trajectory / potential-map graphics.
//
// Usage:
//
//	go run ./cmd/traceonfield > study.json
package main

import (
	"encoding/json"
	"fmt"
	"os"

	"oblikovati.org/api/wire"
	"oblikovati.org/traceon/engine"
)

// captureHost answers the wire methods the study issues and records the client-graphics
// payload instead of forwarding it to a real viewport. Its body sections to a two-aperture
// electrode (a crude einzel-lens-like profile) so the trajectories visibly bend.
type captureHost struct {
	graphics *wire.SetClientGraphicsArgs
}

func (h *captureHost) Call(method string, req []byte) ([]byte, error) {
	switch method {
	case wire.MethodBodyCalculateStrokes:
		return json.Marshal(lensProfileStrokes())
	case wire.MethodClientGraphicsSet:
		var g wire.SetClientGraphicsArgs
		if err := json.Unmarshal(req, &g); err != nil {
			return nil, err
		}
		h.graphics = &g
		return []byte("{}"), nil
	default:
		return []byte("{}"), nil
	}
}

// lensProfileStrokes is a canned axisymmetric electrode: two coaxial cylinder walls at r=1,
// one below and one above the mid-plane, leaving a gap the beam focuses through.
func lensProfileStrokes() wire.StrokeSetResult {
	return wire.StrokeSetResult{
		VertexCount: 6,
		VertexCoordinates: []float64{
			1, -2.0, 0, 1, -0.5, 0, // lower electrode wall
			1, 0.5, 0, 1, 2.0, 0, // upper electrode wall
			1, -0.5, 0, 1, 0.5, 0, // gap bridge (kept thin)
		},
		PolylineCount:   3,
		PolylineLengths: []int{2, 2, 2},
	}
}

func main() {
	host := &captureHost{}
	eng := engine.NewEngine(host)
	res, err := eng.RunStudy(0)
	if err != nil {
		fmt.Fprintln(os.Stderr, "study failed:", err)
		os.Exit(1)
	}
	if host.graphics == nil {
		fmt.Fprintln(os.Stderr, "study pushed no graphics")
		os.Exit(1)
	}
	fmt.Fprintf(os.Stderr, "study: %d elements, %d rays, %d graphics nodes\n",
		res.ElementCount, res.RayCount, len(host.graphics.Nodes))

	out, err := json.MarshalIndent(host.graphics, "", "  ")
	if err != nil {
		fmt.Fprintln(os.Stderr, "encode graphics:", err)
		os.Exit(1)
	}
	fmt.Println(string(out))
}
