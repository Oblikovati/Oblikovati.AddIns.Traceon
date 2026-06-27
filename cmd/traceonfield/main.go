// SPDX-License-Identifier: MPL-2.0

// Command traceonfield runs the Traceon study engine against a built-in capture host and
// writes the exact client-graphics payload the add-in would push to the viewport. It needs
// no live app or render backend, so it is a replayable integration harness: drive the real
// section→solve→trace→render pipeline on a canned electrode and inspect (or replay) the
// resulting electrode / trajectory / potential-map graphics.
//
// Usage:
//
//	go run ./cmd/traceonfield > study.json
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"

	"oblikovati.org/api/wire"
	"oblikovati.org/traceon/engine"
)

// captureHost answers the wire methods the study issues and records the client-graphics
// payload instead of forwarding it to a real viewport. Its single solid body facets to a
// cylindrical electrode about the Y axis.
type captureHost struct {
	graphics *wire.SetClientGraphicsArgs
}

func (h *captureHost) Call(method string, req []byte) ([]byte, error) {
	switch method {
	case wire.MethodBodyList:
		return json.Marshal(wire.BodyListResult{Bodies: []wire.BodyInfo{{Index: 0, Name: "Electrode", Solid: true, Key: "k0"}}})
	case wire.MethodBodyCalculateFacets:
		return json.Marshal(cylinderFacets())
	case wire.MethodDocumentsList:
		return json.Marshal(wire.ListDocumentsResult{Documents: []wire.DocumentInfo{{ID: 1, Active: true}}})
	case wire.MethodAttributesGet:
		return json.Marshal(wire.AttributeResult{Found: false})
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

// cylinderFacets returns facet vertices for a cylinder of radius 1 cm about the Y axis,
// spanning y∈[-2, 2] cm — the meridian extractor turns this into an r=1 electrode profile.
func cylinderFacets() wire.FacetSetResult {
	var coords []float64
	for iy := 0; iy <= 8; iy++ {
		y := -2.0 + 4.0*float64(iy)/8
		for a := 0; a < 12; a++ {
			ang := 2 * math.Pi * float64(a) / 12
			coords = append(coords, math.Cos(ang), y, math.Sin(ang))
		}
	}
	return wire.FacetSetResult{VertexCount: len(coords) / 3, VertexCoordinates: coords}
}

func main() {
	host := &captureHost{}
	res, err := engine.NewEngine(host).RunStudy(0)
	if err != nil {
		fmt.Fprintln(os.Stderr, "study failed:", err)
		os.Exit(1)
	}
	if host.graphics == nil {
		fmt.Fprintln(os.Stderr, "study pushed no graphics")
		os.Exit(1)
	}
	fmt.Fprintf(os.Stderr, "study: %d electrode(s), %d elements, %d rays, focus z = %.3f cm, %d graphics nodes\n",
		res.ElectrodeCount, res.ElementCount, res.RayCount, res.FocusZ, len(host.graphics.Nodes))

	out, err := json.MarshalIndent(host.graphics, "", "  ")
	if err != nil {
		fmt.Fprintln(os.Stderr, "encode graphics:", err)
		os.Exit(1)
	}
	fmt.Println(string(out))
}
