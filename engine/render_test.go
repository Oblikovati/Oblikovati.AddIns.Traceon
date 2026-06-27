// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"testing"

	"oblikovati.org/traceon/core/geom2d"
)

// TestWorldFromEngineAlignsOpticalAxisToY locks the overlay axis convention: the engine's optical
// axis (+Z, the Traceon tracer direction) must map to world Y, because the host bodies are surfaces
// of revolution about Y. A regression here renders the field + rays perpendicular to the lens.
func TestWorldFromEngineAlignsOpticalAxisToY(t *testing.T) {
	// An engine point one unit along the optical axis (z) must land on world +Y.
	if x, y, z := worldFromEngine(0, 0, 1); x != 0 || y != 1 || z != 0 {
		t.Errorf("optical axis: worldFromEngine(0,0,1) = (%g,%g,%g), want (0,1,0)", x, y, z)
	}
	// The radial direction (engine x) stays on world X.
	if x, y, z := worldFromEngine(2, 0, 0); x != 2 || y != 0 || z != 0 {
		t.Errorf("radial: worldFromEngine(2,0,0) = (%g,%g,%g), want (2,0,0)", x, y, z)
	}
}

// TestElectrodeNodePlacesAxialOnY checks an electrode meridian segment along the axis (constant r,
// varying axial) renders as world points whose axial coordinate is world Y, not Z.
func TestElectrodeNodePlacesAxialOnY(t *testing.T) {
	prof := &profile{loops: [][]geom2d.Point2{{{0.5, -1}, {0.5, 1}}}} // r=0.5, axial −1 → 1
	node := electrodeNode([]electrode{{prof: prof}})
	coords := node.Primitives[0].Coordinates
	if len(coords) < 6 {
		t.Fatalf("electrode node has too few coordinates: %d", len(coords))
	}
	// First emitted world point is the engine meridian point (0.5, axial=-1) → (x=0.5, y=-1, z=0).
	x, y, z := coords[0], coords[1], coords[2]
	if x != 0.5 || y != -1 || z != 0 {
		t.Errorf("electrode world point = (%g,%g,%g), want (0.5,-1,0) — axial must be on world Y", x, y, z)
	}
}
