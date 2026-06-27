// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"math"

	"oblikovati.org/api/types"
	"oblikovati.org/api/wire"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/tracing"
)

// graphicsClientID is the single client-graphics group holding the study overlay (electrode
// profiles + trajectories + potential map). Re-running the study replaces the whole group.
const graphicsClientID = "com.oblikovati.traceon.study"

// worldFromEngine maps an engine-frame point to the host world frame. The BEM/tracer work in the
// Traceon convention — radial in the XY plane, optical axis +Z — while the host bodies are surfaces
// of revolution about world Y (CalculateFacets / the section map host Y to the meridian's axial
// coordinate). Without this remap the overlay would render with its optical axis on world Z,
// perpendicular to the lens it describes (the field + trajectories would cross the electrode stack
// sideways). Swapping Y and Z puts the optical axis back on world Y, on top of the geometry.
func worldFromEngine(x, y, z float64) (float64, float64, float64) {
	return x, z, y
}

// appendWorld appends one engine-frame point to a flat xyz coordinate slice, remapped to world axes.
func appendWorld(coords []float64, x, y, z float64) []float64 {
	wx, wy, wz := worldFromEngine(x, y, z)
	return append(coords, wx, wy, wz)
}

// potentialGrid is the (r, z) sampling resolution for the potential heatmap.
const potentialGrid = 40

// renderNodes assembles the study overlay (all coordinates in the host DB unit, cm): the
// potential heatmap (drawn underneath), the electrode and coil profiles, and the traced
// trajectories. The BEM field and the rays are in metres, so positions are scaled by metresToCm.
func renderNodes(electrodes []electrode, coils []coil, magnets []magnet, irons []iron, bem field.FieldRadialBEM, rays [][]tracing.State) []wire.GraphicsNode {
	var nodes []wire.GraphicsNode
	if len(electrodes) > 0 {
		nodes = append(nodes, potentialNode(electrodes, bem), electrodeNode(electrodes))
	}
	if lines, ok := coilNode(coils); ok {
		nodes = append(nodes, outlineNode("traceon.coils", lines, []float32{0.85, 0.45, 0.2, 1})) // copper
	}
	if lines, ok := magnetNode(magnets); ok {
		nodes = append(nodes, outlineNode("traceon.magnets", lines, []float32{0.6, 0.3, 0.8, 1})) // violet
	}
	if lines, ok := ironNode(irons); ok {
		nodes = append(nodes, outlineNode("traceon.iron", lines, []float32{0.5, 0.55, 0.6, 1})) // steel grey
	}
	nodes = append(nodes, trajectoryNodes(rays)...)
	return nodes
}

// outlineNode wraps accumulated (r,z) line segments into an on-top coloured graphics node.
func outlineNode(id string, lines graphicsLines, color []float32) wire.GraphicsNode {
	return wire.GraphicsNode{Id: id, Primitives: []wire.GraphicsPrimitive{{
		Kind: string(types.GraphicsLines), Coordinates: lines.coords, Indices: lines.indices,
		Color: color, OnTop: true,
	}}}
}

// pushGraphics replaces the study client-graphics group with the supplied nodes.
func (e *Engine) pushGraphics(nodes []wire.GraphicsNode) error {
	_, err := e.api.Graphics().Set(wire.SetClientGraphicsArgs{
		ClientId: graphicsClientID,
		Lane:     string(types.GraphicsLanePersistent),
		Nodes:    nodes,
	})
	return err
}

// electrodeNode draws every sectioned electrode profile as line segments in the xz-plane
// (x = r, y = 0, z), in cm — so the overlay tracks the host geometry exactly.
func electrodeNode(electrodes []electrode) wire.GraphicsNode {
	var coords []float64
	var indices []int
	for _, el := range electrodes {
		for _, loop := range el.prof.loops {
			for i := 0; i+1 < len(loop); i++ {
				base := len(coords) / 3
				coords = appendWorld(coords, loop[i][0], 0, loop[i][1])
				coords = appendWorld(coords, loop[i+1][0], 0, loop[i+1][1])
				indices = append(indices, base, base+1)
			}
		}
	}
	return wire.GraphicsNode{Id: "traceon.electrode", Primitives: []wire.GraphicsPrimitive{{
		Kind: string(types.GraphicsLines), Coordinates: coords, Indices: indices,
		Color: []float32{1, 0.85, 0.2, 1}, OnTop: true,
	}}}
}

// trajectoryNodes draws each ray as a connected line strip in the xz-plane, converting the
// traced positions from metres back to the host's cm.
func trajectoryNodes(rays [][]tracing.State) []wire.GraphicsNode {
	nodes := make([]wire.GraphicsNode, 0, len(rays))
	for i, ray := range rays {
		coords := make([]float64, 0, len(ray)*3)
		for _, s := range ray {
			coords = appendWorld(coords, s[0]*metresToCm, s[1]*metresToCm, s[2]*metresToCm)
		}
		nodes = append(nodes, wire.GraphicsNode{
			Id: "traceon.ray." + itoa(i),
			Primitives: []wire.GraphicsPrimitive{{
				Kind: string(types.GraphicsLineStrip), Coordinates: coords,
				Color: []float32{0.2, 1.0, 0.4, 1}, OnTop: true,
			}},
		})
	}
	return nodes
}

// potentialNode samples the electrostatic potential on an (r, z) grid over the study extent
// and draws it as a semi-transparent flood plot colored by a blue→white→red mapper. The grid
// is sampled in metres (the field's units) but its vertices are placed in cm for the viewport.
func potentialNode(electrodes []electrode, bem field.FieldRadialBEM) wire.GraphicsNode {
	rMaxCm, zMinCm, zMaxCm := combinedExtent(electrodes)
	// Extend the grid downstream over the drift region (matching the trace) so the focus
	// crossing is drawn on the same plane as the trajectories.
	r0, r1 := 0.0, driftRadius*rMaxCm
	z0, z1 := zMinCm-boundsMargin*metresToCm, zMaxCm+driftFactor*(zMaxCm-zMinCm)

	var coords, scalars []float64
	for iz := 0; iz < potentialGrid; iz++ {
		for ir := 0; ir < potentialGrid; ir++ {
			rCm := r0 + (r1-r0)*float64(ir)/float64(potentialGrid-1)
			zCm := z0 + (z1-z0)*float64(iz)/float64(potentialGrid-1)
			coords = appendWorld(coords, rCm, 0, zCm)
			scalars = append(scalars, bem.PotentialAtPoint(geom2d.Vertex{rCm * cmToMetres, 0, zCm * cmToMetres}))
		}
	}
	indices := gridTriangles(potentialGrid, potentialGrid)

	mapper := potentialMapper(scalars)
	return wire.GraphicsNode{Id: "traceon.potential", Opacity: 0.5, Primitives: []wire.GraphicsPrimitive{{
		Kind: string(types.GraphicsTriangles), Coordinates: coords, Indices: indices,
		Scalars: scalars, ColorMapper: &mapper, ColorBinding: string(types.GraphicsColorPerVertex),
	}}}
}

// gridTriangles returns the triangle index list for a rows×cols vertex grid (row-major).
func gridTriangles(rows, cols int) []int {
	indices := make([]int, 0, (rows-1)*(cols-1)*6)
	for iz := 0; iz < rows-1; iz++ {
		for ir := 0; ir < cols-1; ir++ {
			a := iz*cols + ir
			b := a + 1
			c := a + cols
			d := c + 1
			indices = append(indices, a, b, c, b, d, c)
		}
	}
	return indices
}

// potentialMapper builds a blue→white→red color mapper spanning the sampled potential range
// (symmetric about zero so ground reads white).
func potentialMapper(scalars []float64) wire.GraphicsColorMapper {
	peak := 0.0
	for _, v := range scalars {
		if a := math.Abs(v); a > peak && !math.IsInf(a, 0) && !math.IsNaN(a) {
			peak = a
		}
	}
	if peak == 0 {
		peak = 1
	}
	return wire.GraphicsColorMapper{
		Values: []float64{-peak, 0, peak},
		Colors: []float32{
			0.1, 0.2, 0.9, 1, // -peak → blue
			0.95, 0.95, 0.95, 1, // 0 → white
			0.9, 0.15, 0.1, 1, // +peak → red
		},
	}
}

// itoa is a tiny integer formatter for node ids (avoids a strconv import here).
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	return string(buf[i:])
}
