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
// profile + trajectories + potential map). Re-running the study replaces the whole group.
const graphicsClientID = "com.oblikovati.traceon.study"

// potentialGrid is the (r, z) sampling resolution for the potential heatmap.
const potentialGrid = 40

// renderNodes assembles the study overlay: the potential heatmap (drawn underneath), the
// electrode profile, and the traced trajectories.
func renderNodes(prof *profile, bem field.FieldRadialBEM, rays [][]tracing.State) []wire.GraphicsNode {
	nodes := []wire.GraphicsNode{potentialNode(prof, bem)}
	nodes = append(nodes, electrodeNode(prof))
	nodes = append(nodes, trajectoryNodes(rays)...)
	return nodes
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

// electrodeNode draws the sectioned electrode profile as line segments in the xz-plane
// (x = r, y = 0, z), so the user sees exactly which boundary was solved.
func electrodeNode(prof *profile) wire.GraphicsNode {
	var coords []float64
	var indices []int
	for _, loop := range prof.loops {
		for i := 0; i+1 < len(loop); i++ {
			base := len(coords) / 3
			coords = append(coords, loop[i][0], 0, loop[i][1], loop[i+1][0], 0, loop[i+1][1])
			indices = append(indices, base, base+1)
		}
	}
	return wire.GraphicsNode{Id: "traceon.electrode", Primitives: []wire.GraphicsPrimitive{{
		Kind: string(types.GraphicsLines), Coordinates: coords, Indices: indices,
		Color: []float32{1, 0.85, 0.2, 1}, OnTop: true,
	}}}
}

// trajectoryNodes draws each ray as a connected line strip in the xz-plane.
func trajectoryNodes(rays [][]tracing.State) []wire.GraphicsNode {
	nodes := make([]wire.GraphicsNode, 0, len(rays))
	for i, ray := range rays {
		coords := make([]float64, 0, len(ray)*3)
		for _, s := range ray {
			coords = append(coords, s[0], s[1], s[2]) // (x, y, z); radial rays lie in y=0
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
// and draws it as a semi-transparent flood plot colored by a blue→red mapper. It is the field
// context the trajectories bend through.
func potentialNode(prof *profile, bem field.FieldRadialBEM) wire.GraphicsNode {
	rMax, zMin, zMax := prof.extent()
	r0, r1 := 0.0, rMax+boundsMargin
	z0, z1 := zMin-boundsMargin, zMax+boundsMargin

	var coords, scalars []float64
	for iz := 0; iz < potentialGrid; iz++ {
		for ir := 0; ir < potentialGrid; ir++ {
			r := r0 + (r1-r0)*float64(ir)/float64(potentialGrid-1)
			z := z0 + (z1-z0)*float64(iz)/float64(potentialGrid-1)
			coords = append(coords, r, 0, z)
			scalars = append(scalars, bem.PotentialAtPoint(geom2d.Vertex{r, 0, z}))
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

// itoa is a tiny non-allocating-ish integer formatter for node ids (avoids strconv import here).
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
