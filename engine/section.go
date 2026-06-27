// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"fmt"
	"math"

	"oblikovati.org/api/wire"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/radial"
)

// sectionTolerance is the chord tolerance (cm, the host DB unit) for the boundary
// tessellation — fine enough that polyline edges track curved profiles.
const sectionTolerance = 0.05

// profile is the axisymmetric electrode cross-section extracted from a host body: the (r, z)
// boundary points of every polyline loop, mapped from the section plane (x→r, the in-plane
// vertical→z). A single body becomes one set of charged BEM elements.
type profile struct {
	loops [][]geom2d.Point2 // each polyline as (r, z) points
}

// extractProfile pulls the body's boundary polylines over api/client and maps them into the
// (r, z) half-plane the radial BEM works in. The section is the body's silhouette in the
// xz-plane: the host returns the strokes' x and the in-plane vertical, which become r and z.
func (e *Engine) extractProfile(bodyIndex int) (*profile, error) {
	strokes, err := e.api.Body().CalculateStrokes(bodyIndex, sectionTolerance)
	if err != nil {
		return nil, fmt.Errorf("calculate strokes: %w", err)
	}
	loops := loopsFromStrokes(strokes)
	if len(loops) == 0 {
		return nil, fmt.Errorf("body %d has no boundary polylines to section", bodyIndex)
	}
	return &profile{loops: loops}, nil
}

// loopsFromStrokes unflattens StrokeSetResult (flat XYZ coords + per-polyline lengths) into
// (r, z) loops: r = |x| (the radial distance), z = the in-plane height (the stroke's y).
func loopsFromStrokes(s wire.StrokeSetResult) [][]geom2d.Point2 {
	loops := make([][]geom2d.Point2, 0, s.PolylineCount)
	at := 0
	for _, n := range s.PolylineLengths {
		loop := make([]geom2d.Point2, 0, n)
		for i := 0; i < n; i++ {
			base := (at + i) * 3
			r := math.Abs(s.VertexCoordinates[base])
			z := s.VertexCoordinates[base+1]
			loop = append(loop, geom2d.Point2{r, z})
		}
		loops = append(loops, loop)
		at += n
	}
	return loops
}

// lineElements builds the radial BEM line elements (GMSH line4 cubics) from the profile: each
// consecutive (r, z) point pair becomes a straight cubic element whose control points are
// collinear (start, end, 1/3, 2/3). Returns the elements and a uniform excitation (every
// element a fixed-voltage electrode at the given voltage).
func (p *profile) lineElements(voltage float64) ([]radial.Line, []radial.ExcitationType, []float64) {
	var lines []radial.Line
	for _, loop := range p.loops {
		for i := 0; i+1 < len(loop); i++ {
			lines = append(lines, straightLine4(loop[i], loop[i+1]))
		}
	}
	types := make([]radial.ExcitationType, len(lines))
	values := make([]float64, len(lines))
	for i := range lines {
		types[i] = radial.VoltageFixed
		values[i] = voltage
	}
	return lines, types, values
}

// straightLine4 builds a GMSH line4 element for the straight (r, z) segment a→b: the two
// interior control points sit at the 1/3 and 2/3 parameter positions. The vertex layout is
// [start, end, 1/3, 2/3] with y = 0 (radial elements use components r=v[0], z=v[2]).
func straightLine4(a, b geom2d.Point2) radial.Line {
	lerp := func(t float64) geom2d.Vertex {
		return geom2d.Vertex{a[0] + (b[0]-a[0])*t, 0, a[1] + (b[1]-a[1])*t}
	}
	return radial.Line{lerp(0), lerp(1), lerp(1.0 / 3), lerp(2.0 / 3)}
}

// extent returns the (r, z) bounding box of the profile (cm), used to place the beam launch
// plane and the field-sampling grid.
func (p *profile) extent() (rMax, zMin, zMax float64) {
	rMax, zMin, zMax = 0, math.Inf(1), math.Inf(-1)
	for _, loop := range p.loops {
		for _, pt := range loop {
			rMax = math.Max(rMax, pt[0])
			zMin = math.Min(zMin, pt[1])
			zMax = math.Max(zMax, pt[1])
		}
	}
	return rMax, zMin, zMax
}
