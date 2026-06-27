// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"fmt"
	"math"

	"oblikovati.org/api/wire"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/radial"
)

// sectionTolerance is the chord tolerance (cm, the host DB unit) for the surface
// tessellation the meridian is extracted from.
const sectionTolerance = 0.1

// meridianBins is the number of axial bands used to extract the outer (r, z) profile.
const meridianBins = 60

// profile is the axisymmetric electrode meridian extracted from a host body: the outer
// (r, z) boundary of its revolved surface. A single body becomes one charged BEM electrode.
type profile struct {
	loops [][]geom2d.Point2 // the meridian as one ordered (r, z) polyline (in loops[0])
}

// extractProfile derives the axisymmetric (r, z) meridian of a host body from its surface
// facets: every facet vertex is mapped to (r = √(x²+z²), z = y) — its distance from the
// optical (Y) axis and its axial position — and the OUTER boundary r(z) is taken per axial
// band. This recovers the electrode's silhouette for any body of revolution (a body's edge
// strokes are only its rims, which collapse to points, so the facet envelope is used instead).
func (e *Engine) extractProfile(bodyIndex int) (*profile, error) {
	facets, err := e.api.Body().CalculateFacets(wire.CalculateFacetsArgs{BodyIndex: bodyIndex, Tolerance: sectionTolerance})
	if err != nil {
		return nil, fmt.Errorf("calculate facets: %w", err)
	}
	if facets.VertexCount == 0 {
		return nil, fmt.Errorf("body %d has no surface facets to section", bodyIndex)
	}
	loop := outerMeridian(facets.VertexCoordinates)
	if len(loop) < 2 {
		return nil, fmt.Errorf("body %d produced a degenerate (r,z) profile", bodyIndex)
	}
	return &profile{loops: [][]geom2d.Point2{loop}}, nil
}

// outerMeridian maps flat xyz vertices to (r = √(x²+z²), z = y) and returns the outer
// boundary r(z): the largest radius seen in each of meridianBins axial bands, as an ordered
// (r, z) polyline. Empty bands are skipped. This is the revolved surface's silhouette.
func outerMeridian(xyz []float64) []geom2d.Point2 {
	n := len(xyz) / 3
	if n == 0 {
		return nil
	}
	zMin, zMax := math.Inf(1), math.Inf(-1)
	for i := 0; i < n; i++ {
		y := xyz[i*3+1]
		zMin, zMax = math.Min(zMin, y), math.Max(zMax, y)
	}
	if zMax <= zMin {
		return nil
	}
	maxR := make([]float64, meridianBins)
	seen := make([]bool, meridianBins)
	for i := 0; i < n; i++ {
		x, y, z := xyz[i*3], xyz[i*3+1], xyz[i*3+2]
		r := math.Hypot(x, z)
		b := int((y - zMin) / (zMax - zMin) * float64(meridianBins-1))
		if b < 0 {
			b = 0
		}
		if b >= meridianBins {
			b = meridianBins - 1
		}
		if !seen[b] || r > maxR[b] {
			maxR[b], seen[b] = r, true
		}
	}
	loop := make([]geom2d.Point2, 0, meridianBins)
	for b := 0; b < meridianBins; b++ {
		if !seen[b] {
			continue
		}
		z := zMin + (zMax-zMin)*float64(b)/float64(meridianBins-1)
		loop = append(loop, geom2d.Point2{maxR[b], z})
	}
	return loop
}

// lineElements builds the radial BEM line elements (GMSH line4 cubics) from the profile: each
// consecutive (r, z) point pair becomes a straight cubic element. Returns the elements and a
// uniform excitation (every element a fixed-voltage electrode at the given voltage).
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
