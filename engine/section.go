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
// optical (Y) axis and its axial position. The full cross-section boundary (outer wall, end
// caps, and any inner bore) is recovered per axial band, so the aperture fields that focus
// the beam are modelled — not just the outer silhouette. (A body's edge strokes are only its
// rims, which collapse to points, so the facet boundary is used instead.)
func (e *Engine) extractProfile(bodyIndex int) (*profile, error) {
	facets, err := e.api.Body().CalculateFacets(wire.CalculateFacetsArgs{BodyIndex: bodyIndex, Tolerance: sectionTolerance})
	if err != nil {
		return nil, fmt.Errorf("calculate facets: %w", err)
	}
	if facets.VertexCount == 0 {
		return nil, fmt.Errorf("body %d has no surface facets to section", bodyIndex)
	}
	loop := fullMeridian(facets.VertexCoordinates)
	if len(loop) < 2 {
		return nil, fmt.Errorf("body %d produced a degenerate (r,z) profile", bodyIndex)
	}
	return &profile{loops: [][]geom2d.Point2{loop}}, nil
}

// fullMeridian maps flat xyz vertices to (r = √(x²+z²), z = y) and returns the closed
// cross-section boundary as an ordered (r, z) loop: the outer wall (max radius per axial
// band, z increasing) followed by the inner wall / bore (min radius per band, z decreasing).
// For a solid electrode the inner wall hugs the axis; for a tube it traces the bore — so the
// aperture that focuses the beam is captured, not only the outer silhouette.
func fullMeridian(xyz []float64) []geom2d.Point2 {
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
	minR := make([]float64, meridianBins)
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
		if !seen[b] {
			minR[b], maxR[b], seen[b] = r, r, true
			continue
		}
		minR[b] = math.Min(minR[b], r)
		maxR[b] = math.Max(maxR[b], r)
	}
	zOf := func(b int) float64 { return zMin + (zMax-zMin)*float64(b)/float64(meridianBins-1) }

	// Outer wall up, then inner wall back down — a closed boundary loop.
	var loop []geom2d.Point2
	for b := 0; b < meridianBins; b++ {
		if seen[b] {
			loop = append(loop, geom2d.Point2{maxR[b], zOf(b)})
		}
	}
	for b := meridianBins - 1; b >= 0; b-- {
		if seen[b] {
			loop = append(loop, geom2d.Point2{minR[b], zOf(b)})
		}
	}
	return dedupeLoop(loop)
}

// coincidentTol is the (r, z) distance (cm) below which two profile points are treated as
// the same — used to collapse the inner wall onto the outer where a body has no bore (a
// solid mid-band or a zero-thickness shell), which would otherwise duplicate BEM elements
// and make the influence matrix singular.
const coincidentTol = 1e-4

// dedupeLoop drops any point that lies within coincidentTol of a point already kept, so
// coincident outer/inner walls collapse to a single wall while a genuine bore (distinct
// inner radius) is preserved.
func dedupeLoop(loop []geom2d.Point2) []geom2d.Point2 {
	out := make([]geom2d.Point2, 0, len(loop))
	for _, p := range loop {
		dup := false
		for _, q := range out {
			if math.Hypot(p[0]-q[0], p[1]-q[1]) < coincidentTol {
				dup = true
				break
			}
		}
		if !dup {
			out = append(out, p)
		}
	}
	return out
}

// segmentDivisions subdivides each meridian segment into this many BEM elements, so a flat
// electrode wall (which tessellates to only its two rims) is still discretized finely enough
// for an accurate boundary-element solve.
const segmentDivisions = 12

// lineElements builds the radial BEM line elements (GMSH line4 cubics) from the profile,
// scaling the (r, z) points by `scale` (e.g. cm→metres) so the solve runs in SI units. Each
// meridian segment is subdivided into segmentDivisions straight cubic elements at the given
// fixed voltage.
func (p *profile) lineElements(voltage, scale float64) ([]radial.Line, []radial.ExcitationType, []float64) {
	var lines []radial.Line
	for _, loop := range p.loops {
		for i := 0; i+1 < len(loop); i++ {
			a, b := loop[i], loop[i+1]
			for s := 0; s < segmentDivisions; s++ {
				t0 := float64(s) / segmentDivisions
				t1 := float64(s+1) / segmentDivisions
				sa := geom2d.Point2{a[0] + (b[0]-a[0])*t0, a[1] + (b[1]-a[1])*t0}
				sb := geom2d.Point2{a[0] + (b[0]-a[0])*t1, a[1] + (b[1]-a[1])*t1}
				lines = append(lines, straightLine4(sa, sb, scale))
			}
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

// straightLine4 builds a GMSH line4 element for the straight (r, z) segment a→b (scaled by
// `scale`): the two interior control points sit at the 1/3 and 2/3 parameter positions. The
// vertex layout is [start, end, 1/3, 2/3] with y = 0 (radial elements use r=v[0], z=v[2]).
func straightLine4(a, b geom2d.Point2, scale float64) radial.Line {
	lerp := func(t float64) geom2d.Vertex {
		return geom2d.Vertex{(a[0] + (b[0]-a[0])*t) * scale, 0, (a[1] + (b[1]-a[1])*t) * scale}
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
