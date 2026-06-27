// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"fmt"
	"math"
)

// discretizePath returns the sorted arc-length sample parameters for meshing a path of the
// given length with the given corner breakpoints. The result always contains 0, every
// breakpoint, and the full length, with at least three elements per smooth segment.
//
// nFactor splices extra interior samples between element endpoints: 1 for flat 2-node
// lines, 3 for curved 4-node "line4" elements (keeping the element count the same while
// adding the two interior nodes). Exactly one of meshSize / meshSizeFactor must be > 0:
// meshSize bounds each element's arc length; meshSizeFactor sets a constant per-segment
// count independent of length. Mirrors geometry.discretize_path.
func discretizePath(length float64, breakpoints []float64, meshSize, meshSizeFactor float64, nFactor int) []float64 {
	if meshSize <= 0 && meshSizeFactor <= 0 {
		panic(fmt.Sprintf("discretizePath: one of meshSize (%g) / meshSizeFactor (%g) must be > 0",
			meshSize, meshSizeFactor))
	}
	points := make([]float64, 0, len(breakpoints)+2)
	points = append(points, 0)
	points = append(points, breakpoints...)
	points = append(points, length)

	var u []float64
	for i := 0; i+1 < len(points); i++ {
		u0, u1 := points[i], points[i+1]
		if u0 == u1 {
			continue
		}
		n := segmentCount(u1-u0, meshSize, meshSizeFactor)
		u = append(u, linspaceOpen(u0, u1, nFactor*n)...)
	}
	return append(u, length)
}

// segmentCount is the number of flat elements spanning one smooth segment of arc length
// span: ceil(span/meshSize) (at least 3) when meshSize is set, else 3·max(factor, 1).
func segmentCount(span, meshSize, meshSizeFactor float64) int {
	if meshSize > 0 {
		n := int(math.Ceil(span / meshSize))
		if n < 3 {
			return 3
		}
		return n
	}
	return int(3 * math.Max(meshSizeFactor, 1))
}

// linspaceOpen returns n evenly-spaced values over [start, end) (endpoint excluded, like
// numpy.linspace(..., endpoint=False)), so the next segment's start is not duplicated.
func linspaceOpen(start, end float64, n int) []float64 {
	out := make([]float64, n)
	step := (end - start) / float64(n)
	for i := 0; i < n; i++ {
		out[i] = start + float64(i)*step
	}
	return out
}
