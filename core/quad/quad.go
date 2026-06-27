// SPDX-License-Identifier: MPL-2.0

// Package quad provides adaptive Gauss-Kronrod (G7/K15) quadrature, a faithful port of
// traceon/backend/kronrod.c. It is the integrator the radial BEM kernels use for the
// singular self-element integrals, so its subdivision behaviour must match the upstream
// to keep the assembled matrix (and hence every field) equivalent.
package quad

import "math"

// Default tolerances match the Python wrapper (and scipy.integrate.quad), so a Go call
// with DefaultTol reproduces the upstream subdivision exactly.
const (
	DefaultAbsTol = 1.49e-08
	DefaultRelTol = 1.49e-08
)

// Gauss-7 weights/nodes (half-set: index 0 is the central node, the rest are mirrored
// about the interval centre). Verbatim from kronrod.c.
var (
	g7Weights = [4]float64{0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870}
	g7Nodes   = [4]float64{0.000000000000000, 0.405845151377397, 0.741531185599394, 0.949107912342759}

	k15Weights = [8]float64{
		0.209482141084728, 0.204432940075299, 0.190350578064785, 0.169004726639268,
		0.140653259715526, 0.104790010322250, 0.063092092629979, 0.022935322010529,
	}
	k15Nodes = [8]float64{
		0.000000000000000, 0.207784955007898, 0.405845151377397, 0.586087235467691,
		0.741531185599394, 0.864864423359769, 0.949107912342759, 0.991455371120813,
	}
)

// Func is the integrand: a scalar function of x. (The upstream passes an opaque args
// pointer; in Go the closure captures any parameters, so no args slot is needed.)
type Func func(x float64) float64

// Adaptive integrates f over [a, b] with the default tolerances, equivalent to the
// upstream kronrod_adaptive called from Python with its defaults.
func Adaptive(f Func, a, b float64) float64 {
	return AdaptiveTol(f, a, b, DefaultAbsTol, DefaultRelTol)
}

// AdaptiveTol integrates f over [a, b] using the G7/K15 pair with greedy bisection: on
// each panel it forms both estimates, and if their difference passes the abs OR rel
// tolerance it accepts K15 and advances; otherwise it bisects the panel's left half.
// This reproduces kronrod.c line-for-line (including its left-to-right sweep), so the
// number and placement of subdivisions match the oracle, not just the final value.
func AdaptiveTol(f Func, a, b, absTol, relTol float64) float64 {
	result := 0.0
	currentStart := a
	currentEnd := b

	for currentStart < b {
		c := (currentStart + currentEnd) / 2.0
		h := (currentEnd - currentStart) / 2.0

		g7 := h * g7Weights[0] * f(c)
		k15 := h * k15Weights[0] * f(c)

		for i := 1; i < len(g7Weights); i++ {
			x := h * g7Nodes[i]
			g7 += h * g7Weights[i] * (f(c-x) + f(c+x))
		}
		for i := 1; i < len(k15Weights); i++ {
			x := h * k15Nodes[i]
			k15 += h * k15Weights[i] * (f(c-x) + f(c+x))
		}

		errEst := math.Abs(k15 - g7)
		if errEst <= absTol || errEst <= relTol*math.Abs(k15) {
			result += k15
			currentStart = currentEnd
			currentEnd = b
		} else {
			currentEnd = c
		}
	}
	return result
}
