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

// gaussKronrod evaluates the G7 and K15 estimates of ∫f over [a, b] on a single panel.
func gaussKronrod(f Func, a, b float64) (g7, k15 float64) {
	c := (a + b) / 2.0
	h := (b - a) / 2.0
	g7 = h * g7Weights[0] * f(c)
	k15 = h * k15Weights[0] * f(c)
	for i := 1; i < len(g7Weights); i++ {
		x := h * g7Nodes[i]
		g7 += h * g7Weights[i] * (f(c-x) + f(c+x))
	}
	for i := 1; i < len(k15Weights); i++ {
		x := h * k15Nodes[i]
		k15 += h * k15Weights[i] * (f(c-x) + f(c+x))
	}
	return g7, k15
}

// maxSubdivisions caps the total number of subintervals a single AdaptiveRecursive call
// may create (matching scipy.integrate.quad's limit=250). A global subdivision budget —
// rather than a per-branch recursion depth — is what keeps an integrand with an interior
// sign change (where |K15|→0 defeats the relative-tolerance check) from blowing up
// exponentially; work is bounded by this many G7/K15 panel evaluations.
const maxSubdivisions = 250

// subinterval is one panel of the global adaptive scheme: its bounds, K15 estimate, and
// the |K15-G7| error estimate that drives which panel is bisected next.
type subinterval struct {
	a, b, value, errEst float64
}

// AdaptiveRecursive integrates f over [a, b] with a global adaptive Gauss-Kronrod scheme:
// it repeatedly bisects the worst (largest-error) panel until the summed error meets the
// abs OR rel tolerance, or the subdivision budget is exhausted. This is the QUADPACK QAG
// strategy (without epsilon extrapolation), which concentrates panels at hard spots (e.g.
// an endpoint singularity) while staying bounded — unlike a depth-limited recursion, which
// branches exponentially near an interior zero of the integrand.
func AdaptiveRecursive(f Func, a, b, absTol, relTol float64) float64 {
	panels := make([]subinterval, 0, maxSubdivisions+1)
	panels = append(panels, newPanel(f, a, b))
	total, totalErr := panels[0].value, panels[0].errEst

	for n := 0; n < maxSubdivisions; n++ {
		if totalErr <= absTol || totalErr <= relTol*math.Abs(total) {
			break
		}
		// Bisect the panel with the largest error estimate.
		worst := 0
		for i := 1; i < len(panels); i++ {
			if panels[i].errEst > panels[worst].errEst {
				worst = i
			}
		}
		p := panels[worst]
		mid := (p.a + p.b) / 2.0
		if mid <= p.a || mid >= p.b { // panel collapsed to roundoff width; cannot refine further
			break
		}
		left := newPanel(f, p.a, mid)
		right := newPanel(f, mid, p.b)
		total += left.value + right.value - p.value
		totalErr += left.errEst + right.errEst - p.errEst
		panels[worst] = left
		panels = append(panels, right)
	}
	return total
}

func newPanel(f Func, a, b float64) subinterval {
	g7, k15 := gaussKronrod(f, a, b)
	return subinterval{a: a, b: b, value: k15, errEst: math.Abs(k15 - g7)}
}

// IntegrateWithSingularities integrates f over [a, b], splitting the interval at each
// interior singular point so the adaptive integrator approaches each singularity from an
// endpoint (where it converges). Mirrors scipy.integrate.quad(..., points=...).
func IntegrateWithSingularities(f Func, a, b float64, points []float64, absTol, relTol float64) float64 {
	bounds := []float64{a}
	for _, p := range points {
		if p > a && p < b {
			bounds = append(bounds, p)
		}
	}
	bounds = append(bounds, b)
	sum := 0.0
	for i := 0; i+1 < len(bounds); i++ {
		sum += AdaptiveRecursive(f, bounds[i], bounds[i+1], absTol, relTol)
	}
	return sum
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
