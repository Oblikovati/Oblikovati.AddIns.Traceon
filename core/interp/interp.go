// SPDX-License-Identifier: MPL-2.0

// Package interp provides the 1-D piecewise-polynomial interpolation the axial field
// series needs, matched to the scipy constructions the upstream Traceon uses:
//
//   - NotAKnotCubic: scipy.interpolate.CubicSpline's default (not-a-knot) cubic spline,
//     used for the two highest axial-derivative orders.
//   - QuinticHermite: scipy BPoly.from_derivatives → PPoly, a per-interval quintic Hermite
//     matching value + 1st + 2nd derivative at both ends, used for the lower orders.
//
// Both return coefficients in DESCENDING powers of the local variable t = x − x[i] on each
// interval, the layout scipy's PPoly.c.T produces and the radial-derivative evaluators expect.
// The axis samples x are assumed equally spaced (the upstream always builds them with
// linspace), which makes the not-a-knot system tridiagonal and exactly solvable.
package interp

import "fmt"

// CubicCoeffs are the per-interval cubic coefficients [c3, c2, c1, c0] (descending powers
// of t = x − x[i]); there is one row per interval (len(x)-1 rows).
type CubicCoeffs [][4]float64

// QuinticCoeffs are the per-interval quintic coefficients [c5, c4, c3, c2, c1, c0].
type QuinticCoeffs [][6]float64

// NotAKnotCubic returns the not-a-knot cubic spline of (x, y) as per-interval power-basis
// coefficients, equivalent to scipy.interpolate.CubicSpline(x, y).c.T for equally-spaced x.
//
// It solves the tridiagonal system for the knot slopes s_i with scipy's not-a-knot end
// rows, then forms each interval's Hermite cubic in power form. Requires len(x) ≥ 3.
func NotAKnotCubic(x, y []float64) (CubicCoeffs, error) {
	n := len(x)
	if n < 3 {
		return nil, fmt.Errorf("interp.NotAKnotCubic: need >=3 points, got %d", n)
	}
	if len(y) != n {
		return nil, fmt.Errorf("interp.NotAKnotCubic: len(y) %d != len(x) %d", len(y), n)
	}
	h := x[1] - x[0]
	slope := make([]float64, n-1) // slope[i] = (y[i+1]-y[i])/h
	for i := 0; i < n-1; i++ {
		slope[i] = (y[i+1] - y[i]) / h
	}

	// Tridiagonal system A·s = d for the knot slopes (equal spacing). Interior rows are the
	// C2-continuity equations; the first/last rows are scipy's not-a-knot conditions.
	lower := make([]float64, n) // sub-diagonal
	diag := make([]float64, n)
	upper := make([]float64, n) // super-diagonal
	rhs := make([]float64, n)

	diag[0], upper[0], rhs[0] = 1, 2, (5*slope[0]+slope[1])/2
	for i := 1; i < n-1; i++ {
		lower[i], diag[i], upper[i] = 1, 4, 1
		rhs[i] = 3 * (slope[i-1] + slope[i])
	}
	lower[n-1], diag[n-1], rhs[n-1] = 2, 1, (slope[n-3]+5*slope[n-2])/2

	s := solveTridiagonal(lower, diag, upper, rhs)

	coeffs := make(CubicCoeffs, n-1)
	for i := 0; i < n-1; i++ {
		// Hermite cubic on [x_i, x_{i+1}] in power basis of t = x - x_i.
		c3 := (s[i] + s[i+1] - 2*slope[i]) / (h * h)
		c2 := (3*slope[i] - 2*s[i] - s[i+1]) / h
		coeffs[i] = [4]float64{c3, c2, s[i], y[i]}
	}
	return coeffs, nil
}

// QuinticHermite returns the per-interval quintic that matches value y, first derivative
// dy, and second derivative d2y at every knot — the piecewise polynomial scipy builds via
// BPoly.from_derivatives([y, dy, d2y]) → PPoly, in power-basis descending coefficients.
// Each interval is an independent local Hermite solve (no global system).
func QuinticHermite(x, y, dy, d2y []float64) (QuinticCoeffs, error) {
	n := len(x)
	if n < 2 {
		return nil, fmt.Errorf("interp.QuinticHermite: need >=2 points, got %d", n)
	}
	if len(y) != n || len(dy) != n || len(d2y) != n {
		return nil, fmt.Errorf("interp.QuinticHermite: y/dy/d2y must all have length %d", n)
	}
	coeffs := make(QuinticCoeffs, n-1)
	for i := 0; i < n-1; i++ {
		coeffs[i] = quinticHermiteInterval(x[i+1]-x[i], y[i], dy[i], d2y[i], y[i+1], dy[i+1], d2y[i+1])
	}
	return coeffs, nil
}

// quinticHermiteInterval solves the unique quintic p on [0, h] (local t) with p(0)=y0,
// p'(0)=v0, p”(0)=a0 and p(h)=y1, p'(h)=v1, p”(h)=a1, returning [c5,c4,c3,c2,c1,c0].
// The closed form is the classic quintic Hermite (b0=y0, b1=v0, b2=a0/2; the cubic..quintic
// terms from the three end conditions at t=h).
func quinticHermiteInterval(h, y0, v0, a0, y1, v1, a1 float64) [6]float64 {
	b0, b1, b2 := y0, v0, a0/2
	// Residual end conditions after subtracting the low-order Taylor part at t=h.
	bigY := y1 - (y0 + v0*h + (a0/2)*h*h)
	bigV := v1 - (v0 + a0*h)
	bigA := a1 - a0
	h2 := h * h
	b3 := (10*bigY - 4*bigV*h + bigA*h2/2) / (h * h2)
	b4 := (-15*bigY + 7*bigV*h - bigA*h2) / (h2 * h2)
	b5 := (6*bigY - 3*bigV*h + bigA*h2/2) / (h2 * h2 * h)
	return [6]float64{b5, b4, b3, b2, b1, b0}
}

// solveTridiagonal solves the tridiagonal system with sub-diagonal lower[1..n-1], diagonal
// diag[0..n-1], super-diagonal upper[0..n-2], and right-hand side rhs, via the Thomas
// algorithm. Returns the solution vector (length n).
func solveTridiagonal(lower, diag, upper, rhs []float64) []float64 {
	n := len(diag)
	cp := make([]float64, n)
	dp := make([]float64, n)
	cp[0] = upper[0] / diag[0]
	dp[0] = rhs[0] / diag[0]
	for i := 1; i < n; i++ {
		m := diag[i] - lower[i]*cp[i-1]
		cp[i] = upper[i] / m
		dp[i] = (rhs[i] - lower[i]*dp[i-1]) / m
	}
	x := make([]float64, n)
	x[n-1] = dp[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = dp[i] - cp[i]*x[i+1]
	}
	return x
}
