// SPDX-License-Identifier: MPL-2.0

package quad

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
)

// TestAdaptive ports tests/test_kronrod.py one-for-one. The expected values for the two
// hard integrands are the upstream Traceon backend's own outputs (the oracle), so this
// asserts the Go subdivision converges to the same number, not merely to the true value.
func TestAdaptive(t *testing.T) {
	cases := []struct {
		name string
		f    Func
		a, b float64
		want float64
	}{
		{"constant", func(x float64) float64 { return 1 }, 0, 1, 1.0},
		{"linear", func(x float64) float64 { return x }, 0, 1, 0.5},
		{"quadratic", func(x float64) float64 { return x * x }, 0, 1, 1.0 / 3.0},
		{"sine", math.Sin, 0, math.Pi, 2.0},
		{"difficult_exponential", func(x float64) float64 { return math.Exp(x * x * math.Cos(10*x)) }, -1.5, 1.5, 4.097655169215941},
		{"almost_singular", func(x float64) float64 { return 1 / (x + 0.001) }, 0, 1, 6.90875477931522},
	}
	for _, c := range cases {
		got := Adaptive(c.f, c.a, c.b)
		oracle.CheckClose(t, c.name, got, c.want)
	}
}
