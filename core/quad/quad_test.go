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

// TestAdaptiveRecursive checks the global adaptive scheme on smooth integrands.
func TestAdaptiveRecursive(t *testing.T) {
	oracle.CheckClose(t, "x^2", AdaptiveRecursive(func(x float64) float64 { return x * x }, 0, 1, DefaultAbsTol, DefaultRelTol), 1.0/3.0)
	oracle.CheckClose(t, "sin", AdaptiveRecursive(math.Sin, 0, math.Pi, DefaultAbsTol, DefaultRelTol), 2.0)
}

// TestAdaptiveRecursiveEndpointSingularity integrates 1/sqrt(x) over (0,1] (= 2) and
// -log(x) over (0,1] (= 1): integrable endpoint singularities that a fixed-depth recursion
// or a left-to-right sweep cannot resolve but the global adaptive scheme can.
func TestAdaptiveRecursiveEndpointSingularity(t *testing.T) {
	invSqrt := AdaptiveRecursive(func(x float64) float64 { return 1 / math.Sqrt(x) }, 0, 1, 1e-9, 1e-9)
	if !oracle.IsClose(invSqrt, 2.0, 1e-4, 1e-4) { // weak singularity; 1e-4 is ample
		t.Errorf("∫1/sqrt(x) = %.10f, want ~2", invSqrt)
	}
	negLog := AdaptiveRecursive(func(x float64) float64 { return -math.Log(x) }, 0, 1, 1e-9, 1e-9)
	if !oracle.IsClose(negLog, 1.0, 1e-6, 1e-6) {
		t.Errorf("∫-log(x) = %.10f, want ~1", negLog)
	}
}

// TestIntegrateWithSingularitiesInterior splits at an interior singularity: ∫_{-1}^{1}
// 1/sqrt(|x|) dx = 4, with a sign-free even singularity at 0. This is the configuration the
// radial self-terms use (split at 0) and the case where a depth-limited recursion blows up.
func TestIntegrateWithSingularitiesInterior(t *testing.T) {
	f := func(x float64) float64 { return 1 / math.Sqrt(math.Abs(x)) }
	got := IntegrateWithSingularities(f, -1, 1, []float64{0}, 1e-9, 1e-9)
	if !oracle.IsClose(got, 4.0, 1e-4, 1e-4) {
		t.Errorf("∫1/sqrt(|x|) over [-1,1] = %.10f, want ~4", got)
	}
}
