// SPDX-License-Identifier: MPL-2.0

package elliptic

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
)

// ellipticCase mirrors one row of core/elliptic/testdata/elliptic.golden.json, produced
// by tools/gen_fixtures.py from the upstream Traceon backend.
type ellipticCase struct {
	M        oracle.F `json:"m"`
	Ellipk   oracle.F `json:"ellipk"`
	Ellipe   oracle.F `json:"ellipe"`
	Ellipkm1 oracle.F `json:"ellipkm1"`
	Ellipem1 oracle.F `json:"ellipem1"`
}

// TestAgainstGolden verifies every elliptic function reproduces the upstream Traceon
// backend value bit-for-bit-to-tolerance across the sampled parameter sweep. This is
// the direct port-equivalence guard (Go vs the exact C the Python wraps).
func TestAgainstGolden(t *testing.T) {
	var fx struct {
		Cases  []ellipticCase `json:"cases"`
		KEOnly []ellipticCase `json:"k_e_only"`
	}
	oracle.LoadGolden(t, "elliptic", &fx)
	if len(fx.Cases) == 0 {
		t.Fatal("no golden cases loaded")
	}
	for _, c := range fx.Cases {
		// The golden file stores ellipkm1/ellipem1 evaluated at p = 1-m (see gen_fixtures).
		m := c.M.Float()
		p := 1 - m
		oracle.CheckClose(t, "Ellipk(m)", Ellipk(m), c.Ellipk.Float())
		oracle.CheckClose(t, "Ellipe(m)", Ellipe(m), c.Ellipe.Float())
		oracle.CheckClose(t, "Ellipkm1(1-m)", Ellipkm1(p), c.Ellipkm1.Float())
		oracle.CheckClose(t, "Ellipem1(1-m)", Ellipem1(p), c.Ellipem1.Float())
	}
	// Reciprocal-modulus branch (m outside [0,1]): the C uses imaginary-modulus transforms
	// whose sqrt(negative) yields NaN for m<=-1 (Ellipk) and m>1 (Ellipe) — these arguments
	// never arise in radial BEM. The faithful port must reproduce the SAME NaN, so assert
	// NaN-equivalence where the oracle is NaN and closeness otherwise.
	for _, c := range fx.KEOnly {
		checkEquivOrNaN(t, "Ellipk(m<0|m>1)", Ellipk(c.M.Float()), c.Ellipk.Float())
		checkEquivOrNaN(t, "Ellipe(m<0|m>1)", Ellipe(c.M.Float()), c.Ellipe.Float())
	}
}

// checkEquivOrNaN asserts got == NaN exactly when want is NaN, and got ≈ want otherwise.
func checkEquivOrNaN(t *testing.T, label string, got, want float64) {
	t.Helper()
	if math.IsNaN(want) {
		if !math.IsNaN(got) {
			t.Errorf("%s: got %v, want NaN (matches upstream out-of-domain)", label, got)
		}
		return
	}
	oracle.CheckClose(t, label, got, want)
}

// TestKnownValues pins the functions to textbook closed forms independent of the
// fixture, so a regenerated-but-wrong fixture cannot mask a real regression.
//
//	K(0) = E(0) = pi/2.
func TestKnownValues(t *testing.T) {
	oracle.CheckClose(t, "Ellipk(0)", Ellipk(0), math.Pi/2)
	oracle.CheckClose(t, "Ellipe(0)", Ellipe(0), math.Pi/2)
}

// TestDomainEdgeMatchesUpstream documents that the parameter m=1 is OUTSIDE the Cody
// approximation's supported domain: the upstream backend returns NaN there too (the
// series evaluates log(1/0)*0). The faithful port must reproduce that, not paper over
// it with the analytic E(1)=1 — Traceon never evaluates exactly at m=1.
func TestDomainEdgeMatchesUpstream(t *testing.T) {
	if !math.IsNaN(Ellipe(1)) {
		t.Errorf("Ellipe(1) = %v, want NaN (matches upstream backend at the domain edge)", Ellipe(1))
	}
}
