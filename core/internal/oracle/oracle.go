// SPDX-License-Identifier: MPL-2.0

// Package oracle holds shared test support for verifying the pure-Go numerical core
// against the upstream Traceon (the oracle): loading golden-fixture JSON and the
// np.isclose-equivalent tolerance check the fixtures are compared with.
//
// It is imported only from _test.go files, but lives in a normal package so every
// core module reuses one definition of "numerically equal to Python".
package oracle

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"strconv"
	"testing"
)

// F is a float64 that decodes from a JSON number OR the non-finite string tokens the
// fixture generator emits ("NaN", "Infinity", "-Infinity"), since standard JSON has no
// way to spell those and radial BEM kernels produce them at singularities.
type F float64

// Float returns the underlying float64.
func (f F) Float() float64 { return float64(f) }

// UnmarshalJSON accepts a JSON number or one of the non-finite string tokens.
func (f *F) UnmarshalJSON(b []byte) error {
	b = bytes.TrimSpace(b)
	if len(b) > 0 && b[0] == '"' {
		switch string(b) {
		case `"NaN"`:
			*f = F(math.NaN())
		case `"Infinity"`:
			*f = F(math.Inf(1))
		case `"-Infinity"`:
			*f = F(math.Inf(-1))
		default:
			return fmt.Errorf("oracle.F: unrecognized non-finite token %s", b)
		}
		return nil
	}
	v, err := strconv.ParseFloat(string(b), 64)
	if err != nil {
		return err
	}
	*f = F(v)
	return nil
}

// Default tolerances mirror numpy.isclose: |a-b| <= Atol + Rtol*|b|.
const (
	DefaultRtol = 1e-5
	DefaultAtol = 1e-8
)

// IsClose reports whether a is within (atol + rtol*|b|) of b, matching numpy.isclose
// (asymmetric in b — b is the reference/oracle value). NaNs are never close; equal
// infinities are close.
func IsClose(a, b, rtol, atol float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	return math.Abs(a-b) <= atol+rtol*math.Abs(b)
}

// Close reports IsClose with the default numpy tolerances.
func Close(a, b float64) bool { return IsClose(a, b, DefaultRtol, DefaultAtol) }

// LoadGolden reads testdata/<name>.golden.json into dest, failing the test if the
// fixture is missing (run `make verify-oracle` to regenerate from upstream Traceon).
func LoadGolden(t *testing.T, name string, dest any) {
	t.Helper()
	path := "testdata/" + name + ".golden.json"
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read golden %s: %v (regenerate with `make verify-oracle`)", path, err)
	}
	if err := json.Unmarshal(raw, dest); err != nil {
		t.Fatalf("decode golden %s: %v", path, err)
	}
}

// CheckClose fails the test (without aborting) when got is not Close to want, reporting
// the case label and the relative error so a drift is diagnosable at a glance.
func CheckClose(t *testing.T, label string, got, want float64) {
	t.Helper()
	if !Close(got, want) {
		rel := math.Abs(got-want) / math.Max(math.Abs(want), 1e-300)
		t.Errorf("%s: got %.17g, want %.17g (rel err %.2e)", label, got, want, rel)
	}
}
