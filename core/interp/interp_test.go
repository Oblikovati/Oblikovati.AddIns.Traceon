// SPDX-License-Identifier: MPL-2.0

package interp

import (
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
)

type interpGolden struct {
	Z       []oracle.F   `json:"z"`
	Y       []oracle.F   `json:"y"`
	Dy      []oracle.F   `json:"dy"`
	D2y     []oracle.F   `json:"d2y"`
	Cubic   [][]oracle.F `json:"cubic"`   // (n-1, 4)
	Quintic [][]oracle.F `json:"quintic"` // (n-1, 6)
}

func floats(fs []oracle.F) []float64 {
	out := make([]float64, len(fs))
	for i, f := range fs {
		out[i] = f.Float()
	}
	return out
}

// TestNotAKnotCubic verifies the cubic spline coefficients match scipy.CubicSpline exactly
// (the default not-a-knot boundary), interval by interval and coefficient by coefficient.
func TestNotAKnotCubic(t *testing.T) {
	var fx interpGolden
	oracle.LoadGolden(t, "interp", &fx)
	z, y := floats(fx.Z), floats(fx.Y)
	got, err := NotAKnotCubic(z, y)
	if err != nil {
		t.Fatalf("NotAKnotCubic: %v", err)
	}
	if len(got) != len(fx.Cubic) {
		t.Fatalf("got %d intervals, want %d", len(got), len(fx.Cubic))
	}
	for i := range got {
		for c := 0; c < 4; c++ {
			oracle.CheckClose(t, "cubic", got[i][c], fx.Cubic[i][c].Float())
		}
	}
}

// TestQuinticHermite verifies the quintic Hermite coefficients match scipy's
// BPoly.from_derivatives → PPoly construction exactly.
func TestQuinticHermite(t *testing.T) {
	var fx interpGolden
	oracle.LoadGolden(t, "interp", &fx)
	z, y, dy, d2y := floats(fx.Z), floats(fx.Y), floats(fx.Dy), floats(fx.D2y)
	got, err := QuinticHermite(z, y, dy, d2y)
	if err != nil {
		t.Fatalf("QuinticHermite: %v", err)
	}
	if len(got) != len(fx.Quintic) {
		t.Fatalf("got %d intervals, want %d", len(got), len(fx.Quintic))
	}
	for i := range got {
		for c := 0; c < 6; c++ {
			oracle.CheckClose(t, "quintic", got[i][c], fx.Quintic[i][c].Float())
		}
	}
}

// TestErrors checks the guards for too-few points and length mismatches.
func TestErrors(t *testing.T) {
	if _, err := NotAKnotCubic([]float64{0, 1}, []float64{0, 1}); err == nil {
		t.Error("NotAKnotCubic: expected error for <3 points")
	}
	if _, err := QuinticHermite([]float64{0}, []float64{0}, []float64{0}, []float64{0}); err == nil {
		t.Error("QuinticHermite: expected error for <2 points")
	}
}
