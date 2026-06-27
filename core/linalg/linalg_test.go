// SPDX-License-Identifier: MPL-2.0

package linalg

import (
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
)

// TestSolveVector solves a small system with a known answer: [[2,1],[1,3]]·x = [3,5] → x=[0.8,1.4].
func TestSolveVector(t *testing.T) {
	a := Matrix{Rows: 2, Cols: 2, Data: []float64{2, 1, 1, 3}}
	x, err := SolveVector(a, []float64{3, 5})
	if err != nil {
		t.Fatalf("SolveVector: %v", err)
	}
	oracle.CheckClose(t, "x0", x[0], 0.8)
	oracle.CheckClose(t, "x1", x[1], 1.4)
}

// TestSolveMultiRHS solves A·X = B with two RHS columns and checks A·X reproduces B.
func TestSolveMultiRHS(t *testing.T) {
	a := Matrix{Rows: 2, Cols: 2, Data: []float64{4, 3, 6, 3}}
	b := Matrix{Rows: 2, Cols: 2, Data: []float64{1, 0, 0, 1}} // identity → X = A^{-1}
	x, err := Solve(a, b)
	if err != nil {
		t.Fatalf("Solve: %v", err)
	}
	// Verify A·X == I.
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			ax := a.At(i, 0)*x.At(0, j) + a.At(i, 1)*x.At(1, j)
			want := 0.0
			if i == j {
				want = 1.0
			}
			oracle.CheckClose(t, "A·X", ax, want)
		}
	}
}

// TestLeastSquares fits y = a + b·t to exactly-linear data and recovers the line.
func TestLeastSquares(t *testing.T) {
	// Rows: design matrix [1, t]; data y = 2 + 3t at t = 0,1,2,3.
	a := Matrix{Rows: 4, Cols: 2, Data: []float64{1, 0, 1, 1, 1, 2, 1, 3}}
	b := []float64{2, 5, 8, 11}
	x, err := LeastSquares(a, b)
	if err != nil {
		t.Fatalf("LeastSquares: %v", err)
	}
	oracle.CheckClose(t, "intercept", x[0], 2.0)
	oracle.CheckClose(t, "slope", x[1], 3.0)
}

// TestSingularReturnsError confirms a singular system errors rather than returning garbage.
func TestSingularReturnsError(t *testing.T) {
	a := Matrix{Rows: 2, Cols: 2, Data: []float64{1, 2, 2, 4}} // rank 1
	if _, err := SolveVector(a, []float64{1, 1}); err == nil {
		t.Error("expected error for singular matrix, got nil")
	}
}
