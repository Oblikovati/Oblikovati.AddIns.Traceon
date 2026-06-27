// SPDX-License-Identifier: MPL-2.0

// Package linalg is the project-owned thin wrapper over the dense linear algebra the
// Traceon core needs: solving the BEM influence system (numpy.linalg.solve in solver.py)
// and least-squares for best-focus fitting (numpy.linalg.lstsq in focus.py). It hides the
// third-party backend (gonum) behind a small Matrix type and two functions, so the rest of
// core depends only on this package, per the project's dependency-wrapping rule.
package linalg

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

// Matrix is a dense row-major matrix owned by this project (no gonum types leak out).
type Matrix struct {
	Rows, Cols int
	Data       []float64 // len == Rows*Cols, row-major
}

// NewMatrix allocates a zeroed Rows×Cols matrix.
func NewMatrix(rows, cols int) Matrix {
	return Matrix{Rows: rows, Cols: cols, Data: make([]float64, rows*cols)}
}

// At returns element (i, j).
func (m Matrix) At(i, j int) float64 { return m.Data[i*m.Cols+j] }

// Set assigns element (i, j).
func (m Matrix) Set(i, j int, v float64) { m.Data[i*m.Cols+j] = v }

func (m Matrix) dense() *mat.Dense {
	// gonum copies/uses the slice as row-major backing — matches our layout exactly.
	return mat.NewDense(m.Rows, m.Cols, append([]float64(nil), m.Data...))
}

// Solve returns X solving A·X = B for a square A. B may have multiple right-hand-side
// columns (as numpy.linalg.solve does); X has B's shape. Equivalent to numpy.linalg.solve.
func Solve(a, b Matrix) (Matrix, error) {
	if a.Rows != a.Cols {
		return Matrix{}, fmt.Errorf("linalg.Solve: A must be square, got %dx%d", a.Rows, a.Cols)
	}
	if b.Rows != a.Rows {
		return Matrix{}, fmt.Errorf("linalg.Solve: B rows %d != A rows %d", b.Rows, a.Rows)
	}
	var x mat.Dense
	if err := x.Solve(a.dense(), b.dense()); err != nil {
		return Matrix{}, fmt.Errorf("linalg.Solve: %w (matrix may be singular)", err)
	}
	return fromDense(&x), nil
}

// SolveVector returns x solving A·x = b for a square A and a single right-hand side.
func SolveVector(a Matrix, b []float64) ([]float64, error) {
	if len(b) != a.Rows {
		return nil, fmt.Errorf("linalg.SolveVector: b len %d != A rows %d", len(b), a.Rows)
	}
	bm := Matrix{Rows: len(b), Cols: 1, Data: append([]float64(nil), b...)}
	x, err := Solve(a, bm)
	if err != nil {
		return nil, err
	}
	return x.Data, nil
}

// LeastSquares returns x minimizing ||A·x - b||₂ for an overdetermined A (Rows ≥ Cols),
// equivalent to numpy.linalg.lstsq's solution vector.
func LeastSquares(a Matrix, b []float64) ([]float64, error) {
	if a.Rows < a.Cols {
		return nil, fmt.Errorf("linalg.LeastSquares: need Rows>=Cols, got %dx%d", a.Rows, a.Cols)
	}
	if len(b) != a.Rows {
		return nil, fmt.Errorf("linalg.LeastSquares: b len %d != A rows %d", len(b), a.Rows)
	}
	bm := mat.NewDense(len(b), 1, append([]float64(nil), b...))
	var x mat.Dense
	if err := x.Solve(a.dense(), bm); err != nil {
		return nil, fmt.Errorf("linalg.LeastSquares: %w", err)
	}
	out := make([]float64, a.Cols)
	for i := 0; i < a.Cols; i++ {
		out[i] = x.At(i, 0)
	}
	return out, nil
}

func fromDense(d *mat.Dense) Matrix {
	r, c := d.Dims()
	m := NewMatrix(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			m.Set(i, j, d.At(i, j))
		}
	}
	return m
}
