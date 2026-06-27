// SPDX-License-Identifier: MPL-2.0

// Package focus finds the focus of a bundle of charged-particle trajectories, porting
// traceon/focus.py. It linearly extends each trajectory's final position backward along
// its final velocity and solves, in a least-squares sense, for the point those lines come
// closest to — the beam focus.
package focus

import (
	"fmt"

	"oblikovati.org/traceon/core/linalg"
	"oblikovati.org/traceon/core/tracing"
)

// FocusPosition returns the (x, y, z) focus of the given trajectories (each a slice of
// 6-states as produced by the tracer). Port of focus_position: it forms the over-determined
// system from the final positions and the paraxial slopes vx/vz, vy/vz and solves it by
// least squares. Requires at least one trajectory, each with at least one state.
func FocusPosition(trajectories [][]tracing.State) ([3]float64, error) {
	n := len(trajectories)
	if n == 0 {
		return [3]float64{}, fmt.Errorf("focus.FocusPosition: no trajectories")
	}
	anglesX := make([]float64, n)
	anglesY := make([]float64, n)
	px := make([]float64, n)
	py := make([]float64, n)
	pz := make([]float64, n)
	for i, traj := range trajectories {
		if len(traj) == 0 {
			return [3]float64{}, fmt.Errorf("focus.FocusPosition: trajectory %d is empty", i)
		}
		last := traj[len(traj)-1]
		anglesX[i] = last[3] / last[5] // vx/vz
		anglesY[i] = last[4] / last[5] // vy/vz
		px[i], py[i], pz[i] = last[0], last[1], last[2]
	}

	// Design matrix A (2N x 3), columns [first, second, third]; unknowns [z, x, y].
	a := linalg.NewMatrix(2*n, 3)
	rhs := make([]float64, 2*n)
	for i := 0; i < n; i++ {
		// Rows 0..N-1 (x equations).
		a.Set(i, 0, -anglesX[i])
		a.Set(i, 1, 1)
		a.Set(i, 2, 0)
		rhs[i] = px[i] - pz[i]*anglesX[i]
		// Rows N..2N-1 (y equations).
		a.Set(n+i, 0, -anglesY[i])
		a.Set(n+i, 1, 0)
		a.Set(n+i, 2, 1)
		rhs[n+i] = py[i] - pz[i]*anglesY[i]
	}

	sol, err := linalg.LeastSquares(a, rhs)
	if err != nil {
		return [3]float64{}, err
	}
	// sol = [z, x, y]; return [x, y, z].
	return [3]float64{sol[1], sol[2], sol[0]}, nil
}
