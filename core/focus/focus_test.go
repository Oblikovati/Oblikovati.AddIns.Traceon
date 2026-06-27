// SPDX-License-Identifier: MPL-2.0

package focus

import (
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/tracing"
)

// TestFocusPosition constructs trajectories whose backward extensions all meet at a chosen
// focus, then checks FocusPosition recovers it exactly. For trajectory i with final z = pz_i
// and paraxial slopes (ax_i, ay_i), the final position consistent with focus (xf,yf,zf) is
// px_i = xf - (zf-pz_i)·ax_i, py_i = yf - (zf-pz_i)·ay_i.
func TestFocusPosition(t *testing.T) {
	xf, yf, zf := 0.2, -0.1, 5.0
	type spec struct{ pz, ax, ay float64 }
	specs := []spec{{0.0, 0.03, 0.01}, {0.5, -0.02, 0.04}, {-0.3, 0.05, -0.02}}

	var trajectories [][]tracing.State
	for _, s := range specs {
		px := xf - (zf-s.pz)*s.ax
		py := yf - (zf-s.pz)*s.ay
		vz := 1.0
		// A two-state trajectory; only the final state is used by FocusPosition.
		traj := []tracing.State{
			{px - s.ax, py - s.ay, s.pz - 1, s.ax * vz, s.ay * vz, vz},
			{px, py, s.pz, s.ax * vz, s.ay * vz, vz},
		}
		trajectories = append(trajectories, traj)
	}

	f, err := FocusPosition(trajectories)
	if err != nil {
		t.Fatalf("FocusPosition: %v", err)
	}
	oracle.CheckClose(t, "focus x", f[0], xf)
	oracle.CheckClose(t, "focus y", f[1], yf)
	oracle.CheckClose(t, "focus z", f[2], zf)
}

// TestFocusErrors checks the guards for no trajectories and an empty trajectory.
func TestFocusErrors(t *testing.T) {
	if _, err := FocusPosition(nil); err == nil {
		t.Error("expected error for no trajectories")
	}
	if _, err := FocusPosition([][]tracing.State{{}}); err == nil {
		t.Error("expected error for empty trajectory")
	}
}
