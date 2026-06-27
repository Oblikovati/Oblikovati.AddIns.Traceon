// SPDX-License-Identifier: MPL-2.0

// Package validation holds end-to-end integration tests that run the whole pure-Go stack
// (solver → field → tracing → focus) on a complete electron-optics problem and assert it
// reproduces upstream Traceon's result. Unlike the per-module oracle tests, these validate
// the composition: a three-electrode einzel lens, solved, traced, and focused, must land on
// the same focal point Traceon computes for the identical geometry.
package validation

import (
	"encoding/json"
	"math"
	"os"
	"testing"

	"oblikovati.org/traceon/core/constants"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/focus"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/solver"
	"oblikovati.org/traceon/core/tracing"
)

type einzelGolden struct {
	Lines            [][][]float64 `json:"lines"`
	Types            []int         `json:"types"`
	Values           []float64     `json:"values"`
	EnergyEV         float64       `json:"energy_eV"`
	LaunchZ          float64       `json:"launch_z"`
	Radii            []float64     `json:"radii"`
	ChargeOverMass   float64       `json:"charge_over_mass"`
	Bounds           [][]float64   `json:"bounds"`
	Focus            []float64     `json:"focus"`
	SampleRadius     float64       `json:"sample_radius"`
	SampleTrajectory [][]float64   `json:"sample_trajectory"`
}

func loadEinzel(t *testing.T) einzelGolden {
	t.Helper()
	raw, err := os.ReadFile("testdata/einzel.golden.json")
	if err != nil {
		t.Fatalf("read golden: %v (regenerate: tools/gen_einzel.py)", err)
	}
	var g einzelGolden
	if err := json.Unmarshal(raw, &g); err != nil {
		t.Fatalf("decode golden: %v", err)
	}
	return g
}

func isClose(a, b, rtol, atol float64) bool {
	return math.Abs(a-b) <= atol+rtol*math.Abs(b)
}

// runEinzel solves the lens, traces every ray, and returns the trajectories.
func runEinzel(t *testing.T, g einzelGolden) [][]tracing.State {
	t.Helper()
	lines := make([]radial.Line, len(g.Lines))
	for i := range g.Lines {
		for v := 0; v < 4; v++ {
			lines[i][v] = geom2d.Vertex{g.Lines[i][v][0], g.Lines[i][v][1], g.Lines[i][v][2]}
		}
	}
	types := make([]radial.ExcitationType, len(g.Types))
	for i, ti := range g.Types {
		types[i] = radial.ExcitationType(ti)
	}

	charges, err := solver.SolveElectrostatic(lines, types, g.Values)
	if err != nil {
		t.Fatalf("solve: %v", err)
	}
	bem := field.NewFieldRadialBEM(charges)
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := bem.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	bounds := tracing.Bounds{
		{g.Bounds[0][0], g.Bounds[0][1]},
		{g.Bounds[1][0], g.Bounds[1][1]},
		{g.Bounds[2][0], g.Bounds[2][1]},
	}
	v0 := tracing.VelocityVec(g.EnergyEV, geom3d.Vec3{0, 0, 1}, constants.ElectronMass)

	var trajectories [][]tracing.State
	for _, r0 := range g.Radii {
		_, states := tracing.TraceParticle(geom3d.Vec3{r0, 0, g.LaunchZ}, v0, g.ChargeOverMass, fieldFn, bounds, 1e-8)
		trajectories = append(trajectories, states)
	}
	return trajectories
}

// TestEinzelLensFocus is the headline end-to-end check: the focal point of a three-electrode
// einzel lens, computed by the full Go stack, must match Traceon's focal point for the same
// geometry. The focus is a sensitive derived quantity, so agreement here exercises the entire
// chain (charge solve, field reconstruction, RKF45 tracing, least-squares focus) at once.
func TestEinzelLensFocus(t *testing.T) {
	g := loadEinzel(t)
	trajectories := runEinzel(t, g)

	got, err := focus.FocusPosition(trajectories)
	if err != nil {
		t.Fatalf("focus: %v", err)
	}
	// The focal z is the headline number; require tight agreement. x,y are ~0 (on-axis).
	if !isClose(got[2], g.Focus[2], 1e-6, 1e-6) {
		t.Errorf("focal z: got %.10g, want %.10g (Traceon)", got[2], g.Focus[2])
	}
	if !isClose(got[0], g.Focus[0], 1e-4, 1e-6) {
		t.Errorf("focal x: got %.10g, want %.10g", got[0], g.Focus[0])
	}
	t.Logf("einzel focus z = %.6f (Traceon %.6f)", got[2], g.Focus[2])
}

// TestEinzelFastTrace checks the fast axial-series interpolation reproduces the einzel focal
// length the direct boundary-integral trace gives — the speed feature must not change the
// physics. The lens is solved once, then traced both ways through the same field.
func TestEinzelFastTrace(t *testing.T) {
	g := loadEinzel(t)

	// Rebuild the solved field (mirrors runEinzel) so we can swap the evaluator.
	lines := make([]radial.Line, len(g.Lines))
	for i := range g.Lines {
		for v := 0; v < 4; v++ {
			lines[i][v] = geom2d.Vertex{g.Lines[i][v][0], g.Lines[i][v][1], g.Lines[i][v][2]}
		}
	}
	types := make([]radial.ExcitationType, len(g.Types))
	for i, ti := range g.Types {
		types[i] = radial.ExcitationType(ti)
	}
	charges, err := solver.SolveElectrostatic(lines, types, g.Values)
	if err != nil {
		t.Fatalf("solve: %v", err)
	}
	bem := field.NewFieldRadialBEM(charges)
	axial, err := field.NewFieldRadialAxial(charges, g.LaunchZ-1, g.Focus[2]+2, 400)
	if err != nil {
		t.Fatalf("axial: %v", err)
	}

	directFocus := traceFocus(t, g, bem.FieldAtPoint)
	fastFocus := traceFocus(t, g, axial.FieldAtPoint)
	if !isClose(fastFocus, directFocus, 1e-3, 1e-3) {
		t.Errorf("fast-trace focus z = %.6f, direct = %.6f (should match)", fastFocus, directFocus)
	}
	t.Logf("einzel focus: direct %.6f, fast-axial %.6f", directFocus, fastFocus)
}

// traceFocus traces the einzel ray bundle through the given electrostatic field evaluator and
// returns the focal z.
func traceFocus(t *testing.T, g einzelGolden, eval func(geom2d.Vertex) geom2d.Vertex) float64 {
	t.Helper()
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := eval(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	bounds := tracing.Bounds{{g.Bounds[0][0], g.Bounds[0][1]}, {g.Bounds[1][0], g.Bounds[1][1]}, {g.Bounds[2][0], g.Bounds[2][1]}}
	v0 := tracing.VelocityVec(g.EnergyEV, geom3d.Vec3{0, 0, 1}, constants.ElectronMass)
	var trajectories [][]tracing.State
	for _, r0 := range g.Radii {
		_, states := tracing.TraceParticle(geom3d.Vec3{r0, 0, g.LaunchZ}, v0, g.ChargeOverMass, fieldFn, bounds, 1e-8)
		trajectories = append(trajectories, states)
	}
	f, err := focus.FocusPosition(trajectories)
	if err != nil {
		t.Fatalf("focus: %v", err)
	}
	return f[2]
}

// TestEinzelTrajectory checks the outermost ray's full trajectory matches Traceon step for
// step (tight early; slightly looser late as -ffast-math rounding accumulates over ~1200 steps).
func TestEinzelTrajectory(t *testing.T) {
	g := loadEinzel(t)
	// Trace just the sample ray.
	gg := g
	gg.Radii = []float64{g.SampleRadius}
	traj := runEinzel(t, gg)[0]

	if len(traj) != len(g.SampleTrajectory) {
		t.Fatalf("trajectory length: got %d, want %d", len(traj), len(g.SampleTrajectory))
	}
	for i := range traj {
		rtol, atol := 1e-9, 1e-9
		if i > 100 {
			rtol, atol = 1e-6, 1e-7
		}
		for c := 0; c < 6; c++ {
			if !isClose(traj[i][c], g.SampleTrajectory[i][c], rtol, atol*math.Max(1, math.Abs(g.SampleTrajectory[i][c]))) {
				t.Fatalf("step %d comp %d: got %.12g, want %.12g", i, c, traj[i][c], g.SampleTrajectory[i][c])
			}
		}
	}
}
