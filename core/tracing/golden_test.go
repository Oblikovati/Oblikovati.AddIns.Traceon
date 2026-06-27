// SPDX-License-Identifier: MPL-2.0

package tracing_test

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/radial"
	"oblikovati.org/traceon/core/tracing"
)

type tracingGolden struct {
	Lines          [][][]oracle.F `json:"lines"`
	Charges        []oracle.F     `json:"charges"`
	P0             []oracle.F     `json:"p0"`
	EnergyEV       oracle.F       `json:"energy_eV"`
	Atol           oracle.F       `json:"atol"`
	ChargeOverMass oracle.F       `json:"charge_over_mass"`
	Bounds         [][]oracle.F   `json:"bounds"`
	Times          []oracle.F     `json:"times"`
	Positions      [][]oracle.F   `json:"positions"` // [N][6]
}

// TestTraceThroughRealField traces an electron through a solved radial BEM field and checks
// the trajectory reproduces the upstream Tracer step for step. This verifies the tracer and
// the field evaluation compose correctly over hundreds of adaptive steps. The first steps
// must match tightly; tiny floating-point differences (the C is built with -ffast-math) may
// accumulate late, so later steps are checked with a looser-but-still-strict tolerance.
func TestTraceThroughRealField(t *testing.T) {
	var fx tracingGolden
	oracle.LoadGolden(t, "tracing", &fx)

	lines := make([]radial.Line, len(fx.Lines))
	for i := range fx.Lines {
		for v := 0; v < 4; v++ {
			lines[i][v] = geom2d.Vertex{fx.Lines[i][v][0].Float(), fx.Lines[i][v][1].Float(), fx.Lines[i][v][2].Float()}
		}
	}
	charges := make([]float64, len(fx.Charges))
	for i := range charges {
		charges[i] = fx.Charges[i].Float()
	}
	jac, pos := radial.FillJacobianBufferRadial(lines)

	// Electrostatic FieldRadialBEM closure: elec = field_radial, mag = 0.
	field := func(p, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := radial.FieldRadial(geom2d.Vertex{p[0], p[1], p[2]}, charges, jac, pos)
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}

	p0 := geom3d.Vec3{fx.P0[0].Float(), fx.P0[1].Float(), fx.P0[2].Float()}
	const mE = 9.1093837139e-31
	v0 := tracing.VelocityVec(fx.EnergyEV.Float(), geom3d.Vec3{0, 0, 1}, mE)
	bounds := tracing.Bounds{
		{fx.Bounds[0][0].Float(), fx.Bounds[0][1].Float()},
		{fx.Bounds[1][0].Float(), fx.Bounds[1][1].Float()},
		{fx.Bounds[2][0].Float(), fx.Bounds[2][1].Float()},
	}
	times, states := tracing.TraceParticle(p0, v0, fx.ChargeOverMass.Float(), field, bounds, fx.Atol.Float())

	if len(states) != len(fx.Positions) {
		t.Fatalf("step count: got %d, want %d", len(states), len(fx.Positions))
	}
	for i := range states {
		// Tight early, slightly looser late (accumulated -ffast-math divergence).
		rtol, atol := 1e-9, 1e-9
		if i > 50 {
			rtol, atol = 1e-6, 1e-7
		}
		if !oracle.IsClose(times[i], fx.Times[i].Float(), rtol, atol) {
			t.Errorf("step %d time: got %.12g, want %.12g", i, times[i], fx.Times[i].Float())
		}
		for c := 0; c < 6; c++ {
			got, want := states[i][c], fx.Positions[i][c].Float()
			// Scale the position tolerance by the velocity magnitude for the velocity slots.
			if !oracle.IsClose(got, want, rtol, atol*math.Max(1, math.Abs(want))) {
				t.Errorf("step %d comp %d: got %.12g, want %.12g", i, c, got, want)
				break
			}
		}
	}
}
