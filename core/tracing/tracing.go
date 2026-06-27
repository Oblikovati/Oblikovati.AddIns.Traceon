// SPDX-License-Identifier: MPL-2.0

// Package tracing integrates charged-particle trajectories through a field, porting
// traceon/backend/tracing.c (the adaptive Runge-Kutta-Fehlberg 4(5) integrator and the
// Lorentz force) and the trajectory helpers in traceon/tracing.py (initial-velocity
// construction and trajectory/plane intersections).
//
// The state is the 6-vector [x, y, z, vx, vy, vz]; the acceleration is
// q/m · (E + μ₀ · v × H). Note the integrator uses the BACKEND's hardcoded μ₀
// (definitions.c), which differs from scipy.constants.mu_0 in its last digits — matching
// it is required for trajectory equivalence.
package tracing

import (
	"math"

	"oblikovati.org/traceon/core/geom3d"
)

// backendMU0 is the vacuum permeability the C backend hardcodes (definitions.c MU_0).
// It is NOT identical to scipy.constants.mu_0 (1.25663706127e-06); the tracer must use
// THIS value to reproduce upstream trajectories exactly.
const backendMU0 = 1.25663706212e-06

// tracingStepMax bounds the step length (TRACING_STEP_MAX in tracing.c): the actual max
// time step is this divided by the speed.
const tracingStepMax = 0.01

// TracingBlockSize caps the number of accepted steps recorded in one call (TRACING_BLOCK_SIZE).
const TracingBlockSize = 100000

// RKF45 Butcher tableau (Runge-Kutta-Fehlberg), verbatim from tracing.c.
var (
	rkB2 = []float64{2. / 9.}
	rkB3 = []float64{1. / 12., 1. / 4.}
	rkB4 = []float64{69. / 128., -243. / 128., 135. / 64.}
	rkB5 = []float64{-17. / 12., 27. / 4., -27. / 5., 16. / 15.}
	rkB6 = []float64{65. / 432., -5. / 16., 13. / 16., 4. / 27., 5. / 144.}
	// rkCoeffs[index] are the weights for building stage `index` (stage 0 has none).
	rkCoeffs = [6][]float64{nil, rkB2, rkB3, rkB4, rkB5, rkB6}
	// rkCH: 5th-order solution weights; rkCT: error (difference) weights.
	rkCH = []float64{47. / 450., 0., 12. / 25., 32. / 225., 1. / 30., 6. / 25.}
	rkCT = []float64{-1. / 150., 0., 3. / 100., -16. / 75., -1. / 20., 6. / 25.}
)

// State is the 6-component phase-space vector [x, y, z, vx, vy, vz].
type State = [6]float64

// Bounds is the tracing box ((xmin,xmax),(ymin,ymax),(zmin,zmax)); tracing stops when the
// particle leaves it.
type Bounds = [3][2]float64

// FieldFunc returns the electrostatic field E and the magnetic field H at the given
// position (the velocity is passed for generality, e.g. velocity-dependent fields). Both
// are returned in the same units the upstream field functions use.
type FieldFunc func(position, velocity geom3d.Vec3) (elec, mag geom3d.Vec3)

// produceNewY builds the intermediate stage state ys[index] from the base state y and the
// already-computed increments ks[0..index-1]. Port of produce_new_y.
func produceNewY(y State, ys *[6]State, ks *[6]State, index int) {
	coeffs := rkCoeffs[index]
	for i := 0; i < 6; i++ {
		ys[index][i] = y[i]
		for j := 0; j < index; j++ {
			ys[index][i] += coeffs[j] * ks[j][i]
		}
	}
}

// produceNewK computes the stage increment ks[index] = h · f(ys[index]), where f is the
// equations of motion: position derivative = velocity, velocity derivative =
// q/m · (E + μ₀ · v × H). Port of produce_new_k.
func produceNewK(ys *[6]State, ks *[6]State, index int, h, chargeOverMass float64, field FieldFunc) {
	pos := geom3d.Vec3{ys[index][0], ys[index][1], ys[index][2]}
	vel := geom3d.Vec3{ys[index][3], ys[index][4], ys[index][5]}
	elec, mag := field(pos, vel)
	cross := geom3d.CrossProduct3D(vel, mag) // v × H
	fx := elec[0] + backendMU0*cross[0]
	fy := elec[1] + backendMU0*cross[1]
	fz := elec[2] + backendMU0*cross[2]
	ks[index][0] = h * ys[index][3]
	ks[index][1] = h * ys[index][4]
	ks[index][2] = h * ys[index][5]
	ks[index][3] = h * chargeOverMass * fx
	ks[index][4] = h * chargeOverMass * fy
	ks[index][5] = h * chargeOverMass * fz
}

// TraceParticle integrates a particle from the initial position/velocity until it leaves
// bounds, returning the accepted (times, states). Adaptive RKF45 with the upstream's
// step-size control. Port of trace_particle. position and velocity are the initial 3-vectors;
// velocity must already be in SI (m/s) — use VelocityVec to build it from an energy in eV.
func TraceParticle(position, velocity geom3d.Vec3, chargeOverMass float64, field FieldFunc, bounds Bounds, atol float64) ([]float64, []State) {
	var y State
	for i := 0; i < 3; i++ {
		y[i] = position[i]
		y[i+3] = velocity[i]
	}
	speed := geom3d.Norm3D(y[3], y[4], y[5])
	hmax := tracingStepMax / speed
	h := hmax

	times := []float64{0}
	states := []State{y}

	inBounds := func(s State) bool {
		return bounds[0][0] <= s[0] && s[0] <= bounds[0][1] &&
			bounds[1][0] <= s[1] && s[1] <= bounds[1][1] &&
			bounds[2][0] <= s[2] && s[2] <= bounds[2][1]
	}

	for inBounds(y) {
		var ks, ys [6]State
		for index := 0; index < 6; index++ {
			produceNewY(y, &ys, &ks, index)
			produceNewK(&ys, &ks, index, h, chargeOverMass, field)
		}

		maxPosErr, maxVelErr := 0.0, 0.0
		for i := 0; i < 3; i++ {
			e := 0.0
			for j := 0; j < 6; j++ {
				e += rkCT[j] * ks[j][i]
			}
			if math.Abs(e) > maxPosErr {
				maxPosErr = math.Abs(e)
			}
		}
		for i := 3; i < 6; i++ {
			e := 0.0
			for j := 0; j < 6; j++ {
				e += rkCT[j] * ks[j][i]
			}
			if math.Abs(e) > maxVelErr {
				maxVelErr = math.Abs(e)
			}
		}
		errEst := maxPosErr + h*maxVelErr

		if errEst <= atol {
			for i := 0; i < 6; i++ {
				y[i] += rkCH[0]*ks[0][i] + rkCH[1]*ks[1][i] + rkCH[2]*ks[2][i] +
					rkCH[3]*ks[3][i] + rkCH[4]*ks[4][i] + rkCH[5]*ks[5][i]
			}
			times = append(times, times[len(times)-1]+h)
			states = append(states, y)
			if len(states) == TracingBlockSize {
				return times, states
			}
		}
		// Step-size update (verbatim from the C). When errEst is 0 — e.g. RKF45 is exact on
		// a constant field — atol/errEst is +Inf and math.Min(..., hmax) returns hmax, exactly
		// as the C's fmin does.
		h = math.Min(0.9*h*math.Pow(atol/errEst, 0.2), hmax)
	}
	return times, states
}
