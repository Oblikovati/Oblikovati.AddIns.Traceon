// SPDX-License-Identifier: MPL-2.0

package tracing

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/internal/oracle"
)

// TestConstantAcceleration ports test_tracing_constant_acceleration: under a constant
// acceleration of 3 along x (and an initial z-velocity of 3), the trajectory is the exact
// parabola x = 3/2·t², z = 3·t. RKF45 is exact on polynomials of its order, so this checks
// the integrator and the recorded times to machine precision.
func TestConstantAcceleration(t *testing.T) {
	const cOverM = -175882001077.2163 // EM = -e/m_e, the value the upstream test uses
	field := func(_, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		// acceleration = cOverM·elec = [3,0,0]  ⟹  elec = [3/cOverM, 0, 0].
		return geom3d.Vec3{3.0 / cOverM, 0, 0}, geom3d.Vec3{}
	}
	bounds := Bounds{{-2, 2}, {-2, 2}, {-2, math.Sqrt(12) + 1}}
	times, states := TraceParticle(geom3d.Vec3{0, 0, 0}, geom3d.Vec3{0, 0, 3}, cOverM, field, bounds, 1e-10)

	if len(states) < 2 {
		t.Fatalf("expected a multi-step trajectory, got %d states", len(states))
	}
	for i, s := range states {
		tt := times[i]
		oracle.CheckClose(t, "x = 3/2 t^2", s[0], 1.5*tt*tt)
		oracle.CheckClose(t, "z = 3 t", s[2], 3*tt)
	}
}

// TestVelocityVec checks the eV→m/s conversion: the resulting speed satisfies
// (1/2)·m·v² = eV·e, and the direction is preserved (normalized).
func TestVelocityVec(t *testing.T) {
	const mE = 9.1093837139e-31
	const eC = 1.602176634e-19
	v := VelocityVec(100, geom3d.Vec3{0, 0, 2}, mE) // direction not unit-length on purpose
	speed := geom3d.Norm3D(v[0], v[1], v[2])
	wantSpeed := math.Sqrt(2 * 100 * eC / mE)
	oracle.CheckClose(t, "speed", speed, wantSpeed)
	oracle.CheckClose(t, "vx", v[0], 0)
	oracle.CheckClose(t, "vy", v[1], 0)
	oracle.CheckClose(t, "vz", v[2], wantSpeed) // +z direction preserved
}

// TestZToBounds checks the three branches of the optical-span margin helper.
func TestZToBounds(t *testing.T) {
	lo, hi := ZToBounds(-3, -1)
	oracle.CheckClose(t, "both-neg lo", lo, -4)
	oracle.CheckClose(t, "both-neg hi", hi, 1)
	lo, hi = ZToBounds(1, 3)
	oracle.CheckClose(t, "both-pos lo", lo, -1)
	oracle.CheckClose(t, "both-pos hi", hi, 4)
	lo, hi = ZToBounds(-2, 2)
	oracle.CheckClose(t, "straddle lo", lo, -3)
	oracle.CheckClose(t, "straddle hi", hi, 3)
}

// TestPlaneIntersection checks a straight trajectory crossing the z=0 plane: a particle
// going from z=-1 to z=+1 at x=0.5 crosses at (0.5, 0, 0).
func TestPlaneIntersection(t *testing.T) {
	positions := []State{
		{0.5, 0, -1, 0, 0, 1},
		{0.5, 0, 1, 0, 0, 1},
	}
	s, ok := XYPlaneIntersection(positions, 0)
	if !ok {
		t.Fatal("expected an intersection")
	}
	oracle.CheckClose(t, "x", s[0], 0.5)
	oracle.CheckClose(t, "z", s[2], 0.0)

	z, ok := AxisIntersection([]State{{-1, 0, 2, 1, 0, 1}, {1, 0, 4, 1, 0, 1}})
	if !ok {
		t.Fatal("expected an axis crossing")
	}
	oracle.CheckClose(t, "axis z", z, 3.0) // crosses x=0 midway, z=3
}
