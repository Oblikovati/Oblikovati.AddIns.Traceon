// SPDX-License-Identifier: MPL-2.0

package validation

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/constants"
	"oblikovati.org/traceon/core/excitation"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/geometry"
	"oblikovati.org/traceon/core/solver"
	"oblikovati.org/traceon/core/tracing"
)

// Spherical deflection analyzer, from validation/spherical_capacitor.py (Cubric, Lencova,
// Read, Zlamal 1999, first benchmark). An electron launched at 1 eV between two concentric
// hemispherical shells traces an arc and re-crosses the axis at an exactly-known position.
const (
	scInnerRadius = 7.5
	scOuterRadius = 12.5
	scAngle       = 0.05 // launch angle (rad) from the r-axis
	scMSF         = 10
)

// scHemisphere builds one concentric shell as two quarter arcs (−z pole → +r equator → +z
// pole), the radial-symmetric meridian of a hemisphere. Port of spherical_capacitor.add_shell.
func scHemisphere(radius float64, name string) geometry.Path {
	return geometry.Arc(geometry.Point{0, 0, 0}, geometry.Point{0, 0, -radius}, geometry.Point{radius, 0, 0}, false).
		ExtendWithArc(geometry.Point{0, 0, 0}, geometry.Point{0, 0, radius}, false).
		WithName(name)
}

// scAnalyticCrossing is the exact axis crossing −10/(2/cos²θ − 1) from the benchmark paper.
func scAnalyticCrossing() float64 { return -10.0 / (2.0/math.Pow(math.Cos(scAngle), 2) - 1.0) }

// TestSphericalCapacitor reproduces the benchmark: solve the two-shell capacitor, trace a 1 eV
// electron through it, and compare its axis re-crossing to the exact analytic value. This
// exercises the arc mesher + electrostatic solve + RKF45 tracer + axis intersection together.
func TestSphericalCapacitor(t *testing.T) {
	m := geometry.MeshGroup(
		[]geometry.Path{scHemisphere(scInnerRadius, "inner"), scHemisphere(scOuterRadius, "outer")},
		geometry.MeshOptions{MeshSizeFactor: scMSF, HigherOrder: true})

	exc := excitation.New(m)
	exc.AddVoltage("inner", 5.0/3.0)
	exc.AddVoltage("outer", 3.0/5.0)
	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve electrostatic: %v", err)
	}
	bem := field.NewFieldRadialBEM(charges)

	// Electrostatic-only Lorentz force: a = q/m · E (no magnetic field).
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := bem.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	v0 := tracing.VelocityVec(1.0, geom3d.Vec3{math.Cos(scAngle), 0, -math.Sin(scAngle)}, constants.ElectronMass)
	qOverM := -constants.ElementaryCharge / constants.ElectronMass
	bounds := tracing.Bounds{{-0.1, 12.5}, {-0.1, 0.1}, {-12.5, 12.5}}

	_, states := tracing.TraceParticle(geom3d.Vec3{0, 0, 10}, v0, qOverM, fieldFn, bounds, 1e-10)
	z, ok := tracing.AxisIntersection(states)
	if !ok {
		t.Fatal("traced electron never re-crossed the axis")
	}

	want := scAnalyticCrossing()
	if rel := math.Abs(z-want) / math.Abs(want); rel > 1e-5 {
		t.Errorf("axis crossing = %.12g, want %.12g (analytic); rel err %.2e > 1e-5", z, want, rel)
	}
}
