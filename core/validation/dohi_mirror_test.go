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

// Dohi electron mirror, from validation/dohi.py (Dohi & Kruit 2018). An electron launched at
// z=15 drifts down, reflects off the −1250 V mirror, and returns; the value of interest is its
// radius back at z=15 — ideally 0 (it should return to the axis). This validates the axial fast
// field through a reflection (a turning point) plus the field-free drift above the field region.
const (
	dohiThickness = 0.15
	dohiRadius    = 0.075
	dohiSpacer    = 0.5
	dohiExtent    = 1.0 - 0.1
	dohiMSF       = 100
	dohiLensVolts = 710.0126605741955 // tuned in dohi.py so the mirror returns the ray to axis
	dohiZ0        = 15.0
	dohiAngle     = 0.5e-3
)

// dohiUpstreamMSF100 is the return radius Traceon computes at MSF=100 (higher-order lines, axial
// field N=500) — the port-equivalence oracle. The ideal value is 0; this near-zero residual is
// the interpolation/trace error, which the Go port must reproduce.
const dohiUpstreamMSF100 = 3.151215467390624e-05

// TestDohiMirror reproduces the Dohi-mirror return radius end to end: solve the mirror/lens/
// ground stack with a boundary, build the axial fast field, reflect a paraxial electron off the
// mirror, and read its radius back at the launch plane.
func TestDohiMirror(t *testing.T) {
	// 'mirror' is an aperture plus a line capping it at z=0 (both named 'mirror' → one electrode).
	mirrorAp := geometry.Aperture(dohiThickness, dohiRadius, dohiExtent, dohiThickness/2).WithName("mirror")
	mirrorLine := geometry.Line(geometry.Point{0, 0, 0}, geometry.Point{dohiRadius, 0, 0}).WithName("mirror")
	lens := geometry.Aperture(dohiThickness, dohiRadius, dohiExtent, dohiThickness+dohiSpacer+dohiThickness/2).WithName("lens")
	ground := geometry.Aperture(dohiThickness, dohiRadius, dohiExtent, 2*dohiThickness+2*dohiSpacer+dohiThickness/2).WithName("ground")
	boundary := geometry.Line(geometry.Point{0, 0, 1.75}, geometry.Point{1.0, 0, 1.75}).
		ExtendWithLine(geometry.Point{1.0, 0, -0.3}).
		ExtendWithLine(geometry.Point{0, 0, -0.3}).WithName("boundary")

	m := geometry.MeshGroup([]geometry.Path{mirrorAp, mirrorLine, lens, ground, boundary},
		geometry.MeshOptions{MeshSizeFactor: dohiMSF, HigherOrder: true, EnsureOutwardNormals: true})

	exc := excitation.New(m)
	exc.AddVoltage("ground", 0)
	exc.AddVoltage("mirror", -1250)
	exc.AddVoltage("lens", dohiLensVolts)
	exc.AddElectrostaticBoundary("boundary")

	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve electrostatic: %v", err)
	}

	// Axial fast field over the electrode region; zero above z=1.7 (field-free drift to z=15).
	fa, err := field.NewFieldRadialAxial(charges, 0.05, 1.7, 500)
	if err != nil {
		t.Fatalf("axial field: %v", err)
	}
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := fa.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	v0 := tracing.VelocityVecXZPlane(1000, dohiAngle, true, constants.ElectronMass)
	qOverM := -constants.ElementaryCharge / constants.ElectronMass
	bounds := tracing.Bounds{{-0.1, 0.1}, {-0.03, 0.03}, {0.05, 19.0}}

	// atol 1e-8 matches Traceon's default tracer tolerance, so the adaptive RKF45 steps align
	// (the mirror reflection makes the trace sensitive to step-size differences).
	_, states := tracing.TraceParticle(geom3d.Vec3{0, 0, dohiZ0}, v0, qOverM, fieldFn, bounds, 1e-8)
	cross, ok := tracing.XYPlaneIntersection(states, dohiZ0)
	if !ok {
		t.Fatal("reflected electron never re-crossed the launch plane z=15")
	}
	r := cross[0]

	t.Logf("return radius = %.10g (upstream %.10g, ideal 0)", r, dohiUpstreamMSF100)
	// Port equivalence: Go must reproduce Traceon's own return radius at this resolution.
	if d := math.Abs(r - dohiUpstreamMSF100); d > 1e-7 {
		t.Errorf("return radius = %.10g, want %.10g (upstream MSF=100); abs diff %.2e > 1e-7", r, dohiUpstreamMSF100, d)
	}
	// Physics: the mirror returns the ray close to the axis.
	if math.Abs(r) > 1e-3 {
		t.Errorf("return radius = %.10g, want |r| < 1e-3 (ray should return near axis)", r)
	}
}
