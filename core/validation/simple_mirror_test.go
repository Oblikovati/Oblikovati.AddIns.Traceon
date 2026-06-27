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

// Simple flat mirror, from validation/simple_mirror.py. A −110 V flat disk mirror inside a 0 V
// boundary reflects an electron launched at z=10; the value of interest is its radius back at
// z=10 after reflection.
const (
	smMSF      = 64
	smAxialN   = 153 // FieldRadialAxial(field, 0.02, 4) auto-samples to this N at MSF=64
	smPaper    = 1.6327355811e-01
	smUpstream = 0.16344628081660606 // Traceon's computed return radius at MSF=64
)

// TestSimpleMirror reproduces the flat-mirror return radius: solve the mirror + boundary, build
// the axial field, reflect an electron, and read its radius back at the launch plane.
func TestSimpleMirror(t *testing.T) {
	boundary := geometry.Line(geometry.Point{0, 0, -1}, geometry.Point{2, 0, -1}).
		ExtendWithLine(geometry.Point{2, 0, 1}).
		ExtendWithLine(geometry.Point{0.3, 0, 1}).WithName("boundary")
	mirror := geometry.Line(geometry.Point{0, 0, 0}, geometry.Point{1, 0, 0}).WithName("mirror")

	m := geometry.MeshGroup([]geometry.Path{boundary, mirror},
		geometry.MeshOptions{MeshSizeFactor: smMSF, HigherOrder: true, EnsureOutwardNormals: true})

	exc := excitation.New(m)
	exc.AddVoltage("mirror", -110)
	exc.AddVoltage("boundary", 0)

	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve electrostatic: %v", err)
	}

	fa, err := field.NewFieldRadialAxial(charges, 0.02, 4, smAxialN)
	if err != nil {
		t.Fatalf("axial field: %v", err)
	}
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := fa.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	v0 := tracing.VelocityVecXZPlane(100, 1e-3, true, constants.ElectronMass)
	qOverM := -constants.ElementaryCharge / constants.ElectronMass
	bounds := tracing.Bounds{{-0.22, 0.22}, {-0.22, 0.22}, {0.02, 11}}

	// atol 1e-8 matches Traceon's default tracer tolerance (mirror reflection is step-sensitive).
	_, states := tracing.TraceParticle(geom3d.Vec3{0, 0, 10}, v0, qOverM, fieldFn, bounds, 1e-8)
	cross, ok := tracing.XYPlaneIntersection(states, 10)
	if !ok {
		t.Fatal("reflected electron never re-crossed z=10")
	}
	r := math.Abs(cross[0])

	t.Logf("return radius = %.12g (upstream %.12g, paper %.12g)", r, smUpstream, smPaper)
	if d := math.Abs(r - smUpstream); d > 1e-6 {
		t.Errorf("return radius = %.12g, want %.12g (upstream MSF=64); abs diff %.2e > 1e-6", r, smUpstream, d)
	}
	if rel := math.Abs(r-smPaper) / smPaper; rel > 2e-3 {
		t.Errorf("return radius = %.12g vs paper %.12g; rel err %.2e > 2e-3", r, smPaper, rel)
	}
}
