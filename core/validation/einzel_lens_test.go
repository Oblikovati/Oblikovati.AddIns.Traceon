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

// Three-aperture einzel lens, from validation/einzel_lens.py. The outer aperture pair is
// grounded, the central one biased to 1000 V; a 1 keV electron launched parallel to the axis
// is focused, and the focal length is read from where it crosses the axis. This case exercises
// same-name electrode merging, an electrostatic boundary, and the axial-series fast tracer.
const (
	elThickness = 0.5
	elSpacing   = 0.5
	elRadius    = 0.15
	elExtent    = 2.0 - 0.1 // aperture outer extent (margin_right = 0.1)
	elMSF       = 20
)

// elPaperFocal is the reference focal length from einzel_lens.py.correct_value_of_interest.
const elPaperFocal = 3.915970140918643

// elUpstreamFocalMSF20 is the focal length Traceon itself computes at MSF=20, higher_order —
// the port-equivalence oracle (regenerate with validation/einzel_lens.py). Traceon's own
// relative error to the reference at this resolution is 1.1e-3.
const elUpstreamFocalMSF20 = 3.920101650400111

// TestEinzelLensFocalLength reproduces the einzel-lens focal length end to end: solve the
// biased three-aperture lens (with a grounded boundary), build the axial-series field, trace a
// paraxial electron, and read its axis crossing. Checks both Traceon's computed value and the
// reference.
func TestEinzelLensFocalLength(t *testing.T) {
	boundary := geometry.Line(geometry.Point{0, 0, 1.75}, geometry.Point{2, 0, 1.75}).
		ExtendWithLine(geometry.Point{2, 0, -1.75}).
		ExtendWithLine(geometry.Point{0, 0, -1.75}).WithName("boundary")
	bottom := geometry.Aperture(elThickness, elRadius, elExtent, -elThickness-elSpacing).WithName("ground")
	middle := geometry.Aperture(elThickness, elRadius, elExtent, 0).WithName("lens")
	top := geometry.Aperture(elThickness, elRadius, elExtent, elThickness+elSpacing).WithName("ground")

	m := geometry.MeshGroup([]geometry.Path{boundary, bottom, middle, top},
		geometry.MeshOptions{MeshSizeFactor: elMSF, HigherOrder: true})

	exc := excitation.New(m)
	exc.AddVoltage("ground", 0)
	exc.AddVoltage("lens", 1000)
	exc.AddElectrostaticBoundary("boundary")

	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		t.Fatalf("solve electrostatic: %v", err)
	}

	// Fast axial-series field over the lens span, then trace a paraxial 1 keV electron.
	fa, err := field.NewFieldRadialAxial(charges, -1.5, 1.5, 600)
	if err != nil {
		t.Fatalf("axial field: %v", err)
	}
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := fa.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	v0 := tracing.VelocityVecXZPlane(1000, 0, true, constants.ElectronMass)
	qOverM := -constants.ElementaryCharge / constants.ElectronMass
	bounds := tracing.Bounds{{-elRadius, elRadius}, {-elRadius, elRadius}, {-5, 3.5}}

	_, states := tracing.TraceParticle(geom3d.Vec3{elRadius / 3, 0, 3}, v0, qOverM, fieldFn, bounds, 1e-10)
	z, ok := tracing.AxisIntersection(states)
	if !ok {
		t.Fatal("focused electron never crossed the axis")
	}
	focal := -z

	t.Logf("focal = %.12g (upstream %.12g, paper %.12g)", focal, elUpstreamFocalMSF20, elPaperFocal)
	if rel := math.Abs(focal-elUpstreamFocalMSF20) / elUpstreamFocalMSF20; rel > 1e-5 {
		t.Errorf("focal = %.12g, want %.12g (upstream MSF=20); rel err %.2e > 1e-5", focal, elUpstreamFocalMSF20, rel)
	}
	if rel := math.Abs(focal-elPaperFocal) / elPaperFocal; rel > 3e-3 {
		t.Errorf("focal = %.12g vs paper %.12g; rel err %.2e > 3e-3", focal, elPaperFocal, rel)
	}
}
