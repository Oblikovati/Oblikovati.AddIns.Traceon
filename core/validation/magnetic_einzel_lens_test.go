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

// Magnetic einzel lens, from validation/magnetic_einzel_lens.py. The same three-aperture
// geometry as the electrostatic einzel, but excited with a magnetostatic scalar potential
// (central pole at 50 A, outer poles grounded, with a magnetostatic boundary). A 1 keV
// electron is focused by the magnetic field; the focal length is its axis crossing.
const (
	meThickness = 0.5
	meSpacing   = 0.5
	meRadius    = 0.15
	meExtent    = 2.0 - 0.1
	meMSF       = 20
)

// mePaperFocal is the reference focal length from magnetic_einzel_lens.py.
const mePaperFocal = 4.08641734

// meUpstreamFocalMSF20 is the focal length Traceon itself computes at MSF=20, higher_order —
// the port-equivalence oracle. Traceon's own relative error to the reference is 1.2e-3.
const meUpstreamFocalMSF20 = 4.09146507143068

// TestMagneticEinzelLensFocalLength reproduces the magnetic einzel focal length end to end:
// solve the magnetostatic BEM (scalar potential + boundary), build the axial-series magnetic
// field, trace a paraxial electron under the v×B force, and read its axis crossing. This is
// the first validation of the magnetostatic solve + magnetic axial fast tracer together.
func TestMagneticEinzelLensFocalLength(t *testing.T) {
	boundary := geometry.Line(geometry.Point{0, 0, 1.75}, geometry.Point{2, 0, 1.75}).
		ExtendWithLine(geometry.Point{2, 0, -1.75}).
		ExtendWithLine(geometry.Point{0, 0, -1.75}).WithName("boundary")
	bottom := geometry.Aperture(meThickness, meRadius, meExtent, -meThickness-meSpacing).WithName("ground")
	middle := geometry.Aperture(meThickness, meRadius, meExtent, 0).WithName("lens")
	top := geometry.Aperture(meThickness, meRadius, meExtent, meThickness+meSpacing).WithName("ground")

	m := geometry.MeshGroup([]geometry.Path{boundary, bottom, middle, top},
		geometry.MeshOptions{MeshSizeFactor: meMSF, HigherOrder: true})

	exc := excitation.New(m)
	exc.AddMagnetostaticPotential("ground", 0)
	exc.AddMagnetostaticPotential("lens", 50)
	exc.AddMagnetostaticBoundary("boundary")

	lines, types, values := exc.Magnetostatic()
	charges, err := solver.SolveMagnetostatic(lines, types, values, nil) // no currents/permanent magnets
	if err != nil {
		t.Fatalf("solve magnetostatic: %v", err)
	}

	// Axial-series interpolation of the magnetic field (H) over the lens span.
	fa, err := field.NewFieldRadialAxial(charges, -1.5, 1.5, 1000)
	if err != nil {
		t.Fatalf("axial field: %v", err)
	}
	// Magnetic-only Lorentz force: a = q/m · μ₀ · v×H (the tracer applies μ₀ internally).
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		h := fa.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{}, geom3d.Vec3{h[0], h[1], h[2]}
	}
	v0 := tracing.VelocityVecXZPlane(1000, 0, true, constants.ElectronMass)
	qOverM := -constants.ElementaryCharge / constants.ElectronMass
	bounds := tracing.Bounds{{-meRadius, meRadius}, {-meRadius, meRadius}, {-5, 3.5}}

	_, states := tracing.TraceParticle(geom3d.Vec3{meRadius / 5, 0, 3}, v0, qOverM, fieldFn, bounds, 1e-10)
	z, ok := tracing.AxisIntersection(states)
	if !ok {
		t.Fatal("focused electron never crossed the axis")
	}
	focal := -z

	t.Logf("focal = %.12g (upstream %.12g, paper %.12g)", focal, meUpstreamFocalMSF20, mePaperFocal)
	if rel := math.Abs(focal-meUpstreamFocalMSF20) / meUpstreamFocalMSF20; rel > 1e-5 {
		t.Errorf("focal = %.12g, want %.12g (upstream MSF=20); rel err %.2e > 1e-5", focal, meUpstreamFocalMSF20, rel)
	}
	if rel := math.Abs(focal-mePaperFocal) / mePaperFocal; rel > 3e-3 {
		t.Errorf("focal = %.12g vs paper %.12g; rel err %.2e > 3e-3", focal, mePaperFocal, rel)
	}
}
