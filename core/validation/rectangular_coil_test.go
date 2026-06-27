// SPDX-License-Identifier: MPL-2.0

package validation

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/excitation"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/geometry"
	"oblikovati.org/traceon/core/solver"
)

// Rectangular current coil, from validation/rectangular_coil.py and its with-circle variant. A
// 1 A rectangular coil cross-section inside a magnetostatic boundary; the value of interest is
// the radial magnetic field at (r=2.5, z=4). The variant adds a magnetizable torus (circle).
const rcMSF = 16

// Upstream oracles (higher-order line mesh, MSF=16): the radial field Hr at (2.5, 0, 4).
const (
	rcUpstreamHr        = 0.09183118623636087  // rectangular_coil
	rcCircleUpstreamHr  = 0.015519708142796379 // rectangular_coil_with_circle
)

// TestRectangularCoil reproduces the rectangular-coil radial field at (2.5, 4).
func TestRectangularCoil(t *testing.T) {
	coilMesh := geometry.RectangleXZSurface(2, 3, 2, 3).WithName("coil").MeshByFactor(rcMSF)
	coils := excitation.NewCoils(coilMesh)
	coils.AddCurrent("coil", 1)

	boundary := geometry.Line(geometry.Point{0, 0, 5}, geometry.Point{5, 0, 5}).
		ExtendWithLine(geometry.Point{5, 0, 0}).
		ExtendWithLine(geometry.Point{0, 0, 0}).WithName("boundary")
	lineMesh := geometry.MeshGroup([]geometry.Path{boundary},
		geometry.MeshOptions{MeshSizeFactor: rcMSF, HigherOrder: true, EnsureOutwardNormals: true})
	lineExc := excitation.New(lineMesh)
	lineExc.AddMagnetostaticBoundary("boundary")

	hr := solveCoilHr(t, coils, lineExc, geom2d.Vertex{2.5, 0, 4})

	t.Logf("Hr = %.12g (upstream %.12g)", hr, rcUpstreamHr)
	if rel := math.Abs(hr-rcUpstreamHr) / math.Abs(rcUpstreamHr); rel > 1e-5 {
		t.Errorf("Hr = %.12g, want %.12g (upstream MSF=16); rel err %.2e > 1e-5", hr, rcUpstreamHr, rel)
	}
}

// TestRectangularCoilWithCircle adds a magnetizable torus (μ=10) the coil field magnetizes.
func TestRectangularCoilWithCircle(t *testing.T) {
	coilMesh := geometry.RectangleXZSurface(2, 3, 2, 3).WithName("coil").Mesh(0.1)
	coils := excitation.NewCoils(coilMesh)
	coils.AddCurrent("coil", 1)

	boundary := geometry.Line(geometry.Point{0, 0, 5}, geometry.Point{5, 0, 5}).
		ExtendWithLine(geometry.Point{5, 0, 0}).
		ExtendWithLine(geometry.Point{0, 0, 0}).WithName("boundary")
	circle := geometry.CircleXZ(2.5, 4, 0.5, 2*math.Pi).WithName("circle")
	lineMesh := geometry.MeshGroup([]geometry.Path{boundary, circle},
		geometry.MeshOptions{MeshSizeFactor: rcMSF, HigherOrder: true, EnsureOutwardNormals: true})
	lineExc := excitation.New(lineMesh)
	lineExc.AddMagnetostaticBoundary("boundary")
	lineExc.AddMagnetizable("circle", 10)

	hr := solveCoilHr(t, coils, lineExc, geom2d.Vertex{2.5, 0, 4})

	t.Logf("Hr = %.12g (upstream %.12g)", hr, rcCircleUpstreamHr)
	if rel := math.Abs(hr-rcCircleUpstreamHr) / math.Abs(rcCircleUpstreamHr); rel > 1e-5 {
		t.Errorf("Hr = %.12g, want %.12g (upstream MSF=16); rel err %.2e > 1e-5", hr, rcCircleUpstreamHr, rel)
	}
}

// solveCoilHr solves the coil current field + the line electrodes' magnetostatic response and
// returns Hr at the sample point.
func solveCoilHr(t *testing.T, coils *excitation.CoilExcitation, lineExc *excitation.Excitation, at geom2d.Vertex) float64 {
	t.Helper()
	currentCharges := coils.Charges()
	lines, types, values := lineExc.Magnetostatic()

	curField := field.NewFieldRadialBEMFull(solver.EffectivePointCharges{}, solver.EffectivePointCharges{}, currentCharges)
	preField := func(p geom3d.Vec3) geom3d.Vec3 {
		h := curField.CurrentFieldAtPoint(geom2d.Vertex{p[0], p[1], p[2]})
		return geom3d.Vec3{h[0], h[1], h[2]}
	}
	magCharges, err := solver.SolveMagnetostatic(lines, types, values, preField)
	if err != nil {
		t.Fatalf("solve magnetostatic: %v", err)
	}
	bem := field.NewFieldRadialBEMFull(solver.EffectivePointCharges{}, magCharges, currentCharges)
	return bem.MagnetostaticFieldAtPoint(at)[0]
}
