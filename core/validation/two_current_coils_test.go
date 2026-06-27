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

// Two opposed current coils with a magnetizable block, from validation/two_current_coils.py.
// Two disk coils carry ±10 A; a rectangular μ=25 block sits between them inside a magnetostatic
// boundary. The value of interest is the radial magnetic field at (r=10 mm, z=0), which lies on
// the block — a slowly-converging point, so the reference match is loose by design.
const (
	tcc2MSF       = 16     // line mesh_size_factor for the boundary + block
	tcc2CoilMesh  = 0.25e-3 // coil cross-section element size
	tcc2BlockPerm = 25.0
	tcc2Current   = 10.0
)

// tccPaperValue is the reference Hr from Lencova's benchmark (two_current_coils.py).
const tccPaperValue = -91.94907464785867

// tccUpstreamMSF16 is the Hr Traceon itself computes at MSF=16 (higher-order line elements,
// 0.25 mm coil mesh) — the port-equivalence oracle. Its relative error to the reference at this
// resolution is ~9.7% (the field point sits on the magnetizable block).
const tccUpstreamMSF16 = -100.83516911750684

// TestTwoCurrentCoils reproduces the two-coil + magnetizable-block magnetic field end to end:
// surface-mesh the coils into current rings, line-mesh the block + boundary, solve the
// magnetostatic BEM with the coil field as the pre-field, and evaluate H. This exercises the
// whole new stack — surface mesher, current coils, magnetizable response, boundary — at once.
func TestTwoCurrentCoils(t *testing.T) {
	// Coils: two disk cross-sections meshed into current rings (the magnetic pre-field source).
	coil1 := geometry.DiskXZ(10e-3, 5e-3, 1e-3).WithName("coil1")
	coil2 := geometry.DiskXZ(10e-3, -5e-3, 1e-3).WithName("coil2")
	coilMesh := geometry.MeshSurfaceGroup([]geometry.Surface{coil1, coil2}, tcc2CoilMesh)

	coils := excitation.NewCoils(coilMesh)
	coils.AddCurrent("coil1", tcc2Current)
	coils.AddCurrent("coil2", -tcc2Current)
	currentCharges := coils.Charges()

	// Boundary + magnetizable block: line elements solved with the coil field as pre-field.
	boundary := geometry.Line(geometry.Point{0, 0, 50e-3}, geometry.Point{100e-3, 0, 50e-3}).
		ExtendWithLine(geometry.Point{100e-3, 0, -50e-3}).
		ExtendWithLine(geometry.Point{0, 0, -50e-3}).WithName("boundary")
	block := geometry.RectangleXZ(5e-3, 15e-3, -1e-3, 1e-3).WithName("block")
	lineMesh := geometry.MeshGroup([]geometry.Path{boundary, block},
		geometry.MeshOptions{MeshSizeFactor: tcc2MSF, HigherOrder: true, EnsureOutwardNormals: true})

	exc := excitation.New(lineMesh)
	exc.AddMagnetizable("block", tcc2BlockPerm)
	exc.AddMagnetostaticBoundary("boundary")
	lines, types, values := exc.Magnetostatic()

	// The pre-existing field the magnetizable block responds to is the coils' magnetic field.
	curField := field.NewFieldRadialBEMFull(solver.EffectivePointCharges{}, solver.EffectivePointCharges{}, currentCharges)
	preField := func(p geom3d.Vec3) geom3d.Vec3 {
		h := curField.CurrentFieldAtPoint(geom2d.Vertex{p[0], p[1], p[2]})
		return geom3d.Vec3{h[0], h[1], h[2]}
	}
	magCharges, err := solver.SolveMagnetostatic(lines, types, values, preField)
	if err != nil {
		t.Fatalf("solve magnetostatic: %v", err)
	}

	// Total magnetic field = coil current field + the magnetizable/boundary surface charges.
	bem := field.NewFieldRadialBEMFull(solver.EffectivePointCharges{}, magCharges, currentCharges)
	hx := bem.MagnetostaticFieldAtPoint(geom2d.Vertex{10e-3, 0, 0})[0]

	t.Logf("Hr = %.12g (upstream %.12g, paper %.12g)", hx, tccUpstreamMSF16, tccPaperValue)
	// Port equivalence: Go must reproduce Traceon's own number at this resolution.
	if rel := math.Abs(hx-tccUpstreamMSF16) / math.Abs(tccUpstreamMSF16); rel > 1e-5 {
		t.Errorf("Hr = %.12g, want %.12g (upstream MSF=16); rel err %.2e > 1e-5", hx, tccUpstreamMSF16, rel)
	}
	// Physics: right sign and order of magnitude (the on-block point converges slowly to paper).
	if rel := math.Abs(hx-tccPaperValue) / math.Abs(tccPaperValue); rel > 0.12 {
		t.Errorf("Hr = %.12g vs paper %.12g; rel err %.2e > 0.12 (unexpectedly far)", hx, tccPaperValue, rel)
	}
}
