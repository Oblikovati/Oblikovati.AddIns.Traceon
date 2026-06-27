// SPDX-License-Identifier: MPL-2.0

// Command einzel-lens reproduces upstream Traceon's einzel-lens example in pure Go: it builds a
// three-aperture einzel lens, solves the electrostatic BEM, traces a bundle of electrons through
// the axial-series field, and renders the electrodes + trajectories over the potential map to
// images/einzel-lens.png.
//
// Run from the repo root:  go run ./examples/einzel-lens
package main

import (
	"fmt"
	"log"
	"math"

	"oblikovati.org/traceon/core/constants"
	"oblikovati.org/traceon/core/excitation"
	"oblikovati.org/traceon/core/field"
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geom3d"
	"oblikovati.org/traceon/core/geometry"
	"oblikovati.org/traceon/core/solver"
	"oblikovati.org/traceon/core/tracing"
	"oblikovati.org/traceon/examples/plot"
)

const (
	thickness  = 0.5
	spacing    = 0.5
	radius     = 0.15
	extent     = 2.0 - 0.1
	lensVolts  = 1800
	numRays    = 7
	beamEnergy = 1000 // eV
)

func main() {
	// Geometry: a grounded boundary box around three apertures — the outer pair grounded, the
	// centre biased to lensVolts (the classic decelerating-then-accelerating einzel lens).
	boundary := geometry.Line(geometry.Point{0, 0, 1.75}, geometry.Point{2, 0, 1.75}).
		ExtendWithLine(geometry.Point{2, 0, -1.75}).
		ExtendWithLine(geometry.Point{0, 0, -1.75}).WithName("boundary")
	bottom := geometry.Aperture(thickness, radius, extent, -thickness-spacing).WithName("ground")
	middle := geometry.Aperture(thickness, radius, extent, 0).WithName("lens")
	top := geometry.Aperture(thickness, radius, extent, thickness+spacing).WithName("ground")
	paths := []geometry.Path{boundary, bottom, middle, top}

	m := geometry.MeshGroup(paths, geometry.MeshOptions{MeshSizeFactor: 30, HigherOrder: true, EnsureOutwardNormals: true})

	exc := excitation.New(m)
	exc.AddVoltage("ground", 0)
	exc.AddVoltage("lens", lensVolts)
	exc.AddElectrostaticBoundary("boundary")

	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		log.Fatalf("solve: %v", err)
	}
	bem := field.NewFieldRadialBEM(charges)
	fa, err := field.NewFieldRadialAxial(charges, -1.5, 1.5, 150)
	if err != nil {
		log.Fatalf("axial field: %v", err)
	}

	rays := traceBeam(fa)
	render(paths, bem, rays)
	fmt.Printf("einzel lens: %d elements, %d rays → images/einzel-lens.png\n", len(lines), len(rays))
}

// traceBeam launches numRays electrons from z=5 heading toward the lens and returns their
// trajectories (r, z) in the meridian plane.
func traceBeam(fa field.FieldRadialAxial) [][2][]float64 {
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := fa.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	v0 := tracing.VelocityVec(beamEnergy, geom3d.Vec3{0, 0, -1}, constants.ElectronMass)
	qOverM := -constants.ElementaryCharge / constants.ElectronMass
	bounds := tracing.Bounds{{-radius, radius}, {-radius, radius}, {-10, 10}}

	var rays [][2][]float64
	for i := 0; i < numRays; i++ {
		frac := float64(i)/float64(numRays-1) - 0.5
		r0 := frac * 2 * radius / 3
		_, states := tracing.TraceParticle(geom3d.Vec3{r0, 0, 5}, v0, qOverM, fieldFn, bounds, 1e-8)
		var zs, rs []float64
		for _, s := range states {
			zs = append(zs, s[2])
			rs = append(rs, s[0])
		}
		rays = append(rays, [2][]float64{zs, rs})
	}
	return rays
}

// render draws the potential map, electrode outlines, and trajectories to the PNG.
func render(paths []geometry.Path, bem field.FieldRadialBEM, rays [][2][]float64) {
	const grid = 220
	zMin, zMax, rMax := -1.8, 1.8, 1.0
	c := plot.New(900, 500, zMin, zMax, -rMax, rMax)

	pot := make([][]float64, grid)
	for iz := 0; iz < grid; iz++ {
		pot[iz] = make([]float64, grid)
		z := zMin + (zMax-zMin)*float64(iz)/float64(grid-1)
		for ir := 0; ir < grid; ir++ {
			r := -rMax + 2*rMax*float64(ir)/float64(grid-1)
			pot[iz][ir] = bem.PotentialAtPoint(geom2d.Vertex{math.Abs(r), 0, z})
		}
	}
	c.Heatmap(pot)
	c.Axis()

	// Electrodes (drawn mirrored about the axis so the figure shows the full bore).
	for _, p := range paths {
		zs, rs := samplePath(p)
		c.Polyline(zs, rs, plot.Electrode, 2)
		c.Polyline(zs, negate(rs), plot.Electrode, 2)
	}
	// Trajectories.
	for _, ray := range rays {
		c.Polyline(ray[0], ray[1], plot.Ray, 2)
	}
	if err := c.Save("images/einzel-lens.png"); err != nil {
		log.Fatalf("save: %v", err)
	}
}

// samplePath samples a path into (z, r) polyline arrays for plotting.
func samplePath(p geometry.Path) (zs, rs []float64) {
	const n = 64
	for i := 0; i < n; i++ {
		pt := p.At(p.Length * float64(i) / float64(n-1))
		zs = append(zs, pt[2])
		rs = append(rs, pt[0])
	}
	return zs, rs
}

func negate(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = -x
	}
	return out
}
