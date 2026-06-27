// SPDX-License-Identifier: MPL-2.0

// Command dohi-mirror reproduces upstream Traceon's mirror example in pure Go: it builds the
// Dohi electron mirror, solves the electrostatic BEM, and traces electrons that fly in, reflect
// off the negatively-biased mirror, and fly back out — rendering the trajectories + electrodes
// to images/dohi-mirror.png.
//
// Run from the repo root:  go run ./examples/dohi-mirror
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
	thickness = 0.15
	radius    = 0.075
	spacer    = 0.5
	extent    = 1.0 - 0.1
	numRays   = 5
)

func main() {
	mirrorAp := geometry.Aperture(thickness, radius, extent, thickness/2).WithName("mirror")
	mirrorLine := geometry.Line(geometry.Point{0, 0, 0}, geometry.Point{radius, 0, 0}).WithName("mirror")
	lens := geometry.Aperture(thickness, radius, extent, thickness+spacer+thickness/2).WithName("lens")
	ground := geometry.Aperture(thickness, radius, extent, 2*thickness+2*spacer+thickness/2).WithName("ground")
	boundary := geometry.Line(geometry.Point{0, 0, 1.75}, geometry.Point{1.0, 0, 1.75}).
		ExtendWithLine(geometry.Point{1.0, 0, -0.3}).
		ExtendWithLine(geometry.Point{0, 0, -0.3}).WithName("boundary")
	paths := []geometry.Path{mirrorAp, mirrorLine, lens, ground, boundary}

	m := geometry.MeshGroup(paths, geometry.MeshOptions{MeshSizeFactor: 50, HigherOrder: true, EnsureOutwardNormals: true})

	exc := excitation.New(m)
	exc.AddVoltage("ground", 0)
	exc.AddVoltage("mirror", -1250)
	exc.AddVoltage("lens", 710.0126605741955)
	exc.AddElectrostaticBoundary("boundary")

	lines, types, values := exc.Electrostatic()
	charges, err := solver.SolveElectrostatic(lines, types, values)
	if err != nil {
		log.Fatalf("solve: %v", err)
	}
	fa, err := field.NewFieldRadialAxial(charges, 0.05, 1.7, 500)
	if err != nil {
		log.Fatalf("axial field: %v", err)
	}

	rays := traceMirror(fa)
	render(paths, field.NewFieldRadialBEM(charges), rays)
	fmt.Printf("dohi mirror: %d elements, %d reflected rays → images/dohi-mirror.png\n", len(lines), len(rays))
}

// traceMirror launches electrons down the axis from z=2.4 at a spread of radii; each decelerates
// in the mirror's retarding field, reverses near z≈0, and returns — so the figure shows the
// reflection as a fan of in-and-out paths rather than one straight axial line.
func traceMirror(fa field.FieldRadialAxial) [][2][]float64 {
	fieldFn := func(pos, _ geom3d.Vec3) (elec, mag geom3d.Vec3) {
		e := fa.FieldAtPoint(geom2d.Vertex{pos[0], pos[1], pos[2]})
		return geom3d.Vec3{e[0], e[1], e[2]}, geom3d.Vec3{}
	}
	v0 := tracing.VelocityVec(1000, geom3d.Vec3{0, 0, -1}, constants.ElectronMass) // straight in (−z)
	qOverM := -constants.ElementaryCharge / constants.ElectronMass
	bounds := tracing.Bounds{{-0.08, 0.08}, {-0.03, 0.03}, {0.04, 3.0}}

	var rays [][2][]float64
	for i := 0; i < numRays; i++ {
		r0 := 0.01 + 0.05*float64(i)/float64(numRays-1) // distinct launch radii 0.01..0.06
		_, states := tracing.TraceParticle(geom3d.Vec3{r0, 0, 2.4}, v0, qOverM, fieldFn, bounds, 1e-9)
		zs := make([]float64, 0, len(states))
		rs := make([]float64, 0, len(states))
		for _, s := range states {
			zs = append(zs, s[2])
			rs = append(rs, s[0])
		}
		if len(zs) > 1 {
			rays = append(rays, [2][]float64{zs, rs})
		}
	}
	return rays
}

func render(paths []geometry.Path, bem field.FieldRadialBEM, rays [][2][]float64) {
	const grid = 220
	zMin, zMax, rMax := -0.2, 2.6, 0.9
	c := plot.New(700, 600, zMin, zMax, -rMax, rMax)

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
	for _, p := range paths {
		zs, rs := samplePath(p)
		c.Polyline(zs, rs, plot.Mirror, 2)
		c.Polyline(zs, negate(rs), plot.Mirror, 2)
	}
	for _, ray := range rays {
		c.Polyline(ray[0], ray[1], plot.Ray, 2)
	}
	if err := c.Save("images/dohi-mirror.png"); err != nil {
		log.Fatalf("save: %v", err)
	}
}

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
