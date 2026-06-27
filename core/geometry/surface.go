// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"math"

	"oblikovati.org/traceon/core/quad"
)

// Surface is a parametric map (u, v) → 3D point over [0, PathLength1] × [0, PathLength2], the
// region a radial current coil's solid cross-section is built from. Breakpoints in each
// parameter direction force section boundaries the mesher respects. Port of geometry.Surface.
type Surface struct {
	Fun          func(u, v float64) Point
	PathLength1  float64
	PathLength2  float64
	Breakpoints1 []float64
	Breakpoints2 []float64
	Name         string
}

// At evaluates the surface at (u, v).
func (s Surface) At(u, v float64) Point { return s.Fun(u, v) }

// WithName tags the surface (its mesh becomes one physical group under this name).
func (s Surface) WithName(name string) Surface {
	s.Name = name
	return s
}

// Move translates the surface by (dx, dy, dz).
func (s Surface) Move(dx, dy, dz float64) Surface {
	f := s.Fun
	s.Fun = func(u, v float64) Point {
		p := f(u, v)
		return Point{p[0] + dx, p[1] + dy, p[2] + dz}
	}
	return s
}

// sections splits the surface at its breakpoints into a grid of sub-surfaces, each spanning one
// breakpoint cell with its own local (0..span) parametrization. Mirrors Surface._sections.
func (s Surface) sections() []Surface {
	b1 := append(append([]float64{0}, s.Breakpoints1...), s.PathLength1)
	b2 := append(append([]float64{0}, s.Breakpoints2...), s.PathLength2)
	f := s.Fun

	var out []Surface
	for i := 0; i+1 < len(b1); i++ {
		for j := 0; j+1 < len(b2); j++ {
			u0, v0 := b1[i], b2[j]
			out = append(out, Surface{
				Fun:         func(u, v float64) Point { return f(u0+u, v0+v) },
				PathLength1: b1[i+1] - u0,
				PathLength2: b2[j+1] - v0,
			})
		}
	}
	return out
}

// quadEpsabs / quadEpsrel are scipy.quad's default tolerances, used so Average matches the
// upstream r_avg integral that sets a revolved surface's circumferential length.
const (
	quadEpsabs = 1.49e-8
	quadEpsrel = 1.49e-8
)

// Average integrates fun(path(s)) over the path's arc length and divides by the length — the
// mean of fun along the path. Port of Path.average (scipy.quad with the breakpoints).
func (p Path) Average(fun func(Point) float64) float64 {
	integrand := func(s float64) float64 { return fun(p.Fun(s)) }
	return quad.IntegrateWithSingularities(integrand, 0, p.Length, p.Breakpoints, quadEpsabs, quadEpsrel) / p.Length
}

// RevolveY revolves the path anticlockwise about the y-axis through the given angle, producing
// a surface. The second parameter spans the circumferential arc length 2π·r_avg (r_avg the
// arc-length-mean radius), so the mesher's element sizing is isotropic. Port of Path.revolve_y.
func (p Path) RevolveY(angle float64) Surface {
	rAvg := p.Average(func(pt Point) float64 { return math.Hypot(pt[0], pt[2]) })
	length2 := 2 * math.Pi * rAvg
	fun := func(u, v float64) Point {
		pt := p.Fun(u)
		theta := math.Atan2(pt[2], pt[0])
		r := math.Hypot(pt[0], pt[2])
		a := theta + v/length2*angle
		return Point{r * math.Cos(a), pt[1], r * math.Sin(a)}
	}
	return Surface{Fun: fun, PathLength1: p.Length, PathLength2: length2, Breakpoints1: p.Breakpoints, Name: p.Name}
}

// DiskXZ builds a flat disk of the given radius centred at (x0, 0, z0) in the xz-plane: the
// radius line revolved a full turn about the y-axis. The meridian of a radial current coil.
// Port of Surface.disk_xz.
func DiskXZ(x0, z0, radius float64) Surface {
	return Line(Point{0, 0, 0}, Point{radius, 0, 0}).RevolveY(2 * math.Pi).Move(x0, 0, z0)
}

// Extrude sweeps the path along vector, producing a surface: the path at v=0, linearly
// translated to path+vector at v=|vector|. Port of Path.extrude.
func (p Path) Extrude(vector Point) Surface {
	length := math.Sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2])
	fun := func(u, v float64) Point {
		pt := p.Fun(u)
		s := v / length
		return Point{pt[0] + s*vector[0], pt[1] + s*vector[1], pt[2] + s*vector[2]}
	}
	return Surface{Fun: fun, PathLength1: p.Length, PathLength2: length, Breakpoints1: p.Breakpoints, Name: p.Name}
}

// RectangleXZSurface builds a filled rectangle region in the xz-plane (a coil's solid
// cross-section): the left edge line [xmin,zmin]→[xmin,zmax] extruded across in +x by
// (xmax−xmin). Port of Surface.rectangle_xz.
func RectangleXZSurface(xmin, xmax, zmin, zmax float64) Surface {
	return Line(Point{xmin, 0, zmin}, Point{xmin, 0, zmax}).Extrude(Point{xmax - xmin, 0, 0})
}
