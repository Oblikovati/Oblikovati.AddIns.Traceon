// SPDX-License-Identifier: MPL-2.0

// Package geometry is the pure-Go port of Traceon's radially-symmetric parametric
// geometry: a Path is an arc-length-parametrized curve in the (x=r, y=0, z) meridian
// plane, built from lines and arcs and discretized into the line elements the radial
// BEM solver consumes (see Path.Mesh in mesh.go). The 3D Surface/revolve machinery in
// upstream geometry.py is closed-source "Traceon Pro" and is intentionally not ported.
//
// The parameter of a Path is true arc length: fun is unit-speed for the analytic
// builders, so uniform sampling in the parameter is uniform in distance along the curve.
package geometry

import (
	"fmt"
	"math"

	"oblikovati.org/traceon/core/geom3d"
)

// Point is a point on a path (x=r, y, z=axial); radial paths keep y == 0.
type Point = geom3d.Vec3

// Path is an arc-length-parametrized curve. Fun maps t in [0, Length] to a point;
// Breakpoints are arc-length positions of corners that must land on a mesh node; Name is
// the optional electrode/physical-group name carried into the mesh.
type Path struct {
	Fun         func(t float64) Point
	Length      float64
	Breakpoints []float64
	Name        string
}

// At evaluates the path at arc length t.
func (p Path) At(t float64) Point { return p.Fun(t) }

// StartingPoint, MiddlePoint and Endpoint sample the path's ends and centre.
func (p Path) StartingPoint() Point { return p.Fun(0) }
func (p Path) MiddlePoint() Point   { return p.Fun(p.Length / 2) }
func (p Path) Endpoint() Point      { return p.Fun(p.Length) }

// Line is the straight segment from `from` to `to`, parametrized by arc length.
func Line(from, to Point) Path {
	length := geom3d.Distance3D(from, to)
	fun := func(pl float64) Point {
		f := pl / length
		return add(scale(from, 1-f), scale(to, f))
	}
	return Path{Fun: fun, Length: length}
}

// Arc is the circular arc on the circle through `start` about `center`, sweeping to the
// angular position of `end` (projected onto the circle's plane). With reverse it takes the
// long way round (theta_max -= 2π). Mirrors geometry.Path.arc.
func Arc(center, start, end Point, reverse bool) Path {
	xUnit := geom3d.Normalize3D(sub(start, center))
	vector := sub(end, center)
	yUnit := geom3d.Normalize3D(sub(vector, scale(xUnit, geom3d.Dot3D(vector, xUnit))))
	radius := geom3d.Distance3D(start, center)
	thetaMax := math.Atan2(geom3d.Dot3D(vector, yUnit), geom3d.Dot3D(vector, xUnit))
	if reverse {
		thetaMax -= 2 * math.Pi
	}
	length := math.Abs(thetaMax * radius)
	fun := func(l float64) Point {
		theta := (l / length) * thetaMax
		return add(center, add(scale(xUnit, radius*math.Cos(theta)), scale(yUnit, radius*math.Sin(theta))))
	}
	return Path{Fun: fun, Length: length}
}

// CircleXZ is a full (or partial) circle of the given radius in the meridian (x, z) plane,
// centred at (x0, z0), starting on the +x side. Radial (y == 0). Mirrors circle_xz.
func CircleXZ(x0, z0, radius, angle float64) Path {
	fun := func(u float64) Point {
		theta := u / radius
		return Point{x0 + radius*math.Cos(theta), 0, z0 + radius*math.Sin(theta)}
	}
	return Path{Fun: fun, Length: angle * radius}
}

// Then concatenates two paths (the `>>` operator): the end of p must coincide with the
// start of other. The join arc-length becomes a breakpoint so it lands on a mesh node.
func (p Path) Then(other Path) Path {
	if d := geom3d.Distance3D(p.Endpoint(), other.StartingPoint()); d > joinTolerance {
		panic(fmt.Sprintf("cannot join paths: endpoint %v and start %v differ by %g (> %g)",
			p.Endpoint(), other.StartingPoint(), d, joinTolerance))
	}
	l1 := p.Length
	fun := func(t float64) Point {
		if t <= l1 {
			return p.Fun(t)
		}
		return other.Fun(t - l1)
	}
	bps := append([]float64{}, p.Breakpoints...)
	bps = append(bps, l1)
	for _, b := range other.Breakpoints {
		bps = append(bps, b+l1)
	}
	return Path{Fun: fun, Length: l1 + other.Length, Breakpoints: bps, Name: p.Name}
}

// joinTolerance is the maximum endpoint gap allowed when concatenating paths.
const joinTolerance = 1e-9

// ExtendWithLine appends a straight segment from the current endpoint to point.
func (p Path) ExtendWithLine(point Point) Path {
	return p.Then(Line(p.Endpoint(), point))
}

// ExtendWithArc appends an arc from the current endpoint about center to end.
func (p Path) ExtendWithArc(center, end Point, reverse bool) Path {
	return p.Then(Arc(center, p.Endpoint(), end, reverse))
}

// Close appends a straight segment back to the starting point.
func (p Path) Close() Path { return p.ExtendWithLine(p.StartingPoint()) }

// Reverse traverses the path from end to start, mirroring the breakpoints.
func (p Path) Reverse() Path {
	l := p.Length
	fun := func(t float64) Point { return p.Fun(l - t) }
	bps := make([]float64, len(p.Breakpoints))
	for i, b := range p.Breakpoints {
		bps[len(p.Breakpoints)-1-i] = l - b
	}
	return Path{Fun: fun, Length: l, Breakpoints: bps, Name: p.Name}
}

// Move translates every point on the path by (dx, dy, dz).
func (p Path) Move(dx, dy, dz float64) Path {
	delta := Point{dx, dy, dz}
	fun := func(t float64) Point { return add(p.Fun(t), delta) }
	return Path{Fun: fun, Length: p.Length, Breakpoints: p.Breakpoints, Name: p.Name}
}

// WithName tags the path with a physical-group name carried into the mesh.
func (p Path) WithName(name string) Path {
	p.Name = name
	return p
}

// RectangleXZ is the closed counter-clockwise rectangle in the meridian (x, z) plane.
// Its four corners become breakpoints. Mirrors rectangle_xz.
func RectangleXZ(xmin, xmax, zmin, zmax float64) Path {
	return Line(Point{xmin, 0, zmin}, Point{xmax, 0, zmin}).
		ExtendWithLine(Point{xmax, 0, zmax}).
		ExtendWithLine(Point{xmin, 0, zmax}).
		Close()
}

// Aperture is an open three-sided slab in the meridian plane: an electrode with a bore of
// the given radius, total height, and outer extent, centred at axial position z. Mirrors
// geometry.Path.aperture.
func Aperture(height, radius, extent, z float64) Path {
	return Line(Point{extent, 0, -height / 2}, Point{radius, 0, -height / 2}).
		ExtendWithLine(Point{radius, 0, height / 2}).
		ExtendWithLine(Point{extent, 0, height / 2}).
		Move(0, 0, z)
}

// --- small Vec3 helpers (geom3d exposes dot/cross/norm but not add/sub/scale) ---

func sub(a, b Point) Point { return Point{a[0] - b[0], a[1] - b[1], a[2] - b[2]} }
func add(a, b Point) Point { return Point{a[0] + b[0], a[1] + b[1], a[2] + b[2]} }
func scale(a Point, s float64) Point {
	return Point{a[0] * s, a[1] * s, a[2] * s}
}
