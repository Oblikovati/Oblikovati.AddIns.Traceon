// SPDX-License-Identifier: MPL-2.0

// Package radial ports the electrostatic core of traceon/backend/radial.c: the radial
// BEM matrix assembly and field evaluation from effective point charges. A line element
// is four control points (GMSH "line4" ordering) in the (r, z) half-plane; the upstream
// remaps them to (v1,v2,v3,v4) = (p0, p2, p3, p1) everywhere, which this package mirrors.
//
// Scope of this file is the electrostatic path (voltage + dielectric). Current-coil
// (magnetostatic) assembly and the axial-series interpolation evaluators are ported in
// their own PBIs.
package radial

import (
	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/ring"
)

// NQuad2D is the Gauss-Legendre order used per line element (== geom2d.N_QUAD_2D).
const NQuad2D = geom2d.N_QUAD_2D

// Line is a radial line element: four control points (each a 3D vertex with y == 0).
type Line = [4]geom2d.Vertex

// ExcitationType mirrors the C enum (definitions.c) used to select the assembly kernel.
type ExcitationType uint8

const (
	VoltageFixed     ExcitationType = 1
	VoltageFun       ExcitationType = 2
	Dielectric       ExcitationType = 3
	Current          ExcitationType = 4
	MagnetostaticPot ExcitationType = 5
	Magnetizable     ExcitationType = 6
)

// JacobianBuffer holds the per-element quadrature Jacobians (weight*dl) at the NQuad2D
// Gauss points; PositionBuffer the corresponding (r, z) positions. Both are indexed
// [lineIndex][quadIndex].
type (
	JacobianBuffer [][NQuad2D]float64
	PositionBuffer [][NQuad2D]geom2d.Point2
)

// reorder returns the GMSH line4 control-point remap (v1,v2,v3,v4) = (p0,p2,p3,p1).
func reorder(line Line) (v1, v2, v3, v4 geom2d.Vertex) {
	return line[0], line[2], line[3], line[1]
}

// ElementCenter returns the (r, z) point at the element's parametric midpoint (α=0), used
// as the collocation point for matrix rows and the field-sampling point for the
// magnetostatic right-hand side. Mirrors get_center_of_element (higher-order radial).
func ElementCenter(line Line) geom2d.Point2 {
	v1, v2, v3, v4 := reorder(line)
	_, pos := geom2d.PositionAndJacobianRadial(0, v1, v2, v3, v4)
	return pos
}

// ElementNormal returns the unit (r, z) normal at the element's parametric midpoint (α=0).
func ElementNormal(line Line) geom2d.Point2 {
	v1, v2, v3, v4 := reorder(line)
	return geom2d.HigherOrderNormalRadial(0, v1, v2, v3, v4)
}

// FillJacobianBufferRadial precomputes, for every line element, the quadrature Jacobians
// (weight*dl) and (r, z) positions at the NQuad2D Gauss points. Port of
// fill_jacobian_buffer_radial.
func FillJacobianBufferRadial(lines []Line) (JacobianBuffer, PositionBuffer) {
	jac := make(JacobianBuffer, len(lines))
	pos := make(PositionBuffer, len(lines))
	for i, line := range lines {
		v1, v2, v3, v4 := reorder(line)
		for k := 0; k < NQuad2D; k++ {
			j, p := geom2d.PositionAndJacobianRadial(geom2d.GaussQuadPoints[k], v1, v2, v3, v4)
			jac[i][k] = geom2d.GaussQuadWeights[k] * j
			pos[i][k] = p
		}
	}
	return jac, pos
}

// ChargeRadial returns the total charge on a line element carrying uniform surface charge
// density `charge`, integrating 2π r over the element. Port of charge_radial.
func ChargeRadial(line Line, charge float64) float64 {
	v1, v2, v3, v4 := reorder(line)
	sum := 0.0
	for k := 0; k < NQuad2D; k++ {
		jac, pos := geom2d.PositionAndJacobianRadial(geom2d.GaussQuadPoints[k], v1, v2, v3, v4)
		sum += 2 * piConst * pos[0] * geom2d.GaussQuadWeights[k] * jac * charge
	}
	return sum
}

// piConst is M_PI as used in the C (== math.Pi).
const piConst = 3.14159265358979323846

// PotentialRadial returns the potential at point (x, y, z) due to effective point charges
// `charges` with precomputed buffers. r0 = sqrt(x^2+y^2). Port of potential_radial.
func PotentialRadial(point geom2d.Vertex, charges []float64, jac JacobianBuffer, pos PositionBuffer) float64 {
	r0 := geom2d.Norm2D(point[0], point[1])
	z0 := point[2]
	sum := 0.0
	for i := range charges {
		for k := 0; k < NQuad2D; k++ {
			p := pos[i][k]
			sum += charges[i] * jac[i][k] * ring.PotentialRadialRing(r0, z0, p[0]-r0, p[1]-z0)
		}
	}
	return sum
}

// FieldRadial returns the electric field (Ex, Ey, Ez) at point due to effective point
// charges. On the axis (r < MinDistanceAxis) the radial components are zeroed. Port of
// field_radial.
func FieldRadial(point geom2d.Vertex, charges []float64, jac JacobianBuffer, pos PositionBuffer) geom2d.Vertex {
	r := geom2d.Norm2D(point[0], point[1])
	er, ez := 0.0, 0.0
	for i := range charges {
		for k := 0; k < NQuad2D; k++ {
			p := pos[i][k]
			er -= charges[i] * jac[i][k] * ring.Dr1PotentialRadialRing(r, point[2], p[0]-r, p[1]-point[2])
			ez -= charges[i] * jac[i][k] * ring.Dz1PotentialRadialRing(r, point[2], p[0]-r, p[1]-point[2])
		}
	}
	var result geom2d.Vertex
	if r >= ring.MinDistanceAxis {
		result[0] = point[0] / r * er
		result[1] = point[1] / r * er
	}
	result[2] = ez
	return result
}

// SelfPotentialRadialIntegrand is the α-integrand of the singular self-potential of a line
// element (the diagonal of a voltage row). Port of self_potential_radial; the diagonal is
// the integral of this over α ∈ [-1, 1] (computed by the solver with a singularity-aware
// quadrature, since the integrand has a log singularity at α = 0).
func SelfPotentialRadialIntegrand(alpha float64, line Line) float64 {
	v1, v2, v3, v4 := reorder(line)
	jac, deltaPos := geom2d.DeltaPositionAndJacobianRadial(alpha, v1, v2, v3, v4)
	_, target := geom2d.PositionAndJacobianRadial(0, v1, v2, v3, v4)
	return jac * ring.PotentialRadialRing(target[0], target[1], deltaPos[0], deltaPos[1])
}

// SelfFieldDotNormalRadialIntegrand is the α-integrand of the singular self-field-dot-normal
// of a dielectric line element (the diagonal of a dielectric row, before the −1). K is the
// relative permittivity. Port of self_field_dot_normal_radial.
func SelfFieldDotNormalRadialIntegrand(alpha float64, line Line, k float64) float64 {
	v1, v2, v3, v4 := reorder(line)
	jac, deltaPos := geom2d.DeltaPositionAndJacobianRadial(alpha, v1, v2, v3, v4)
	_, target := geom2d.PositionAndJacobianRadial(0, v1, v2, v3, v4)
	normal := geom2d.HigherOrderNormalRadial(0.0, v1, v2, v3, v4)
	return jac * ring.FieldDotNormalRadial(target[0], target[1], deltaPos[0], deltaPos[1], normal, k)
}

// FillMatrixRadial fills rows [start, end] of the dense BEM influence matrix (row-major,
// N x N) for the given line elements, excitation types and values, and precomputed
// buffers. Voltage/magnetostatic-potential rows accumulate the ring potential; dielectric/
// magnetizable rows accumulate the dielectric-weighted normal field. Port of
// fill_matrix_radial. The singular diagonal is left as the quadrature approximation here
// and overwritten by the solver with the accurate self-integral.
func FillMatrixRadial(matrix []float64, lines []Line, types []ExcitationType, values []float64,
	jac JacobianBuffer, pos PositionBuffer, start, end int) {
	n := len(lines)
	for i := start; i <= end; i++ {
		tv1, tv2, tv3, tv4 := reorder(lines[i])
		_, target := geom2d.PositionAndJacobianRadial(0.0, tv1, tv2, tv3, tv4)

		switch types[i] {
		case VoltageFixed, VoltageFun, MagnetostaticPot:
			for j := 0; j < n; j++ {
				acc := 0.0
				for k := 0; k < NQuad2D; k++ {
					p := pos[j][k]
					acc += jac[j][k] * ring.PotentialRadialRing(target[0], target[1], p[0]-target[0], p[1]-target[1])
				}
				matrix[i*n+j] += acc
			}
		case Dielectric, Magnetizable:
			normal := geom2d.HigherOrderNormalRadial(0.0, tv1, tv2, tv3, tv4)
			for j := 0; j < n; j++ {
				acc := 0.0
				for k := 0; k < NQuad2D; k++ {
					p := pos[j][k]
					acc += jac[j][k] * ring.FieldDotNormalRadial(target[0], target[1], p[0]-target[0], p[1]-target[1], normal, values[i])
				}
				matrix[i*n+j] += acc
			}
		default:
			panic("radial.FillMatrixRadial: unknown excitation type")
		}
	}
}
