// SPDX-License-Identifier: MPL-2.0

package geom2d

import (
	"math"
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/quad"
)

// TestNormal2D ports test_normal_2d: the normal is unit-length and orthogonal to p2-p1.
func TestNormal2D(t *testing.T) {
	p1, p2 := Point2{1.0, -1.0}, Point2{2.0, 3.0}
	n := Normal2D(p1, p2)
	dot := n[0]*(p2[0]-p1[0]) + n[1]*(p2[1]-p1[1])
	oracle.CheckClose(t, "normal·tangent", dot, 0.0)
	oracle.CheckClose(t, "|normal|", math.Hypot(n[0], n[1]), 1.0)
}

// vertexLine builds the four radial-element vertices from r values (z and y zero).
func vertexLine(rs [4]float64) (a, b, c, d Vertex) {
	return Vertex{rs[0], 0, 0}, Vertex{rs[1], 0, 0}, Vertex{rs[2], 0, 0}, Vertex{rs[3], 0, 0}
}

// TestPositionRadial ports test_position_and_jacobian_radial: the cubic element passes
// through its four control points at α = -1, -1+2/3, -1+4/3, 1.
func TestPositionRadial(t *testing.T) {
	v1, v2, v3, v4 := vertexLine([4]float64{0.0, 0.25, 0.5, 1.0})
	alphas := []float64{-1, -1 + 2.0/3, -1 + 4.0/3, 1}
	wantR := []float64{0.0, 0.25, 0.5, 1.0}
	for i, alpha := range alphas {
		_, pos := PositionAndJacobianRadial(alpha, v1, v2, v3, v4)
		oracle.CheckClose(t, "pos.r", pos[0], wantR[i])
		oracle.CheckClose(t, "pos.z", pos[1], 0.0)
	}
}

// TestDeltaRadial ports test_delta_position_and_jacobian_radial: delta == pos - pos(α=0).
func TestDeltaRadial(t *testing.T) {
	v1, v2, v3, v4 := vertexLine([4]float64{0.0, 0.25, 0.5, 1.0})
	_, mid := PositionAndJacobianRadial(0, v1, v2, v3, v4)
	for _, alpha := range []float64{-1, -1 + 2.0/3, -1 + 4.0/3, 1} {
		_, pos := PositionAndJacobianRadial(alpha, v1, v2, v3, v4)
		_, delta := DeltaPositionAndJacobianRadial(alpha, v1, v2, v3, v4)
		oracle.CheckClose(t, "delta.r", delta[0], pos[0]-mid[0])
		oracle.CheckClose(t, "delta.z", delta[1], pos[1]-mid[1])
	}
}

// TestArcLengthViaQuad ports the first half of test_position_and_jacobian_radial_length:
// integrating the Jacobian over α ∈ [-1, 1] recovers the element's arc length (1.5). This
// also cross-exercises core/quad against the same integrand the upstream uses.
func TestArcLengthViaQuad(t *testing.T) {
	v1, v2, v3, v4 := vertexLine([4]float64{0.0, 0.5, 1.0, 1.5})
	length := quad.Adaptive(func(x float64) float64 {
		jac, _ := PositionAndJacobianRadial(x, v1, v2, v3, v4)
		return jac
	}, -1, 1)
	oracle.CheckClose(t, "arc length", length, 1.5)
}
