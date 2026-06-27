// SPDX-License-Identifier: MPL-2.0

// Package solver assembles and solves the radial BEM linear system, porting the matrix
// build + solve of traceon/solver.py (SolverRadial / ElectrostaticSolverRadial). It fills
// the influence matrix with radial.FillMatrixRadial, overwrites the diagonal with the
// accurate singular self-term integrals, builds the right-hand side from the excitations,
// and solves for the effective point charges via core/linalg.
//
// This file covers the electrostatic path (voltage + dielectric). The magnetostatic /
// current-coil right-hand side (which needs a pre-existing field) is added in its PBI.
package solver

import (
	"oblikovati.org/traceon/core/linalg"
	"oblikovati.org/traceon/core/quad"
	"oblikovati.org/traceon/core/radial"
)

// selfTermAbsTol / selfTermRelTol match the scipy.integrate.quad tolerances the upstream
// uses for the singular self-term integrals (quad(..., epsabs=1e-9, epsrel=1e-9)).
const (
	selfTermAbsTol = 1e-9
	selfTermRelTol = 1e-9
)

// EffectivePointCharges is the solved surface-charge distribution: one charge per line
// element plus the precomputed quadrature buffers needed to evaluate the field it produces.
type EffectivePointCharges struct {
	Charges []float64
	Jac     radial.JacobianBuffer
	Pos     radial.PositionBuffer
}

// SelfPotentialRadial returns the diagonal entry for a voltage element: the singular
// self-potential, ∫_{-1}^{1} of the self-potential integrand over α, with the log
// singularity at α=0 handled by splitting there. Port of the Python self_potential_radial
// wrapper (quad with points=(0,)).
func SelfPotentialRadial(line radial.Line) float64 {
	return quad.IntegrateWithSingularities(
		func(alpha float64) float64 { return radial.SelfPotentialRadialIntegrand(alpha, line) },
		-1, 1, []float64{0}, selfTermAbsTol, selfTermRelTol)
}

// SelfFieldDotNormalRadial returns the integrated self-field-dot-normal for a dielectric/
// magnetizable element with relative permittivity K. Port of the Python
// self_field_dot_normal_radial wrapper.
func SelfFieldDotNormalRadial(line radial.Line, k float64) float64 {
	return quad.IntegrateWithSingularities(
		func(alpha float64) float64 { return radial.SelfFieldDotNormalRadialIntegrand(alpha, line, k) },
		-1, 1, []float64{0}, selfTermAbsTol, selfTermRelTol)
}

// AssembleMatrix builds the N×N influence matrix: the off-diagonal (and a rough diagonal)
// from radial.FillMatrixRadial, then the accurate singular self-term on the diagonal —
// self_potential for voltage rows, (self_field_dot_normal − 1) for dielectric rows. Mirrors
// SolverRadial.get_matrix.
func AssembleMatrix(lines []radial.Line, types []radial.ExcitationType, values []float64,
	jac radial.JacobianBuffer, pos radial.PositionBuffer) linalg.Matrix {
	n := len(lines)
	m := linalg.NewMatrix(n, n)
	radial.FillMatrixRadial(m.Data, lines, types, values, jac, pos, 0, n-1)
	for i := 0; i < n; i++ {
		switch types[i] {
		case radial.Dielectric, radial.Magnetizable:
			m.Set(i, i, SelfFieldDotNormalRadial(lines[i], values[i])-1) // −1 from the matrix equation
		default:
			m.Set(i, i, SelfPotentialRadial(lines[i]))
		}
	}
	return m
}

// rightHandSideElectrostatic builds F for the electrostatic system: the prescribed voltage
// on voltage rows, 0 on dielectric rows. Mirrors Solver.get_right_hand_side (electrostatic).
func rightHandSideElectrostatic(types []radial.ExcitationType, values []float64) []float64 {
	f := make([]float64, len(types))
	for i, t := range types {
		switch t {
		case radial.VoltageFixed, radial.VoltageFun, radial.MagnetostaticPot:
			f[i] = values[i]
		case radial.Dielectric:
			f[i] = 0
		default:
			f[i] = 0
		}
	}
	return f
}

// SolveElectrostatic assembles and solves the electrostatic radial BEM system, returning
// the effective point charges. Mirrors ElectrostaticSolverRadial.solve_matrix for a single
// right-hand side. types/values are per-element (VoltageFixed/VoltageFun/Dielectric).
func SolveElectrostatic(lines []radial.Line, types []radial.ExcitationType, values []float64) (EffectivePointCharges, error) {
	jac, pos := radial.FillJacobianBufferRadial(lines)
	if len(lines) == 0 {
		return EffectivePointCharges{Charges: nil, Jac: jac, Pos: pos}, nil
	}
	matrix := AssembleMatrix(lines, types, values, jac, pos)
	f := rightHandSideElectrostatic(types, values)
	charges, err := linalg.SolveVector(matrix, f)
	if err != nil {
		return EffectivePointCharges{}, err
	}
	return EffectivePointCharges{Charges: charges, Jac: jac, Pos: pos}, nil
}
