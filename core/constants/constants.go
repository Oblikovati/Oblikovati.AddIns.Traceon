// SPDX-License-Identifier: MPL-2.0

// Package constants holds the physical constants the Traceon numerical core uses,
// matched bit-for-bit to the CODATA values in scipy.constants so the Go port stays
// numerically equivalent to the upstream Python (which imports them from
// scipy.constants — see solver.py, tracing.py, excitation.py).
package constants

// CODATA values, identical to scipy.constants (verified against the oracle venv).
const (
	// ElementaryCharge is the elementary charge e in coulombs (scipy.constants.e).
	ElementaryCharge = 1.602176634e-19
	// ElectronMass is the electron rest mass m_e in kilograms (scipy.constants.m_e).
	ElectronMass = 9.1093837139e-31
	// VacuumPermeability is mu_0 in henry/metre (scipy.constants.mu_0).
	VacuumPermeability = 1.25663706127e-06
	// VacuumPermittivity is epsilon_0 in farad/metre (scipy.constants.epsilon_0).
	VacuumPermittivity = 8.8541878188e-12
	// SpeedOfLight is c in metres/second (scipy.constants.c).
	SpeedOfLight = 299792458.0
)
