// SPDX-License-Identifier: MPL-2.0

package constants

import "testing"

// TestConstantsMatchScipy pins each constant to the scipy.constants value the upstream
// Traceon imports. These are exact (CODATA) literals, so equality must be bit-exact —
// any drift here silently shifts every downstream field/trajectory result.
func TestConstantsMatchScipy(t *testing.T) {
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"ElementaryCharge", ElementaryCharge, 1.602176634e-19},
		{"ElectronMass", ElectronMass, 9.1093837139e-31},
		{"VacuumPermeability", VacuumPermeability, 1.25663706127e-06},
		{"VacuumPermittivity", VacuumPermittivity, 8.8541878188e-12},
		{"SpeedOfLight", SpeedOfLight, 299792458.0},
	}
	for _, c := range cases {
		if c.got != c.want {
			t.Errorf("%s = %g, want %g (scipy.constants)", c.name, c.got, c.want)
		}
	}
}
