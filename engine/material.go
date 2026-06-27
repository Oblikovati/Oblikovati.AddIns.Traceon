// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"oblikovati.org/api/types"
	"oblikovati.org/api/wire"

	"oblikovati.org/traceon/core/constants"
)

// materialRole classifies a body by its assigned material's magnetic constitutive class, so a host
// model can drive the simulation by MATERIAL rather than by electrode naming: assign a soft-iron
// material to a pole piece and it becomes magnetizable iron; assign a NdFeB material and it becomes
// a permanent magnet. A soft-magnetic material yields ("iron", its relative permeability); a
// hard-magnetic material yields ("magnet", its magnetisation M = Br/μ0 in A/m); a non-magnetic or
// unassigned material yields no role, so the body falls through to an electrode.
func (e *Engine) materialRole(materialID string) (role string, value float64, ok bool) {
	if materialID == "" {
		return "", 0, false
	}
	mat, err := e.api.Materials().Get(materialID)
	if err != nil {
		return "", 0, false
	}
	switch mat.Magnetic.Class {
	case types.SoftMagnetic:
		if mu := mat.Magnetic.RelativePermeability; mu > 1 {
			return "iron", mu, true
		}
	case types.HardMagnetic:
		if br := mat.Magnetic.Remanence; br != 0 {
			return "magnet", br / constants.VacuumPermeability, true
		}
	}
	return "", 0, false
}

// classifyByMaterial sorts a body into magnetizable iron or a permanent magnet by its assigned
// material's magnetic class, returning whether it matched (so a material-classified body is
// excluded from the electrode default). Name/attribute conventions take precedence; this is the
// fallback that reads the host material.
func (e *Engine) classifyByMaterial(b wire.BodyInfo, magnets *[]magnet, irons *[]iron, prof *profile) bool {
	role, value, ok := e.materialRole(b.MaterialID)
	if !ok {
		return false
	}
	switch role {
	case "iron":
		*irons = append(*irons, iron{prof: prof, permeability: value})
	case "magnet":
		*magnets = append(*magnets, magnet{prof: prof, magnetisation: value})
	}
	return true
}
