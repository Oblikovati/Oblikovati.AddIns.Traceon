// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"fmt"

	"oblikovati.org/traceon/core/geom2d"
	"oblikovati.org/traceon/core/geometry"
)

// paramLens selects how a study's electrodes are defined. lensHost sections the host part
// (the default); the others build the electrodes parametrically via core/geometry, so a study
// runs with no host model at all.
type paramLens string

const (
	lensHost     paramLens = "host"     // section the live host geometry
	lensEinzel   paramLens = "einzel"   // three-aperture einzel lens (outer grounded, centre biased)
	lensCylinder paramLens = "cylinder" // two coaxial cylinders at 0 V and the bias voltage
)

// parseLens maps a panel string to a lens mode, falling back to host on anything unrecognised.
func parseLens(s string) paramLens {
	switch paramLens(s) {
	case lensEinzel:
		return lensEinzel
	case lensCylinder:
		return lensCylinder
	default:
		return lensHost
	}
}

// lensMeshFactor is the mesh-size factor used to discretize parametric electrodes into BEM
// line elements — enough elements for an accurate paraxial solve without a heavy matrix.
const lensMeshFactor = 8

// buildParametricLens builds the electrodes for the selected parametric lens template from the
// panel dimensions (cm). The central/biased electrode carries params.voltage; ground electrodes
// are at 0 V. Each electrode's BEM elements come from the core/geometry parametric mesher.
func buildParametricLens(params studyParams) ([]electrode, error) {
	switch params.lens {
	case lensEinzel:
		return einzelElectrodes(params), nil
	case lensCylinder:
		return cylinderElectrodes(params), nil
	default:
		return nil, fmt.Errorf("unknown parametric lens %q", params.lens)
	}
}

// einzelElectrodes builds the classic three-aperture einzel lens: two grounded outer apertures
// straddling a central aperture biased to params.voltage, spaced by lensThickness+lensSpacing.
func einzelElectrodes(params studyParams) []electrode {
	t, sp, r := params.lensThickness, params.lensSpacing, params.lensRadius
	extent := r + 3*t // outer wall well clear of the bore so the field is the aperture's
	step := t + sp

	return []electrode{
		paramElectrode(geometry.Aperture(t, r, extent, -step), 0),
		paramElectrode(geometry.Aperture(t, r, extent, 0), params.voltage),
		paramElectrode(geometry.Aperture(t, r, extent, step), 0),
	}
}

// cylinderElectrodes builds a two-cylinder immersion lens: two coaxial tubes of bore radius r
// separated by a small gap, the lower at 0 V and the upper at params.voltage.
func cylinderElectrodes(params studyParams) []electrode {
	t, sp, r := params.lensThickness, params.lensSpacing, params.lensRadius
	gap := sp
	length := 4 * t // each cylinder's axial length

	lower := geometry.Line(geometry.Point{r, 0, -gap/2 - length}, geometry.Point{r, 0, -gap / 2})
	upper := geometry.Line(geometry.Point{r, 0, gap / 2}, geometry.Point{r, 0, gap/2 + length})

	return []electrode{
		paramElectrode(lower, 0),
		paramElectrode(upper, params.voltage),
	}
}

// renderSamples is the number of points sampled along a parametric path for the render profile.
const renderSamples = 48

// paramElectrode meshes a parametric path into BEM line elements (the core/geometry mesher) and
// samples it into a render profile, returning a study electrode at the given voltage. All in cm.
func paramElectrode(path geometry.Path, voltage float64) electrode {
	mesh := path.Mesh(geometry.MeshOptions{MeshSizeFactor: lensMeshFactor, HigherOrder: true})
	return electrode{
		prof:    pathProfile(path),
		voltage: voltage,
		lines:   geometry.RadialLines(mesh),
	}
}

// pathProfile samples a path into an (r, z) polyline (cm) for the viewport overlay and the
// study extent. The path lies in the meridian plane, so r = x and z = z of each sample.
func pathProfile(path geometry.Path) *profile {
	loop := make([]geom2d.Point2, renderSamples)
	for i := 0; i < renderSamples; i++ {
		t := path.Length * float64(i) / float64(renderSamples-1)
		p := path.At(t)
		loop[i] = geom2d.Point2{p[0], p[2]}
	}
	return &profile{loops: [][]geom2d.Point2{loop}}
}
