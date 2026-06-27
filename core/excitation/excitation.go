// SPDX-License-Identifier: MPL-2.0

// Package excitation assigns excitation types (fixed/functional voltage, dielectric,
// electrostatic boundary) to the named electrodes of a radial mesh and emits the per-element
// solver inputs. It is the radial-symmetric port of traceon/excitation.py — the bridge
// between a meshed geometry (core/geometry, core/mesh) and the BEM solve (core/solver).
//
// Scope: the electrostatic excitations (voltage and dielectric/boundary). The magnetostatic
// excitations (current, magnetizable, permanent magnet) are applied by the engine directly
// from host geometry today; they are added here when that path is unified.
package excitation

import (
	"fmt"

	"oblikovati.org/traceon/core/geometry"
	"oblikovati.org/traceon/core/mesh"
	"oblikovati.org/traceon/core/radial"
)

// VoltageFunc gives the prescribed voltage at a point on an electrode. In radial symmetry it
// is evaluated at the element centre (x = r, y = 0, z), matching traceon's get_center_of_element
// convention — so a z-dependent ramp (e.g. a gap voltage) is sampled per element.
type VoltageFunc func(x, y, z float64) float64

// entry is one electrode's assigned excitation: its type plus the fixed value (volts for a
// voltage, relative permittivity for a dielectric) or, for VoltageFun, the position function.
type entry struct {
	typ   radial.ExcitationType
	value float64
	fn    VoltageFunc
}

// Excitation maps electrode names (the mesh's physical groups) to their excitation. Build it
// from a meshed geometry, assign voltages/boundaries by name, then call Electrostatic to get
// the solver inputs. Mirrors traceon.excitation.Excitation for RADIAL symmetry.
type Excitation struct {
	m      *mesh.Mesh
	byName map[string]entry
}

// New creates an excitation over the radial mesh m. Electrode names assigned later must be
// physical groups present in m (m.PhysicalToLines).
func New(m *mesh.Mesh) *Excitation {
	return &Excitation{m: m, byName: map[string]entry{}}
}

// AddVoltage assigns a fixed voltage (volts) to the named electrode. Port of
// add_voltage with a scalar value. Panics if name is not an electrode in the mesh.
func (e *Excitation) AddVoltage(name string, volts float64) {
	e.assign(name, entry{typ: radial.VoltageFixed, value: volts})
}

// AddVoltageFunc assigns a position-dependent voltage to the named electrode; fn is evaluated
// at each element's centre. Port of add_voltage with a callable value. Panics if name is not
// an electrode in the mesh, or if fn is nil.
func (e *Excitation) AddVoltageFunc(name string, fn VoltageFunc) {
	if fn == nil {
		panic(fmt.Sprintf("excitation: nil voltage function for electrode %q", name))
	}
	e.assign(name, entry{typ: radial.VoltageFun, fn: fn})
}

// AddDielectric assigns a relative permittivity to the named electrode. Port of add_dielectric.
// Panics if name is not an electrode in the mesh.
func (e *Excitation) AddDielectric(name string, permittivity float64) {
	e.assign(name, entry{typ: radial.Dielectric, value: permittivity})
}

// AddElectrostaticBoundary marks the named electrodes as electrostatic boundaries (E·n = 0),
// implemented — as upstream — as a dielectric with relative permittivity zero. Each boundary's
// normals are first made inward-pointing (the orientation the boundary kernel assumes), then
// the dielectric is assigned. Placing a boundary around the modelled region greatly improves
// the solve's conditioning. Port of add_electrostatic_boundary (ensure_inward_normals=True).
// Panics if any name is not an electrode in the mesh.
func (e *Excitation) AddElectrostaticBoundary(names ...string) {
	for _, name := range names {
		e.m.EnsureInwardNormals(name)
		e.AddDielectric(name, 0)
	}
}

// Electrodes returns the electrode names present in the mesh (its physical groups).
func (e *Excitation) Electrodes() []string {
	names := make([]string, 0, len(e.m.PhysicalToLines))
	for name := range e.m.PhysicalToLines {
		names = append(names, name)
	}
	return names
}

// Electrostatic returns the active electrostatic elements as solver inputs (lines, per-element
// excitation types, per-element values), in mesh-line order. Functional voltages are evaluated
// at each element's centre. Elements whose electrode carries no electrostatic excitation are
// omitted. Feed the result straight into solver.SolveElectrostatic. Mirrors
// ElectrostaticSolverRadial's active-element + right-hand-side construction.
func (e *Excitation) Electrostatic() ([]radial.Line, []radial.ExcitationType, []float64) {
	all := geometry.RadialLines(e.m) // 1:1 with e.m.Lines, preserving order
	owner := e.lineOwners()

	var (
		lines  []radial.Line
		types  []radial.ExcitationType
		values []float64
	)
	for i, line := range all {
		ent, ok := e.byName[owner[i]]
		if !ok || !isElectrostatic(ent.typ) {
			continue
		}
		v := ent.value
		if ent.typ == radial.VoltageFun {
			c := radial.ElementCenter(line) // (r, z)
			v = ent.fn(c[0], 0, c[1])
		}
		lines = append(lines, line)
		types = append(types, ent.typ)
		values = append(values, v)
	}
	return lines, types, values
}

// assign records an electrode's excitation, validating the name against the mesh first.
func (e *Excitation) assign(name string, ent entry) {
	if _, ok := e.m.PhysicalToLines[name]; !ok {
		panic(fmt.Sprintf("excitation: %q is not an electrode in the mesh; have %v", name, e.Electrodes()))
	}
	e.byName[name] = ent
}

// lineOwners maps each mesh-line index to the electrode (physical group) it belongs to.
func (e *Excitation) lineOwners() map[int]string {
	owner := make(map[int]string, len(e.m.Lines))
	for name, idxs := range e.m.PhysicalToLines {
		for _, i := range idxs {
			owner[i] = name
		}
	}
	return owner
}

// isElectrostatic reports whether an excitation type contributes to the electrostatic solve.
// Mirrors ExcitationType.is_electrostatic.
func isElectrostatic(t radial.ExcitationType) bool {
	return t == radial.VoltageFixed || t == radial.VoltageFun || t == radial.Dielectric
}
