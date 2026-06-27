// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"strconv"
	"strings"

	"oblikovati.org/api/client"
	"oblikovati.org/api/types"
	"oblikovati.org/api/wire"
)

// TraceonPanelID is the stable dockable-window id the Traceon add-in owns.
const TraceonPanelID = "com.oblikovati.traceon.panel"

// ShowPanel creates (or replaces) the Traceon simulation-parameters dockable window: the
// editable study settings plus a Run button. Edits arrive as panel.valueChanged events
// (applyPanelEdit).
func (e *Engine) ShowPanel() (wire.OKResult, error) {
	e.mu.Lock()
	p := e.params
	e.mu.Unlock()
	return e.api.DockableWindows().Set(wire.DockableWindowSpec{
		ID:      TraceonPanelID,
		Title:   "Traceon Electron Optics",
		Dock:    types.DockRight,
		Visible: true,
		Controls: []wire.PanelControlSpec{
			client.PanelLabel("hdr", "— Simulation parameters —"),
			client.PanelValueEditor("voltage", "Central electrode (V)", strconv.FormatFloat(p.voltage, 'g', -1, 64)),
			client.PanelValueEditor("energy", "Beam energy (eV)", strconv.FormatFloat(p.energyEV, 'g', -1, 64)),
			client.PanelValueEditor("rays", "Number of rays", strconv.Itoa(p.numRays)),
			client.PanelValueEditor("coil_current", "Coil current (A)", strconv.FormatFloat(p.coilCurrent, 'g', -1, 64)),
			client.PanelValueEditor("magnetisation", "Magnet M (A/m)", strconv.FormatFloat(p.magnetisation, 'g', -1, 64)),
			client.PanelValueEditor("permeability", "Iron μr", strconv.FormatFloat(p.permeability, 'g', -1, 64)),
			client.PanelValueEditor("fast_trace", "Fast trace (axial) 0/1", boolField(p.fastTrace)),
			client.PanelSeparator(),
			client.PanelLabel("lens_hdr", "— Parametric lens (no host geometry) —"),
			client.PanelValueEditor("lens", "Lens (host/einzel/cylinder)", string(p.lens)),
			client.PanelValueEditor("lens_radius", "Bore radius (cm)", strconv.FormatFloat(p.lensRadius, 'g', -1, 64)),
			client.PanelValueEditor("lens_thickness", "Electrode thickness (cm)", strconv.FormatFloat(p.lensThickness, 'g', -1, 64)),
			client.PanelValueEditor("lens_spacing", "Electrode spacing (cm)", strconv.FormatFloat(p.lensSpacing, 'g', -1, 64)),
			client.PanelSeparator(),
			client.PanelButton("run", "Run Electron-Optics Study", RunStudyCommandID),
		},
	})
}

// applyPanelEdit writes one edited simulation parameter back into the engine, keyed by control id.
func (e *Engine) applyPanelEdit(controlID, value string) {
	e.mu.Lock()
	defer e.mu.Unlock()
	switch controlID {
	case "voltage":
		e.params.voltage = simNum(value, e.params.voltage)
	case "energy":
		e.params.energyEV = simNum(value, e.params.energyEV)
	case "rays":
		if n := int(simNum(value, float64(e.params.numRays))); n > 0 {
			e.params.numRays = n
		}
	case "coil_current":
		e.params.coilCurrent = simNum(value, e.params.coilCurrent)
	case "magnetisation":
		e.params.magnetisation = simNum(value, e.params.magnetisation)
	case "permeability":
		e.params.permeability = simNum(value, e.params.permeability)
	case "fast_trace":
		e.params.fastTrace = simNum(value, 0) != 0
	case "lens":
		e.params.lens = parseLens(strings.TrimSpace(value))
	case "lens_radius":
		if r := simNum(value, e.params.lensRadius); r > 0 {
			e.params.lensRadius = r
		}
	case "lens_thickness":
		if t := simNum(value, e.params.lensThickness); t > 0 {
			e.params.lensThickness = t
		}
	case "lens_spacing":
		if s := simNum(value, e.params.lensSpacing); s > 0 {
			e.params.lensSpacing = s
		}
	}
}

// boolField formats a boolean as the "0"/"1" the panel shows for a toggle.
func boolField(b bool) string {
	if b {
		return "1"
	}
	return "0"
}

// simNum reads the leading number from a form value (e.g. "1000 V" → 1000), keeping fallback
// when the field is empty or half-typed.
func simNum(value string, fallback float64) float64 {
	fields := strings.Fields(value)
	if len(fields) == 0 {
		return fallback
	}
	v, err := strconv.ParseFloat(fields[0], 64)
	if err != nil {
		return fallback
	}
	return v
}
