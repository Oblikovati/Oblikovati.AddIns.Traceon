// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"fmt"
	"math"
	"strconv"

	"oblikovati.org/api/wire"
)

// sweepClientID is the client-graphics group holding the focus-vs-parameter plot. It is separate
// from the study overlay (graphicsClientID) so a sweep leaves the last study's field + trajectories
// in place and adds its plot alongside, and re-running the sweep replaces only the plot.
const sweepClientID = "com.oblikovati.traceon.sweep"

// minSweepSteps / maxSweepSteps bound a sweep: at least two samples define a curve, and the upper
// cap keeps an interactive sweep (each step recomputes the model and solves a BEM matrix) bounded.
const (
	minSweepSteps = 2
	maxSweepSteps = 64
)

// sweepConfig is a validated parameter-sweep request: sample the host parameter `param` at `steps`
// values evenly spaced over [start, stop], in the parameter's own display unit `unit`.
type sweepConfig struct {
	param string
	unit  string
	start float64
	stop  float64
	steps int
}

// valueAt returns the parameter value at sample i ∈ [0, steps).
func (c sweepConfig) valueAt(i int) float64 {
	if c.steps < 2 {
		return c.start
	}
	return c.start + (c.stop-c.start)*float64(i)/float64(c.steps-1)
}

// expression formats sample value v as a host parameter expression in the parameter's unit
// (e.g. "0.45 cm"); a unit-less parameter gets a bare number.
func (c sweepConfig) expression(v float64) string {
	s := strconv.FormatFloat(v, 'g', -1, 64)
	if c.unit == "" {
		return s
	}
	return s + " " + c.unit
}

// sweepPoint is one sampled (parameter value, axial focus) pair; focusZ is NaN when the beam does
// not cross the optical axis at that parameter value (no focus forms).
type sweepPoint struct {
	value  float64
	focusZ float64
}

// runSweep varies a host parameter over the configured range, recomputing the model and re-solving
// the study at each step, collects the axial focus, renders a focus-vs-parameter plot, and reports
// a one-line summary. The swept parameter is ALWAYS restored to its original expression afterward
// (deferred), so a sweep never leaves the user's model mutated. MUST run off the host session
// goroutine (it issues many host calls).
func (e *Engine) runSweep() (string, error) {
	e.mu.Lock()
	params := e.params
	e.mu.Unlock()
	if params.lens != lensHost {
		return "Traceon: parameter sweep needs host geometry (set the lens field to 'host')", nil
	}
	cfg, err := e.sweepConfig(params)
	if err != nil {
		return "Traceon: " + err.Error(), nil
	}

	orig, err := e.api.Parameters().GetDetail(cfg.param)
	if err != nil {
		return "", fmt.Errorf("read parameter %q: %w", cfg.param, err)
	}
	defer e.restoreParam(cfg.param, orig.Expression)

	points, err := e.collectSweep(cfg, params)
	if err != nil {
		return "", err
	}
	if err := e.pushSweepPlot(cfg, points); err != nil {
		return "", err
	}
	return sweepSummary(cfg, points), nil
}

// defaultSweepFraction is the half-width of the default sweep range, as a fraction of the
// parameter's current value: with no explicit range the sweep spans ±50 % around it. This lets
// "Run Parameter Sweep" work with zero configuration — sweep the model's design parameter around
// where it sits — while the panel range still overrides it.
const defaultSweepFraction = 0.5

// sweepConfig turns the panel sweep settings into a runnable config. The parameter name defaults
// to the model's sole/first user parameter; the range defaults to ±defaultSweepFraction around the
// parameter's current value. Values are in the parameter's DISPLAY unit (what the user types and
// sees), so each sample is set as an unambiguous unit-qualified expression.
func (e *Engine) sweepConfig(params studyParams) (sweepConfig, error) {
	steps := params.sweepSteps
	if steps < minSweepSteps || steps > maxSweepSteps {
		return sweepConfig{}, fmt.Errorf("sweep steps %d out of range [%d, %d]", steps, minSweepSteps, maxSweepSteps)
	}
	name, err := e.resolveSweepParam(params.sweepParam)
	if err != nil {
		return sweepConfig{}, err
	}
	detail, err := e.api.Parameters().GetDetail(name)
	if err != nil {
		return sweepConfig{}, fmt.Errorf("parameter %q not found: %w", name, err)
	}
	start, stop, err := sweepRange(params, detail)
	if err != nil {
		return sweepConfig{}, err
	}
	return sweepConfig{param: name, unit: detail.Units, start: start, stop: stop, steps: steps}, nil
}

// resolveSweepParam returns the panel parameter name, or — when none is set — the model's first
// user parameter, so a one-parameter design sweeps with no configuration.
func (e *Engine) resolveSweepParam(panelName string) (string, error) {
	if panelName != "" {
		return panelName, nil
	}
	list, err := e.api.Parameters().List()
	if err != nil {
		return "", fmt.Errorf("list parameters: %w", err)
	}
	for _, p := range list.Parameters {
		if p.Kind == "user" {
			return p.Name, nil
		}
	}
	return "", fmt.Errorf("no user parameter to sweep — add a model parameter or name one in the panel")
}

// sweepRange returns the configured [start, stop] in the parameter's display unit, or — when the
// panel range is empty — a default ±defaultSweepFraction band around the parameter's current value.
func sweepRange(params studyParams, detail wire.ParameterDetail) (start, stop float64, err error) {
	if params.sweepStop != params.sweepStart {
		return params.sweepStart, params.sweepStop, nil
	}
	v0 := leadingFloat(detail.Value)
	if v0 == 0 {
		return 0, 0, fmt.Errorf("set a sweep range: parameter %q has no nonzero default to sweep around", detail.Name)
	}
	return (1 - defaultSweepFraction) * v0, (1 + defaultSweepFraction) * v0, nil
}

// leadingFloat parses the leading number of a value string like "3 mm" → 3 (0 when absent), so the
// default sweep range can be built from a parameter's display value.
func leadingFloat(value string) float64 {
	return simNum(value, 0)
}

// collectSweep runs the study at each sampled parameter value and returns the (value, focus) points.
func (e *Engine) collectSweep(cfg sweepConfig, params studyParams) ([]sweepPoint, error) {
	points := make([]sweepPoint, 0, cfg.steps)
	for i := 0; i < cfg.steps; i++ {
		v := cfg.valueAt(i)
		focus, err := e.focusAtParam(cfg, v, params)
		if err != nil {
			return nil, err
		}
		points = append(points, sweepPoint{value: v, focusZ: focus})
	}
	return points, nil
}

// focusAtParam sets the swept parameter to v, recomputes the model, re-solves the study, and
// returns the axial focus (cm, NaN if the beam does not cross the axis). It does not push the
// study overlay — the sweep renders its own plot.
func (e *Engine) focusAtParam(cfg sweepConfig, v float64, params studyParams) (float64, error) {
	if _, err := e.api.Parameters().Set(wire.ParameterSetArgs{Name: cfg.param, Expression: cfg.expression(v)}); err != nil {
		return 0, fmt.Errorf("set %s = %s: %w", cfg.param, cfg.expression(v), err)
	}
	if _, err := e.api.Documents().Update(true); err != nil {
		return 0, fmt.Errorf("recompute at %s = %s: %w", cfg.param, cfg.expression(v), err)
	}
	res, _, err := e.computeStudy(params)
	if err != nil {
		return 0, fmt.Errorf("study at %s = %s: %w", cfg.param, cfg.expression(v), err)
	}
	return res.FocusZ, nil
}

// restoreParam puts a swept parameter back to its original expression and recomputes, so the sweep
// leaves the model exactly as it found it. Errors are swallowed: this runs deferred on the cleanup
// path, and a failure to restore must not mask the sweep's own outcome.
func (e *Engine) restoreParam(name, expression string) {
	if _, err := e.api.Parameters().Set(wire.ParameterSetArgs{Name: name, Expression: expression}); err != nil {
		return
	}
	_, _ = e.api.Documents().Update(true)
}

// sweepSummary reports the swept range and the focus span over it (ignoring values where no focus
// formed), so the status bar shows the sweep's outcome at a glance.
func sweepSummary(cfg sweepConfig, points []sweepPoint) string {
	lo, hi, found := focusSpan(points)
	if !found {
		return fmt.Sprintf("Traceon: swept %s over [%g, %g] %s in %d steps — beam never focused",
			cfg.param, cfg.start, cfg.stop, cfg.unit, cfg.steps)
	}
	return fmt.Sprintf("Traceon: swept %s over [%g, %g] %s in %d steps — focus z = %.3f … %.3f cm",
		cfg.param, cfg.start, cfg.stop, cfg.unit, cfg.steps, lo, hi)
}

// focusSpan returns the min and max finite focus across the sweep, and whether any was finite.
func focusSpan(points []sweepPoint) (lo, hi float64, found bool) {
	lo, hi = math.Inf(1), math.Inf(-1)
	for _, p := range points {
		if math.IsNaN(p.focusZ) || math.IsInf(p.focusZ, 0) {
			continue
		}
		lo, hi, found = math.Min(lo, p.focusZ), math.Max(hi, p.focusZ), true
	}
	return lo, hi, found
}
