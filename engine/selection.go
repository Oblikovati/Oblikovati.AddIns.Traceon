// SPDX-License-Identifier: MPL-2.0

package engine

import (
	"encoding/base64"
	"fmt"
	"strings"

	"oblikovati.org/api/types"
)

// attrRole / attrValue name the per-body excitation a viewport pick assigns: a role
// (electrode/coil/magnet/iron) and its value (volts / amperes / A·m⁻¹ / μr), anchored to the
// body's reference key so the tag survives recompute.
const (
	attrRole  = "role"
	attrValue = "value"
)

// assignRoles are the excitation roles the panel offers; any other text falls back to electrode.
var assignRoles = map[string]bool{"electrode": true, "coil": true, "magnet": true, "iron": true}

// normalizeRef converts a viewport selection ref to the raw reference-key form that
// model.referenceKeys returns. The selection canonicalises a picked entity to "<kind>/<base64url
// of the reference key>" (e.g. "face/A1Jl…"), whereas referenceKeys returns the raw key bytes; so
// strip the kind prefix and base64-decode to match. A ref without that shape passes through.
func normalizeRef(ref string) string {
	i := strings.IndexByte(ref, '/')
	if i < 0 {
		return ref
	}
	raw, err := base64.RawURLEncoding.DecodeString(ref[i+1:])
	if err != nil {
		return ref
	}
	return string(raw)
}

// faceToBody maps every face/edge/vertex reference key (and each body's own key) to the reference
// key of the body that owns it, so a viewport pick — which selects faces — resolves to the body
// an excitation is assigned to. ReferenceKeys() lists topology per body in body.list order.
func (e *Engine) faceToBody() (map[string]string, error) {
	bodies, err := e.api.Body().List()
	if err != nil {
		return nil, fmt.Errorf("list bodies: %w", err)
	}
	refs, err := e.api.Model().ReferenceKeys()
	if err != nil {
		return nil, fmt.Errorf("reference keys: %w", err)
	}
	owner := map[string]string{}
	for i, bt := range refs.Bodies {
		if i >= len(bodies.Bodies) {
			break
		}
		key := bodies.Bodies[i].Key
		owner[key] = key // selecting the body itself resolves to the body
		for _, f := range bt.Faces {
			owner[f.Key] = key
		}
		for _, ed := range bt.Edges {
			owner[ed.Key] = key
		}
		for _, vx := range bt.Vertices {
			owner[vx.Key] = key
		}
	}
	return owner, nil
}

// assignToSelection tags the body owning each currently-selected face with the panel's role +
// value, stored as per-body attributes (anchored to the body's reference key). Picking any face
// of an electrode and clicking Assign therefore sets that whole electrode's excitation. Returns a
// one-line status for the host status bar. MUST run off the host session goroutine (it calls the
// host).
func (e *Engine) assignToSelection() (string, error) {
	docID, ok := e.activeDocID()
	if !ok {
		return "Traceon: no active document", nil
	}
	sel, err := e.api.Model().Selection()
	if err != nil {
		return "", fmt.Errorf("read selection: %w", err)
	}
	if len(sel.Refs) == 0 {
		return "Traceon: select a face (or body), then Assign", nil
	}
	owner, err := e.faceToBody()
	if err != nil {
		return "", err
	}

	e.mu.Lock()
	role, value := e.params.assignRole, e.params.assignValue
	e.mu.Unlock()
	if !assignRoles[role] {
		role = "electrode"
	}

	bodies := map[string]bool{}
	for _, ref := range sel.Refs {
		if key, ok := owner[normalizeRef(ref)]; ok {
			bodies[key] = true
		}
	}
	for key := range bodies {
		if _, err := e.api.Attributes().SetOn(docID, key, attrSet, attrRole, types.StringVariant(role)); err != nil {
			return "", fmt.Errorf("assign role: %w", err)
		}
		if _, err := e.api.Attributes().SetOn(docID, key, attrSet, attrValue, types.DoubleVariant(value)); err != nil {
			return "", fmt.Errorf("assign value: %w", err)
		}
	}
	return fmt.Sprintf("Traceon: assigned %s = %g to %d body(ies)", role, value, len(bodies)), nil
}

// bodyAssignment reads a body's viewport-assigned excitation (role + value) from its per-body
// attributes, or ok=false when the body was never assigned. This is how a viewport pick drives
// classification, overriding the name/material conventions.
func (e *Engine) bodyAssignment(docID uint64, bodyKey string) (role string, value float64, ok bool) {
	r, err := e.api.Attributes().GetOn(docID, bodyKey, attrSet, attrRole)
	if err != nil || !r.Found {
		return "", 0, false
	}
	role, _ = r.Attribute.Value.Str()
	if role == "" {
		return "", 0, false
	}
	v, err := e.api.Attributes().GetOn(docID, bodyKey, attrSet, attrValue)
	if err == nil && v.Found {
		value, _ = v.Attribute.Value.Double()
	}
	return role, value, true
}

// applyAssignment sorts an assigned body into the right excitation set by its role, returning
// whether it was assigned (so it bypasses the name/material/default classification).
func applyAssignment(role string, value float64, prof *profile, coils *[]coil, magnets *[]magnet, irons *[]iron, electrodes *[]electrode) {
	switch role {
	case "coil":
		*coils = append(*coils, coil{prof: prof, current: value})
	case "magnet":
		*magnets = append(*magnets, magnet{prof: prof, magnetisation: value})
	case "iron":
		*irons = append(*irons, iron{prof: prof, permeability: value})
	default: // "electrode"
		*electrodes = append(*electrodes, electrode{prof: prof, voltage: value})
	}
}
