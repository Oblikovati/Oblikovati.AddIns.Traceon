// SPDX-License-Identifier: MPL-2.0

package geometry

import (
	"testing"

	"oblikovati.org/traceon/core/internal/oracle"
	"oblikovati.org/traceon/core/mesh"
)

// mesherGolden mirrors the JSON the fixture generator emits for the parametric mesher.
type mesherGolden struct {
	Cases []struct {
		Name   string       `json:"name"`
		Points [][]oracle.F `json:"points"`
		Lines  [][]int      `json:"lines"`
	} `json:"cases"`
	Discretize struct {
		Basic        []oracle.F `json:"basic"`
		Line4Factor3 []oracle.F `json:"line4_factor3"`
	} `json:"discretize"`
}

// buildCase reconstructs the Go Path + mesh for a fixture case by name, exactly mirroring
// the upstream construction in tools/gen_fixtures.py::_mesher.
func buildCase(name string) *mesh.Mesh {
	switch name {
	case "line_flat":
		return Line(Point{1, 0, 0}, Point{1, 0, 2}).Mesh(MeshOptions{MeshSize: 0.5, EnsureOutwardNormals: true})
	case "line_line4":
		return Line(Point{1, 0, 0}, Point{1, 0, 2}).Mesh(MeshOptions{MeshSize: 2.0, HigherOrder: true, EnsureOutwardNormals: true})
	case "arc_flat":
		return Arc(Point{0, 0, 0}, Point{2, 0, 0}, Point{0, 0, 2}, false).Mesh(MeshOptions{MeshSize: 1.0, EnsureOutwardNormals: true})
	case "rect_named_line4":
		return RectangleXZ(0.5, 1.0, -0.5, 0.5).Mesh(MeshOptions{MeshSize: 0.5, HigherOrder: true, Name: "rect", EnsureOutwardNormals: true})
	case "aperture_flat":
		return Aperture(0.5, 0.3, 1.5, 0.0).Mesh(MeshOptions{MeshSize: 0.4, EnsureOutwardNormals: true})
	default:
		return nil
	}
}

// TestMesherGolden verifies the parametric mesher reproduces upstream Path.mesh: the same
// post-deduplication point coordinates (ordering included) and the exact integer line
// connectivity, across straight/curved lines, arcs, a named rectangle (with outward-normal
// orientation), and an aperture.
func TestMesherGolden(t *testing.T) {
	var fx mesherGolden
	oracle.LoadGolden(t, "mesher", &fx)

	for _, c := range fx.Cases {
		t.Run(c.Name, func(t *testing.T) {
			m := buildCase(c.Name)
			if m == nil {
				t.Fatalf("no Go builder for case %q", c.Name)
			}
			if len(m.Points) != len(c.Points) {
				t.Fatalf("point count = %d, want %d", len(m.Points), len(c.Points))
			}
			for i, p := range c.Points {
				for k := 0; k < 3; k++ {
					oracle.CheckClose(t, c.Name+" point", m.Points[i][k], p[k].Float())
				}
			}
			if len(m.Lines) != len(c.Lines) {
				t.Fatalf("line count = %d, want %d", len(m.Lines), len(c.Lines))
			}
			for i, want := range c.Lines {
				if len(m.Lines[i]) != len(want) {
					t.Fatalf("line %d arity = %d, want %d", i, len(m.Lines[i]), len(want))
				}
				for j, idx := range want {
					if m.Lines[i][j] != idx {
						t.Errorf("line %d = %v, want %v", i, m.Lines[i], want)
						break
					}
				}
			}
		})
	}
}

// TestDiscretizeGolden verifies discretizePath reproduces upstream sample parameters.
func TestDiscretizeGolden(t *testing.T) {
	var fx mesherGolden
	oracle.LoadGolden(t, "mesher", &fx)

	basic := discretizePath(10, []float64{3.33, 5, 9}, 1.0, 0, 1)
	checkParams(t, "basic", basic, fx.Discretize.Basic)

	line4 := discretizePath(2, nil, 2.0, 0, 3)
	checkParams(t, "line4_factor3", line4, fx.Discretize.Line4Factor3)
}

func checkParams(t *testing.T, label string, got []float64, want []oracle.F) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: len = %d, want %d", label, len(got), len(want))
	}
	for i := range want {
		oracle.CheckClose(t, label, got[i], want[i].Float())
	}
}
