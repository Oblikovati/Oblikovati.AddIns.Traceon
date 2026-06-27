# Traceon examples

Runnable Go programs that drive the pure-Go `core/` packages the same way upstream Traceon's
Python examples drive the library — build a lens, solve the boundary-element field, trace
electrons, and render the result to `../images/`. They double as living documentation of the
public `core/` API.

Run from the repository root:

```sh
go run ./examples/einzel-lens   # -> images/einzel-lens.png
go run ./examples/dohi-mirror   # -> images/dohi-mirror.png
```

| Example | What it shows | Output |
|---|---|---|
| `einzel-lens` | A three-aperture einzel lens (centre at 1.8 kV) focusing a 1 keV electron bundle. | `images/einzel-lens.png` |
| `dohi-mirror` | The Dohi electron mirror (−1.25 kV) reflecting electrons back along the axis. | `images/dohi-mirror.png` |

Each program is self-contained: geometry via `core/geometry`, excitation via `core/excitation`,
the BEM solve via `core/solver`, the fast axial field via `core/field`, and tracing via
`core/tracing`. Figures are drawn by the dependency-free `examples/plot` helper (the optical
axis runs left→right; the potential is a blue→white→red map; electrodes are outlined and the
electron trajectories are green).

The same lenses can be built and run **interactively in the host** via the add-in's parametric
lens panel (no host geometry needed) — see the engine's `Parametric lens` panel section.
