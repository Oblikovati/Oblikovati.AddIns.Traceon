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

## Live test (MCP bridge)

The parametric lenses are also exposed as commands, so the whole study runs live in the app
driven over the MCP bridge — no host model required:

```sh
# 1. build + install the add-in next to the MCP bridge, then launch the head
make -C .. install ADDINS_DIR=../Oblikovati/head/addins
( cd ../Oblikovati/head && make run )            # GUI on $DISPLAY; MCP bridge on 127.0.0.1:7800

# 2. drive it (from the MCPBridge repo)
go run ./cmd/mcpdrive call create_document  '{"type":"part","name":"Lens"}'
go run ./cmd/mcpdrive call execute_command  '{"id":"Traceon.RunEinzelLens"}'
go run ./cmd/mcpdrive call status_get_text  '{}'      # -> "3 electrode(s), 216 elements, 7 rays — focus z = 4.460 cm"
go run ./cmd/mcpdrive call capture_viewport '{"path":"shot.png"}'
```

`images/live-traceon-panel.png` shows the parametric-lens panel in the host, and
`images/live-einzel-mcp.png` shows the resulting overlay — the three einzel apertures, the
potential map (red at the biased centre electrode), and the seven electrons converging to the
focus — all built parametrically and rendered live.
