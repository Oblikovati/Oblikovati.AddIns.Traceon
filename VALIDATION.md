# Validation

Two layers of validation back the Traceon port: per-module oracle tests (every `core/`
package, checked against the upstream Traceon backend at `np.isclose` tolerances) and the
**end-to-end** integration test below, which runs the whole pure-Go stack on a complete
electron-optics problem and asserts it reproduces Traceon's result.

## End-to-end: einzel-lens focal length

`validation/einzel_test.go` builds a three-electrode einzel lens (two grounded outer tubes,
a 5 kV centre tube) directly as (r, z) BEM elements, then runs the **entire** pipeline —
charge solve → field reconstruction → RKF45 ray tracing → least-squares focus — and compares
the focal point to Traceon computing the identical case.

Result: the Go stack lands the beam on the axis at

```
focus z = 8.190615   (Traceon: 8.190615)
```

— agreement to ~1e-6, and the full 1240-step trajectory of the outermost ray matches Traceon
step for step. Because the focal point is a sensitive derived quantity, this single number
exercises the solver, field, tracer, and focus together.

Regenerate the golden + plot from the upstream clone:

```
../Traceon/.venv/bin/python tools/gen_einzel.py --traceon-dir ../Traceon
go test ./validation/...
```

### Visual confirmation

`validation/testdata/einzel.png` (regenerated above) shows the parallel beam entering from
the left, passing the three electrodes (gold), and converging on the axis — the lens focuses,
exactly as the trajectory numbers assert.

![einzel lens](validation/testdata/einzel.png)

## Live in-application validation

The add-in itself is driven from the running application; the study is `Traceon.RunStudy`
(a ribbon command, also invokable over the MCP bridge's `execute_command`). To validate live:

1. `make install` — builds the c-shared library and copies it + `manifest.json` into the
   host's `addins/` directory (the host scans it at startup).
2. Launch Oblikovati with the MCP bridge add-in loaded.
3. Over the MCP bridge: `close_all_documents(force)`, build or import an axisymmetric electrode
   (a revolved profile), then `execute_command("Traceon.RunStudy")`.
4. `viewport.capture` a screenshot and confirm the electrode, the potential heatmap, and the
   focusing trajectories render.

`cmd/traceonfield` is the offline stand-in: it runs the real section→solve→trace→render
pipeline against a canned electrode profile and dumps the exact client-graphics payload the
add-in would push, with no live app or render backend:

```
go run ./cmd/traceonfield > study.json   # 3 elements, 7 rays, 9 graphics nodes
```

### Known calibration items (tracked for a follow-up)

- **Units.** Host geometry arrives in the DB unit (cm) and the study traces in those units.
  The trajectory *shape* is faithful; absolute beam dynamics need an explicit cm→m conversion
  at the engine boundary (the `core/` physics is unit-consistent, so the validation above —
  run entirely in `core/` units — is unaffected).
- **Per-electrode voltages.** The engine v1 applies one voltage to the whole sectioned body;
  a true multi-electrode lens (as in the `core/`-level validation) needs per-face voltage
  assignment, keyed by face reference key.
