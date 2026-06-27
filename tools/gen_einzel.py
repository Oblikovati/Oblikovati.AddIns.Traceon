#!/usr/bin/env python3
# SPDX-License-Identifier: MPL-2.0
"""Generate the einzel-lens end-to-end validation golden from upstream Traceon.

Builds a three-electrode einzel lens directly as (r, z) BEM line elements (no mesher),
solves the electrostatic field, traces a paraxial beam, and records the focal length plus
one full trajectory. The Go validation runs the same case through core/ and asserts it
reproduces these numbers — exercising solver → field → tracing → focus end to end against
the upstream as oracle.

    ../Traceon/.venv/bin/python tools/gen_einzel.py --traceon-dir ../Traceon
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys


def line4(r0, z0, r1, z1):
    return [[r0, 0.0, z0], [r1, 0.0, z1],
            [r0 + (r1 - r0) / 3, 0.0, z0 + (z1 - z0) / 3],
            [r0 + 2 * (r1 - r0) / 3, 0.0, z0 + 2 * (z1 - z0) / 3]]


def electrode(r, z_lo, z_hi, n):
    """A cylindrical tube electrode at radius r over [z_lo, z_hi], as n line4 elements."""
    out = []
    for i in range(n):
        z0 = z_lo + (z_hi - z_lo) * i / n
        z1 = z_lo + (z_hi - z_lo) * (i + 1) / n
        out.append(line4(r, z0, r, z1))
    return out


def build():
    import numpy as np
    import traceon.backend as B
    import traceon.field as F
    import traceon.tracing as T
    import traceon.focus as FO

    R = 1.0          # electrode radius
    Vc = 5000.0      # centre-electrode voltage (outer two grounded)
    per = 8          # elements per electrode

    lines = []
    types = []
    values = []
    for (zlo, zhi, volt) in [(-3.0, -1.0, 0.0), (-0.5, 0.5, Vc), (1.0, 3.0, 0.0)]:
        els = electrode(R, zlo, zhi, per)
        lines.extend(els)
        types.extend([1] * len(els))      # VOLTAGE_FIXED
        values.extend([volt] * len(els))

    lines = np.array(lines)
    types = np.array(types, dtype=np.uint8)
    values = np.array(values)
    n = len(lines)

    jac, pos = B.fill_jacobian_buffer_radial(lines)
    matrix = np.zeros((n, n))
    B.fill_matrix_radial(matrix, lines, types, values, jac, pos, 0, n - 1)
    for i in range(n):
        matrix[i, i] = B.self_potential_radial(lines[i])
    charges = np.linalg.solve(matrix, values)

    field = F.FieldRadialBEM(electrostatic_point_charges=F.EffectivePointCharges(charges, jac, pos))
    bounds = ((-R, R), (-0.01, 0.01), (-6.0, 8.0))
    tracer = T.Tracer(field, bounds)

    energy_eV = 3000.0
    launch_z = -5.0
    radii = [0.02, 0.04, 0.06, 0.08, 0.10]
    v0 = T.velocity_vec(energy_eV, [0, 0, 1])

    trajectories = []
    for r0 in radii:
        _, p = tracer(np.array([r0, 0.0, launch_z]), v0)
        trajectories.append(p)

    focus = FO.focus_position(trajectories)

    # Record one full trajectory (the outermost ray) for a step-wise trajectory check.
    sample_traj = trajectories[-1].tolist()

    _plot(lines, trajectories, focus, Vc)

    return {
        "lines": lines.tolist(),
        "types": types.tolist(),
        "values": values.tolist(),
        "charges": charges.tolist(),
        "energy_eV": energy_eV,
        "launch_z": launch_z,
        "radii": radii,
        "charge_over_mass": float(-1.602176634e-19 / 9.1093837139e-31),
        "bounds": [list(b) for b in bounds],
        "focus": focus.tolist(),
        "sample_radius": radii[-1],
        "sample_trajectory": sample_traj,
    }


def _plot(lines, trajectories, focus, Vc):
    """Save an (r, z) plot of the electrode + trajectories converging on the focus — a visual
    confirmation that the lens focuses (written next to the golden as einzel.png)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, ax = plt.subplots(figsize=(9, 4))
    for el in lines:
        ax.plot([el[0][2], el[1][2]], [el[0][0], el[1][0]], color="goldenrod", lw=3)
        ax.plot([el[0][2], el[1][2]], [-el[0][0], -el[1][0]], color="goldenrod", lw=3)
    for p in trajectories:
        ax.plot(p[:, 2], p[:, 0], color="tab:green", lw=1)
        ax.plot(p[:, 2], -p[:, 0], color="tab:green", lw=1)
    ax.axvline(focus[2], color="tab:red", ls="--", lw=1, label=f"focus z={focus[2]:.3f}")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlabel("z"); ax.set_ylabel("r")
    ax.set_title(f"Traceon einzel lens (centre V={Vc:g}) — beam focuses on axis")
    ax.legend(loc="upper left")
    out = "validation/testdata/einzel.png"
    fig.tight_layout(); fig.savefig(out, dpi=110)
    print(f"wrote {out}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--traceon-dir", required=True)
    ap.add_argument("--out", default="validation/testdata/einzel.golden.json")
    args = ap.parse_args()
    sys.path.insert(0, os.path.abspath(args.traceon_dir))

    data = build()

    def safe(o):
        if isinstance(o, float):
            if math.isnan(o):
                return "NaN"
            if math.isinf(o):
                return "Infinity" if o > 0 else "-Infinity"
        if isinstance(o, dict):
            return {k: safe(v) for k, v in o.items()}
        if isinstance(o, list):
            return [safe(v) for v in o]
        return o

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(safe(data), f, indent=1, sort_keys=True, allow_nan=False)
        f.write("\n")
    print(f"wrote {args.out}: focus={data['focus']}, {len(data['lines'])} elements, "
          f"{len(data['sample_trajectory'])}-step sample trajectory")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
