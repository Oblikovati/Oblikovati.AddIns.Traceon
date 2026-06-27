#!/usr/bin/env python3
# SPDX-License-Identifier: MPL-2.0
"""Generate golden-fixture JSON from the upstream Traceon backend.

Each fixture drives a Traceon backend/solver function over representative + sampled
inputs and writes (inputs, outputs) to core/<pkg>/testdata/<name>.golden.json. The Go
port loads these and asserts numerical equivalence to np.isclose tolerances, with NO
Python at test time. This is the oracle that keeps the pure-Go core honest.

Run via `make verify-oracle` (which points --traceon-dir at the sibling upstream clone
whose .venv has the C backend built), or directly:

    ../Traceon/.venv/bin/python tools/gen_fixtures.py --traceon-dir ../Traceon --out core

Add a fixture: write a function decorated with @fixture("<pkg>", "<name>") returning a
{"cases": [...]} dict, then it is emitted automatically. Keep inputs deterministic
(seed any RNG) so regeneration is reproducible and diffs are meaningful.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Any, Callable

# Registry of (pkg, name, generator) populated by the @fixture decorator.
_FIXTURES: list[tuple[str, str, Callable[[], dict]]] = []


def _json_safe(obj: Any) -> Any:
    """Recursively replace non-finite floats with strings so the output is valid JSON.

    Standard JSON has no NaN/Infinity; Python's json emits the bare tokens NaN/Infinity
    which Go's encoding/json rejects. Radial BEM kernels legitimately produce these at
    singularities, so we encode them as "NaN"/"Infinity"/"-Infinity" strings and decode
    them on the Go side via oracle.F. Finite floats stay as JSON numbers.
    """
    if isinstance(obj, float):
        if math.isnan(obj):
            return "NaN"
        if math.isinf(obj):
            return "Infinity" if obj > 0 else "-Infinity"
        return obj
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def fixture(pkg: str, name: str):
    """Register a golden-fixture generator emitting to core/<pkg>/testdata/<name>.golden.json."""
    def deco(fn: Callable[[], dict]) -> Callable[[], dict]:
        _FIXTURES.append((pkg, name, fn))
        return fn
    return deco


# --------------------------------------------------------------------------------------
# Fixtures. Import the backend lazily inside generators so --help works without the venv.
# --------------------------------------------------------------------------------------

@fixture("elliptic", "elliptic")
def _elliptic() -> dict:
    """Complete elliptic integrals K and E (and the m-1 variants) over m in (0, 1)."""
    import numpy as np
    import traceon.backend as B

    # Dense deterministic sweep plus values near the endpoints where the Chebyshev
    # approximations (Cody 1965) are most sensitive.
    ms = sorted(set(
        list(np.linspace(1e-6, 1.0 - 1e-6, 50))
        + [1e-9, 1e-6, 1e-3, 0.5, 0.9, 0.99, 1.0 - 1e-9]
    ))
    cases = []
    for m in ms:
        cases.append({
            "m": m,
            "ellipk": float(B.ellipk(m)),
            "ellipe": float(B.ellipe(m)),
            "ellipkm1": float(B.ellipkm1(1.0 - m)),
            "ellipem1": float(B.ellipem1(1.0 - m)),
        })

    # Reciprocal-modulus branch: Ellipk/Ellipe at m outside [0, 1] exercise the
    # imaginary-modulus transforms in the C (which neither test_elliptic.py nor the
    # km1/em1 columns above reach). km1/em1 are undefined there, so omit them.
    k_e_only = []
    for m in [-5.0, -2.0, -1.5, 2.0, 5.0, 10.0]:
        k_e_only.append({
            "m": m,
            "ellipk": float(B.ellipk(m)),
            "ellipe": float(B.ellipe(m)),
        })

    return {"cases": cases, "k_e_only": k_e_only}


# --------------------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--traceon-dir", required=True,
                    help="path to the upstream Traceon checkout (its backend .so must be built)")
    ap.add_argument("--out", default="core",
                    help="output root; fixtures land in <out>/<pkg>/testdata/")
    ap.add_argument("--only", default=None,
                    help="comma-separated pkg/name filters; default = all")
    args = ap.parse_args()

    sys.path.insert(0, os.path.abspath(args.traceon_dir))

    only = set(args.only.split(",")) if args.only else None
    wrote = 0
    for pkg, name, gen in _FIXTURES:
        if only is not None and pkg not in only and name not in only:
            continue
        data = gen()
        dest_dir = os.path.join(args.out, pkg, "testdata")
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, f"{name}.golden.json")
        with open(dest, "w") as f:
            json.dump(_json_safe(data), f, indent=1, sort_keys=True, allow_nan=False)
            f.write("\n")
        n = len(data.get("cases", []))
        print(f"wrote {dest} ({n} cases)")
        wrote += 1

    if wrote == 0:
        print("no fixtures matched", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
