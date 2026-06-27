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


@fixture("ring", "ring")
def _ring() -> dict:
    """Single-ring kernels: potential, r/z derivatives, current potential/field, and the
    on-axis derivative recurrences. Deterministic samples (seeded) so regen is stable."""
    import ctypes as C
    import numpy as np
    import traceon.backend as B

    # axial_derivatives_radial_ring / current_axial_derivatives_radial_ring have no Python
    # wrapper (used internally), so call the C symbols directly. void fn(z0, r, z, out[9]).
    def _ring_derivs(symbol: str, z0: float, r: float, z: float) -> list[float]:
        fn = getattr(B.backend_lib, symbol)
        fn.restype = None
        fn.argtypes = [C.c_double, C.c_double, C.c_double, C.POINTER(C.c_double)]
        out = (C.c_double * 9)()
        fn(z0, r, z, out)
        return [float(v) for v in out]

    rng = np.random.default_rng(20240617)
    # Generic off-axis samples in (0,1]^4 plus a few structured near-axis / near-singular
    # points that stress the guards (r0→0, delta→0).
    samples = rng.uniform(0.01, 1.0, size=(60, 4)).tolist()
    samples += [
        [1.0, 0.0, 0.0, 0.5],    # pure axial offset
        [1.0, 0.0, 0.5, 0.0],    # pure radial offset
        [2.0, 1.0, 0.3, -0.4],
        [1e-11, 0.0, 0.5, 0.5],  # r0 BELOW MinDistanceAxis (1e-10) → dr1 guard returns 0
        [55.0, 0.0, 0.0, 1.0],
    ]

    potential, dr1, dz1, axial_derivs = [], [], [], []
    for r0, z0, a, b in samples:
        # potential/dr1/dz1 take (r0, z0, delta_r, delta_z); use a,b as the deltas.
        potential.append(float(B.potential_radial_ring(r0, z0, a, b)))
        dr1.append(float(B.dr1_potential_radial_ring(r0, z0, a, b)))
        dz1.append(float(B.dz1_potential_radial_ring(r0, z0, a, b)))
        axial_derivs.append(_ring_derivs("axial_derivatives_radial_ring", z0, max(a, 1e-3), z0 + b))

    # Current ring: potential + 2-vector field + axial derivative recurrence.
    cur_samples = rng.uniform(0.01, 3.0, size=(40, 4)).tolist()
    cur_samples += [[0.0, 0.0, 2.0, 0.0], [0.0, 1.5, 2.0, 0.0], [1e-9, 0.0, 1.0, 0.0]]
    cur_potential, cur_field, cur_axial_derivs = [], [], []
    for x0, y0, x, y in cur_samples:
        cur_potential.append(float(B.current_potential_axial_radial_ring(y0, x, y)))
        f = B.current_field_radial_ring(x0, y0, x, y)
        cur_field.append([float(f[0]), float(f[1])])
        cur_axial_derivs.append(_ring_derivs("current_axial_derivatives_radial_ring", y0, x, y))

    return {
        "samples": samples,
        "potential": potential,
        "dr1": dr1,
        "dz1": dz1,
        "axial_derivs": axial_derivs,
        "cur_samples": cur_samples,
        "cur_potential": cur_potential,
        "cur_field": cur_field,
        "cur_axial_derivs": cur_axial_derivs,
    }


@fixture("radial", "radial")
def _radial() -> dict:
    """Electrostatic radial BEM: jacobian buffers, charge, potential/field evaluation,
    the singular self-term integrands, and dense matrix assembly on a small line set."""
    import ctypes as C
    import numpy as np
    import traceon.backend as B

    def line4(r0, z0, r1, z1):
        # GMSH line4 ordering: [start, end, 1/3-point, 2/3-point].
        return [[r0, 0.0, z0], [r1, 0.0, z1],
                [r0 + (r1 - r0) / 3, 0.0, z0 + (z1 - z0) / 3],
                [r0 + 2 * (r1 - r0) / 3, 0.0, z0 + 2 * (z1 - z0) / 3]]

    lines = np.array([
        line4(1.0, 0.0, 1.0, 1.0),    # vertical at r=1 (charge_radial → 2π)
        line4(2.0, 0.0, 2.0, 0.5),    # vertical at r=2
        line4(1.0, 0.0, 1.5, 0.5),    # slanted
    ])
    n = len(lines)

    jac, pos = B.fill_jacobian_buffer_radial(lines)
    charges = np.array([1.0, 0.7, -0.4])
    charge_radial = [float(B.charge_radial(lines[i], 1.0)) for i in range(n)]

    eval_points = np.array([
        [0.0, 0.0, 0.5], [1.2, 0.0, 0.3], [3.0, 0.0, -1.0], [0.5, 0.0, 2.0],
    ])
    potential = [float(B.potential_radial(p, charges, jac, pos)) for p in eval_points]
    field = [[float(c) for c in B.field_radial(p, charges, jac, pos)] for p in eval_points]

    # Singular self-term integrands at sampled α (raw C functions; the diagonal is their
    # integral over α — verified in the solver PBI). self_potential: double fn(double, double[4][3]).
    alphas = [-0.9, -0.5, -0.123, 0.25, 0.5, 0.9]
    sp = B.backend_lib.self_potential_radial
    sp.restype = C.c_double
    sp.argtypes = [C.c_double, C.POINTER(C.c_double)]

    class _SFArgs(C.Structure):
        _fields_ = [("line_points", C.POINTER(C.c_double)), ("K", C.c_double)]

    sf = B.backend_lib.self_field_dot_normal_radial
    sf.restype = C.c_double
    sf.argtypes = [C.c_double, C.POINTER(_SFArgs)]

    self_potential, self_field = [], []
    K = 2.0
    for i in range(n):
        lp = np.ascontiguousarray(lines[i], dtype=np.float64)
        lp_ptr = lp.ctypes.data_as(C.POINTER(C.c_double))
        self_potential.append([float(sp(C.c_double(a), lp_ptr)) for a in alphas])
        args = _SFArgs(line_points=lp_ptr, K=K)
        self_field.append([float(sf(C.c_double(a), C.byref(args))) for a in alphas])

    # Dense matrix assembly: 2 voltage rows + 1 dielectric row.
    exc_types = np.array([1, 1, 3], dtype=np.uint8)        # VOLTAGE_FIXED, VOLTAGE_FIXED, DIELECTRIC
    exc_values = np.array([1.0, 0.5, K])
    matrix = np.zeros((n, n))
    B.fill_matrix_radial(matrix, lines, exc_types, exc_values, jac, pos, 0, n - 1)

    # On-axis accumulated derivatives (AxialDerivativesRadial) for unit charges.
    z_axis = np.linspace(-1.0, 2.0, 8)
    unit_charges = np.array([1.0, 0.7, -0.4])
    axial_derivs = B.axial_derivatives_radial(z_axis, unit_charges, jac, pos).tolist()

    return {
        "z_axis": z_axis.tolist(),
        "unit_charges": unit_charges.tolist(),
        "axial_derivs": axial_derivs,
        "lines": lines.tolist(),
        "jac": jac.tolist(),
        "pos": pos.tolist(),
        "charges": charges.tolist(),
        "charge_radial": charge_radial,
        "eval_points": eval_points.tolist(),
        "potential": potential,
        "field": field,
        "alphas": alphas,
        "K": K,
        "self_potential": self_potential,
        "self_field": self_field,
        "exc_types": exc_types.tolist(),
        "exc_values": exc_values.tolist(),
        "matrix": matrix.tolist(),
    }


@fixture("solver", "solver")
def _solver() -> dict:
    """Electrostatic radial solve: the integrated singular self-terms, the assembled
    matrix (with overwritten diagonal), the right-hand side, and the solved charges —
    reproducing traceon.solver.ElectrostaticSolverRadial's get_matrix + solve_matrix on a
    small explicit line set (no mesher needed)."""
    import numpy as np
    import traceon.backend as B

    def line4(r0, z0, r1, z1):
        return [[r0, 0.0, z0], [r1, 0.0, z1],
                [r0 + (r1 - r0) / 3, 0.0, z0 + (z1 - z0) / 3],
                [r0 + 2 * (r1 - r0) / 3, 0.0, z0 + 2 * (z1 - z0) / 3]]

    lines = np.array([line4(1.0, 0.0, 1.0, 1.0), line4(2.0, 0.0, 2.0, 0.5), line4(1.0, 0.0, 1.5, 0.5)])
    n = len(lines)
    types = np.array([1, 1, 3], dtype=np.uint8)   # VOLTAGE_FIXED, VOLTAGE_FIXED, DIELECTRIC
    values = np.array([1.0, 0.5, 2.0])            # volts, volts, relative permittivity

    jac, pos = B.fill_jacobian_buffer_radial(lines)
    matrix = np.zeros((n, n))
    B.fill_matrix_radial(matrix, lines, types, values, jac, pos, 0, n - 1)

    self_potential = [float(B.self_potential_radial(lines[i])) for i in range(n)]
    self_field = [float(B.self_field_dot_normal_radial(lines[i], values[i])) for i in range(n)]

    # Overwrite the diagonal exactly as SolverRadial.get_matrix does.
    for i in range(n):
        if types[i] == 3:  # DIELECTRIC
            matrix[i, i] = self_field[i] - 1
        else:
            matrix[i, i] = self_potential[i]

    rhs = np.array([1.0, 0.5, 0.0])  # voltage→value, dielectric→0
    charges = np.linalg.solve(matrix, rhs)

    return {
        "lines": lines.tolist(),
        "types": types.tolist(),
        "values": values.tolist(),
        "self_potential": self_potential,
        "self_field": self_field,
        "matrix": matrix.tolist(),
        "rhs": rhs.tolist(),
        "charges": charges.tolist(),
    }


@fixture("interp", "interp")
def _interp() -> dict:
    """scipy not-a-knot cubic and BPoly-derived quintic Hermite coefficients, the two
    interpolations the axial field series uses, over equally-spaced sample data."""
    import numpy as np
    from scipy.interpolate import CubicSpline, BPoly, PPoly

    rng = np.random.default_rng(424242)
    z = np.linspace(-2.0, 3.0, 12)
    y = rng.uniform(-1, 1, size=z.size)
    dy = rng.uniform(-1, 1, size=z.size)
    d2y = rng.uniform(-1, 1, size=z.size)

    cubic = CubicSpline(z, y).c.T.tolist()  # (n-1, 4) descending: [c3,c2,c1,c0]

    # _get_one_dimensional_high_order_ppoly: quintic Hermite via Bernstein → power basis.
    bpoly = BPoly.from_derivatives(z, np.array([y, dy, d2y]).T)
    quintic = PPoly.from_bernstein_basis(bpoly).c.T.tolist()  # (n-1, 6): [c5..c0]

    return {
        "z": z.tolist(),
        "y": y.tolist(),
        "dy": dy.tolist(),
        "d2y": d2y.tolist(),
        "cubic": cubic,
        "quintic": quintic,
    }


@fixture("field", "field")
def _field() -> dict:
    """Axial-series field interpolation: per-z axial derivatives, the assembled quintic
    coefficients, and the interpolated potential/field — plus the direct-BEM field/potential
    the FieldRadialBEM wrappers expose. Driven by a solved charge distribution."""
    import numpy as np
    import traceon.backend as B
    from traceon.field import _quintic_spline_coefficients

    def line4(r0, z0, r1, z1):
        return [[r0, 0.0, z0], [r1, 0.0, z1],
                [r0 + (r1 - r0) / 3, 0.0, z0 + (z1 - z0) / 3],
                [r0 + 2 * (r1 - r0) / 3, 0.0, z0 + 2 * (z1 - z0) / 3]]

    # A small ring of charged elements around the axis (a crude electrode), then solve.
    lines = np.array([
        line4(1.0, -0.5, 1.0, 0.5),
        line4(1.0, 0.5, 1.0, 1.5),
        line4(1.0, -1.5, 1.0, -0.5),
    ])
    n = len(lines)
    types = np.array([1, 1, 1], dtype=np.uint8)
    values = np.array([1.0, 0.5, 0.5])
    jac, pos = B.fill_jacobian_buffer_radial(lines)
    matrix = np.zeros((n, n))
    B.fill_matrix_radial(matrix, lines, types, values, jac, pos, 0, n - 1)
    for i in range(n):
        matrix[i, i] = B.self_potential_radial(lines[i])
    charges = np.linalg.solve(matrix, values)

    # Axial interpolation over [zmin, zmax].
    zmin, zmax, N = -2.0, 2.0, 20
    z = np.linspace(zmin, zmax, N)
    derivs = B.axial_derivatives_radial(z, charges, jac, pos)           # (N, 9)
    coeffs = _quintic_spline_coefficients(z, derivs.T)                  # (N-1, 9, 6)

    eval_points = np.array([
        [0.0, 0.0, 0.0], [0.1, 0.0, 0.3], [0.2, 0.0, -0.7], [0.05, 0.0, 1.1],
        [0.0, 0.0, 5.0],  # outside [zmin, zmax] → zero
    ])
    pot_direct = [float(B.potential_radial(p, charges, jac, pos)) for p in eval_points]
    field_direct = [[float(c) for c in B.field_radial(p, charges, jac, pos)] for p in eval_points]
    pot_interp = [float(B.potential_radial_derivs(p, z, coeffs)) for p in eval_points]
    field_interp = [[float(c) for c in B.field_radial_derivs(p, z, coeffs)] for p in eval_points]

    return {
        "lines": lines.tolist(),
        "charges": charges.tolist(),
        "z": z.tolist(),
        "derivs": derivs.tolist(),
        "coeffs": coeffs.tolist(),
        "eval_points": eval_points.tolist(),
        "pot_direct": pot_direct,
        "field_direct": field_direct,
        "pot_interp": pot_interp,
        "field_interp": field_interp,
    }


@fixture("tracing", "tracing")
def _tracing() -> dict:
    """A full RKF45 trajectory through a real radial BEM field, to verify the tracer and
    the field evaluation compose correctly against upstream over hundreds of steps."""
    import numpy as np
    import traceon.backend as B
    import traceon.field as F
    import traceon.tracing as T

    def line4(r0, z0, r1, z1):
        return [[r0, 0.0, z0], [r1, 0.0, z1],
                [r0 + (r1 - r0) / 3, 0.0, z0 + (z1 - z0) / 3],
                [r0 + 2 * (r1 - r0) / 3, 0.0, z0 + 2 * (z1 - z0) / 3]]

    lines = np.array([line4(1.0, -0.5, 1.0, 0.5), line4(1.0, 0.5, 1.0, 1.5), line4(1.0, -1.5, 1.0, -0.5)])
    n = len(lines)
    types = np.array([1, 1, 1], dtype=np.uint8)
    values = np.array([1.0, 0.5, 0.5])
    jac, pos = B.fill_jacobian_buffer_radial(lines)
    matrix = np.zeros((n, n))
    B.fill_matrix_radial(matrix, lines, types, values, jac, pos, 0, n - 1)
    for i in range(n):
        matrix[i, i] = B.self_potential_radial(lines[i])
    charges = np.linalg.solve(matrix, values)

    field = F.FieldRadialBEM(electrostatic_point_charges=F.EffectivePointCharges(charges, jac, pos))
    bounds = ((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0))
    tracer = T.Tracer(field, bounds)

    energy_eV = 100.0
    p0 = [0.05, 0.0, -1.9]
    v0 = T.velocity_vec(energy_eV, [0, 0, 1])
    times, positions = tracer(np.array(p0), v0)

    return {
        "lines": lines.tolist(),
        "charges": charges.tolist(),
        "energy_eV": energy_eV,
        "p0": p0,
        "atol": 1e-8,
        "charge_over_mass": float(-1.602176634e-19 / 9.1093837139e-31),
        "bounds": [list(b) for b in bounds],
        "times": times.tolist(),
        "positions": positions.tolist(),
    }


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
