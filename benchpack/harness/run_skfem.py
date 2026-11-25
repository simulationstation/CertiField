"""Run scikit-fem on the benchmark cases and write benchpack/results/skfem_<case>.csv.

    python benchpack/harness/run_skfem.py P1 [P2 P3c ...]   (all cases if none given)

Element variants:
    P1           linear triangles / tetrahedra, straight-sided (polygonal / polyhedral hole)
    P2-straight  quadratic elements on the same straight-sided mesh
    P2-curved    quadratic isoparametric elements, mid-edge nodes on the circle / sphere (gmsh order 2)
All values are FLOAT (floating point, no error control); checked by the agent only.
Resumable: each run is appended to the CSV when it finishes; runs with a row already are skipped.  Runs whose
estimated memory exceeds the cap of runlog.py (4 GB) are not done and get a row saying so.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runlog  # noqa: E402,F401  (caps the threads before numpy is loaded)

import numpy as np
import skfem
from skfem import (BilinearForm, CellBasis, ElementTetP1, ElementTetP2, ElementTriP1, ElementTriP2, ElementVector,
                   Functional, MeshTet1, MeshTet2, MeshTri1, MeshTri2, asm, condense)
from skfem.helpers import ddot, sym_grad, trace

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases  # noqa: E402

try:
    import pypardiso

    def spsolve(A, b):
        return pypardiso.spsolve(A.tocsr(), b)
    SOLVER = "pypardiso"
except ImportError:  # pragma: no cover
    from scipy.sparse.linalg import spsolve
    SOLVER = "scipy-superlu"

RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
VARIANTS = ["P1", "P2-straight", "P2-curved"]


def lame(dim):
    E, nu = cases.E, cases.NU
    mu = E / (2 * (1 + nu))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)) if dim == 3 else E * nu / (1 - nu ** 2)   # plane stress in 2D
    return lam, mu


def build(case, h, variant):
    c = cases.CASES[case]
    dim = c["dim"]
    if variant == "P1":
        m = cases.make_mesh(case, h, 1, False)
        mesh = (MeshTri1 if dim == 2 else MeshTet1)(m["points"].T.copy(), m["cells"].T.copy())
        elem = ElementTriP1() if dim == 2 else ElementTetP1()
    elif variant == "P2-straight":
        m = cases.make_mesh(case, h, 1, False)
        mesh = (MeshTri1 if dim == 2 else MeshTet1)(m["points"].T.copy(), m["cells"].T.copy())
        elem = ElementTriP2() if dim == 2 else ElementTetP2()
    elif variant == "P2-curved":
        m = cases.make_mesh(case, h, 2, True)
        mesh = (MeshTri2 if dim == 2 else MeshTet2)(m["points"].T.copy(), m["cells"].T.copy())
        elem = ElementTriP2() if dim == 2 else ElementTetP2()
    else:
        raise KeyError(variant)
    return mesh, ElementVector(elem)


def stiffness(basis, dim):
    lam, mu = lame(dim)

    @BilinearForm
    def a(u, v, w):
        eu, ev = sym_grad(u), sym_grad(v)
        return 2 * mu * ddot(eu, ev) + lam * trace(eu) * trace(ev)

    return asm(a, basis)


def patch_average(case, mesh, e, u):
    """(Q2, patch area): the average of sigma_xx of the displacement vector u over the patch elements."""
    c = cases.CASES[case]
    lam, mu = lame(c["dim"])
    cent = mesh.p[:, mesh.t].mean(axis=1)
    sel = cases.in_p1_patch(cent) if c["q2"] == "p1" else cases.in_p2_patch(cent)
    idx = np.nonzero(sel)[0]
    pb = CellBasis(mesh, e, elements=idx, intorder=max(4, 2 * e.maxdeg + 2))

    @Functional
    def sxx(w):
        eps = sym_grad(w["u"])
        return 2 * mu * eps[0, 0] + lam * trace(eps)

    @Functional
    def one(w):
        return 1.0 + 0.0 * w.x[0]

    area = one.assemble(pb)
    return sxx.assemble(pb, u=pb.interpolate(u)) / area, area


def memory_gb(dim, ndof, quadratic):
    """estimated peak memory of a run (GB), from measurements with pypardiso and 2 threads on the harness machine:
    2D 0.68 GB at 148 000 unknowns (quadratic), 3D 3.67 GB at 145 000 (quadratic) and 1.14 GB at 55 000"""
    per_million = 4.6 if dim == 2 else (26.0 if quadratic else 12.0)
    return 0.3 + per_million * ndof / 1e6


def solve(case, h, variant):
    c = cases.CASES[case]
    dim = c["dim"]
    t0 = time.perf_counter()
    mesh, e = build(case, h, variant)
    t_mesh = time.perf_counter() - t0
    t0 = time.perf_counter()
    basis = CellBasis(mesh, e)
    est = memory_gb(dim, basis.N, variant != "P1")
    if est > runlog.MAX_GB:
        return dict(code="scikit-fem", case=case, element=variant, h=h, ndof=basis.N, nelem=mesh.nelements,
                    note="not run: estimated %.1f GB, above the %.0f GB cap" % (est, runlog.MAX_GB))
    A = stiffness(basis, dim)
    x = np.zeros(basis.N)
    D = []
    for (axis, val), comp, g in c["bc"]:
        facets = mesh.facets_satisfying(lambda X, axis=axis, val=val: np.abs(X[axis] - val) < 1e-9,
                                        boundaries_only=True)
        dofs = basis.get_dofs(facets).all("u^%d" % (comp + 1))
        x[dofs] = g
        D.append(dofs)
    D = np.unique(np.concatenate(D))
    AII, bI, xx, I = condense(A, x=x, D=D)
    xx[I] = spsolve(AII, bI)
    u = xx
    r = A @ u
    (gaxis, gval), gcomp = c["grip"]
    gf = mesh.facets_satisfying(lambda X: np.abs(X[gaxis] - gval) < 1e-9, boundaries_only=True)
    gdofs = basis.get_dofs(gf).all("u^%d" % (gcomp + 1))
    F = float(np.sum(r[gdofs]))
    K = c["kscale"] * F * np.sign(c["d_model"])
    Ken = c["energy_scale"] * float(u @ r)
    q2, area = float("nan"), float("nan")
    if c["q2"]:
        q2, area = patch_average(case, mesh, e, u)
    t_solve = time.perf_counter() - t0
    return dict(code="scikit-fem", case=case, element=variant, h=h, ndof=basis.N, nfree=len(I),
                nelem=mesh.nelements, time_s=round(t_solve, 3), mesh_time_s=round(t_mesh, 3),
                K=K, K_energy=Ken, Q2=q2, patch_area=area)


FIELDS = ["code", "case", "element", "h", "ndof", "nfree", "nelem", "time_s", "mesh_time_s", "K", "K_energy", "Q2",
          "patch_area", "note"]


def main(argv):
    sys.stdout.reconfigure(encoding="utf-8")
    todo = argv or list(cases.CASES)
    os.makedirs(RESULTS, exist_ok=True)
    for case in todo:
        log = runlog.CsvLog(os.path.join(RESULTS, "skfem_%s.csv" % case), FIELDS)
        for variant in VARIANTS:
            if variant == "P2-curved" and cases.CASES[case]["mesh"] == "lbracket":
                continue          # no curved boundary: identical to P2-straight
            levels = cases.LEVELS[case] if variant == "P1" else cases.LEVELS_QUAD.get(case, cases.LEVELS[case])
            for h in levels:
                if log.has(element=variant, h="%.12g" % h):
                    continue
                row = solve(case, h, variant)
                log.append(row)
                if "K" in row:
                    print("%-8s %-12s h=%-8g ndof=%-8d t=%7.2fs  K=%.10f  (energy %.10f)  Q2=%.8f" % (
                        case, variant, h, row["ndof"], row["time_s"], row["K"], row["K_energy"], row["Q2"]),
                        flush=True)
                else:
                    print("%-8s %-12s h=%-8g %s" % (case, variant, h, row["note"]), flush=True)
    print("solver:", SOLVER, " scikit-fem", skfem.__version__)


if __name__ == "__main__":
    main(sys.argv[1:])
