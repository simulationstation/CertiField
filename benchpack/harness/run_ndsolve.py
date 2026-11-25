"""Run Wolfram NDSolve (finite elements) on the harness meshes and write benchpack/results/ndsolve_<case>.csv.

    python benchpack/harness/run_ndsolve.py P1 [P2 P3c ...]

One kernel per run, one run at a time: $WOLFRAMSCRIPT -file benchpack/harness/ndsolve_fem.wls <job.wl> (default
wolframscript on PATH; set WOLFRAMSCRIPT to a wrapper command, e.g. '"C:/Program Files/Git/bin/bash.exe" /path/wl.sh',
where the kernel seats are shared).  A run whose kernel stops twice gets a row saying so.
Resumable: each run is appended to the CSV when it finishes; runs with a row already are skipped.  Runs estimated
above the memory cap of runlog.py get a row saying so.
Meshes are the gmsh meshes of run_skfem.py / run_calculix.py; element variants:
    linear (order-1 mesh), quadratic-straight (order 2, mid-edge nodes at the chord midpoints),
    quadratic-curved (order 2, mid-edge nodes on the circle / sphere).
Columns: K = K_system of ndsolve_fem.wls (a(u_h, u_h) from NDSolve's assembled system);
K_nintegrate, Q2_nintegrate = NIntegrate of the InterpolatingFunction (the usual NDSolve post-processing);
Q2 and K_check = the quantities of NDSolve's nodal displacement field, evaluated here with the element shape
functions (scikit-fem basis on the same mesh, same element space).
All values are FLOAT; checked by the agent only.
"""
import csv
import os
import shlex
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runlog  # noqa: E402  (caps the threads before numpy is loaded)

import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from skfem import CellBasis  # noqa: E402

import cases  # noqa: E402
import run_skfem  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")
VARIANTS = [("linear", 1, False, "P1"), ("quadratic-straight", 2, False, "P2-straight"),
            ("quadratic-curved", 2, True, "P2-curved")]
# NDSolve runs: the first levels of cases.LEVELS (NIntegrate over large meshes is slow)
NLEVELS = {"P1solid": 3, "P1": 6, "P1c": 6, "P2": 6, "P3c": 5, "P3f": 5}
NLEVELS_QUAD_3D = 4
# node order of a 10-node TetrahedronElement in Mathematica: vertices 1-4, then the mid-edge nodes of the edges
# (1,2) (2,3) (3,1) (2,4) (3,4) (1,4) -- checked against ToElementMesh's own second-order tetrahedra and the volume of
# a unit tetrahedron; the VTK / Abaqus order used by the harness ends with (1,4) (2,4) (3,4)
MMA_TET10 = [0, 1, 2, 3, 4, 5, 6, 8, 9, 7]


def wpath(p):
    return p.replace("\\", "/")


def job_entry(name, case, h, order, curved, workdir):
    c = cases.CASES[case]
    m = cases.make_mesh(case, h, order, curved)
    P, C = m["points"], m["cells"]
    stem = os.path.join(workdir, name)
    np.savetxt(stem + "_coords.csv", P, delimiter=",", fmt="%.17g")
    CW = C[:, MMA_TET10] if C.shape[1] == 10 else C
    np.savetxt(stem + "_cells.csv", CW + 1, delimiter=",", fmt="%d")
    pc = pcells = "None"
    if c["q2"]:
        cent = P[C[:, :3]].mean(axis=1).T
        sel = cases.in_p1_patch(cent) if c["q2"] == "p1" else cases.in_p2_patch(cent)
        sub = C[sel]
        used = np.unique(sub)
        remap = np.full(len(P), -1)
        remap[used] = np.arange(len(used))
        np.savetxt(stem + "_pcoords.csv", P[used], delimiter=",", fmt="%.17g")
        np.savetxt(stem + "_pcells.csv", remap[sub] + 1, delimiter=",", fmt="%d")
        pc = '"%s"' % wpath(stem + "_pcoords.csv")
        pcells = '"%s"' % wpath(stem + "_pcells.csv")
    bc = "{" + ", ".join("{%d, %r, %d, %r}" % (ax, float(val), comp, float(g)) for (ax, val), comp, g in c["bc"]) + "}"
    bc = bc.replace("e-", "*^-").replace("e+", "*^")
    out, nodal = wpath(stem + "_out.csv"), wpath(stem + "_nodal.csv")
    entry = ('<|"name" -> "%s", "dim" -> %d, "coords" -> "%s", "cells" -> "%s", "patchcoords" -> %s, '
             '"patchcells" -> %s, "bc" -> %s, "energyScale" -> %r, "out" -> "%s", "nodal" -> "%s"|>') % (
        name, c["dim"], wpath(stem + "_coords.csv"), wpath(stem + "_cells.csv"), pc, pcells, bc,
        float(c["energy_scale"]), out, nodal)
    return entry, out, nodal, P


def evaluate_nodal(case, h, skvariant, P, nodal_csv):
    """K and Q2 of NDSolve's nodal displacements, with the scikit-fem basis of the same mesh and element space."""
    c = cases.CASES[case]
    U = np.loadtxt(nodal_csv, delimiter=",", ndmin=2)
    mesh, e = run_skfem.build(case, h, skvariant)
    basis = CellBasis(mesh, e)
    tree = cKDTree(P)
    u = np.zeros(basis.N)
    for k, dofs in enumerate(basis.split_indices()):
        dist, idx = tree.query(basis.doflocs[:, dofs].T)
        if dist.max() > 1e-9:
            raise RuntimeError("degree of freedom without a mesh node (%g)" % dist.max())
        u[dofs] = U[idx, k]
    A = run_skfem.stiffness(basis, c["dim"])
    kc = c["energy_scale"] * float(u @ (A @ u))
    q2 = float("nan")
    if c["q2"]:
        q2, _ = run_skfem.patch_average(case, mesh, e, u)
    return kc, q2


FIELDS = ["code", "case", "element", "h", "ndof", "nelem", "time_s", "integrate_time_s", "K", "K_check",
          "K_nintegrate", "Q2", "Q2_nintegrate", "patch_area_nintegrate", "mesh_measure", "note"]


def memory_gb(dim, ndof, quadratic):
    """estimated peak memory of the kernel for one run (GB), same scale as run_skfem.memory_gb"""
    return run_skfem.memory_gb(dim, ndof, quadratic)


def run_kernel(jobfile):
    cmd = shlex.split(os.environ.get("WOLFRAMSCRIPT", "wolframscript"), posix=True)
    cmd += ["-file", wpath(os.path.join(HERE, "ndsolve_fem.wls")), wpath(jobfile)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.stdout[-3000:] + p.stderr[-1000:]


def main(argv):
    sys.stdout.reconfigure(encoding="utf-8")
    todo = argv or list(cases.CASES)
    # work directory for the mesh and result files of the single runs; kept between invocations, so a run whose
    # kernel result exists there is only post-processed again
    workdir = os.environ.get("ND_WORKDIR") or os.path.join(tempfile.gettempdir(), "benchpack_ndsolve_work")
    os.makedirs(workdir, exist_ok=True)
    for case in todo:
        c = cases.CASES[case]
        log = runlog.CsvLog(os.path.join(RESULTS, "ndsolve_%s.csv" % case), FIELDS)
        for variant, order, curved, skv in VARIANTS:
            if curved and c["mesh"] == "lbracket":
                continue
            n = NLEVELS[case] if not (c["dim"] == 3 and order == 2) else NLEVELS_QUAD_3D
            n = min(n, int(os.environ.get("ND_MAXLEVELS", "99")))
            for h in cases.LEVELS[case][:n]:
                if log.has(element=variant, h="%.12g" % h):
                    continue
                name = "%s_%s_%g" % (case, variant, h)
                entry, out, nodal, P = job_entry(name, case, h, order, curved, workdir)
                row = dict(code="NDSolve", case=case, element=variant, h=h, ndof=c["dim"] * len(P))
                est = memory_gb(c["dim"], c["dim"] * len(P), order == 2)
                if est > runlog.MAX_GB:
                    row["note"] = "not run: estimated %.1f GB, above the %.0f GB cap" % (est, runlog.MAX_GB)
                    log.append(row)
                    print(case, variant, h, row["note"], flush=True)
                    continue
                jobfile = os.path.join(workdir, name + "_job.wl")
                with open(jobfile, "w") as f:
                    f.write("{\n" + entry + "\n}\n")
                for attempt in range(2):
                    if os.path.exists(out):
                        break
                    text = run_kernel(jobfile)
                    if not os.path.exists(out):
                        print(name, "kernel stopped (attempt %d):" % (attempt + 1), text[-500:], flush=True)
                if not os.path.exists(out):
                    row["note"] = "kernel stopped twice"
                    log.append(row)
                    continue
                with open(out) as f:
                    r = list(csv.reader(f))
                row.update(dict(zip(r[0], r[1])))
                row["K"] = row.pop("K_system")
                kc, q2 = evaluate_nodal(case, h, skv, P, nodal)
                row["K_check"] = "%.12g" % kc
                row["Q2"] = "%.12g" % q2
                log.append(row)
                print(case, variant, h, "K", row.get("K"), "K_check", row.get("K_check"), "K_nintegrate",
                      row.get("K_nintegrate"), "Q2", row.get("Q2"), "Q2_nintegrate", row.get("Q2_nintegrate"),
                      flush=True)
    print("work directory:", workdir)


if __name__ == "__main__":
    main(sys.argv[1:])
