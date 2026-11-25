"""P3 certificates with higher-degree splines (block with a spherical cavity; see p3float.py for the model).

usage:
    python p3cert.py build STEM N ps pu L_solve L_cert Lu_solve Lu_cert [frictionless|clamped]
To resume after a stop rerun the same command: the solved fields and the partial sums of the certification are kept in
benchpack/work/STEM and reused.  The stored curved3d certificates (Reference/Curved3D/certificates) are not touched.
"""
import limits  # noqa: F401  (thread and process limits; must come before numpy)
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np

import p3float as pf
from p3float import Setup3, solve_dense
import p3rig
import certio as cio


def _fields(S, Ls, Lus, verbose, centre=False):
    t0 = time.time()
    A, f, red, nr = S.assemble_dual(Ls, verbose=verbose, centre=centre)
    pidx, npar = S.param_map()
    z = solve_dense(A, f, protect=npar)[:, 0]
    if verbose:
        print(f"  dual solved: {A.shape[0]} dofs, float K >= {2 * (f @ z) ** 2 / (z @ A @ z):.8f}  "
              f"({time.time() - t0:.0f}s)", flush=True)
    del A
    Psi = np.zeros(red.shape)
    m = red >= 0
    Psi[m] = z[npar + red[m]]
    gam = z[pidx].copy()
    t1 = time.time()
    K, b, k0, redu, nru = S.assemble_primal(Lus, centre=centre)
    v = solve_dense(K, -b)[:, 0]
    if verbose:
        print(f"  primal solved: {K.shape[0]} dofs, float K <= {8 * (v @ K @ v + 2 * b @ v + k0):.8f}  "
              f"({time.time() - t1:.0f}s)", flush=True)
    V = np.zeros(redu.shape)
    mu = redu >= 0
    V[mu] = v[redu[mu]]
    return Psi.reshape(6, *S.ns), gam, V.reshape(3, *S.nu_)


def check_ties(S, Psi, gam, V):
    """exact checks of the coefficient conditions (zeros and equalities of stored floats)"""
    red, _ = S.stress_map()
    flat = Psi.reshape(-1)
    assert np.all(flat[red < 0] == 0)
    # coefficients that share a reduced index must be equal
    order = np.argsort(red)
    r, vals = red[order], flat[order]
    same = (r[1:] == r[:-1]) & (r[1:] >= 0)
    assert np.array_equal(vals[1:][same], vals[:-1][same])
    redu, _ = S.disp_map()
    assert np.all(V.reshape(-1)[redu < 0] == 0)
    assert np.all(np.isfinite(Psi)) and np.all(np.isfinite(gam)) and np.all(np.isfinite(V))


def force(S, gam):
    """exact F_oct = int over the grip face of tau_xx (the Beltrami part carries none)"""
    ng = len(S.gexp)
    return sum((Fr(float(gam[i])) + Fr(float(gam[ng + i]))) * Fr(1, 2) ** (e + 1) / (e + 1)
               for i, e in enumerate(S.gexp))


def build(stem, N, ps, pu, Ls, Lc, Lus, Luc, grip='frictionless', verbose=True, mode='superset'):
    """mode 'superset': curved3d's octree boxes (depths Lc, Luc); mode 'exact': full cells by Gauss rules and the
    cells cut by the sphere exactly (p3exact); Lc, Luc are then unused"""
    t0 = time.time()
    S = Setup3(N, ps, pu, grip)
    wd = os.path.join(limits.WORK_DIR, stem)
    os.makedirs(wd, exist_ok=True)
    fp = json.dumps([N, ps, pu, Ls, Lus, grip, list(S.gexp)] + ([] if mode == "superset" else [mode]))
    fn = os.path.join(wd, 'fields.npz')
    if verbose:
        print(f"  P3 {grip}: N={N}, degrees {ps}/{pu}; to resume after a stop rerun the same command "
              f"(checkpoints in benchpack/work/{stem})", flush=True)
    if os.path.exists(fn) and str(np.load(fn)['fp']) == fp:
        d = np.load(fn)
        Psi, gam, V = d['Psi'], d['gam'], d['V']
        if verbose:
            print("  solved fields loaded from the checkpoint", flush=True)
    else:
        for old in os.listdir(wd):
            os.remove(os.path.join(wd, old))
        Psi, gam, V = _fields(S, Ls, Lus, verbose, centre=(mode == 'exact'))
        np.savez(fn + '.tmp.npz', Psi=Psi, gam=gam, V=V, fp=fp)
        os.replace(fn + '.tmp.npz', fn)
    check_ties(S, Psi, gam, V)
    if mode == 'exact':
        import p3exact
        full = p3exact.full_cell_boxes(N)
        a_up = Fr(p3rig.primal_upper(S, V, 0, verbose=verbose, checkpoint=os.path.join(wd, 'primal_full.json'),
                                     boxes=full))
        a_up += p3exact.primal_cut_upper(S, V, checkpoint=os.path.join(wd, 'primal_cut.json'), verbose=verbose)
    else:
        a_up = Fr(p3rig.primal_upper(S, V, Luc, verbose=verbose, checkpoint=os.path.join(wd, f'primal_L{Luc}.json')))
    if verbose:
        print(f"  certified upper bound K <= {float(8 * a_up):.8f}  ({time.time() - t0:.0f}s)", flush=True)
    if mode == 'exact':
        Q_up = Fr(p3rig.dual_upper(S, Psi, gam, 0, verbose=verbose, checkpoint=os.path.join(wd, 'dual_full.json'),
                                   boxes=full))
        Q_up += p3exact.dual_cut_upper(S, Psi, gam, checkpoint=os.path.join(wd, 'dual_cut.json'), verbose=verbose)
    else:
        Q_up = Fr(p3rig.dual_upper(S, Psi, gam, Lc, verbose=verbose, checkpoint=os.path.join(wd, f'dual_L{Lc}.json')))
    F = force(S, gam)
    K_lo = 2 * F * F / Q_up
    K_hi = 8 * a_up
    assert K_lo <= K_hi
    cK = cio.claim_strings(K_lo, K_hi)
    meta = dict(
        format='benchpack-certificate-1', problem='P3', variant=grip, model='octant', N=N,
        displacement_space=dict(degree=pu, knots=[[cio.fstr(t) for t in S.U[k].t] for k in range(3)]),
        stress_space=dict(degree=ps, knots=[[cio.fstr(t) for t in S.S[k].t] for k in range(3)],
                          g_exponents=list(S.gexp), g_base='1/2', rho='eta(x) + y^2 + z^2 - 1/4',
                          eta_flat_from='3/4'),
        integration=(dict(mode='superset', stress_depth=Lc, displacement_depth=Luc) if mode == 'superset'
                     else dict(mode='exact')),
        claims=dict(Q1=dict(lower=cK[0], upper=cK[1])),
        computed=dict(Q_oct_upper=float(Q_up), F_oct=cio.fstr(F), a_oct_upper=float(a_up),
                      K=[float(K_lo), float(K_hi)]),
        label=cio.LABEL,
    )
    cio.save(stem, meta, dict(Psi=Psi, gam=gam, V=V))
    if verbose:
        print(f"  CERTIFIED  {float(K_lo):.8f} <= K <= {float(K_hi):.8f}   rel width "
              f"{float((K_hi - K_lo) / K_hi):.2e};  claim [{cK[0]}, {cK[1]}]  ({time.time() - t0:.0f}s)", flush=True)
    return meta


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    if sys.argv[1] == 'build':
        a = sys.argv[2:]
        build(a[0], *[int(v) for v in a[1:8]], grip=(a[8] if len(a) > 8 else 'frictionless'),
              mode=(a[9] if len(a) > 9 else 'superset'))
