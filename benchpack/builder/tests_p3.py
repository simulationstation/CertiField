"""Tests of the P3 builder (python tests_p3.py [STEM]).

1. Known answer of the exact cut-cell integration: the lift alone, v = (x/2, 0, 0), has the constant energy density
   (lambda + 2 mu)/4, so a_o = (lambda + 2 mu)(1 - pi/48)/4 exactly (volume of the octant minus the ball).
2. Exact identities of a random admissible stress field (N = 4, degree 3, both grip variants, exact rationals):
   equilibrium, tau grad(rho) divisible by rho on every cut cell, zero traction on y = 1 and z = 1, zero shear on the
   mirror planes and (frictionless) on the grip, traction continuity across cell faces, no grip force of the
   Beltrami part.
3. The rigorous evaluator (interval Gauss) agrees with the float assembly on a solved field.
4. With STEM: falsification against the stored certificate (perturbed fields never cross the claim) and its
   re-evaluation.
"""
import limits  # noqa: F401  (thread and process limits; must come before numpy)
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from flint import arb, fmpq

import p3float as pf
from p3float import Setup3, solve_dense, COMPS, G, eta_poly
import p3rig
import p3exact
from p3exact import Pieces3, tensor_cell, cut_cells
from p3cut import X3, Y3, Z3, ZERO3, uni3, outside_integral
from exact import q
import certio as cio

rng = np.random.default_rng(5)
T0 = time.time()


def head(s):
    print(f"[{time.time() - T0:5.0f}s] {s}", flush=True)


def stress_cell(S, Psi, gam, cell):
    ex, ey, ez = cell
    px, py, pz = Pieces3(S.S[0], X3), Pieces3(S.S[1], Y3), Pieces3(S.S[2], Z3)
    fx, fy, fz = [int(S.S[k].first[e]) for k, e in enumerate(cell)]
    rho = uni3(eta_poly(S.S[0].a[ex]), X3) + Y3 * Y3 + Z3 * Z3 - q(Fr(1, 4))
    Phi = [rho * rho * tensor_cell(px(ex), py(ey), pz(ez), Psi[c], fx, fy, fz) for c in range(6)]
    H = [[Phi[c].derivative('xyz'[k]).derivative('xyz'[m]) for c in range(6)] for (k, m) in COMPS]
    tau = []
    for s in range(6):
        t = ZERO3
        for p in range(6):
            for c in range(6):
                if G[s, p, c]:
                    t += H[c][p] * int(G[s, p, c])
        tau.append(t)
    ng = len(S.gexp)
    if Fr(ey, S.N) >= Fr(1, 2):
        tau[0] += sum((q(float(gam[i])) * (Y3 - q(Fr(1, 2))) ** e for i, e in enumerate(S.gexp)), ZERO3)
    if Fr(ez, S.N) >= Fr(1, 2):
        tau[0] += sum((q(float(gam[ng + i])) * (Z3 - q(Fr(1, 2))) ** e for i, e in enumerate(S.gexp)), ZERO3)
    return tau, rho


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    head("1. known answer: energy of the lift (x/2, 0, 0) over the octant minus the ball (exact cut cells)")
    S = Setup3(4, 2, 2, 'frictionless')
    V = np.zeros((3, *S.nu_))
    full = p3exact.full_cell_boxes(S)
    a_full = p3rig.primal_upper(S, V, 0, boxes=full)
    a_cut = p3exact.primal_cut_upper(S, V)
    lam, mu = pf.LAM, pf.MU
    exact = (float(lam) + 2 * float(mu)) / 4 * (1 - np.pi / 48)
    tot = float(a_full) + float(a_cut)
    print(f"   computed upper bound {tot!r}, exact value {exact!r}, difference {tot - exact:.2e}")
    assert 0 <= tot - exact < 1e-12
    head("   ok")

    head("2. exact identities of random admissible stress fields (N = 4, degree 3)")
    for grip in ('frictionless', 'clamped'):
        S = Setup3(4, 3, 2, grip)
        red, nr = S.stress_map()
        pidx, npar = S.param_map()
        z = rng.standard_normal(nr)
        Psi = np.zeros(red.shape)
        Psi[red >= 0] = z[red[red >= 0]]
        Psi = Psi.reshape(6, *S.ns)
        gam = rng.standard_normal(npar)[pidx]
        cells = {}
        cuts = set(cut_cells(S))
        for ex in range(S.N):
            for ey in range(S.N):
                for ez in range(S.N):
                    hi = Fr((ex + 1) ** 2 + (ey + 1) ** 2 + (ez + 1) ** 2, S.N ** 2)
                    if hi > Fr(1, 4):
                        cells[(ex, ey, ez)] = stress_cell(S, Psi, gam, (ex, ey, ez))
        nchk = 0
        for cell, (tau, rho) in cells.items():
            t = {(0, 0): tau[0], (1, 1): tau[1], (2, 2): tau[2], (0, 1): tau[3], (1, 2): tau[4], (0, 2): tau[5]}
            T = lambda i, j: t[(min(i, j), max(i, j))]
            for i in range(3):
                assert sum((T(i, j).derivative('xyz'[j]) for j in range(3)), ZERO3).is_zero(), "equilibrium"
            if cell in cuts:
                grad = [rho.derivative(v) for v in 'xyz']
                for i in range(3):
                    comp = sum((T(i, j) * grad[j] for j in range(3)), ZERO3)
                    assert divmod(comp, rho)[1].is_zero(), "traction on the sphere"
                nchk += 1
            lo = [Fr(c, S.N) for c in cell]
            hi_ = [Fr(c + 1, S.N) for c in cell]
            for k, v in enumerate('xyz'):
                if lo[k] == 0:              # mirror plane: shear zero
                    for j in range(3):
                        if j != k:
                            assert T(k, j).subs({v: q(0)}).is_zero(), ("mirror", k, j)
                if hi_[k] == 1 and k > 0:   # free faces y = 1, z = 1
                    for j in range(3):
                        assert T(k, j).subs({v: q(1)}).is_zero(), ("free", k, j)
                if hi_[k] == 1 and k == 0 and grip == 'frictionless':
                    assert T(0, 1).subs({'x': q(1)}).is_zero() and T(0, 2).subs({'x': q(1)}).is_zero(), "grip"
                nb = list(cell)
                nb[k] += 1
                nb = tuple(nb)
                if nb in cells:             # traction continuity across x_k = hi
                    tn = cells[nb][0]
                    tt = {(0, 0): tn[0], (1, 1): tn[1], (2, 2): tn[2], (0, 1): tn[3], (1, 2): tn[4], (0, 2): tn[5]}
                    for j in range(3):
                        key = (min(k, j), max(k, j))
                        assert (t[key] - tt[key]).subs({v: q(hi_[k])}).is_zero(), ("continuity", cell, k, j)
        # grip force of the whole field equals the closed form
        F = fmpq(0)
        for ey in range(S.N):
            for ez in range(S.N):
                txx = cells[(S.N - 1, ey, ez)][0][0].subs({'x': q(1)})
                p = txx.integral('y')
                p = p.subs({'y': q(Fr(ey + 1, S.N))}) - p.subs({'y': q(Fr(ey, S.N))})
                p = p.integral('z')
                p = p.subs({'z': q(Fr(ez + 1, S.N))}) - p.subs({'z': q(Fr(ez, S.N))})
                F += p.coefficient(0) if not p.is_zero() else 0
        import p3cert
        assert Fr(int(F.p), int(F.q)) == p3cert.force(S, gam), "force"
        print(f"   [{grip}] all identities hold exactly ({len(cells)} cells, {nchk} cut by the sphere)")
    head("   ok")

    head("3. rigorous evaluator against the float assembly (N = 4, degree 3)")
    S = Setup3(4, 3, 2, 'frictionless')
    A, f, red, nr = S.assemble_dual(1)
    pidx, npar = S.param_map()
    z = solve_dense(A, f, protect=npar)[:, 0]
    Psi = np.zeros(red.shape)
    Psi[red >= 0] = z[npar + red[red >= 0]]
    Qr = p3rig.dual_upper(S, Psi.reshape(6, *S.ns), z[pidx], 1)
    Qf = z @ A @ z
    print(f"   Q float {Qf!r}, rigorous upper bound {Qr!r}")
    assert Qf <= Qr * (1 + 1e-12) and Qr <= Qf * (1 + 1e-9)
    head("   ok")

    if len(sys.argv) > 1:
        stem = sys.argv[1]
        head(f"4. falsification against {stem}")
        meta, arr = cio.load(stem)
        S = Setup3(meta['N'], meta['stress_space']['degree'], meta['displacement_space']['degree'], meta['variant'],
                   grade=meta.get('grade'))
        Klo, Khi = Fr(meta['claims']['Q1']['lower']), Fr(meta['claims']['Q1']['upper'])
        integ = meta['integration']
        for eps in (0.0, 1e-3):
            Psi = arr['Psi'] * (1 + eps * rng.standard_normal(arr['Psi'].shape) * (arr['Psi'] != 0))
            red, _ = S.stress_map()
            # restore the ties after the perturbation: average over each tie class
            flat = Psi.reshape(-1)
            vals = np.zeros(red.max() + 1)
            cnt = np.zeros(red.max() + 1)
            np.add.at(vals, red[red >= 0], flat[red >= 0])
            np.add.at(cnt, red[red >= 0], 1)
            flat[red >= 0] = (vals / np.maximum(cnt, 1))[red[red >= 0]]
            flat[red < 0] = 0
            Psi = flat.reshape(arr['Psi'].shape)
            gam = arr['gam']
            if integ['mode'] == 'exact':
                Q = Fr(p3rig.dual_upper(S, Psi, gam, 0, boxes=p3exact.full_cell_boxes(S))) + \
                    p3exact.dual_cut_upper(S, Psi, gam)
            else:
                Q = Fr(p3rig.dual_upper(S, Psi, gam, integ['stress_depth']))
            F = p3cert.force(S, gam)
            Kl = 2 * F * F / Q
            ok = Kl <= Khi
            print(f"   eps={eps:g}: lower bound {float(Kl):.8f} vs claim upper {float(Khi):.8f}: {'no crossing' if ok else 'CROSSING'}")
            assert ok
        head("   ok")
        head(f"5. re-evaluation of {stem}")
        assert p3cert.verify(stem)
    head("ALL P3 TESTS PASSED")


if __name__ == '__main__':
    main()
