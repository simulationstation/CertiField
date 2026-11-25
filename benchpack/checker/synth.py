"""Synthetic certificates for the checker's own tests (random admissible fields in the format of FORMAT.md, and the
exact fields of the plate without the hole).  Test use only."""
from fractions import Fraction as Fr

import numpy as np

import exact as X


def open_knots(p, breaks, mult=None):
    """open knot vector of degree p over the breakpoints (first/last repeated p+1 times), interior multiplicity mult"""
    mult = mult or {}
    kn = [breaks[0]] * (p + 1)
    for b in breaks[1:-1]:
        kn += [b] * mult.get(b, 1)
    kn += [breaks[-1]] * (p + 1)
    return kn


def greville(kn, p):
    n = len(kn) - p - 1
    return [sum(kn[j + 1:j + p + 1], Fr(0)) / p for j in range(n)]


def rnd(rng, shape, scale=1.0):
    """random float64 array, as exact Fractions"""
    return X.farr(rng.standard_normal(shape) * scale)


def space(p, kx, ky, **extra):
    d = dict(degree=p, knots_x=[str(v) for v in kx], knots_y=[str(v) for v in ky])
    d.update(extra)
    return d


def make_P1(rng, breaks_x, breaks_y, pu=2, ps=3, gexp=(0, 1, 2), variant='frictionless', mult_s=None, scale=0.3):
    """P1 quarter model, random admissible fields v, tau (and the same for Q2 fields)"""
    kxu, kyu = open_knots(pu, breaks_x), open_knots(pu, breaks_y)
    kxs, kys = open_knots(ps, breaks_x, mult_s), open_knots(ps, breaks_y, mult_s)
    assert Fr(1, 2) in breaks_y
    meta = dict(format='benchpack-certificate-1', problem='P1', variant=variant, model='quarter',
                displacement_space=space(pu, kxu, kyu),
                stress_space=space(ps, kxs, kys, g_exponents=list(gexp), g_base='1/2'),
                fields=dict(Q1=dict(v='v', tau='tau'),
                            Q2=dict(kappa='1/1000', v_plus='vp', v_minus='vm', tau_plus='tp', tau_minus='tm')),
                label='test')
    arrays = {}
    nxu, nyu = len(kxu) - pu - 1, len(kyu) - pu - 1
    nxs, nys = len(kxs) - ps - 1, len(kys) - ps - 1
    for nm in ('v', 'vp', 'vm'):
        A = rnd(rng, (nxu, nyu), 0.1)
        B = rnd(rng, (nxu, nyu), 0.1)
        A[0, :] = Fr(0)
        A[-1, :] = Fr(0)
        B[:, 0] = Fr(0)
        if variant == 'clamped':
            B[-1, :] = Fr(0)
        arrays[nm + '_A'], arrays[nm + '_B'] = A, B
    for nm in ('tau', 'tp', 'tm'):
        C = rnd(rng, (nxs, nys), scale)
        C[0, :] = C[1, :]
        if variant != 'clamped':
            C[-1, :] = C[-2, :]
        C[:, 0] = C[:, 1]
        C[:, -1] = Fr(0)
        C[:, -2] = Fr(0)
        arrays[nm + '_C'] = C
        arrays[nm + '_par'] = rnd(rng, (len(gexp),), 1.0)
    return meta, arrays


def make_nohole(kappa='1/1000'):
    """the plate without the hole (test-only variant): v = (x/2, -3y/20), tau = (1/2) e_x e_x  (K = 1, Q2 = 1/2)"""
    pu, ps = 2, 3
    br = [Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(1)]
    kxu = kyu = open_knots(pu, br)
    kxs = kys = open_knots(ps, br)
    meta = dict(format='benchpack-certificate-1', problem='P1', variant='nohole-test', model='quarter',
                displacement_space=space(pu, kxu, kyu),
                stress_space=space(ps, kxs, kys, g_exponents=[0], g_base='0'),
                fields=dict(Q1=dict(v='v', tau='tau'),
                            Q2=dict(kappa=kappa, v_plus='v', v_minus='v', tau_plus='tau', tau_minus='tau')),
                label='test')
    n = len(kxu) - pu - 1
    A = X.zeros((n, n))
    gy = greville(kyu, pu)
    B = X.zeros((n, n))
    for i in range(n):
        for j in range(n):
            B[i, j] = Fr(-3, 20) * gy[j]
    ns = len(kxs) - ps - 1
    arrays = dict(v_A=A, v_B=B, tau_C=X.zeros((ns, ns)), tau_par=np.array([Fr(1, 2)], dtype=object))
    return meta, arrays


def p2_classes(kn, p):
    n = len(kn) - p - 1
    L, R, S = [], [], []
    for i in range(n):
        a, b = kn[i], kn[i + p + 1]
        (L if b <= 0 else R if a >= 0 else S).append(i)
    return L, R, S


def make_P2(rng, breaks, pu=2, ps=3, scale=0.3):
    """P2 full model, random admissible fields"""
    kxu = kyu = open_knots(pu, breaks)
    kxs = kys = open_knots(ps, breaks, {Fr(0): ps - 1})
    meta = dict(format='benchpack-certificate-1', problem='P2', variant='standard', model='full',
                displacement_space=space(pu, kxu, kyu), stress_space=space(ps, kxs, kys),
                fields=dict(Q1=dict(v='v', tau='tau'),
                            Q2=dict(kappa='1/100', v_plus='vp', v_minus='vm', tau_plus='tp', tau_minus='tm')),
                label='test')
    arrays = {}
    nu = len(kxu) - pu - 1
    ns = len(kxs) - ps - 1
    for nm in ('v', 'vp', 'vm'):
        A, B = rnd(rng, (nu, nu), 0.1), rnd(rng, (nu, nu), 0.1)
        for M in (A, B):
            M[0, :] = Fr(0)
            M[-1, :] = Fr(0)
        arrays[nm + '_A'], arrays[nm + '_B'] = A, B
    Lx, Rx, Sx = p2_classes(kxs, ps)
    Ly, Ry, Sy = p2_classes(kys, ps)
    for nm in ('tau', 'tp', 'tm'):
        C = rnd(rng, (ns, ns), scale)
        for i in range(ns):
            for j in range(ns):
                if (j in (ns - 2, ns - 1) or ((i in Lx or i in Sx) and j in (0, 1)) or (i in Sx and (j in Ly or j in Sy))
                        or (i in Rx and (j in Sy or j in Ly))):
                    C[i, j] = Fr(0)
        arrays[nm + '_C'] = C
        arrays[nm + '_par'] = rnd(rng, (3,), 1.0)
    return meta, arrays


# ------------------------------------------------------------------------------------------------- P3 (FORMAT.md 9)
P3_COMP = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]


def p3_breaks(N, grade=None):
    """cell edges of FORMAT.md 9: k/N, and with grade (L, m) also 1 - (1/N) 2^-(i+1) (1 + j/m)"""
    br = {Fr(i, N) for i in range(N + 1)}
    if grade:
        L, m = grade
        br |= {1 - Fr(1, N) * Fr(1, 2 ** (i + 1)) * (1 + Fr(j, m)) for i in range(L) for j in range(m)}
    return sorted(br)


def p3_space(p, N, extra_x=None, grade=None, **kw):
    br = p3_breaks(N, grade)
    kx = open_knots(p, br, extra_x)
    ky = open_knots(p, br)
    d = dict(degree=p, knots=[[str(v) for v in kx], [str(v) for v in ky], [str(v) for v in ky]])
    d.update(kw)
    return d, [len(kx) - p - 1, len(ky) - p - 1, len(ky) - p - 1]


def p3_conditions(Psi, frictionless):
    """impose the parity, grip and free-face conditions of FORMAT.md 9 on a (6, nx, ny, nz) array"""
    for c, (a, b) in enumerate(P3_COMP):
        A = Psi[c]
        for ax in range(3):
            even = ((a == ax) + (b == ax)) % 2 == 0
            s0 = [slice(None)] * 3
            s1 = [slice(None)] * 3
            s0[ax], s1[ax] = 0, 1
            A[tuple(s0)] = A[tuple(s1)] if even else Fr(0)
        if frictionless:
            even = ((a == 0) + (b == 0)) % 2 == 0
            A[-1] = A[-2] if even else Fr(0)
        for ax in (1, 2):
            for idx in (-1, -2):
                s = [slice(None)] * 3
                s[ax] = idx
                A[tuple(s)] = Fr(0)
    return Psi


def make_P3(rng, N=4, pu=2, ps=3, gexp=(0, 1, 2), variant='frictionless', cavity=True, L=2, Lu=2, scale=0.3,
            grade=None):
    """random admissible P3 fields (octant model) in the format of FORMAT.md 9"""
    fr = variant in ('frictionless', 'nocavity-test')
    dsp, dn = p3_space(pu, N, grade=grade)
    ssp, sn = p3_space(ps, N, {Fr(3, 4): ps - 1} if cavity else None, grade=grade, g_exponents=list(gexp),
                       g_base='1/2' if cavity else '0', eta_flat_from='3/4')
    meta = dict(format='benchpack-certificate-1', problem='P3',
                variant=variant if cavity else ('nocavity-test' if fr else 'nocavity-clamped-test'),
                model='octant', N=N, displacement_space=dsp, stress_space=ssp,
                integration=dict(mode='superset', stress_depth=L, displacement_depth=Lu), label='test')
    if grade:
        meta['grade'] = list(grade)
    V = rnd(rng, (3,) + tuple(dn), 0.1)
    V[0][0] = Fr(0)
    V[0][-1] = Fr(0)
    V[1][:, 0, :] = Fr(0)
    V[2][:, :, 0] = Fr(0)
    if not fr:
        V[1][-1] = Fr(0)
        V[2][-1] = Fr(0)
    Psi = p3_conditions(rnd(rng, (6,) + tuple(sn), scale), fr)
    gam = rnd(rng, (2 * len(gexp),), 1.0)
    return meta, dict(V=V, Psi=Psi, gam=gam)


def make_P3_solid(N=4):
    """the solid cube with frictionless grips (test-only variant): tau = (1/2) e_x e_x, v = (x/2, -nu y/2, -nu z/2);
    K = 2 exactly"""
    pu, ps = 2, 3
    dsp, dn = p3_space(pu, N)
    ssp, sn = p3_space(ps, N, None, g_exponents=[0], g_base='0', eta_flat_from='3/4')
    meta = dict(format='benchpack-certificate-1', problem='P3', variant='nocavity-test', model='octant', N=N,
                displacement_space=dsp, stress_space=ssp, label='test')
    ky = [Fr(v) for v in dsp['knots'][1]]
    g = greville(ky, pu)
    V = X.zeros((3,) + tuple(dn))
    for i in range(dn[0]):
        for j in range(dn[1]):
            for k in range(dn[2]):
                V[1][i, j, k] = Fr(-3, 20) * g[j]
                V[2][i, j, k] = Fr(-3, 20) * g[k]
    return meta, dict(V=V, Psi=X.zeros((6,) + tuple(sn)), gam=np.array([Fr(1, 4), Fr(1, 4)], dtype=object))
