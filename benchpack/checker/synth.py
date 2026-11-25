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
