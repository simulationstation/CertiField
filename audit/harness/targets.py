"""Targets of the audit: the library call, its documented domain, samplers, a grid and the
reference.  Library modules are imported lazily, so worker processes (which only compute
references) never import scipy or torch.
"""
import math

import numpy as np
from mpmath import mp, mpf

import refs

# =======================================================================================
# samplers: each generator is g(rng, n) -> float64 array
# =======================================================================================


def LU(lo, hi):
    return lambda rng, n: np.exp(rng.uniform(math.log(lo), math.log(hi), n))


def SLU(lo, hi):
    return lambda rng, n: LU(lo, hi)(rng, n) * rng.choice([-1.0, 1.0], n)


def U(lo, hi):
    return lambda rng, n: rng.uniform(lo, hi, n)


def C(*vals):
    return lambda rng, n: rng.choice(np.array(vals, dtype=float), n)


def ONEM(lo, hi):
    """1 - LU(lo, hi): values just below 1."""
    return lambda rng, n: 1.0 - LU(lo, hi)(rng, n)


def NEAR(c, lo=1e-17, hi=1e-2):
    return lambda rng, n: c * (1.0 + SLU(lo, hi)(rng, n))


def NEARINT(kmin, kmax, sign=1.0, lo=1e-16, hi=0.5):
    def g(rng, n):
        k = np.floor(LU(max(kmin, 1), kmax + 1)(rng, n))
        if kmin == 0:
            k = np.where(rng.random(n) < 0.1, 0.0, k)
        return sign * k + SLU(lo, hi)(rng, n)
    return g


def INT(lo, hi):
    return lambda rng, n: np.floor(LU(lo, hi + 1)(rng, n))


def M(*pairs):
    w = np.array([p[0] for p in pairs], dtype=float)
    w /= w.sum()

    def g(rng, n):
        idx = rng.choice(len(pairs), n, p=w)
        out = np.empty(n)
        for i, (_, gen) in enumerate(pairs):
            m = idx == i
            if m.any():
                out[m] = gen(rng, int(m.sum()))
        return out
    return g


UNIT = M((4, LU(1e-300, 0.5)), (3, ONEM(1.2e-16, 0.5)), (3, U(0, 1)))
PROB = UNIT
PAR = M((1.5, LU(1e-6, 1e-2)), (3.5, LU(1e-2, 10)), (3, LU(10, 1e4)), (2, LU(1e4, 1e8)))
PAR_HALF = M((8, PAR), (2, lambda rng, n: np.floor(LU(1, 1e4)(rng, n)) / 2))  # half-integers
DF = M((1.5, LU(1e-3, 1)), (3.5, LU(1, 100)), (3, LU(100, 1e8)), (2, LU(1e8, 1e300)))
REAL_WIDE = M((4, SLU(1e-300, 1e300)), (3, U(-40, 40)), (3, SLU(1e-3, 1e3)))

G_PAR = [1e-6, 1e-3, 0.1, 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7]
G_UNIT = [1e-300, 1e-100, 1e-30, 1e-10, 1e-5, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99,
          1 - 1e-5, 1 - 1e-10, 1 - 2.0 ** -52]
G_REAL = sorted({s * v for s in (-1.0, 1.0) for v in
                 [1e-300, 1e-100, 1e-20, 1e-8, 1e-3, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 10.0,
                  15.0, 20.0, 26.0, 27.0, 30.0, 37.0, 38.0, 40.0, 50.0, 100.0, 1e3, 1e5, 1e10,
                  1e20, 1e100, 1e300]})
G_DF = [1e-3, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 1e3, 1e5, 1e7, 1e10, 1e15, 1e20,
        1e50, 1e100, 1e300]


def _prod(*lists):
    out = [()]
    for L in lists:
        out = [o + (v,) for o in out for v in L]
    return out


# =======================================================================================
# library call wrappers
# =======================================================================================


def sc(name):
    def f(*a):
        import scipy.special as s
        return np.asarray(getattr(s, name)(*a), dtype=float)
    f.__name__ = name
    return f


def st(dist, method):
    def f(x, *shape):
        import scipy.stats as s
        return np.asarray(getattr(getattr(s, dist), method)(x, *shape), dtype=float)
    return f


def th(name, dtype):
    def f(*a):
        import torch
        tdt = torch.float32 if dtype == 'float32' else torch.float64
        ndt = np.float32 if dtype == 'float32' else np.float64
        fn = getattr(torch.special, name)
        if name == 'polygamma':
            n, x = a
            out = np.empty(len(x))
            for k in np.unique(n):
                m = n == k
                out[m] = fn(int(k), torch.from_numpy(np.ascontiguousarray(x[m], dtype=ndt))).double().numpy()
            return out
        ts = [torch.from_numpy(np.ascontiguousarray(v, dtype=ndt)) for v in a]
        return fn(*ts).to(torch.float64).numpy()
    return f


# =======================================================================================
# Target
# =======================================================================================


class Target:
    def __init__(self, name, lib, dtype, call, kinds, gens=None, sampler=None, grid=(),
                 valid=None, mref=None, aref=None, malt=None, inv=None, out='value',
                 osc=None, note=''):
        self.name, self.lib, self.dtype, self.call = name, lib, dtype, call
        self.kinds = kinds
        self.nargs = len(kinds)
        self.gens, self._sampler, self.grid = gens, sampler, list(grid)
        self.valid = valid or (lambda *a: True)
        self.mref, self.aref, self.malt, self.inv, self.out = mref, aref, malt, inv, out
        self.osc = osc
        self.note = note

    def sample(self, rng, n):
        if self._sampler is not None:
            return self._sampler(rng, n)
        return np.stack([g(rng, n) for g in self.gens], axis=1)


REG = {}


def add(t):
    REG[t.name] = t
    return t


def get(name):
    if not REG:
        build()
    return REG[name]


def names():
    if not REG:
        build()
    return list(REG)


# ---------------------------------------------------------------------------------------
# domain checks
# ---------------------------------------------------------------------------------------
fin = math.isfinite


def pos(*v):
    return all(fin(x) and x > 0 for x in v)


def unit_open(x):
    return fin(x) and 0 < x < 1


def is_nonpos_int(x):
    return x <= 0 and x == math.floor(x)


# ---------------------------------------------------------------------------------------
# special samplers
# ---------------------------------------------------------------------------------------


def _digamma_roots(kmax=40):
    mp.dps = 25
    roots = [1.4616321449683623]
    for k in range(kmax):
        lo, hi = mpf(-k - 1) + mpf('1e-20'), mpf(-k) - mpf('1e-20')
        r = refs.solve_monotone(mp.digamma, mpf(0), lo, hi, True)
        roots.append(float(r))
    return roots


def _lgamma_zeros(roots, kmax=40):
    """Zeros of log|Gamma| on the negative axis, two per interval (-k-1, -k), k >= 2."""
    mp.dps = 25
    zs = []
    f = lambda x: mp.re(mp.loggamma(x))
    for k in range(2, kmax):
        r = mpf(roots[k + 1])  # minimum of log|Gamma| on (-k-1, -k)
        if f(r) > 0:
            continue
        a = refs.solve_monotone(f, mpf(0), mpf(-k - 1) + mpf(2) ** -60, r, False)
        b = refs.solve_monotone(f, mpf(0), r, mpf(-k) - mpf(2) ** -60, True)
        zs += [float(z) for z in (a, b) if z is not None]
    return zs


_CACHE = {}


def special_points(kind):
    if kind not in _CACHE:
        if kind == 'digamma_roots':
            _CACHE[kind] = _digamma_roots()
        elif kind == 'lgamma_zeros':
            _CACHE[kind] = _lgamma_zeros(special_points('digamma_roots'))
    return _CACHE[kind]


def NEAR_SET(kind, lo=1e-17, hi=1e-3):
    def g(rng, n):
        pts = np.array(special_points(kind))
        c = rng.choice(pts, n)
        return c * (1.0 + SLU(lo, hi)(rng, n))
    return g


def tail_x(dist_name, shape_gens, scale_free=None, side='both'):
    """Sampler (x, *shape) with x placed in the tails by the library's own isf/ppf at a
    log-uniform tail probability in [1e-300, 0.5]; the library's inaccuracy only moves
    the sample point."""
    def g(rng, n):
        import scipy.stats as s
        d = getattr(s, dist_name)
        shapes = [gen(rng, n) for gen in shape_gens]
        q = LU(1e-300, 0.5)(rng, n)
        with np.errstate(all='ignore'):
            xr = np.asarray(d.isf(q, *shapes), dtype=float)
            xl = np.asarray(d.ppf(q, *shapes), dtype=float)
        sd = rng.random(n) < 0.5 if side == 'both' else np.full(n, side == 'right')
        x = np.where(sd, xr, xl)
        bad = ~np.isfinite(x)
        if scale_free is not None:
            alt = scale_free(rng, n)
            x = np.where(bad | (rng.random(n) < 0.2), alt, x)
        return np.stack([x] + shapes, axis=1)
    return g


# =======================================================================================
# registry
# =======================================================================================


def build():
    if REG:
        return
    _build_scipy_special()
    _build_scipy_stats()
    for dt in ('float64', 'float32'):
        _build_torch(dt)


def _beta_sampler(rng, n):
    a, b = PAR(rng, n), PAR(rng, n)
    x = UNIT(rng, n)
    # 40%: x in the transition region around the mean a/(a+b)
    m = rng.random(n) < 0.4
    mu = a / (a + b)
    sd = np.sqrt(a * b / ((a + b) ** 2 * (a + b + 1)))
    z = rng.normal(0, 1, n) * LU(0.1, 40)(rng, n)
    xt = mu + z * sd
    ok = m & (xt > 0) & (xt < 1)
    x = np.where(ok, xt, x)
    return np.stack([a, b, x], axis=1)


def _gamma_sampler(rng, n):
    a = M((1.5, LU(1e-6, 1e-2)), (3.5, LU(1e-2, 10)), (3, LU(10, 1e5)), (2, LU(1e5, 1e12)))(rng, n)
    x = LU(1e-300, 1e300)(rng, n)
    m = rng.random(n) < 0.6
    z = rng.normal(0, 1, n) * LU(0.1, 40)(rng, n)
    xt = a + z * np.sqrt(a)
    x = np.where(m & (xt > 0), xt, x)
    m2 = rng.random(n) < 0.2
    x = np.where(m2, a * LU(1e-4, 1e4)(rng, n), x)
    return np.stack([a, x], axis=1)


def _bdtr_sampler(rng, n):
    N = INT(1, 1e9)(rng, n)
    p = M((4, LU(1e-12, 0.5)), (3, ONEM(1e-12, 0.5)), (3, U(0, 1)))(rng, n)
    k1 = np.floor(rng.random(n) * N)
    sd = np.sqrt(N * p * (1 - p))
    k2 = np.floor(N * p + U(-40, 40)(rng, n) * sd)
    k = np.where(rng.random(n) < 0.5, k1, k2)
    k = np.clip(k, 0, N - 1)
    return np.stack([k, N, p], axis=1)


def _inv_dict(solve, fwd_a, increasing, yi=None):
    d = dict(solve=solve, fwd_a=fwd_a, increasing=increasing)
    if yi is not None:
        d['yi'] = yi
    return d


def _build_scipy_special():
    S = 'scipy'
    D = 'float64'
    sp = 'scipy.special.'
    grid3 = _prod(G_PAR, G_PAR, G_UNIT)
    add(Target(sp + 'betainc', S, D, sc('betainc'), ('pos', 'pos', 'unit'), sampler=_beta_sampler,
               grid=grid3, valid=lambda a, b, x: pos(a, b) and 0 <= x <= 1,
               mref=refs.m_betainc, aref=refs.A_betainc, malt=refs.mb_betainc))
    add(Target(sp + 'betaincc', S, D, sc('betaincc'), ('pos', 'pos', 'unit'), sampler=_beta_sampler,
               grid=grid3, valid=lambda a, b, x: pos(a, b) and 0 <= x <= 1,
               mref=refs.m_betaincc, aref=refs.A_betaincc, malt=refs.mb_betaincc))
    add(Target(sp + 'betaincinv', S, D, sc('betaincinv'), ('pos', 'pos', 'unit'),
               gens=(PAR, PAR, PROB), grid=grid3,
               valid=lambda a, b, y: pos(a, b) and unit_open(y),
               inv=_inv_dict(refs.inv_betaincinv, refs.A_betainc, True)))
    grid2 = _prod(G_PAR + [1e8, 1e10, 1e12], [1e-300, 1e-100, 1e-20, 1e-5, 0.01, 0.1, 0.5, 1.0, 2.0,
                                               5.0, 10.0, 30.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e8, 1e10,
                                               1e12, 1e15, 1e100, 1e300])
    add(Target(sp + 'gammainc', S, D, sc('gammainc'), ('pos', 'pos'), sampler=_gamma_sampler,
               grid=grid2, valid=lambda a, x: pos(a) and fin(x) and x >= 0,
               mref=refs.m_gammainc, aref=refs.A_gammainc, malt=refs.mb_gammainc))
    add(Target(sp + 'gammaincc', S, D, sc('gammaincc'), ('pos', 'pos'), sampler=_gamma_sampler,
               grid=grid2, valid=lambda a, x: pos(a) and fin(x) and x >= 0,
               mref=refs.m_gammaincc, aref=refs.A_gammaincc, malt=refs.mb_gammaincc))
    add(Target(sp + 'gammaincinv', S, D, sc('gammaincinv'), ('pos', 'unit'),
               gens=(M((1.5, LU(1e-6, 1e-2)), (3.5, LU(1e-2, 10)), (3, LU(10, 1e5)), (2, LU(1e5, 1e12))), PROB),
               grid=_prod(G_PAR + [1e8, 1e10, 1e12], G_UNIT[:-1] + [1 - 1e-15]),
               valid=lambda a, y: pos(a) and unit_open(y),
               inv=_inv_dict(refs.inv_gammaincinv, refs.A_gammainc, True)))
    LG = M((3, REAL_WIDE), (2, NEARINT(0, 170, -1.0, 1e-300, 0.5)), (1, NEAR(1.0)), (1, NEAR(2.0)),
           (2, NEAR_SET('lgamma_zeros', 1e-17, 1e-3)), (1, LU(1e-300, 1e300)))
    gx = [(v,) for v in G_REAL] + [(v,) for v in (1.0 + 2.0 ** -52, 1.0 - 2.0 ** -53, 2.0 + 2.0 ** -51, -0.5, -1.5, -2.5,
                                                    -2.4570247382208006, -2.7476826467274127, -170.5, 171.62, 2.5e305)]
    add(Target(sp + 'gammaln', S, D, sc('gammaln'), ('real',), gens=(LG,), grid=gx,
               valid=lambda x: fin(x) and not is_nonpos_int(x),
               mref=refs.m_gammaln, aref=refs.A_gammaln))
    add(Target(sp + 'loggamma', S, D, sc('loggamma'), ('pos',),
               gens=(M((3, LU(1e-300, 1e300)), (1, NEAR(1.0)), (1, NEAR(2.0)), (2, LU(1e-3, 30))),),
               grid=[g for g in gx if g[0] > 0], valid=lambda x: pos(x),
               mref=refs.m_gammaln, aref=refs.A_gammaln, note='real x > 0'))
    DG = M((3, REAL_WIDE), (2, NEARINT(0, 1e6, -1.0, 1e-300, 0.5)), (3, NEAR_SET('digamma_roots', 1e-17, 1e-3)),
           (1, LU(1e-300, 1e300)))
    add(Target(sp + 'digamma', S, D, sc('digamma'), ('real',), gens=(DG,),
               grid=[(v,) for v in G_REAL] + [(1.4616321449683623,), (1.4616321449683622,), (-0.5040830082644554,)],
               valid=lambda x: fin(x) and not is_nonpos_int(x),
               mref=refs.m_digamma, aref=refs.A_digamma))
    EX = M((4, SLU(1e-300, 1e300)), (4, U(-30, 30)), (2, SLU(0.1, 1e4)))
    gr = [(v,) for v in G_REAL]
    for nm, mr, ar in (('erf', refs.m_erf, refs.A_erf), ('erfc', refs.m_erfc, refs.A_erfc),
                       ('erfcx', refs.m_erfcx, refs.A_erfcx), ('ndtr', refs.m_ndtr, refs.A_ndtr),
                       ('log_ndtr', refs.m_log_ndtr, refs.A_log_ndtr)):
        add(Target(sp + nm, S, D, sc(nm), ('real',), gens=(EX,), grid=gr, valid=lambda x: fin(x),
                   mref=mr, aref=ar))
    P2 = M((4, LU(1e-300, 1.0)), (3, lambda rng, n: 2.0 - LU(1.2e-16, 1.0)(rng, n)), (3, U(0, 2)))
    add(Target(sp + 'erfcinv', S, D, sc('erfcinv'), ('p2',), gens=(P2,),
               grid=[(v,) for v in [1e-300, 1e-100, 1e-20, 1e-5, 0.1, 0.5, 1.0, 1.5, 1.9, 2 - 1e-10, 2 - 2.0 ** -52]],
               valid=lambda y: fin(y) and 0 < y < 2,
               inv=_inv_dict(refs.inv_erfcinv, refs.A_erfc, False)))
    add(Target(sp + 'ndtri', S, D, sc('ndtri'), ('unit',), gens=(PROB,), grid=[(v,) for v in G_UNIT],
               valid=lambda p: unit_open(p), inv=_inv_dict(refs.inv_ndtri, refs.A_ndtr, True)))

    def st_sampler(rng, n):
        df = DF(rng, n)
        import scipy.special as s
        q = LU(1e-300, 0.5)(rng, n)
        with np.errstate(all='ignore'):
            t = -np.asarray(s.stdtrit(df, q), dtype=float)
        t = np.where(rng.random(n) < 0.5, t, -t)
        alt = SLU(1e-10, 1e300)(rng, n)
        t = np.where(~np.isfinite(t) | (rng.random(n) < 0.25), alt, t)
        return np.stack([df, t], axis=1)
    add(Target(sp + 'stdtr', S, D, sc('stdtr'), ('pos', 'real'), sampler=st_sampler,
               grid=_prod(G_DF, G_REAL), valid=lambda df, t: pos(df) and fin(t),
               mref=refs.m_stdtr, aref=refs.A_stdtr))
    add(Target(sp + 'stdtrit', S, D, sc('stdtrit'), ('pos', 'unit'), gens=(DF, PROB),
               grid=_prod(G_DF, G_UNIT), valid=lambda df, p: pos(df) and unit_open(p),
               inv=_inv_dict(refs.inv_stdtrit, refs.A_stdtr, True)))
    FD = M((1.5, LU(1e-3, 1)), (4, LU(1, 1e3)), (3, LU(1e3, 1e8)), (1, LU(1e8, 1e15)))

    def fd_sampler(rng, n):
        import scipy.special as s
        d1, d2 = FD(rng, n), FD(rng, n)
        q = PROB(rng, n)
        with np.errstate(all='ignore'):
            x = np.asarray(s.fdtri(d1, d2, q), dtype=float)
        alt = LU(1e-300, 1e300)(rng, n)
        x = np.where(~np.isfinite(x) | (x <= 0) | (rng.random(n) < 0.25), alt, x)
        return np.stack([d1, d2, x], axis=1)
    gfd = [1e-3, 0.5, 1.0, 2.0, 5.0, 30.0, 1e3, 1e6, 1e10]
    add(Target(sp + 'fdtr', S, D, sc('fdtr'), ('pos', 'pos', 'pos'), sampler=fd_sampler,
               grid=_prod(gfd, gfd, [1e-300, 1e-50, 1e-10, 1e-3, 0.1, 0.5, 1.0, 2.0, 10.0, 1e3, 1e10, 1e50, 1e300]),
               valid=lambda a, b, x: pos(a, b) and fin(x) and x >= 0,
               mref=refs.m_fdtr, aref=refs.A_fdtr))

    def ch_sampler(rng, n):
        import scipy.special as s
        v = M((1.5, LU(1e-6, 1e-2)), (3.5, LU(1e-2, 20)), (3, LU(20, 1e5)), (2, LU(1e5, 1e12)))(rng, n)
        q = PROB(rng, n)
        with np.errstate(all='ignore'):
            x = np.asarray(s.chdtri(v, q), dtype=float)
        alt = LU(1e-300, 1e300)(rng, n)
        x = np.where(~np.isfinite(x) | (x <= 0) | (rng.random(n) < 0.25), alt, x)
        return np.stack([v, x], axis=1)
    add(Target(sp + 'chdtr', S, D, sc('chdtr'), ('pos', 'pos'), sampler=ch_sampler,
               grid=[(2 * a, 2 * x) for a, x in grid2], valid=lambda v, x: pos(v) and fin(x) and x >= 0,
               mref=refs.m_chdtr, aref=refs.A_chdtr))
    add(Target(sp + 'bdtr', S, D, sc('bdtr'), ('int', 'int', 'unit'), sampler=_bdtr_sampler,
               grid=[(k, n, p) for n in (1.0, 2.0, 10.0, 100.0, 1e4, 1e6, 1e9) for k in (0.0, 1.0, n // 2, n - 1)
                     for p in G_UNIT],
               valid=lambda k, n, p: fin(k) and fin(n) and 0 <= k < n and n == math.floor(n) and 0 <= p <= 1,
               mref=refs.m_bdtr, aref=refs.A_bdtr))
    HP = M((4, U(-20, 20)), (3, SLU(1e-3, 1e3)), (2, NEARINT(0, 60, -1.0, 1e-15, 1e-2)),
           (1, lambda rng, n: -np.floor(LU(1, 60)(rng, n))))

    def h2_valid(a, b, c, x):
        if not all(fin(v) for v in (a, b, c, x)):
            return False
        if x > 1 or is_nonpos_int(c):
            return False
        if x == 1 and not (c - a - b > 0):
            return False
        return True
    add(Target(sp + 'hyp2f1', S, D, sc('hyp2f1'), ('real', 'real', 'real', 'real'),
               gens=(HP, HP, HP, M((4, U(-1, 1)), (2, ONEM(1e-16, 0.5)), (2, lambda rng, n: -LU(1, 1e5)(rng, n)),
                                   (1, SLU(1e-10, 1e-2)))),
               grid=_prod([-10.0, -2.0, -0.5, 0.5, 1.0, 2.5, 10.0, 100.0], [-3.0, -0.5, 0.25, 1.5, 10.0, 100.0],
                          [-2.5, 0.5, 1.0, 3.0, 20.0, 150.0], [-100.0, -1.0, -0.5, 0.1, 0.5, 0.9, 0.999, 1 - 1e-10]),
               valid=h2_valid, mref=refs.m_hyp2f1, aref=refs.A_hyp2f1))
    H1 = M((4, U(-50, 50)), (3, SLU(1e-3, 1e3)), (2, NEARINT(0, 60, -1.0, 1e-15, 1e-2)))
    add(Target(sp + 'hyp1f1', S, D, sc('hyp1f1'), ('real', 'real', 'real'),
               gens=(H1, H1, M((5, U(-60, 60)), (5, SLU(1e-5, 1e3)))),
               grid=_prod([-50.5, -10.0, -1.5, -0.5, 0.5, 1.0, 3.0, 20.0, 200.0], [-9.5, -0.5, 0.5, 1.0, 2.0, 20.0, 200.0],
                          [-500.0, -50.0, -5.0, -0.1, 0.1, 5.0, 50.0, 500.0]),
               valid=lambda a, b, x: all(fin(v) for v in (a, b, x)) and not is_nonpos_int(b),
               mref=refs.m_hyp1f1, aref=refs.A_hyp1f1))
    BV = M((3, U(0, 5)), (3, LU(1e-3, 1e3)), (2, lambda rng, n: np.floor(U(0, 100)(rng, n))), (1, U(-50, 0)))
    BX = M((3, LU(1e-3, 1e3)), (5, LU(10, 1e4)), (2, LU(1e4, 1e300)))
    gb = _prod([0.0, 0.5, 1.0, 2.5, 10.0, 50.5, 100.0, 1e3], [0.1, 1.0, 10.0, 100.0, 700.0, 1e3, 1e4, 1e6, 1e9, 1e15, 1e100, 1e300])
    add(Target(sp + 'iv', S, D, sc('iv'), ('real', 'pos'), gens=(BV, BX), grid=gb,
               valid=lambda v, x: fin(v) and pos(x), mref=refs.m_iv, aref=refs.A_iv))
    add(Target(sp + 'kv', S, D, sc('kv'), ('real', 'pos'), gens=(BV, BX), grid=gb,
               valid=lambda v, x: fin(v) and pos(x), mref=refs.m_kv, aref=refs.A_kv))

    def jv_sampler(rng, n):
        v = M((4, BV), (2, LU(1e2, 1e4)))(rng, n)
        x = BX(rng, n)
        m = rng.random(n) < 0.3
        x = np.where(m, np.abs(v) * (1 + SLU(1e-6, 0.5)(rng, n)) + 1e-3, x)
        return np.stack([v, x], axis=1)
    add(Target(sp + 'jv', S, D, sc('jv'), ('real', 'pos'), sampler=jv_sampler, grid=gb,
               valid=lambda v, x: fin(v) and pos(x), mref=refs.m_jv, aref=refs.A_jv, osc='bessel'))
    ZX = M((3, U(-170, 60)), (2, NEAR(1.0, 1e-16, 0.5)), (2, NEARINT(1, 85, -2.0, 1e-16, 0.2)),
           (2, LU(1, 1e3)), (1, SLU(1e-300, 1e-3)))
    add(Target(sp + 'zeta', S, D, sc('zeta'), ('real',), gens=(ZX,),
               grid=[(v,) for v in [-170.5, -100.0, -50.5, -10.0, -2.0 - 1e-10, -1.0, -0.5, 1e-300, 0.5, 1 - 1e-10,
                                    1 + 1e-10, 1.5, 2.0, 3.0, 10.0, 54.0, 60.0, 1e3, 1e300]],
               valid=lambda x: fin(x) and x != 1, mref=refs.m_zeta, aref=refs.A_zeta, note='Riemann zeta'))
    add(Target(sp + 'zeta2', S, D, sc('zeta'), ('pos', 'pos'),
               gens=(M((3, lambda rng, n: 1.0 + LU(1e-16, 1)(rng, n)), (5, LU(1, 1e3))),
                     M((4, LU(1e-300, 1e300)), (4, LU(1e-3, 1e3)), (2, LU(1e-20, 1e-3)))),
               grid=_prod([1 + 1e-15, 1 + 1e-8, 1.5, 2.0, 3.0, 10.0, 100.0, 1e3], [1e-300, 1e-20, 1e-3, 0.5, 1.0, 10.0, 1e3, 1e10, 1e100, 1e300]),
               valid=lambda s, q: fin(s) and s > 1 and pos(q), mref=refs.m_zeta2, aref=refs.A_zeta2,
               note='Hurwitz zeta(x, q)'))


def _build_scipy_stats():
    S, D = 'scipy', 'float64'
    sx = 'scipy.stats.'
    shapes = {
        'norm': (),
        't': (DF,),
        'chi2': (M((1.5, LU(1e-6, 1e-2)), (3.5, LU(1e-2, 20)), (3, LU(20, 1e5)), (2, LU(1e5, 1e12))),),
        'f': (M((1.5, LU(1e-3, 1)), (4, LU(1, 1e3)), (3, LU(1e3, 1e8)), (1, LU(1e8, 1e15))),) * 2,
        'beta': (PAR, PAR),
        'gamma': (M((1.5, LU(1e-6, 1e-2)), (3.5, LU(1e-2, 10)), (3, LU(10, 1e5)), (2, LU(1e5, 1e12))),),
    }
    xgen = {'norm': M((5, SLU(1e-5, 1e300)), (5, U(-40, 40))), 't': SLU(1e-10, 1e300), 'chi2': LU(1e-300, 1e300),
            'f': LU(1e-300, 1e300), 'beta': UNIT, 'gamma': LU(1e-300, 1e300)}
    xkind = {'norm': 'real', 't': 'real', 'chi2': 'pos', 'f': 'pos', 'beta': 'unit', 'gamma': 'pos'}
    gridx = {'norm': [v for v in G_REAL], 't': G_REAL, 'chi2': [1e-300, 1e-20, 1e-3, 0.5, 1.0, 10.0, 100.0, 1e3, 1e5, 1e10, 1e100, 1e300],
             'f': [1e-300, 1e-20, 1e-3, 0.5, 1.0, 10.0, 100.0, 1e3, 1e5, 1e10, 1e100, 1e300], 'beta': G_UNIT,
             'gamma': [1e-300, 1e-20, 1e-3, 0.5, 1.0, 10.0, 100.0, 1e3, 1e5, 1e10, 1e100, 1e300]}
    gridshape = {'norm': [()], 't': [(v,) for v in G_DF], 'chi2': [(v,) for v in [1e-3, 0.5, 1.0, 2.0, 5.0, 30.0, 1e3, 1e6, 1e10]],
                 'f': _prod([0.5, 1.0, 5.0, 30.0, 1e3, 1e6], [0.5, 1.0, 5.0, 30.0, 1e3, 1e6]),
                 'beta': _prod([1e-3, 0.5, 1.0, 5.0, 100.0, 1e5], [1e-3, 0.5, 1.0, 5.0, 100.0, 1e5]),
                 'gamma': [(v,) for v in [1e-3, 0.5, 1.0, 2.0, 5.0, 30.0, 1e3, 1e6, 1e10]]}

    def valid_for(dist):
        def v(x, *s):
            if not all(fin(t) and t > 0 for t in s):
                return False
            if dist == 'beta':
                return 0 <= x <= 1
            if dist in ('chi2', 'f', 'gamma'):
                return fin(x) and x >= 0
            return fin(x)
        return v

    for dist in ('norm', 't', 'chi2', 'f', 'beta', 'gamma'):
        kinds = (xkind[dist],) + ('pos',) * len(shapes[dist])
        grid = [(x,) + s for s in gridshape[dist] for x in gridx[dist]]
        for meth in ('sf', 'logsf', 'logcdf'):
            add(Target(sx + dist + '.' + meth, S, D, st(dist, meth), kinds,
                       sampler=tail_x(dist, shapes[dist], xgen[dist]), grid=grid, valid=valid_for(dist),
                       mref=refs.m_dist(dist, meth), aref=refs.a_dist(dist, meth)))
        mf = refs.DISTS[dist][1]
        add(Target(sx + dist + '.isf', S, D, st(dist, 'isf'), ('unit',) + ('pos',) * len(shapes[dist]),
                   gens=(PROB,) + shapes[dist], grid=[(q,) + s for s in gridshape[dist] for q in G_UNIT],
                   valid=lambda q, *s: unit_open(q) and all(fin(t) and t > 0 for t in s),
                   inv=_inv_dict(refs.inv_isf(dist), [(lambda mf: lambda x, *s: mf(x, *s)[1])(mf)], False, yi=0),
                   note='isf(q, *shape)'))
    # binomial: (k, n, p)

    def bsampler(rng, n):
        k, N, p = _bdtr_sampler(rng, n).T
        return np.stack([k, N, p], axis=1)
    bgrid = [(k, n, p) for n in (1.0, 10.0, 100.0, 1e4, 1e6, 1e9) for k in (0.0, 1.0, n // 2, n - 1)
             for p in G_UNIT]
    bvalid = lambda k, n, p: fin(k) and fin(n) and 0 <= k < n and n == math.floor(n) and n >= 1 and 0 <= p <= 1
    for meth in ('sf', 'logsf', 'logcdf'):
        add(Target(sx + 'binom.' + meth, S, D, st('binom', meth), ('int', 'int', 'unit'), sampler=bsampler,
                   grid=bgrid, valid=bvalid, mref=refs.m_dist('binom', meth), aref=refs.a_dist('binom', meth)))

    def bisf_sampler(rng, n):
        k, N, p = _bdtr_sampler(rng, n).T
        q = PROB(rng, n)
        return np.stack([q, N, p], axis=1)
    add(Target(sx + 'binom.isf', S, D, st('binom', 'isf'), ('unit', 'int', 'unit'), sampler=bisf_sampler,
               grid=[(q, n, p) for n in (1.0, 10.0, 1e4, 1e9) for p in (1e-10, 0.1, 0.5, 0.9, 1 - 1e-10) for q in G_UNIT],
               valid=lambda q, n, p: unit_open(q) and fin(n) and n >= 1 and n == math.floor(n) and unit_open(p),
               out='discrete', note='discrete: smallest k with sf(k) <= q'))


def _build_torch(dt):
    T = 'torch'
    tp = 'torch.special.'
    sfx = '' if dt == 'float64' else '[f32]'
    LG = M((3, REAL_WIDE), (2, NEARINT(0, 170, -1.0, 1e-300, 0.5)), (1, NEAR(1.0)), (1, NEAR(2.0)),
           (2, NEAR_SET('lgamma_zeros', 1e-17, 1e-3)))
    add(Target(tp + 'gammaln' + sfx, T, dt, th('gammaln', dt), ('real',), gens=(LG,),
               grid=[(v,) for v in G_REAL], valid=lambda x: fin(x) and not is_nonpos_int(x),
               mref=refs.m_gammaln, aref=refs.A_gammaln))
    DG = M((3, REAL_WIDE), (2, NEARINT(0, 1e6, -1.0, 1e-300, 0.5)), (3, NEAR_SET('digamma_roots', 1e-17, 1e-3)))
    add(Target(tp + 'digamma' + sfx, T, dt, th('digamma', dt), ('real',), gens=(DG,),
               grid=[(v,) for v in G_REAL] + [(1.4616321449683623,)], valid=lambda x: fin(x) and not is_nonpos_int(x),
               mref=refs.m_digamma, aref=refs.A_digamma))
    PN = M((5, C(1.0, 2.0, 3.0)), (3, C(4.0, 5.0, 6.0, 8.0, 10.0)), (2, C(15.0, 20.0, 30.0)))
    add(Target(tp + 'polygamma' + sfx, T, dt, th('polygamma', dt), ('fixed', 'real'),
               gens=(PN, M((4, SLU(1e-3, 1e6)), (3, U(-30, 30)), (2, NEARINT(0, 1e4, -1.0, 1e-12, 0.5)), (1, SLU(1e-30, 1e300)))),
               grid=_prod([1.0, 2.0, 3.0, 5.0, 10.0], [-10.5, -2.5, -0.5, 1e-8, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1e4, 1e8]),
               valid=lambda n, x: fin(x) and not is_nonpos_int(x), mref=refs.m_polygamma, aref=refs.A_polygamma))
    EX = M((4, SLU(1e-300, 1e300)), (4, U(-30, 30)), (2, SLU(0.1, 1e4)))
    for nm, mr, ar in (('erfcx', refs.m_erfcx, refs.A_erfcx), ('log_ndtr', refs.m_log_ndtr, refs.A_log_ndtr),
                       ('i0e', refs.m_i0e, refs.A_i0e), ('i1e', refs.m_i1e, refs.A_i1e)):
        add(Target(tp + nm + sfx, T, dt, th(nm, dt), ('real',), gens=(EX,), grid=[(v,) for v in G_REAL],
                   valid=lambda x: fin(x), mref=mr, aref=ar))
    add(Target(tp + 'ndtri' + sfx, T, dt, th('ndtri', dt), ('unit',), gens=(PROB,), grid=[(v,) for v in G_UNIT],
               valid=lambda p: unit_open(p), inv=_inv_dict(refs.inv_ndtri, refs.A_ndtr, True)))
    add(Target(tp + 'zeta' + sfx, T, dt, th('zeta', dt), ('pos', 'pos'),
               gens=(M((3, lambda rng, n: 1.0 + LU(1e-7, 1)(rng, n)), (5, LU(1, 100))),
                     M((4, LU(1e-30, 1e30)), (4, LU(1e-3, 1e3)), (2, LU(1e-10, 1e-3)))),
               grid=_prod([1 + 1e-6, 1.5, 2.0, 3.0, 10.0, 100.0], [1e-30, 1e-3, 0.5, 1.0, 10.0, 1e3, 1e10, 1e30]),
               valid=lambda s, q: fin(s) and s > 1 and pos(q), mref=refs.m_zeta2, aref=refs.A_zeta2,
               note='Hurwitz zeta(x, q)'))
