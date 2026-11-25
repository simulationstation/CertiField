"""Per-point work done in worker processes: reference value, error metrics, condition
number; and the confirmation of a candidate (mpmath at 50 and 100 digits, Arb ball)."""
import math
import sys

import numpy as np
import mpmath
from mpmath import mp, mpf

import refs
import targets

DBL_MAX = sys.float_info.max
DBL_MIN = sys.float_info.min
FLT_MAX = float(np.finfo(np.float32).max)
FLT_MIN = float(np.finfo(np.float32).tiny)
U64 = 2.0 ** -53
U32 = 2.0 ** -24
H = 2.0 ** -30  # relative perturbation for condition numbers


# ---------------------------------------------------------------------------------------
# error metrics
# ---------------------------------------------------------------------------------------


def _ord32(v):
    b = np.array([v], dtype=np.float32).view(np.int32)[0]
    b = int(b)
    return b if b >= 0 else -(b & 0x7FFFFFFF)


def metrics(dtype, y, r):
    """Compare the library value y (float) with the reference r (mpf).

    Returns dict(relerr, ulps, flag).  flag is None, or one of nan / inf / zero /
    finite_for_overflow / nonzero_for_zero."""
    out = {'relerr': None, 'ulps': None, 'flag': None}
    if r is None:
        out['flag'] = 'noref'
        return out
    big, small = (DBL_MAX, DBL_MIN) if dtype == 'float64' else (FLT_MAX, FLT_MIN)
    ar = abs(r)
    if math.isnan(y):
        if ar == 0 or ar <= big * 1.0000001:
            out['flag'] = 'nan'
            out['relerr'] = math.inf
            return out
        out['flag'] = 'nan'
        out['relerr'] = math.inf
        return out
    # expected value in the target format
    if dtype == 'float32':
        r32 = float(np.float32(float(r))) if ar < 2 * FLT_MAX else math.copysign(math.inf, float(r))
        rd = r32
    else:
        rd = float(r) if ar < 2 * DBL_MAX else math.copysign(math.inf, float(r))
    if math.isinf(rd):
        if y == rd:
            out['relerr'] = 0.0
            out['ulps'] = 0
        else:
            out['flag'] = 'finite_for_overflow' if math.isfinite(y) else 'inf'
            out['relerr'] = abs(mpf(y) - r) / ar if math.isfinite(y) else math.inf
        return out
    if math.isinf(y):
        out['flag'] = 'inf'
        out['relerr'] = math.inf
        return out
    if ar == 0:
        out['relerr'] = 0.0 if y == 0 else math.inf
        out['ulps'] = 0 if y == 0 else None
        if y != 0:
            out['flag'] = 'nonzero_for_zero'
        return out
    if y == 0 and rd != 0:
        out['flag'] = 'zero'
    out['relerr'] = float(abs(mpf(y) - r) / max(ar, mpf(small)))
    if dtype == 'float32':
        out['ulps'] = abs(_ord32(y) - _ord32(rd))
    return out


def is_candidate(dtype, m):
    if m['flag'] in ('nan', 'inf', 'zero', 'finite_for_overflow', 'nonzero_for_zero', 'wrong_k'):
        return True
    if dtype == 'float32':
        return m['ulps'] is not None and m['ulps'] > 64
    return m['relerr'] is not None and m['relerr'] >= 1e-10


def score(dtype, m, kappa):
    """Objective for the adversarial search: error in excess of what conditioning explains."""
    if m['flag'] in ('nan', 'inf', 'zero', 'finite_for_overflow', 'nonzero_for_zero', 'wrong_k'):
        return 1e300
    if m['relerr'] is None:
        return 0.0
    k = kappa if (kappa is not None and math.isfinite(kappa)) else 1.0
    if dtype == 'float32':
        return float(m['ulps'] or 0) / max(1.0, 32 * k)
    return m['relerr'] / max(16 * k * U64, U64)


# ---------------------------------------------------------------------------------------
# references
# ---------------------------------------------------------------------------------------


def _mref(t, args, dps):
    with mp.workdps(dps):
        try:
            v = t.mref(*[mpf(a) for a in args])
            return +v if v is not None else None
        except (mpmath.libmp.NoConvergence, ZeroDivisionError, ValueError, OverflowError):
            return None


def ref_screen(t, args, y):
    """Reference for screening: Arb (>= 64 bits) with an mpmath (30 digits) fallback."""
    if t.out == 'discrete':
        with mp.workdps(30):
            k, sfm1, sfk = refs.binom_isf_ref(*args)
        return mpf(k), 'mp30', (sfm1, sfk)
    if t.inv is not None:
        with mp.workdps(30):
            try:
                x = t.inv['solve'](*args, x0=y if math.isfinite(y) else None)
            except (mpmath.libmp.NoConvergence, ZeroDivisionError, ValueError, OverflowError):
                x = None
        return x, 'mp30-solve', None
    v, ball, i = refs.arb_eval(t.aref, args) if t.aref else (None, None, None)
    if v is not None:
        return v, 'arb%d' % i, None
    v = _mref(t, args, 30)
    return v, 'mp30', None


def kappa(t, args, r, pts=None):
    """Condition number sum_i |x_i d f/d x_i / f| (continuous arguments only)."""
    if r is None or r == 0 or not mpmath.isfinite(r):
        return None
    cont = [i for i, k in enumerate(t.kinds) if k not in ('int', 'fixed')]
    if t.out == 'discrete':
        return None
    if t.inv is None:
        tot = 0.0
        for i in cont:
            if args[i] == 0:
                continue
            v, _, _ = refs.arb_eval(t.aref, args, pert=(i, H)) if t.aref else (None, None, None)
            if v is None:
                with mp.workdps(30):
                    try:
                        a = [mpf(x) for x in args]
                        a[i] = a[i] * (1 + mpf(H))
                        v = t.mref(*a)
                    except Exception:
                        return None
            tot += float(abs((v - r) / (H * r)))
        return tot
    # inverse: x solves F(fwd args) = y, the forward args being args with y replaced by x
    yi = inv_yi(t, args)
    y = args[yi]
    fa = t.inv['fwd_a']
    base = list(args)
    base[yi] = r
    F0, _, _ = refs.arb_eval(fa, base)
    Fx, _, _ = refs.arb_eval(fa, base, pert=(yi, H))
    if F0 is None or Fx is None:
        return None
    dFx = (Fx - F0) / (H * r)  # dF/dx
    if dFx == 0:
        return None
    tot = float(abs(mpf(y) / (r * dFx)))
    for i in cont:
        if i == yi or args[i] == 0:
            continue
        Fi, _, _ = refs.arb_eval(fa, base, pert=(i, H))
        if Fi is None:
            continue
        dFi = (Fi - F0) / (H * mpf(args[i]))
        tot += float(abs(mpf(args[i]) * dFi / (r * dFx)))
    return tot


def inv_yi(t, args):
    """Index of the probability argument of an inverse function (isf: first; else last)."""
    return t.inv.get('yi', len(args) - 1)


def screen(name, args, y):
    t = targets.get(name)
    try:
        r, src, extra = ref_screen(t, args, y)
    except Exception as e:  # reference failure: recorded, never counted
        return {'flag': 'noref', 'err': repr(e)[:200]}
    if t.out == 'discrete':
        return _discrete_metrics(t, args, y, r, extra, src)
    m = metrics(t.dtype, y, r)
    rec = dict(m)
    rec['r'] = mpmath.nstr(r, 22) if r is not None else None
    rec['src'] = src
    rec['kappa'] = None
    near = (m['flag'] not in (None, 'noref')) or (m['relerr'] is not None and (
        (t.dtype == 'float64' and m['relerr'] >= 1e-12) or (t.dtype == 'float32' and (m['ulps'] or 0) >= 16)))
    if near and m['flag'] != 'noref':
        try:
            rec['kappa'] = kappa(t, args, r)
        except Exception:
            rec['kappa'] = None
    if t.osc == 'bessel' and r is not None and m['relerr'] is not None and m['flag'] is None:
        v, x = args
        if x > abs(v) and x > 1:
            amp = mpf(2 / (math.pi * x)) ** 0.5
            rec['abserr_env'] = float(abs(mpf(y) - r) / amp)
    rec['score'] = score(t.dtype, m, rec['kappa'])
    rec['cand'] = is_candidate(t.dtype, m)
    return rec


def _discrete_metrics(t, args, y, k, extra, src):
    sfm1, sfk = extra
    q = mpf(args[0])
    rec = {'r': mpmath.nstr(k, 22), 'src': src, 'kappa': None, 'ulps': None}
    if math.isnan(y):
        rec.update(flag='nan', relerr=math.inf)
    elif y == float(k):
        rec.update(flag=None, relerr=0.0)
    else:
        border = any(v is not None and v != 0 and abs(v / q - 1) < 1e-9 for v in (sfm1, sfk))
        rec.update(flag=None if border else 'wrong_k', relerr=float(abs(mpf(y) - k) / max(1, abs(k))))
        rec['border'] = border
    rec['score'] = 1e300 if rec['flag'] else 0.0
    rec['cand'] = rec['flag'] is not None
    return rec


# ---------------------------------------------------------------------------------------
# confirmation
# ---------------------------------------------------------------------------------------


class _null:
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _agree(a, b, digits=20):
    if a is None or b is None:
        return False
    if a == b:
        return True
    return abs(a - b) <= mpf(10) ** (-digits) * abs(b)


def confirm(name, args, y):
    t = targets.get(name)
    res = {'args': list(args), 'y': y}
    if t.out == 'discrete':
        ks = []
        for dps in (50, 100):
            with mp.workdps(dps):
                ks.append(refs.binom_isf_ref(*args))
        k = ks[1][0]
        res['r50'], res['r100'] = str(ks[0][0]), str(k)
        res['mp_agree'] = ks[0][0] == k
        # Arb check: sf(k-1) > q >= sf(k)
        q, n, p = args
        ok = True
        for kk, rel in ((k - 1, 'gt'), (k, 'le')):
            if kk < 0:
                continue
            v, ball, _ = refs.arb_eval(refs.a_dist('binom', 'sf'), (float(kk), n, p), min_bits=80)
            if ball is None:
                ok = None
                break
            if rel == 'gt' and not (ball > q):
                ok = False
            if rel == 'le' and not (ball <= q):
                ok = False
        res['arb_check'] = ok
        rec = _discrete_metrics(t, args, y, mpf(k), (ks[1][1], ks[1][2]), 'mp100')
        res.update(rec)
        res['status'] = 'verified' if (res['mp_agree'] and ok) else ('mp-only' if res['mp_agree'] and ok is None else 'unverified')
        return res
    if t.inv is not None:
        x0 = y if math.isfinite(y) else None
        with mp.workdps(50):
            r50 = t.inv['solve'](*args, x0=x0)
        with mp.workdps(100):
            r100 = t.inv['solve'](*args, x0=float(r50) if r50 is not None else x0)
        res['r50'] = mpmath.nstr(r50, 30) if r50 is not None else None
        res['r100'] = mpmath.nstr(r100, 40) if r100 is not None else None
        res['mp_agree'] = _agree(r50, r100)
        # Arb bracket: F(x(1-d)) and F(x(1+d)) on opposite sides of y
        ok = None
        if r100 is not None and r100 != 0:
            d = mpf(10) ** -25
            yi = inv_yi(t, args)
            lo_, hi_ = (r100 * (1 - d), r100 * (1 + d)) if r100 > 0 else (r100 * (1 + d), r100 * (1 - d))
            alo, ahi = list(args), list(args)
            alo[yi], ahi[yi] = lo_, hi_
            _, blo, _ = refs.arb_eval(t.inv['fwd_a'], alo, min_bits=100, precs=(256, 512, 1024, 2048, 4096))
            _, bhi, _ = refs.arb_eval(t.inv['fwd_a'], ahi, min_bits=100, precs=(256, 512, 1024, 2048, 4096))
            if blo is not None and bhi is not None:
                yy = refs.arb(args[yi])
                if t.inv['increasing']:
                    ok = bool(blo < yy) and bool(bhi > yy)
                else:
                    ok = bool(blo > yy) and bool(bhi < yy)
        res['arb_check'] = ok
        r = r100
    else:
        r50 = _mref(t, args, 50)
        r100 = _mref(t, args, 100)
        res['r50'] = mpmath.nstr(r50, 30) if r50 is not None else None
        res['r100'] = mpmath.nstr(r100, 40) if r100 is not None else None
        res['mp_agree'] = _agree(r50, r100)
        ok = None
        if t.aref:
            va, ball, i = refs.arb_eval(t.aref, args, min_bits=110, precs=(256, 512, 1024, 2048, 4096))
            if ball is not None and r100 is not None:
                ok = _agree(va, r100, 25)
                res['arb'] = mpmath.nstr(va, 30)
                with refs.flint.ctx.workprec(512) if hasattr(refs.flint.ctx, 'workprec') else _null():
                    res['arb_contains'] = bool(ball.overlaps(refs.mpf_to_arb(r100) * (1 + refs.arb(0, 1e-40))))
        res['arb_check'] = ok
        if t.malt is not None:
            try:
                with mp.workdps(50):
                    alt = t.malt(*[mpf(a) for a in args])
                res['malt_agree'] = _agree(alt, r100)
            except Exception:
                res['malt_agree'] = None
        r = r100
    if res['mp_agree'] and ok:
        res['status'] = 'verified'
    elif res['mp_agree'] and ok is None:
        res['status'] = 'mp-only'
    else:
        res['status'] = 'unverified'
    m = metrics(t.dtype, y, r)
    res.update(m)
    try:
        res['kappa'] = kappa(t, args, r)
    except Exception:
        res['kappa'] = None
    k = res['kappa']
    u = U64 if t.dtype == 'float64' else U32
    res['cond_limited'] = bool(m['flag'] is None and m['relerr'] is not None and k is not None and m['relerr'] <= 16 * k * u)
    res['cand'] = is_candidate(t.dtype, m)
    return res


def run(task):
    kind, name, args, y = task
    if kind == 'screen':
        return screen(name, args, y)
    if kind == 'confirm':
        return confirm(name, args, y)
    raise ValueError(kind)
