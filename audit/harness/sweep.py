"""Screening sweep: grid + random sampling + adversarial search, per target.

    python sweep.py --targets scipy.special.betainc,torch.special  [--n 3000] [--workers 12]

A target is selected when its name starts with one of the given prefixes.  Every
evaluated point is written to audit/results/<target>.jsonl (git-ignored) with the
library value, the screening reference, the relative error (and float32 ULPs), the
condition number for points near the threshold, and flags.
"""
import argparse
import json
import math
import os
import sys
import time
import warnings
import zlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), 'results')


def cast(t, pts):
    """Make the points exactly representable in the target dtype; drop invalid ones."""
    if t.dtype == 'float32':
        with np.errstate(all='ignore'):
            p32 = pts.astype(np.float32)
        ok = np.all(np.isfinite(p32), axis=1) & np.all((p32 != 0) | (pts == 0), axis=1)
        pts = p32[ok].astype(np.float64)
    keep = np.array([bool(t.valid(*p)) for p in pts], dtype=bool) if len(pts) else np.zeros(0, bool)
    return pts[keep]


def call(t, pts):
    out = np.empty(len(pts))
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('ignore')
        for i in range(0, len(pts), 4096):
            c = pts[i:i + 4096]
            try:
                out[i:i + 4096] = t.call(*[c[:, j] for j in range(t.nargs)])
            except Exception:
                for k in range(len(c)):
                    try:
                        out[i + k] = t.call(*[c[k:k + 1, j] for j in range(t.nargs)])[0]
                    except Exception:
                        out[i + k] = np.nan
    return out


def perturb(rng, t, x, sigma):
    x = list(x)
    idx = [i for i, k in enumerate(t.kinds) if k != 'fixed']
    which = idx if rng.random() < 0.5 else [idx[rng.integers(len(idx))]]
    for i in which:
        k, v = t.kinds[i], x[i]
        g = rng.normal() * sigma
        if k == 'pos':
            v = v * math.exp(g) if v > 0 else 1e-300
        elif k == 'real':
            if v == 0:
                v = rng.normal() * 1e-300
            elif rng.random() < 0.05:
                v = v + g * (abs(v) + 1)
            else:
                v = v * math.exp(g)
        elif k == 'unit':
            v = v * math.exp(g) if v < 0.5 else 1 - (1 - v) * math.exp(g)
        elif k == 'p2':
            v = v * math.exp(g) if v < 1 else 2 - (2 - v) * math.exp(g)
        elif k == 'int':
            if rng.random() < 0.3:
                v = v + rng.choice([-1.0, 1.0])
            else:
                v = float(round(v * math.exp(g)))
        x[i] = v
    return tuple(x)


def run_target(pool, t, args):
    rng = np.random.default_rng(args.seed * 1000003 + zlib.crc32(t.name.encode()))
    t0 = time.time()
    grid = np.array(t.grid, dtype=float).reshape(-1, t.nargs) if t.grid else np.zeros((0, t.nargs))
    with np.errstate(all='ignore'):
        rnd = t.sample(rng, args.n)
    grid, rnd = cast(t, grid), cast(t, rnd)
    pts = np.concatenate([grid, rnd])
    phase = ['grid'] * len(grid) + ['random'] * len(rnd)
    y = call(t, pts)
    tasks_ = [('screen', t.name, tuple(float(v) for v in p), float(yy)) for p, yy in zip(pts, y)]
    res = pool.map(tasks_, progress=lambda d, n: print('   %s %d/%d' % (t.name, d, n), flush=True))
    recs = []
    for p, yy, r, ph in zip(pts, y, res, phase):
        r = dict(r or {})
        r.update(x=[float(v) for v in p], y=float(yy), phase=ph)
        recs.append(r)
    # adversarial search from the worst points (by excess score) and some random points
    if args.adv_steps > 0 and recs:
        order = sorted(range(len(recs)), key=lambda i: -(recs[i].get('score') or 0))
        starts, seen = [], set()
        for i in order:
            key = tuple(np.round(np.log10(np.abs(np.array(recs[i]['x'])) + 1e-300), 1))
            if key in seen or recs[i].get('flag') in ('noref', 'timeout'):
                continue
            seen.add(key)
            starts.append(i)
            if len(starts) >= args.adv_starts * 2 // 3:
                break
        extra = rng.choice(len(recs), min(len(recs), args.adv_starts - len(starts)), replace=False)
        starts += [int(i) for i in extra if recs[int(i)].get('flag') not in ('noref', 'timeout')]
        cur = [(tuple(recs[i]['x']), recs[i].get('score') or 0.0) for i in starts]
        sig = [0.3] * len(cur)
        floor = 2.0 ** -22 if t.dtype == 'float32' else 1e-14
        for step in range(args.adv_steps):
            props = []
            for (x, s), sg in zip(cur, sig):
                for _ in range(30):
                    q = perturb(rng, t, x, sg)
                    qa = cast(t, np.array([q], dtype=float))
                    if len(qa):
                        props.append(tuple(float(v) for v in qa[0]))
                        break
                else:
                    props.append(x)
            pa = np.array(props, dtype=float)
            ya = call(t, pa)
            rr = pool.map([('screen', t.name, p, float(v)) for p, v in zip(props, ya)])
            for j, (p, v, r) in enumerate(zip(props, ya, rr)):
                r = dict(r or {})
                r.update(x=list(p), y=float(v), phase='adv')
                recs.append(r)
                sc = r.get('score') or 0.0
                if sc > cur[j][1]:
                    cur[j] = (p, sc)
                    sig[j] = min(sig[j] * 1.6, 3.0)
                else:
                    sig[j] = max(sig[j] * 0.8, floor)
                    if sig[j] <= floor * 1.0001:
                        sig[j] = 0.3  # restart the step size
    os.makedirs(RESULTS, exist_ok=True)
    fn = os.path.join(RESULTS, t.name.replace('[', '_').replace(']', '') + '.jsonl')
    with open(fn, 'w', encoding='utf-8') as f:
        for r in recs:
            f.write(json.dumps(r, default=float) + '\n')
    summarize(t, recs, time.time() - t0)


def summarize(t, recs, dt):
    n = len(recs)
    flags = {}
    for r in recs:
        if r.get('flag'):
            flags[r['flag']] = flags.get(r['flag'], 0) + 1
    good = [r for r in recs if r.get('flag') is None and r.get('relerr') is not None]
    mx = max(good, key=lambda r: r['relerr'], default=None)
    mxs = max(good, key=lambda r: r.get('score') or 0, default=None)
    ncand = sum(1 for r in recs if r.get('cand'))
    line = '%-40s n=%6d cand=%5d  max_rel=%.2e at %s  max_excess=%.2e at %s  flags=%s  %.0fs' % (
        t.name, n, ncand, mx['relerr'] if mx else float('nan'), mx['x'] if mx else None,
        (mxs.get('score') or 0) if mxs else float('nan'), mxs['x'] if mxs else None, flags, dt)
    if t.dtype == 'float32':
        mu = max(good, key=lambda r: r.get('ulps') or 0, default=None)
        line += '  max_ulps=%s at %s' % (mu.get('ulps') if mu else None, mu['x'] if mu else None)
    print(line, flush=True)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    ap = argparse.ArgumentParser()
    ap.add_argument('--targets', default='')
    ap.add_argument('--exclude', default='')
    ap.add_argument('--n', type=int, default=3000)
    ap.add_argument('--adv-starts', type=int, default=24)
    ap.add_argument('--adv-steps', type=int, default=60)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--timeout', type=float, default=30)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--list', action='store_true')
    args = ap.parse_args()
    import targets
    names = targets.names()
    pre = [p for p in args.targets.split(',') if p]
    exc = [p for p in args.exclude.split(',') if p]
    sel = [n for n in names if (not pre or any(n.startswith(p) for p in pre)) and not any(n.startswith(p) for p in exc)]
    if args.list:
        print('\n'.join(sel))
        return
    from pool import Pool
    pool = Pool(args.workers, args.timeout)
    try:
        for nm in sel:
            run_target(pool, targets.get(nm), args)
    finally:
        pool.close()


if __name__ == '__main__':
    main()
