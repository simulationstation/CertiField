"""Classify and cluster the screening candidates of finished targets; for SciPy functions
that call Boost.Math, evaluate each cluster's worst point with the Boost commit SciPy
bundles and with Boost develop (boostcheck.py), to see whether upstream has fixed it.

    python triage.py [--targets prefix,...]      -> audit/results/triage/<target>.json

Kinds: special (NaN/inf/0/finite-for-overflow/wrong k where a finite non-zero value is
correct), accuracy (relative error above the threshold and above 16 kappa u),
conditioning (within 16 kappa u), underflow (true value below the normal range: flushed
to zero or clamped to the smallest normal number).
"""
import argparse
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import runner  # noqa: E402

RESULTS = os.environ.get('AUDIT_RESULTS') or os.path.join(os.path.dirname(HERE), 'results')
OUT = os.path.join(RESULTS, 'triage')
TINY = {'float64': 2.2250738585072014e-308, 'float32': 1.1754943508222875e-38}
U = {'float64': 2.0 ** -53, 'float32': 2.0 ** -24}

# SciPy target -> (Boost call, argument map, value transform); None: not Boost-backed
NEG = lambda v: -v
LOG = lambda v: math.log(v) if v > 0 else (-math.inf if v == 0 else math.nan)
BOOST = {
    'scipy.special.betainc': ('ibeta', lambda a, b, x: (a, b, x), None),
    'scipy.special.betaincc': ('ibetac', lambda a, b, x: (a, b, x), None),
    'scipy.special.betaincinv': ('ibeta_inv', lambda a, b, y: (a, b, y), None),
    'scipy.special.stdtr': ('t_cdf', lambda df, t: (df, t), None),
    'scipy.special.stdtrit': ('t_ppf', lambda df, p: (df, p), None),
    'scipy.special.hyp1f1': ('hyp1f1', lambda a, b, x: (a, b, x), None),
    'scipy.special.fdtr': ('f_cdf', lambda dfn, dfd, x: (dfn, dfd, x), None),
    'scipy.stats.f.sf': ('f_sf', lambda x, dfn, dfd: (dfn, dfd, x), None),
    'scipy.stats.f.logsf': ('f_sf', lambda x, dfn, dfd: (dfn, dfd, x), LOG),
    'scipy.stats.f.logcdf': ('f_cdf', lambda x, dfn, dfd: (dfn, dfd, x), LOG),
    'scipy.stats.t.sf': ('t_cdf', lambda x, df: (df, -x), None),
    'scipy.stats.t.logsf': ('t_cdf', lambda x, df: (df, -x), LOG),
    'scipy.stats.t.logcdf': ('t_cdf', lambda x, df: (df, x), LOG),
    'scipy.stats.t.isf': ('t_ppf', lambda q, df: (df, q), NEG),
    'scipy.stats.beta.sf': ('ibetac', lambda x, a, b: (a, b, x), None),
    'scipy.stats.beta.logsf': ('ibetac', lambda x, a, b: (a, b, x), LOG),
    'scipy.stats.beta.logcdf': ('ibeta', lambda x, a, b: (a, b, x), LOG),
    'scipy.stats.beta.isf': ('ibetac_inv', lambda q, a, b: (a, b, q), None),
    'scipy.stats.binom.sf': ('binom_sf', lambda k, n, p: (k, n, p), None),
    'scipy.stats.binom.logsf': ('binom_sf', lambda k, n, p: (k, n, p), LOG),
    'scipy.stats.binom.logcdf': ('binom_cdf', lambda k, n, p: (k, n, p), LOG),
    'scipy.stats.binom.isf': ('binom_isf', lambda q, n, p: (q, n, p), None),
}


def fname(name):
    return name.replace('[', '_').replace(']', '')


def load(path):
    out = []
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            for line in f:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    pass
    return out


def rval(r):
    try:
        return float(r.get('r'))
    except (TypeError, ValueError):
        return None


def kind(t, r):
    f = r.get('flag')
    if f == 'lib_error':
        return 'special'  # the library call failed or did not return
    rv = rval(r)
    tiny = TINY[t.dtype]
    if f == 'zero' and rv is not None and abs(rv) < tiny:
        return 'underflow'
    if f is None and rv is not None and abs(rv) < tiny and abs(r.get('y', 0)) <= tiny * (1 + 1e-12):
        return 'underflow'
    k = r.get('kappa')
    if f == 'zero' and k is not None and 16 * k * U[t.dtype] >= 1:
        return 'conditioning'  # 0 returned at (the double nearest to) a zero of the function
    if f in ('nan', 'inf', 'zero', 'finite_for_overflow', 'nonzero_for_zero', 'wrong_k'):
        return 'special'
    if k is not None and r.get('relerr') is not None and r['relerr'] <= 16 * k * U[t.dtype]:
        return 'conditioning'
    return 'accuracy'


def bucket(t, x):
    key = []
    for kd, v in zip(t.kinds, x):
        if v == 0 or not math.isfinite(v):
            key.append(str(v))
        elif kd == 'unit' and v > 0.5:
            key.append('1-1e%d' % math.floor(math.log10(max(1 - v, 1e-300))))
        elif kd == 'fixed':
            key.append('%g' % v)
        else:
            key.append(('' if v > 0 else '-') + '1e%d' % math.floor(math.log10(abs(v))))
    return ','.join(key)


def boost_eval(calls):
    """calls: list of (func, args) -> list of (old, new) strings."""
    if not calls:
        return []
    d = os.path.join(runner.state_dir(), 'boost')
    lines = ['%s %s' % (f, ' '.join(repr(float(a)) for a in args)) for f, args in calls]
    res = {}
    for tag in ('old', 'new'):
        exe = os.path.join(d, 'check_%s.exe' % tag)
        if not os.path.exists(exe):
            return [(None, None)] * len(calls)
        p = subprocess.run([exe], input='\n'.join(lines) + '\n', capture_output=True, text=True, timeout=600)
        res[tag] = [l.split(' -> ', 1)[1] if ' -> ' in l else l for l in p.stdout.splitlines()]
    return list(zip(res['old'], res['new']))


def to_float(s):
    try:
        return float(s.split()[0])
    except (ValueError, AttributeError, IndexError):
        return None


def triage(t, recs):
    cands = [r for r in recs if r.get('cand')]
    cl = {}
    for r in cands:
        k = kind(t, r)
        cl.setdefault((k, bucket(t, r['x'])), []).append(r)
    rows = []
    for (k, b), rs in cl.items():
        fin = [r for r in rs if isinstance(r.get('relerr'), float) and math.isfinite(r['relerr'])]
        worst = max(rs, key=lambda r: ((r.get('score') or 0), (r.get('relerr') if r in fin else 0)))
        xs = [r['x'] for r in rs]
        rows.append({'kind': k, 'bucket': b, 'n': len(rs),
                     'xmin': [min(v[i] for v in xs) for i in range(t.nargs)],
                     'xmax': [max(v[i] for v in xs) for i in range(t.nargs)],
                     'max_relerr': max((r['relerr'] for r in fin), default=None),
                     'max_ulps': max((r.get('ulps') or 0 for r in rs), default=None) if t.dtype == 'float32' else None,
                     'flags': sorted({r.get('flag') for r in rs if r.get('flag')}),
                     'worst': {k2: worst.get(k2) for k2 in ('x', 'y', 'r', 'relerr', 'ulps', 'kappa', 'flag')}})
    if t.name in BOOST:
        f, amap, tr = BOOST[t.name]
        calls = [(f, amap(*row['worst']['x'])) for row in rows]
        for row, (o, n) in zip(rows, boost_eval(calls)):
            row['boost_old'], row['boost_new'] = o, n
            rv = rval(row['worst'])
            nv = to_float(n) if n else None
            if nv is not None and tr is not None:
                nv = tr(nv)
            if rv is not None and nv is not None and math.isfinite(nv) and rv != 0:
                row['boost_new_relerr'] = abs(nv - rv) / abs(rv)
            elif rv is not None and nv is not None:
                row['boost_new_relerr'] = 0.0 if (nv == rv) else math.inf
    order = {'special': 0, 'accuracy': 1, 'conditioning': 2, 'underflow': 3}
    rows.sort(key=lambda r: (order[r['kind']], -(r['max_relerr'] or 0)))
    return rows


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    ap = argparse.ArgumentParser()
    ap.add_argument('--targets', default='')
    ap.add_argument('--show', type=int, default=6)
    args = ap.parse_args()
    import targets
    pre = [p for p in args.targets.split(',') if p]
    os.makedirs(OUT, exist_ok=True)
    for nm in targets.names():
        if pre and not any(nm.startswith(p) for p in pre):
            continue
        base = os.path.join(RESULTS, fname(nm))
        if not os.path.exists(base + '.done'):
            continue
        t = targets.get(nm)
        recs = load(base + '.jsonl')
        rows = triage(t, recs)
        with open(os.path.join(OUT, fname(nm) + '.json'), 'w', encoding='utf-8') as f:
            json.dump({'target': nm, 'points': len(recs),
                       'timeouts': sum(1 for r in recs if r.get('flag') == 'timeout'),
                       'noref': sum(1 for r in recs if r.get('flag') == 'noref'), 'clusters': rows}, f, indent=1, default=float)
        cnt = {}
        for r in rows:
            cnt[r['kind']] = cnt.get(r['kind'], 0) + r['n']
        print('%s  points=%d  %s' % (nm, len(recs), cnt))
        for r in [r for r in rows if r['kind'] in ('special', 'accuracy')][:args.show]:
            w = r['worst']
            extra = ''
            if 'boost_new' in r:
                extra = '  boost-old=%s  boost-new=%s' % ((r['boost_old'] or '')[:40], (r['boost_new'] or '')[:40])
            print('   %-8s %-34s n=%-4d max_rel=%-9s x=%s y=%r r=%s kappa=%s%s' % (
                r['kind'], r['bucket'][:34], r['n'], '%.2e' % r['max_relerr'] if r['max_relerr'] is not None else None,
                w['x'], w['y'], w['r'], '%.3g' % w['kappa'] if w['kappa'] else None, extra))


if __name__ == '__main__':
    main()
