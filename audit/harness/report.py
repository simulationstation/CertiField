"""Build audit/summary.md (git-ignored) from the triage and confirmation results and the
local notes file audit/results/notes.json, which records the duplicate checks:

    {"<target>": [{"kind": "accuracy", "match": "<regex on cluster bucket>",
                   "group": "<short name>", "status": "new|duplicate|fixed-upstream|...",
                   "dup": "<issues / release notes checked>", "mechanism": "...",
                   "draft": "<drafts/... file or empty>"}]}

    python report.py
"""
import json
import math
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
AUDIT = os.path.dirname(HERE)
RESULTS = os.environ.get('AUDIT_RESULTS') or os.path.join(AUDIT, 'results')

ARGS = {
    'betainc': 'a,b,x', 'betaincc': 'a,b,x', 'betaincinv': 'a,b,y', 'gammainc': 'a,x', 'gammaincc': 'a,x',
    'gammaincinv': 'a,y', 'gammaln': 'x', 'loggamma': 'x', 'digamma': 'x', 'erf': 'x', 'erfc': 'x',
    'erfcx': 'x', 'ndtr': 'x', 'log_ndtr': 'x', 'erfcinv': 'y', 'ndtri': 'p', 'stdtr': 'df,t',
    'stdtrit': 'df,p', 'fdtr': 'dfn,dfd,x', 'chdtr': 'v,x', 'bdtr': 'k,n,p', 'hyp2f1': 'a,b,c,x',
    'hyp1f1': 'a,b,x', 'iv': 'v,x', 'kv': 'v,x', 'jv': 'v,x', 'zeta': 'x', 'zeta2': 'x,q',
    'polygamma': 'n,x', 'i0e': 'x', 'i1e': 'x',
}
STATS_ARGS = {'norm': '', 't': 'df', 'chi2': 'df', 'f': 'dfn,dfd', 'beta': 'a,b', 'gamma': 'a', 'binom': 'n,p'}


OTHER = {
    'R.pt': 'x,df', 'R.pf': 'x,df1,df2', 'R.pchisq': 'x,df', 'R.pbeta': 'x,a,b', 'R.pgamma': 'x,shape',
    'R.pbinom': 'k,n,p', 'R.qt': 'lp,df', 'R.qbeta': 'lp,a,b',
    'julia.SpecialFunctions.gamma_inc': 'a,x', 'julia.SpecialFunctions.gamma_inc_inv': 'a,p',
    'julia.SpecialFunctions.beta_inc': 'a,b,x', 'julia.SpecialFunctions.beta_inc_inv': 'a,b,p',
    'julia.SpecialFunctions.logabsgamma': 'x', 'julia.SpecialFunctions.digamma': 'x',
    'julia.SpecialFunctions.trigamma': 'x', 'julia.SpecialFunctions.erfcinv': 'y',
    'julia.SpecialFunctions.erfcx': 'x', 'julia.SpecialFunctions.zeta(s,q)': 's,q',
    'julia.SpecialFunctions.besselj': 'v,x', 'julia.SpecialFunctions.besselix': 'v,x',
    'julia.SpecialFunctions.besselkx': 'v,x',
    'julia.Distributions.logccdf(normal)': 'x', 'julia.Distributions.logccdf(t)': 'x,df',
    'julia.Distributions.logccdf(chisq)': 'x,k', 'julia.Distributions.logccdf(f)': 'x,d1,d2',
    'julia.Distributions.logccdf(beta)': 'x,a,b', 'julia.Distributions.logccdf(gamma)': 'x,a',
    'julia.Distributions.logccdf(binom)': 'k,n,p', 'julia.Distributions.invlogccdf(t)': 'lp,df',
    'julia.Distributions.invlogccdf(beta)': 'lp,a,b',
}


def argnames(target):
    for k, v in OTHER.items():
        if target.startswith(k) and (len(target) == len(k) or target[len(k)] in '([' or k.endswith(')')):
            return v.split(',')
    base = target.split('[')[0]
    lib, rest = base.split('.', 1)
    if rest.startswith('stats.'):
        dist, meth = rest.split('.')[1:3]
        first = {'isf': 'q', 'sf': 'x', 'logsf': 'x', 'logcdf': 'x'}[meth]
        if dist == 'binom' and meth != 'isf':
            first = 'k'
        return [first] + [s for s in STATS_ARGS[dist].split(',') if s]
    return ARGS[rest.split('.')[-1]].split(',')


def fname(name):
    return name.replace('[', '_').replace(']', '')


def load_jsonl(path):
    out = []
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            for line in f:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    pass
    return out


def fmt(v):
    if v is None:
        return '-'
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return str(v)
    if v == 0:
        return '0'
    if abs(v) >= 1e4 or abs(v) < 1e-3:
        return '%.3g' % v
    return '%.4g' % v


def rng(names, lo, hi):
    out = []
    for n, a, b in zip(names, lo, hi):
        out.append('%s=%s' % (n, fmt(a)) if a == b else '%s in [%s, %s]' % (n, fmt(a), fmt(b)))
    return '; '.join(out)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    import targets
    notes = {}
    np_ = os.path.join(RESULTS, 'notes.json')
    if os.path.exists(np_):
        with open(np_, encoding='utf-8') as f:
            notes = json.load(f)
    rows, quiet, other, pending = [], [], [], []
    for nm in targets.names():
        tj = os.path.join(RESULTS, 'triage', fname(nm) + '.json')
        if not os.path.exists(tj):
            pending.append(nm)
            continue
        with open(tj, encoding='utf-8') as f:
            tri = json.load(f)
        t = targets.get(nm)
        names = argnames(nm)
        conf = load_jsonl(os.path.join(RESULTS, 'confirm', fname(nm) + '.jsonl'))
        groups = {}
        for c in tri['clusters']:
            if c['kind'] not in ('special', 'accuracy'):
                other.append((nm, c['kind'], c['n']))
                continue
            note = None
            for n in notes.get(nm, []):
                if n.get('kind', c['kind']) == c['kind'] and re.search(n.get('match', ''), c['bucket']):
                    note = n
                    break
            key = note['group'] if note else 'unclassified %s' % c['kind']
            groups.setdefault(key, {'note': note, 'clusters': [], 'kind': c['kind']})['clusters'].append(c)
        if not groups:
            quiet.append((nm, tri['points'], tri.get('timeouts', 0), tri.get('noref', 0)))
        for key, g in groups.items():
            cl = g['clusters']
            lo = [min(c['xmin'][i] for c in cl) for i in range(len(names))]
            hi = [max(c['xmax'][i] for c in cl) for i in range(len(names))]
            rel = [c['max_relerr'] for c in cl if c['max_relerr'] is not None]
            flags = sorted({f for c in cl for f in c['flags']})
            if t.dtype == 'float32':
                err = 'max %s ulps' % max((c.get('max_ulps') or 0) for c in cl)
            else:
                err = 'max rel. err. %s' % fmt(max(rel)) if rel else ''
            if flags:
                err = (err + '; ' if err else '') + 'returns ' + '/'.join(flags)
            # confirmation status of the confirmed points that fall in these clusters
            bks = {c['bucket'] for c in cl}
            import triage
            st = [c.get('status') for c in conf if triage.bucket(t, c['x']) in bks]
            ver, mpo, unv = st.count('verified'), st.count('mp-only'), st.count('unverified') + st.count('timeout')
            if ver and not unv:
                ref = 'CERTIFIED: Arb ball/bracket + mpmath 50/100 digits (%d pts)' % ver
            elif mpo and not unv:
                ref = 'FLOAT: mpmath 50/100 digits agree, no Arb (%d pts)' % (mpo + ver)
            elif st:
                ref = 'reference unverified at %d of %d pts' % (unv, len(st))
            else:
                ref = 'DIAGNOSTIC: screening only'
            if g['kind'] == 'special':
                sev = 'wrong special value'
            elif t.dtype == 'float32':
                sev = 'reportable (> 64 ulps)'
            else:
                sev = 'severe (>= 1e-6)' if rel and max(rel) >= 1e-6 else 'reportable (>= 1e-10)'
            n = g['note'] or {}
            rows.append({'lib': t.lib, 'fn': nm.split('.', 1)[1], 'group': key, 'range': rng(names, lo, hi),
                         'err': err, 'ref': ref, 'sev': sev, 'dup': n.get('dup', 'not yet checked'),
                         'status': n.get('status', 'open'), 'mech': n.get('mechanism', ''), 'draft': n.get('draft', ''),
                         'npts': sum(c['n'] for c in cl)})
    order = {'new': 0, 'open': 1, 'fixed-upstream': 2, 'duplicate': 3}
    rows.sort(key=lambda r: (order.get(r['status'], 4), r['fn']))
    from importlib.metadata import version
    vers = ', '.join('%s %s' % (p, version(p)) for p in ('numpy', 'scipy', 'torch', 'mpmath', 'python-flint'))
    L = ['# Numerical accuracy audit: summary', '',
         'Local working file (git-ignored); checked by the agent only. Generated %s.' % time.strftime('%Y-%m-%d %H:%M'), '',
         'Tested (latest releases): Python %s; %s (CPU wheel).' % (sys.version.split()[0], vers),
         'Boost.Math: SciPy 1.18.1 bundles commit c114cf4a6a (2026-04-23); compared with develop ef1019b87a '
         '(2026-11-24) through audit/harness/boostcheck.py.',
         'References: mpmath at 50 and 100 digits (agreement to >= 20 digits) and an Arb ball (python-flint), or an '
         'Arb bracket of the root for inverse functions. Errors are relative to the exact value at the exact binary input.', '']
    L += ['Status counts: ' + ', '.join('%s %d' % (s, sum(1 for r in rows if r['status'] == s))
                                         for s in ('new', 'open', 'fixed-upstream', 'duplicate', 'not-a-defect')), '']
    L += ['| library | function | input range | maximum error | reference | severity | duplicate check | status |',
          '|---|---|---|---|---|---|---|---|']
    for r in rows:
        L.append('| %s | %s (%s; %d pts) | %s | %s | %s | %s | %s | %s |' % (
            r['lib'], r['fn'], r['group'], r['npts'], r['range'], r['err'], r['ref'], r['sev'], r['dup'],
            r['status'] + (' ([draft](%s))' % r['draft'] if r['draft'] else '')))
    L += ['', '## Mechanisms', '']
    for r in rows:
        if r['mech']:
            L.append('- %s (%s): %s' % (r['fn'], r['group'], r['mech']))
    L += ['', '## Not counted: conditioning-limited or underflow handling (points)', '']
    agg = {}
    for nm, k, n in other:
        agg[(nm, k)] = agg.get((nm, k), 0) + n
    for (nm, k), n in sorted(agg.items()):
        L.append('- %s: %s %d' % (nm, k, n))
    L += ['', '## No candidates (points screened, timeouts, no reference)', '']
    for nm, n, to, nr in quiet:
        L.append('- %s: %d, %d, %d' % (nm, n, to, nr))
    if pending:
        L += ['', '## Pending', ''] + ['- ' + p for p in pending]
    with open(os.path.join(AUDIT, 'summary.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')
    print('rows=%d new=%d open=%d pending=%d' % (len(rows), sum(r['status'] == 'new' for r in rows),
                                                 sum(r['status'] == 'open' for r in rows), len(pending)))


if __name__ == '__main__':
    main()
