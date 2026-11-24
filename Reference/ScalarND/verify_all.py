#!/usr/bin/env python3
"""Recompute both bounds of every stored certificate in results/ from its exact data (problem, candidate, field) and
compare with the stored claims (exact decimals, see scalarnd.claim_strings).
usage: python verify_all.py [--update]
--update rewrites every claim from the recomputed bounds; use it only after changing the rigorous code.
A certificate whose .npz file is absent (the large ones may be distributed separately, see README.md) is reported
as skipped; it does not count as a failure."""
import glob, json, math, os, sys, time
import scalarnd as S

update = '--update' in sys.argv
bad = 0
skipped = 0
rows = []
here = os.path.dirname(os.path.abspath(__file__))
for js in sorted(glob.glob(os.path.join(here, 'results', '*.json'))):
    stem = js[:-5]
    t0 = time.time()
    meta = json.load(open(js))
    if meta.get('format') != 'CertiField-ND-certificate-1':
        continue
    if not os.path.exists(stem + '.npz'):
        skipped += 1
        print(f"{os.path.basename(stem):15s} SKIPPED: {os.path.basename(stem)}.npz is not here   (claim {meta['lower']} <= inf E <= {meta['upper']})", flush=True)
        continue
    ok, lb, ub = S.verify_certificate(stem, verbose=False)
    if update:
        meta['lower'], meta['upper'] = S.claim_strings(lb, ub)
        json.dump(meta, open(js, 'w'), indent=1)
        ok = lb <= ub
        tag = 'UPDATED'
    else:
        tag = 'OK' if ok else 'FAILED'
    bad += (not ok)
    ref = (meta.get('info') or {}).get('reference_float')
    inside = '' if ref is None else ('   reference %+.8f %s' % (ref, 'inside' if lb <= ref <= ub else 'OUTSIDE'))
    if ref is not None and not (lb <= ref <= ub):
        bad += 1
    P = meta['problem']
    lo, hi = float(meta['lower']), float(meta['upper'])
    rows.append(f"{os.path.basename(stem):15s} N={'x'.join(map(str, P['N'])):9s} K={meta['K']:2d}  claim {meta['lower']} <= inf E <= {meta['upper']}"
                f"   gap {math.ceil(1e5 * (hi - lo) / abs(hi)) / 1e3:6.3f} %   {tag}  ({time.time() - t0:.0f} s){inside}")
    print(rows[-1], flush=True)
if skipped:
    print(f'{skipped} certificate(s) skipped (no .npz file)')
print('ALL CERTIFICATES VERIFIED' if bad == 0 else f'{bad} PROBLEM(S)')
sys.exit(1 if bad else 0)
