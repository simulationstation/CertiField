#!/usr/bin/env python3
"""results table (Markdown) from the stored certificates; bounds rounded in the safe direction.
Run it from this directory; it writes results/table.md (UTF-8)."""
import json, math, sys
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')            # the table has non-ASCII characters (Windows consoles)

order = [('screened2d', 'u² − u', 'u = 0 on ∂(0,1)²'), ('screened2d_128', 'u² − u', 'the same'),
         ('screened3d', 'u² − u', 'u = 0 on ∂(0,1)³'), ('screened3d_48', 'u² − u', 'the same'),
         ('quartic2d', 'u⁴/4 + u²/2 − u', 'u = 0 on ∂(0,1)²'), ('quartic3d', 'u⁴/4 + u²/2 − u', 'u = 0 on ∂(0,1)³'),
         ('wall2d', '(u² − 1)²/4', '(−3,3)×(0,1): u = −1 and +1 on the ends, sides free'),
         ('wall3d', '(u² − 1)²/4', '(−3,3)×(0,1)²: the same'),
         ('ramp2d', '4(u² − 1)²', 'u = x on ∂(−1,1)²'), ('saddle2d', '4(u² − 1)²', 'u = xy on ∂(−1,1)²'),
         ('ramp3d_32', '(u² − 1)²', 'u = x on ∂(−1,1)³')]


def dn(x, q='0.000001'):                                # x: a stored claim (exact decimal string)
    return str(Decimal(str(x)).quantize(Decimal(q), rounding=ROUND_FLOOR))


def up(x, q='0.000001'):
    return str(Decimal(str(x)).quantize(Decimal(q), rounding=ROUND_CEILING))


lines = ['| Potential W | Domain and data | Mode | Grid | K | Certified: lower ≤ inf E ≤ upper | Gap | Time | Reference value |',
         '|---|---|---|---|---|---|---|---|---|']
for name, W, data in order:
    m = json.load(open(f'results/{name}.json'))
    lo, hi = m['lower'], m['upper']
    gap = (float(hi) - float(lo)) / abs(float(hi))
    info = m.get('info') or {}
    sec = info['levels'][-1]['seconds']
    ref = info.get('reference_float')
    q = '0.00000001' if abs(float(hi)) < 0.1 else '0.000001'
    t = f'{sec:.0f} s' if sec < 120 else f'{sec / 60:.0f} min'
    if sec < 1.5:
        t = '1 s'
    gp = math.ceil(gap * 1e5) / 1e3
    lines.append(f"| {W} | {data} | {'convex' if m['K'] == 0 else 'calibration'} | {'×'.join(map(str, m['problem']['N']))} | "
                 f"{m['K'] if m['K'] else '–'} | {dn(lo, q)} ≤ inf E ≤ {up(hi, q)} | {gp:.3f} % | {t} | "
                 f"{('%.8f' % ref) if ref is not None else 'none known'} |")
open('results/table.md', 'w', encoding='utf-8', newline='\n').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
