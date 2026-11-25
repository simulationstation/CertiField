"""Certificate files (JSON + NPZ, see benchpack/certificates/FORMAT.md), claims and their weakening."""
import json
import math
import os
from fractions import Fraction as Fr

import numpy as np
from flint import arb

from exact import arb_bounds

CERT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificates')
LABEL = 'CERTIFIED (checked by the agent only)'


def fstr(x):
    return str(Fr(x))


def dec_round(x, ex, up):
    """x (Fraction) rounded down/up to a multiple of 10^ex, as an exact decimal string"""
    qn = Fr(10) ** ex
    m = math.ceil(x / qn) if up else math.floor(x / qn)
    if ex >= 0:
        return str(m * 10 ** ex)
    sgn = '-' if m < 0 else ''
    digits = str(abs(m)).rjust(-ex + 1, '0')
    return sgn + digits[:ex] + '.' + digits[ex:]


def claim_strings(lo, hi, digits=10):
    """weaken [lo, hi] (Fractions) by w = max(1e-9 s, 1e-3 (hi - lo)), s = max(|lo|, |hi|), then round outward to
    `digits` significant digits (FORMAT.md section 6)"""
    lo, hi = Fr(lo), Fr(hi)
    assert lo <= hi
    s = max(abs(lo), abs(hi))
    w = max(s / 10 ** 9, (hi - lo) / 1000)
    a, b = lo - w, hi + w
    ex = math.floor(math.log10(float(max(abs(a), abs(b))))) - digits + 1
    return dec_round(a, ex, False), dec_round(b, ex, True)


def ball_lo(b):
    return arb_bounds(b)[0]


def ball_hi(b):
    return arb_bounds(b)[1]


def encl(b):
    """float enclosure [lo, hi] of a ball (outward)"""
    lo, hi = arb_bounds(b)
    a, c = float(lo), float(hi)
    while Fr(a) > lo:
        a = math.nextafter(a, -math.inf)
    while Fr(c) < hi:
        c = math.nextafter(c, math.inf)
    return [a, c]


def save(stem, meta, arrays):
    os.makedirs(CERT_DIR, exist_ok=True)
    with open(os.path.join(CERT_DIR, stem + '.json'), 'w') as f:
        json.dump(meta, f, indent=1)
    np.savez_compressed(os.path.join(CERT_DIR, stem + '.npz'), **arrays)


def load(stem):
    with open(os.path.join(CERT_DIR, stem + '.json')) as f:
        meta = json.load(f)
    z = np.load(os.path.join(CERT_DIR, stem + '.npz'))
    return meta, {k: np.array(z[k], dtype=float) for k in z.files}


def kappa_dyadic(k, bits=4):
    """a dyadic rational close to k > 0 with `bits` significant bits"""
    e = math.floor(math.log2(k)) - bits + 1
    m = max(1, round(k / 2.0 ** e))
    return Fr(m) * Fr(2) ** e
