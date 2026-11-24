#!/usr/bin/env python3
"""Checks of the functions that docs/SCALAR_ND_2.md adds to ia.py: the logarithm, s ln s, the verified square root and
the lower bound for pi, each against exact rational arithmetic that shares no code with them.
usage: python check_ia.py        (a few seconds)"""
import math
import sys
from fractions import Fraction as Fr

import numpy as np

from ia import I, ln, xlnx, sqrt_up, pi_lower


def ln_bounds(q, n=60):
    """rational lower and upper bounds for ln q, q a positive rational: q = r 2^e with r in [2/3, 2), then the artanh
    series for ln r and for ln 2, each with its geometric remainder"""
    q = Fr(q)
    e = 0
    while q >= 2:
        q /= 2
        e += 1
    while q < Fr(2, 3):
        q *= 2
        e -= 1
    z = (q - 1) / (q + 1)
    s = 2 * sum(z ** (2 * j + 1) / (2 * j + 1) for j in range(n + 1))
    rem = 2 * abs(z) ** (2 * n + 3) / ((2 * n + 3) * (1 - z * z))
    z2 = Fr(1, 3)
    s2 = 2 * sum(z2 ** (2 * j + 1) / (2 * j + 1) for j in range(n + 1))
    rem2 = 2 * z2 ** (2 * n + 3) / ((2 * n + 3) * (1 - z2 * z2))
    return (s - rem + e * (s2 if e >= 0 else s2 + rem2)), (s + rem + e * (s2 + rem2 if e >= 0 else s2))


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    rng = np.random.default_rng(1)
    ok = True

    xs = np.concatenate([rng.random(300), 10.0 ** rng.uniform(-300, 300, 300),
                         [1.0, 0.5, 2.0, 5e-324, 1.7976931348623157e308, 0.7071067811865476, 0.7071067811865475]])
    L = ln(I(xs))
    good, worst = True, 0.0
    for x, lo, hi in zip(xs, L.lo, L.hi):
        a, b = ln_bounds(Fr(float(x)))
        good &= Fr(float(lo)) <= a and b <= Fr(float(hi))
        worst = max(worst, (hi - lo) / max(abs(math.log(x)), 1e-3))
    ok &= good and worst < 1e-13
    print(f"ln: the exact value is enclosed at {len(xs)} points: {good};  largest width relative to max(|ln x|, 1e-3): {worst:.1e}")

    s = np.concatenate([rng.random(200), 10.0 ** rng.uniform(-200, 2, 100), [1.0, 0.25, 0.3]])
    X = xlnx(I(s))
    good = True
    for t, lo, hi in zip(s, X.lo, X.hi):
        a, b = ln_bounds(Fr(float(t)))
        v = sorted([Fr(float(t)) * a, Fr(float(t)) * b])
        good &= Fr(float(lo)) <= v[0] and v[1] <= Fr(float(hi))
    Z = xlnx(I(np.zeros(4), np.array([0.0, 0.2, 0.3, 2.0])))
    # on [0, b]: the range of s ln s is [b ln b, 0] for b <= 1/e, and contains -1/e and max(0, b ln b) otherwise
    good &= (Z.lo[0] == 0.0 == Z.hi[0]) and Z.lo[1] <= 0.2 * math.log(0.2) and Z.hi[1] >= 0.0
    good &= Z.lo[2] <= -0.3679 and Z.hi[2] >= 0.0 and Z.lo[3] <= -0.3679 and Z.hi[3] >= 2.0 * math.log(2.0)
    ok &= good
    print(f"s ln s: enclosed at {len(s)} points and on four intervals [0, b]: {good}")

    good = True
    for _ in range(3000):
        x = float(10.0 ** rng.uniform(-200, 200))
        r = sqrt_up(x)
        good &= Fr(r) ** 2 >= Fr(x) and r <= math.sqrt(x) * (1 + 1e-15)
    good &= sqrt_up(0.0) == 0.0 and sqrt_up(4.0) >= 2.0
    ok &= good
    print(f"sqrt_up: s^2 >= x exactly and s within 1e-15 of the square root, 3000 points: {good}")

    p = pi_lower()
    # pi = 3.14159265358979323846264338327950288...  (50 digits are in every table; the test uses 35)
    good = Fr(31415926535897932384626433832795028, 10 ** 34) < p < Fr(31415926535897932384626433832795029, 10 ** 34)
    ok &= good
    print(f"pi_lower: {float(p)!r}, between 3.1415926535897932384626433832795028 and ...029: {good}")

    print('ALL CHECKS PASSED' if ok else 'A CHECK FAILED')
    sys.exit(0 if ok else 1)
