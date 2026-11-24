#!/usr/bin/env python3
"""FLOAT reference values (not rigorous): exact minimum of the 1D problem  int_{-L}^{L} k u'^2 + W(u),  u(-L) = -1,
u(L) = 1, for the double well W = (u^2-1)^2/4, from the first integral  k u'^2 = W(u) + c."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

def wall_energy(L, k=0.5, lam=1.0):
    W = lambda u: lam * (u * u - 1.0) ** 2 / 4.0
    def length(c):
        return quad(lambda u: np.sqrt(k / (W(u) + c)), -1, 1, epsabs=1e-13, epsrel=1e-13, limit=400)[0] - 2 * L
    c = brentq(length, 1e-16, 10.0, xtol=1e-18, rtol=1e-14)
    E = quad(lambda u: 2 * np.sqrt(k * (W(u) + c)), -1, 1, epsabs=1e-13, epsrel=1e-13, limit=400)[0] - 2 * L * c
    return E, c

if __name__ == '__main__':
    for L in (2, 3, 4, 6):
        print(L, wall_energy(L))
    print('infinite line', 2 * np.sqrt(2) / 3)
