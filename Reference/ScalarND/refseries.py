#!/usr/bin/env python3
"""FLOAT reference: exact minimum of  int_{(0,1)^n} k |grad u|^2 + c u^2 - f u,  u = 0 on the boundary  (sine series)."""
import numpy as np

def screened_min(n, k, c, f, M=4001):
    odd = np.arange(1, M + 1, 2, dtype=float)
    if n == 2:
        m, p = np.meshgrid(odd, odd, indexing='ij')
        lam = 2 * k * np.pi ** 2 * (m * m + p * p) + 2 * c
        return -32 * f * f / np.pi ** 4 * np.sum(1.0 / (m * m * p * p * lam))
    if n == 3:
        odd = np.arange(1, 401, 2, dtype=float)
        m, p, q = np.meshgrid(odd, odd, odd, indexing='ij')
        lam = 2 * k * np.pi ** 2 * (m * m + p * p + q * q) + 2 * c
        return -256 * f * f / np.pi ** 6 * np.sum(1.0 / (m * m * p * p * q * q * lam))

if __name__ == '__main__':
    print('2D k=1/2 c=1 f=1:', screened_min(2, 0.5, 1.0, 1.0))
    print('3D k=1/2 c=1 f=1:', screened_min(3, 0.5, 1.0, 1.0))
    print('2D k=1/2 c=0 f=2 (torsion: -J/2, J = 0.1406):', screened_min(2, 0.5, 0.0, 2.0))
