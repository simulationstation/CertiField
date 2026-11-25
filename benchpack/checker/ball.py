"""Closed-form integral of a polynomial over a box intersected with a ball (P3 checker, mode "exact").

For a box B = [x0,x1] x [y0,y1] x [z0,z1] in the first octant, the ball |x| <= R (R^2 = R2 rational, R rational) and a
polynomial f (global coordinates), with Q(x,y,z) = int_{z0}^z f dz',  r_i^2 = R2 - z_i^2  and  Z = sqrt(R2 - x^2 - y^2):

    int_{B cap ball} f = int_{rect cap D(r1)} Q(x,y,z1) + int_{rect cap (D(r0) minus D(r1))} Q(x,y,Z) dx dy

(for (x,y) in the rectangle the vertical segment meets the ball in z0 <= z < min(z1, Z)).  Q(x,y,Z) = P(x,y) +
W(x,y) Z.  The polynomial parts are integrated by region.integrate_region (exact closed forms over rectangles cut by
circles).  For the part with Z the y-integral is exact,

    int y^n sqrt(c^2 - y^2) dy = P_n(y, c^2) sqrt(c^2 - y^2) + kappa_n c^(n+2) asin(y/c),     c^2 = R2 - x^2,

and the x-integral of what remains is elementary: polynomials times sqrt(a^2 - x^2) (region.sqrt_antider), and
terms T(x) asin(g(x)) integrated by parts, which leave  int A(x) x / ((R2 - x^2) sqrt(a^2 - x^2)) dx  with
a^2 + b^2 = R2:  polynomial division by R2 - x^2, then

    int x^m / sqrt(a^2 - x^2) dx = U_m(x) sqrt(a^2 - x^2) + nu_m asin(x/a),
    int x / ((R2 - x^2) sqrt(a^2 - x^2)) dx = -(1/b) asin( sqrt(a^2 - x^2) / sqrt(R2 - x^2) ),
    int 1 / ((R2 - x^2) sqrt(a^2 - x^2)) dx = (1/(R b)) asin( b x / (a sqrt(R2 - x^2)) ).

Every closed form is verified symbolically (an exact polynomial identity) when it is used.  Only the values at the
break points (rationals or square roots of rationals) are evaluated, with mpmath.iv (outward rounding).
"""
from fractions import Fraction as Fr
from functools import lru_cache, cmp_to_key

import numpy as np
from mpmath import iv

import exact as X
import region as RG
from region import Num, num_cmp, ivq


# ------------------------------------------------------------------------------------------------- closed forms
@lru_cache(maxsize=None)
def Pn(n):
    """P_n(y, u) as {j: coefficients in u} (u = c^2) and kappa_n:
    int y^n sqrt(u - y^2) dy = P_n sqrt(u - y^2) + kappa_n u^((n+2)/2) asin(y/sqrt(u))"""
    if n == 0:
        return {1: (Fr(1, 2),)}, Fr(1, 2)
    if n == 1:
        return {0: (Fr(0), Fr(-1, 3)), 2: (Fr(1, 3),)}, Fr(0)
    prev, kprev = Pn(n - 2)
    out = {}

    def add(j, poly):
        out[j] = tuple(X.padd(list(out.get(j, ())), list(poly)))
    add(n - 1, (Fr(0), Fr(-1, n + 2)))                 # -y^(n-1) u/(n+2)
    add(n + 1, (Fr(1, n + 2),))                        # +y^(n+1)/(n+2)
    f = Fr(n - 1, n + 2)
    for j, poly in prev.items():
        add(j, tuple(X.pscale([Fr(0)] + list(poly), f)))
    return out, f * kprev


def check_Pn(n):
    """exact identity d/dy [P sqrt(u-y^2) + kappa u^((n+2)/2) asin(y/sqrt u)] = y^n sqrt(u - y^2), as polynomials:
    P_y (u - y^2) - y P + kappa u^((n+2)/2) = y^n (u - y^2)"""
    P, kap = Pn(n)
    lhs = {}

    def add(j, poly):
        lhs[j] = X.padd(lhs.get(j, []), list(poly))
    for j, pu in P.items():
        if j >= 1:
            add(j - 1, X.pscale(X.pmul(list(pu), [Fr(0), Fr(1)]), j))      # j y^(j-1) p(u) u
            add(j + 1, X.pscale(list(pu), -j))                              # -j y^(j+1) p(u)
        add(j + 1, X.pscale(list(pu), -1))                                  # -y * y^j p(u)
    if kap:
        assert n % 2 == 0
        add(0, [Fr(0)] * ((n + 2) // 2) + [kap])
    rhs = {n: [Fr(0), Fr(1)], n + 2: [Fr(-1)]}
    keys = set(lhs) | set(rhs)
    return all(X.ptrim(lhs.get(j, [])) == X.ptrim(rhs.get(j, [])) for j in keys)


def inv_sqrt_antider(q, a2):
    """U, nu with d/dx [U sqrt(a2 - x^2) + nu asin(x/a)] = q / sqrt(a2 - x^2), verified exactly"""
    U, nu = [], Fr(0)
    Us = {0: ([], Fr(1)), 1: ([Fr(-1)], Fr(0))}
    for m in range(2, len(q)):
        Um2, num2 = Us[m - 2]
        Us[m] = (X.padd([Fr(0)] * (m - 1) + [Fr(-1, m)], X.pscale(Um2, (m - 1) * a2 / m)), (m - 1) * a2 / m * num2)
    for m, c in enumerate(q):
        if c:
            Um, num = Us[m]
            U = X.padd(U, X.pscale(Um, c))
            nu += c * num
    # check: U'(a2 - x^2) - x U + nu == q
    lhs = X.padd(X.padd(X.pmul(X.pder(U), [a2, Fr(0), Fr(-1)]), X.pscale([Fr(0)] + U, -1)), [nu])
    assert X.same_poly(lhs, q), 'closed form of int x^m / sqrt(a^2 - x^2) failed'
    return U, nu


def divmod_R2(P, R2):
    """P = q (R2 - x^2) + lam x + mu, exactly"""
    p = list(P) + [Fr(0)] * 2
    q = [Fr(0)] * max(len(P) - 2, 0)
    for d in range(len(P) - 1, 1, -1):
        c = p[d]
        if c:
            q[d - 2] -= c                  # -c x^(d-2) (R2 - x^2) = c x^d - c R2 x^(d-2)
            p[d - 2] += c * R2
            p[d] = Fr(0)
    lam, mu = p[1], p[0]
    assert X.same_poly(X.padd(X.pmul(q, [R2, Fr(0), Fr(-1)]), [mu, lam]), list(P))
    return q, lam, mu


# ------------------------------------------------------------------------------------------------- end-point values
def _x(xn):
    return xn.iv()


def _sq(xn):
    s = xn.sq()
    assert s is not None
    return s


def _sqrt_q(q):
    assert q >= 0
    return iv.sqrt(ivq(q)) if q else iv.mpf(0)


def _asin01(Z):
    return RG.iv_asin(RG._clip01(Z))


def H_int(A, a2, b, R, R2, xa, xb):
    """int_{xa}^{xb} A(x) x / ((R2 - x^2) sqrt(a2 - x^2)) dx,  a2 + b^2 = R2, b > 0"""
    assert a2 + b * b == R2 and b > 0 and R > 0
    q, lam, mu = divmod_R2(X.pmul(A, [Fr(0), Fr(1)]), R2)

    def F(xn):
        x2 = _sq(xn)
        sig = _sqrt_q(a2 - x2)
        sR = _sqrt_q(R2 - x2)
        val = iv.mpf(0)
        if X.ptrim(q):
            U, nu = inv_sqrt_antider(q, a2)
            val += RG.iv_poly(U, _x(xn)) * sig
            if nu:
                val += ivq(nu) * _asin01(_x(xn) / _sqrt_q(a2))
        if lam:
            val += ivq(-lam / b) * _asin01(sig / sR)
        if mu:
            val += ivq(mu / (R * b)) * _asin01(ivq(b) * _x(xn) / (_sqrt_q(a2) * sR))
        return val
    return F(xb) - F(xa)


def asin_term(T, kind, par, R, R2, xa, xb):
    """int_{xa}^{xb} T(x) asin(g(x)) dx, by parts:
    const y_c:   g = y_c / sqrt(R2 - x^2),                 = [A asin g] - y_c H(A, R2 - y_c^2, y_c)
    circle z_i:  g = sqrt(rho2 - x^2)/sqrt(R2 - x^2),     = [A asin g] + z_i H(A, rho2, z_i)   (rho2 = R2 - z_i^2)"""
    A = X.pint(T)
    if kind == 'const':
        yc = par
        if yc == 0:
            return iv.mpf(0)

        def g(xn):
            return ivq(yc) / _sqrt_q(R2 - _sq(xn))
        bt = RG.iv_poly(A, _x(xb)) * _asin01(g(xb)) - RG.iv_poly(A, _x(xa)) * _asin01(g(xa))
        return bt - ivq(yc) * H_int(A, R2 - yc * yc, yc, R, R2, xa, xb)
    rho2, zi = par
    if zi == 0:
        return (iv.pi / 2) * (RG.iv_poly(A, _x(xb)) - RG.iv_poly(A, _x(xa)))

    def g(xn):
        return _sqrt_q(rho2 - _sq(xn)) / _sqrt_q(R2 - _sq(xn))
    bt = RG.iv_poly(A, _x(xb)) * _asin01(g(xb)) - RG.iv_poly(A, _x(xa)) * _asin01(g(xa))
    return bt + ivq(zi) * H_int(A, rho2, zi, R, R2, xa, xb)


def sqrt_poly_int(S, a2, xa, xb):
    """int_{xa}^{xb} S(x) sqrt(a2 - x^2) dx (closed form of region.sqrt_antider)"""
    if not X.ptrim(S):
        return iv.mpf(0)
    P, c = RG.sqrt_antider(S, a2)

    def F(xn):
        v = RG.iv_poly(P, _x(xn)) * _sqrt_q(a2 - _sq(xn))
        if c:
            v += ivq(c) * _asin01(_x(xn) / _sqrt_q(a2))
        return v
    return F(xb) - F(xa)


def poly_int(E, xa, xb):
    if not X.ptrim(E):
        return iv.mpf(0)
    A = X.pint(E)
    return RG.iv_poly(A, _x(xb)) - RG.iv_poly(A, _x(xa))


# ------------------------------------------------------------------------------------------------- the Z part
def compose_c2(pu, R2):
    """p(R2 - x^2) for p given in powers of u"""
    out = []
    for c in reversed(list(pu)):
        out = X.padd(X.pmul(out, [R2, Fr(0), Fr(-1)]), [c])
    return out


def S_and_T(Wn, R2):
    """for W(x,y) = sum_n Wn[n](x) y^n:  int W sqrt(c^2 - y^2) dy = S(x,y) sqrt(c^2 - y^2) + T(x) asin(y/c),
    c^2 = R2 - x^2;  returns S as {j: x-poly} and T"""
    S, T = {}, []
    c2 = [R2, Fr(0), Fr(-1)]
    for n, w in enumerate(Wn):
        if not X.ptrim(w):
            continue
        assert check_Pn(n)
        P, kap = Pn(n)
        for j, pu in P.items():
            S[j] = X.padd(S.get(j, []), X.pmul(w, compose_c2(pu, R2)))
        if kap:
            T = X.padd(T, X.pscale(X.pmul(w, X.ppow(c2, (n + 2) // 2)), kap))
    return S, T


def Y_int(S, T, fn, R, R2, xa, xb):
    """int_{xa}^{xb} [S(x, y) sqrt(R2 - x^2 - y^2) + T(x) asin(y / sqrt(R2 - x^2))]_{y = f(x)} dx"""
    kind = fn[0]
    if kind == 'const':
        yc = fn[1]
        Sc = []
        p = Fr(1)
        for j in range(max(S) + 1 if S else 0):
            if j in S:
                Sc = X.padd(Sc, X.pscale(S[j], p))
            p *= yc
        val = sqrt_poly_int(Sc, R2 - yc * yc, xa, xb)
        if X.ptrim(T):
            val += asin_term(T, 'const', yc, R, R2, xa, xb)
        return val
    rho2, zi = fn[1], fn[2]
    # y = sqrt(rho2 - x^2): sqrt(R2 - x^2 - y^2) = zi;  S(x, y) = E(x) + O(x) sqrt(rho2 - x^2)
    E, O = [], []
    w = [rho2, Fr(0), Fr(-1)]
    for j, sj in S.items():
        part = X.pmul(sj, X.ppow(w, j // 2))
        if j % 2 == 0:
            E = X.padd(E, part)
        else:
            O = X.padd(O, part)
    val = ivq(zi) * (poly_int(E, xa, xb) + sqrt_poly_int(O, rho2, xa, xb)) if zi else iv.mpf(0)
    if X.ptrim(T):
        val += asin_term(T, 'circle', (rho2, zi), R, R2, xa, xb)
    return val


def Z_part(Wn, x0, x1, y0, y1, z0, z1, R, R2):
    """int over {(x,y) in rect : r1^2 < x^2 + y^2 < r0^2} of W(x,y) sqrt(R2 - x^2 - y^2), r_i^2 = R2 - z_i^2"""
    r02, r12 = R2 - z0 * z0, R2 - z1 * z1
    S, T = S_and_T(Wn, R2)
    if not S and not X.ptrim(T):
        return iv.mpf(0)
    lows = [('const', y0)] + ([('circle', r12, z1)] if r12 > 0 else [])
    ups = [('const', y1), ('circle', r02, z0)]
    bps = [Num('rat', x0), Num('rat', x1)]
    for (kind, *par) in lows + ups:
        if kind == 'circle':
            rho2 = par[0]
            bps.append(Num('sqrt', rho2))
            for yc in (y0, y1):
                if rho2 - yc * yc >= 0:
                    bps.append(Num('sqrt', rho2 - yc * yc))
    lo, hi = bps[0], bps[1]
    pts = sorted([b for b in bps if num_cmp(b, lo) >= 0 and num_cmp(b, hi) <= 0], key=cmp_to_key(num_cmp))
    uniq = []
    for b in pts:
        if not uniq or num_cmp(uniq[-1], b) != 0:
            uniq.append(b)
    total = iv.mpf(0)
    for a, b in zip(uniq[:-1], uniq[1:]):
        A_, B_ = a.iv(), b.iv()
        assert A_.b < B_.a
        xq = (RG.hi_q(A_) + RG.lo_q(B_)) / 2
        # values squared at the test point (all functions are >= 0): exact comparisons
        def val2(fn):
            if fn[0] == 'const':
                return fn[1] * fn[1]
            return fn[1] - xq * xq if fn[1] - xq * xq > 0 else None
        lo_fn, lo_v = None, None
        for fn in lows:
            v = val2(fn)
            if v is None:
                continue
            if lo_v is None or v > lo_v:
                lo_fn, lo_v = fn, v
            else:
                assert v != lo_v
        up_fn, up_v, empty = None, None, False
        for fn in ups:
            v = val2(fn)
            if v is None:
                empty = True
                break
            if up_v is None or v < up_v:
                up_fn, up_v = fn, v
            else:
                assert v != up_v
        if empty or not lo_v < up_v:
            continue
        total += Y_int(S, T, up_fn, R, R2, a, b) - Y_int(S, T, lo_fn, R, R2, a, b)
    return total


# ------------------------------------------------------------------------------------------------- the whole box
def _local2(P2, x0, h, y0, k):
    """global bivariate array P2[a, b] (x^a y^b) -> local coefficients for x = x0 + h s, y = y0 + k t"""
    na, nb = P2.shape
    cols = [X.pcompose_lin(list(P2[:, b]), x0, h) for b in range(nb)]
    tmp = X.zeros((max(len(c) for c in cols), nb))
    for b, c in enumerate(cols):
        for a, v in enumerate(c):
            tmp[a, b] = v
    out = X.zeros(tmp.shape)
    for a in range(tmp.shape[0]):
        row = X.pcompose_lin(list(tmp[a, :]), y0, k)
        for b, v in enumerate(row):
            out[a, b] = v
    return out


def _b_pow(base, k):
    out = X.zeros((1, 1)) + 1
    for _ in range(k):
        out = X.b_mul(out, base)
    return out


def ball_part(f, box, R2):
    """interval enclosure of the integral of the global polynomial f (array f[a,b,c] of x^a y^b z^c) over
    box cap {x^2 + y^2 + z^2 <= R2}; box = ((x0,x1),(y0,y1),(z0,z1)) in the first octant"""
    (x0, x1), (y0, y1), (z0, z1) = box
    assert min(x0, y0, z0) >= 0
    R = Fr(1)
    R = Fr(int(round(float(R2) ** 0.5 * 10 ** 6)), 10 ** 6)
    assert R * R == R2, 'R must be rational'
    if z0 * z0 >= R2 or x0 * x0 + y0 * y0 + z0 * z0 >= R2:
        return iv.mpf(0)
    na, nb, nc = f.shape
    G = X.zeros((na, nb, nc + 1))
    for c in range(nc):
        G[:, :, c + 1] = f[:, :, c] / (c + 1)
    Q = G.copy()
    at0 = X.zeros((na, nb))
    for c in range(nc + 1):
        at0 = at0 + G[:, :, c] * z0 ** c
    Q[:, :, 0] = Q[:, :, 0] - at0
    r02, r12 = R2 - z0 * z0, R2 - z1 * z1
    h, k = x1 - x0, y1 - y0
    total = iv.mpf(0)
    # part A: rect cap D(r1), Q(x, y, z1)
    if r12 > 0:
        Qz1 = X.zeros((na, nb))
        for c in range(nc + 1):
            Qz1 = Qz1 + Q[:, :, c] * z1 ** c
        total += RG.integrate_region(_local2(Qz1, x0, h, y0, k), x0, h, y0, k, [], [('circle', r12)])
    # part B: rect cap (D(r0) minus D(r1)), Q(x, y, Z) = P + W Z
    base = X.zeros((3, 3))
    base[0, 0], base[2, 0], base[0, 2] = R2, Fr(-1), Fr(-1)          # R2 - x^2 - y^2
    P = X.zeros((na + nc + 1, nb + nc + 1))
    W = X.zeros((na + nc + 1, nb + nc + 1))
    for c in range(nc + 1):
        Qc = Q[:, :, c]
        if all(v == 0 for v in Qc.flat):
            continue
        term = X.b_mul(Qc, _b_pow(base, c // 2))
        if c % 2 == 0:
            P[:term.shape[0], :term.shape[1]] += term
        else:
            W[:term.shape[0], :term.shape[1]] += term
    lows = [('circle', r12)] if r12 > 0 else []
    total += RG.integrate_region(_local2(P, x0, h, y0, k), x0, h, y0, k, lows, [('circle', r02)])
    Wn = [list(W[:, n]) for n in range(W.shape[1])]
    total += Z_part(Wn, x0, x1, y0, y1, z0, z1, R, R2)
    return total


def to_global(f_loc, e, N):
    """local trivariate array on cell e of the grid 1/N (s = N x - e_x, ...) -> global array (exact)"""
    from math import comb
    out = f_loc
    for d in range(3):
        n = out.shape[d]
        Tm = np.zeros((n, n), dtype=object)           # (N x - e)^a = sum_a2 C(a, a2) N^a2 (-e)^(a - a2) x^a2, integers
        for a in range(n):
            for a2 in range(a + 1):
                Tm[a, a2] = comb(a, a2) * N ** a2 * (-e[d]) ** (a - a2)
        out = np.moveaxis(np.tensordot(Tm.T, np.moveaxis(out, d, 0), axes=([1], [0])), 0, d)
    return out
