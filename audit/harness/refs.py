"""Reference values for the accuracy audit.

Each mathematical function has
  * Arb formulations (python-flint): rigorous balls; the precision is raised until the
    ball has >= `min_bits` bits of relative accuracy (several formulations are tried, for
    example P = 1 - Q when the direct series of P is slow);
  * an mpmath formulation, evaluated at the current mp precision.  For the incomplete
    gamma and beta functions with large parameters, where mpmath's built-ins raise
    NoConvergence or stall, the formulation is a continued fraction / power series in
    mpmath arithmetic (DLMF 8.7.1, 8.9.2, 8.17.22) with a convergence test below the
    working precision;
  * optionally a second mpmath formulation (the built-in) used as an extra check.

Inverse functions are a safeguarded Newton/bisection solve of the forward function.
All entry points take Python floats (exact binary inputs).
"""
import math

import flint
import mpmath
from mpmath import mp, mpf

arb = flint.arb
acb = flint.acb

# =======================================================================================
# Arb helpers
# =======================================================================================


def _arb_to_mpf(v):
    man, exp = v.mid().man_exp()
    return mpf((int(man), int(exp)))


def mpf_to_arb(x):
    """Exact conversion of an mpf to an arb (the context precision must exceed the
    mantissa length; callers use >= 512 bits)."""
    x = mpf(x)
    if x == 0:
        return arb(0)
    s, mant, e, bc = x._mpf_
    v = arb(flint.fmpz(int(mant))) * arb(2) ** e
    return -v if s else v


def _to_arb(x):
    if isinstance(x, (float, int)):
        return arb(x)
    return mpf_to_arb(x)


def arb_eval(forms, args, min_bits=64, precs=(128, 256, 512, 1024, 2048), pert=None):
    """Try each Arb formulation with rising precision until the ball is tight.

    args are floats or mpf (converted exactly); pert = (i, h) multiplies args[i] by
    (1 + h) exactly (h a power of two) for derivative estimates.
    Returns (mpf midpoint, arb ball, formulation index) or (None, None, None)."""
    old = flint.ctx.prec
    try:
        for i, fn in enumerate(forms):
            for p in precs:
                flint.ctx.prec = p
                a = [_to_arb(x) for x in args]
                if pert is not None:
                    a[pert[0]] = a[pert[0]] * (1 + arb(pert[1]))
                try:
                    v = fn(*a)
                except (ValueError, ZeroDivisionError, AssertionError, NotImplementedError):
                    break
                if v is None:
                    break
                if isinstance(v, acb):
                    if not v.imag.contains(0):
                        break
                    v = v.real
                if v.is_zero():
                    return mpf(0), v, i
                if v.is_finite() and v.rel_accuracy_bits() >= min_bits:
                    return _arb_to_mpf(v), v, i
        return None, None, None
    finally:
        flint.ctx.prec = old


def lt_half(v):
    """Decide v < 1/2 for a ball (by its midpoint; only used to pick a formulation)."""
    return float(v.mid()) < 0.5


# =======================================================================================
# mpmath: incomplete gamma and beta in mpmath arithmetic (series / continued fractions)
# =======================================================================================
MAXIT = 2_000_000


def _tol():
    return mpf(2) ** (-mp.prec - 6)


def _lentz(an, bn, b0, tol, maxit=MAXIT):
    """Modified Lentz for b0 + a1/(b1 + a2/(b2 + ...)); an(n), bn(n) for n>=1."""
    tiny = mpf(2) ** (-mp.prec * 4)
    f = b0 if b0 != 0 else tiny
    C, D = f, mpf(0)
    for n in range(1, maxit):
        a, b = an(n), bn(n)
        D = b + a * D
        D = tiny if D == 0 else D
        C = b + a / C
        C = tiny if C == 0 else C
        D = 1 / D
        delta = C * D
        f *= delta
        if abs(delta - 1) < tol:
            return f
    raise mpmath.libmp.NoConvergence("continued fraction did not converge")


def _extra_bits(*mags):
    m = max([abs(float(mpmath.log(abs(v) + 2, 2))) for v in mags] + [1.0])
    return int(m) + 20


def gamma_P_series(a, x):
    """P(a,x) = x^a e^-x / Gamma(a+1) * sum_k x^k / ((a+1)...(a+k))   (DLMF 8.7.1)."""
    with mp.extraprec(_extra_bits(a * mp.log(x), x, mp.loggamma(a + 1)) + 10):
        lpre = a * mp.log(x) - x - mp.loggamma(a + 1)
        s = term = mpf(1)
        tol = _tol()
        k = 0
        while True:
            k += 1
            term = term * x / (a + k)
            s += term
            if term < tol * s and k > x - a:
                break
            if k > MAXIT:
                raise mpmath.libmp.NoConvergence("series")
        return mp.exp(lpre + mp.log(s))


def gamma_Q_cf(a, x):
    """Q(a,x) by the Legendre continued fraction (DLMF 8.9.2), for x > a - 1."""
    with mp.extraprec(_extra_bits(a * mp.log(x), x, mp.loggamma(a)) + 10):
        lpre = a * mp.log(x) - x - mp.loggamma(a)
        # Gamma(a,x) = e^-x x^a / (x+1-a - 1(1-a)/(x+3-a - 2(2-a)/(x+5-a - ...)))
        cf = _lentz(lambda n: -n * (n - a), lambda n: x + 2 * n + 1 - a, x + 1 - a, _tol())
        return mp.exp(lpre) / cf


def gamma_PQ(a, x):
    """Return (P, Q) at the current precision, each with relative accuracy."""
    a, x = mpf(a), mpf(x)
    if x == 0:
        return mpf(0), mpf(1)
    if x == mp.inf:
        return mpf(1), mpf(0)
    if a < 30:
        # mpmath built-ins are reliable and fast for small a
        P = mp.gammainc(a, 0, x, regularized=True)
        Q = mp.gammainc(a, x, mp.inf, regularized=True)
        return P, Q
    if x < a + 1:
        P = gamma_P_series(a, x)
        with mp.extraprec(10):
            Q = 1 - P
        return P, Q
    Q = gamma_Q_cf(a, x)
    with mp.extraprec(10):
        P = 1 - Q
    return P, Q


def _beta_cf(a, b, x, y):
    """I_x(a,b) by the continued fraction DLMF 8.17.22; y = 1 - x exactly.  Best for
    x < (a+1)/(a+b+2)."""
    with mp.extraprec(_extra_bits(a * mp.log(x), b * mp.log(y), mp.loggamma(a + b)) + 10):
        lpre = a * mp.log(x) + b * mp.log(y) - mp.log(a) - (mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b))

        def dn(n):
            m = n // 2
            if n % 2 == 0:
                return m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
            return -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))

        cf = _lentz(dn, lambda n: mpf(1), mpf(1), _tol())
        return mp.exp(lpre) / cf


def beta_IJ(a, b, x, y=None):
    """Return (I_x(a,b), I_y(b,a)) with y = 1-x (exact if not given), both relative-accurate."""
    a, b, x = mpf(a), mpf(b), mpf(x)
    if y is None:
        y = mp.fsub(1, x, exact=True)
    else:
        y = mpf(y)
    if x == 0:
        return mpf(0), mpf(1)
    if y == 0:
        return mpf(1), mpf(0)
    if x * (a + b + 2) < a + 1:
        I = _beta_cf(a, b, x, y)
        with mp.extraprec(10):
            return I, 1 - I
    J = _beta_cf(b, a, y, x)
    with mp.extraprec(10):
        return 1 - J, J


# =======================================================================================
# Forward functions: mpmath (args mpf) and Arb (args arb balls)
# =======================================================================================

# ---- incomplete beta / gamma and the distributions built on them ----------------------


def m_betainc(a, b, x):
    return beta_IJ(a, b, x)[0]


def m_betaincc(a, b, x):
    return beta_IJ(a, b, x)[1]


def mb_betainc(a, b, x):
    return mp.betainc(a, b, 0, x, regularized=True)


def mb_betaincc(a, b, x):
    return mp.betainc(b, a, 0, mp.fsub(1, x, exact=True), regularized=True)


A_betainc = [lambda a, b, x: x.beta_lower(a, b, 1),
             lambda a, b, x: 1 - (1 - x).beta_lower(b, a, 1)]
A_betaincc = [lambda a, b, x: (1 - x).beta_lower(b, a, 1),
              lambda a, b, x: 1 - x.beta_lower(a, b, 1)]


def m_gammainc(a, x):
    return gamma_PQ(a, x)[0]


def m_gammaincc(a, x):
    return gamma_PQ(a, x)[1]


def mb_gammainc(a, x):
    return mp.gammainc(a, 0, x, regularized=True)


def mb_gammaincc(a, x):
    return mp.gammainc(a, x, mp.inf, regularized=True)


A_gammainc = [lambda a, x: x.gamma_lower(a, 1) if x < a + 1 else 1 - x.gamma_upper(a, 1),
              lambda a, x: 1 - x.gamma_upper(a, 1),
              lambda a, x: x.gamma_lower(a, 1)]
A_gammaincc = [lambda a, x: x.gamma_upper(a, 1),
               lambda a, x: 1 - x.gamma_lower(a, 1)]


# Student t: tails via I_w(df/2, 1/2), w = df/(df+t^2), 1-w = t^2/(df+t^2)
def _t_tails_m(df, t):
    """Return (P(T<=t), P(T>t))."""
    if t == 0:
        return mpf(1) / 2, mpf(1) / 2
    with mp.extraprec(1200):
        t2 = t * t
        w = df / (df + t2)
        y = t2 / (df + t2)
    I, J = beta_IJ(df / 2, mpf(1) / 2, w, y)  # I = I_w(df/2,1/2) = 2*tail, J = 1 - I
    tail = I / 2
    with mp.extraprec(10):
        body = 1 - tail
    return (tail, body) if t < 0 else (body, tail)


def _t_tails_a(df, t):
    if t.is_zero():
        return arb(1) / 2, arb(1) / 2
    t2 = t * t
    w = df / (df + t2)
    y = t2 / (df + t2)
    tail = w.beta_lower(df / 2, arb(1) / 2, 1) / 2
    if not (tail.rel_accuracy_bits() > 60):
        tail = (1 - y.beta_lower(arb(1) / 2, df / 2, 1)) / 2
    body = 1 - tail
    return (tail, body) if t < 0 else (body, tail)


def _t_tails_a2(df, t):
    if t.is_zero():
        return arb(1) / 2, arb(1) / 2
    t2 = t * t
    y = t2 / (df + t2)
    tail = (1 - y.beta_lower(arb(1) / 2, df / 2, 1)) / 2
    body = 1 - tail
    return (tail, body) if t < 0 else (body, tail)


def m_stdtr(df, t):
    return _t_tails_m(df, t)[0]


A_stdtr = [lambda df, t: _t_tails_a(df, t)[0], lambda df, t: _t_tails_a2(df, t)[0]]


def m_fdtr(dfn, dfd, x):
    with mp.extraprec(1200):
        w = dfn * x / (dfn * x + dfd)
        y = dfd / (dfn * x + dfd)
    return beta_IJ(dfn / 2, dfd / 2, w, y)[0]


A_fdtr = [lambda dfn, dfd, x: (dfn * x / (dfn * x + dfd)).beta_lower(dfn / 2, dfd / 2, 1),
          lambda dfn, dfd, x: 1 - (dfd / (dfn * x + dfd)).beta_lower(dfd / 2, dfn / 2, 1)]


def m_chdtr(v, x):
    return gamma_PQ(v / 2, x / 2)[0]


A_chdtr = [lambda v, x: A_gammainc[0](v / 2, x / 2), lambda v, x: A_gammainc[1](v / 2, x / 2),
           lambda v, x: A_gammainc[2](v / 2, x / 2)]


def m_bdtr(k, n, p):
    k = mp.floor(k)
    if k < 0:
        return mpf(0)
    if k >= n:
        return mpf(1)
    # P(X <= k) = I_{1-p}(n-k, k+1) = 1 - I_p(k+1, n-k)
    return beta_IJ(k + 1, n - k, p)[1]


def _a_bdtr(k, n, p, flip):
    k = k.floor()
    if k < 0:
        return arb(0)
    if k >= n:
        return arb(1)
    if flip:
        return 1 - p.beta_lower(k + 1, n - k, 1)
    return (1 - p).beta_lower(n - k, k + 1, 1)


A_bdtr = [lambda k, n, p: _a_bdtr(k, n, p, False), lambda k, n, p: _a_bdtr(k, n, p, True)]

# ---- gamma-function family -------------------------------------------------------------


def m_gammaln(x):
    return mp.re(mp.loggamma(x))


def _a_lgamma(x):
    if x > 0:
        return x.lgamma()
    return acb(x).lgamma().real


A_gammaln = [_a_lgamma]


def m_digamma(x):
    return mp.digamma(x)


A_digamma = [lambda x: x.digamma()]


def m_polygamma(n, x):
    return mp.psi(int(n), x)


A_polygamma = [lambda n, x: acb(x).polygamma(acb(n)).real]

# ---- error function family ---------------------------------------------------------------


def m_erf(x):
    return mp.erf(x)


def m_erfc(x):
    return mp.erfc(x)


def m_erfcx(x):
    return mp.erfc(x) * mp.exp(x * x)


A_erf = [lambda x: x.erf()]
A_erfc = [lambda x: x.erfc()]
A_erfcx = [lambda x: x.erfc() * (x * x).exp()]


def m_ndtr(x):
    return mp.erfc(-x / mp.sqrt(2)) / 2


def m_log_ndtr(x):
    if x > 0:
        return mp.log1p(-m_ndtr(-x))
    return mp.log(m_ndtr(x))


def _a_ndtr(x):
    return (-x / arb(2).sqrt()).erfc() / 2


A_ndtr = [_a_ndtr]
A_log_ndtr = [lambda x: (-_a_ndtr(-x)).log1p() if x > 0 else _a_ndtr(x).log()]

# ---- Bessel, hypergeometric, zeta -----------------------------------------------------


def m_jv(v, x):
    return mp.besselj(v, x)


def m_iv(v, x):
    return mp.besseli(v, x)


def m_kv(v, x):
    return mp.besselk(v, x)


def m_i0e(x):
    return mp.besseli(0, x) * mp.exp(-abs(x))


def m_i1e(x):
    return mp.besseli(1, x) * mp.exp(-abs(x))


A_jv = [lambda v, x: x.bessel_j(v)]
A_iv = [lambda v, x: x.bessel_i(v)]
A_kv = [lambda v, x: x.bessel_k(v)]
A_i0e = [lambda x: x.bessel_i(0, scaled=True)]
A_i1e = [lambda x: x.bessel_i(1, scaled=True)]


def m_hyp2f1(a, b, c, x):
    return mp.hyp2f1(a, b, c, x)


def m_hyp1f1(a, b, x):
    return mp.hyp1f1(a, b, x)


A_hyp2f1 = [lambda a, b, c, x: x.hypgeom_2f1(a, b, c)]
A_hyp1f1 = [lambda a, b, x: x.hypgeom_1f1(a, b)]


def m_zeta(x):
    return mp.zeta(x)


def m_zeta2(s, q):
    return mp.zeta(s, q)


A_zeta = [lambda x: x.zeta()]
A_zeta2 = [lambda s, q: s.zeta(q)]

# =======================================================================================
# Distributions (scipy.stats): (lower, upper) tail probabilities, args (x, *shape)
# =======================================================================================


def _m_norm(x):
    return m_ndtr(x), m_ndtr(-x)


def _a_norm(x):
    return _a_ndtr(x), _a_ndtr(-x)


def _m_t(x, df):
    return _t_tails_m(df, x)


def _a_t(x, df):
    return _t_tails_a(df, x)


def _m_chi2(x, df):
    return gamma_PQ(df / 2, x / 2)


def _a_chi2(x, df):
    return _a_gamma(x / 2, df / 2)


def _m_gamma(x, a):
    return gamma_PQ(a, x)


def _a_gamma(x, a):
    Q = x.gamma_upper(a, 1)
    if lt_half(Q) or not (x < a + 1):
        P = x.gamma_lower(a, 1) if x < a + 1 else 1 - Q
        if not (P.rel_accuracy_bits() > 60):
            P = 1 - Q
        return P, Q
    P = x.gamma_lower(a, 1)
    return P, 1 - P


def _m_f(x, dfn, dfd):
    with mp.extraprec(1200):
        w = dfn * x / (dfn * x + dfd)
        y = dfd / (dfn * x + dfd)
    return beta_IJ(dfn / 2, dfd / 2, w, y)


def _a_f(x, dfn, dfd):
    w = dfn * x / (dfn * x + dfd)
    y = dfd / (dfn * x + dfd)
    return w.beta_lower(dfn / 2, dfd / 2, 1), y.beta_lower(dfd / 2, dfn / 2, 1)


def _m_beta(x, a, b):
    return beta_IJ(a, b, x)


def _a_beta(x, a, b):
    return x.beta_lower(a, b, 1), (1 - x).beta_lower(b, a, 1)


def _m_binom(k, n, p):
    k = mp.floor(k)
    if k < 0:
        return mpf(0), mpf(1)
    if k >= n:
        return mpf(1), mpf(0)
    I, J = beta_IJ(k + 1, n - k, p)  # I = I_p(k+1, n-k) = P(X > k)
    return J, I


def _a_binom(k, n, p):
    k = k.floor()
    if k < 0:
        return arb(0), arb(1)
    if k >= n:
        return arb(1), arb(0)
    return (1 - p).beta_lower(n - k, k + 1, 1), p.beta_lower(k + 1, n - k, 1)


DISTS = {
    'norm': (_m_norm, _a_norm),
    't': (_m_t, _a_t),
    'chi2': (_m_chi2, _a_chi2),
    'gamma': (_m_gamma, _a_gamma),
    'f': (_m_f, _a_f),
    'beta': (_m_beta, _a_beta),
    'binom': (_m_binom, _a_binom),
}


def m_dist(dist, method):
    mf = DISTS[dist][0]
    if method == 'sf':
        return lambda *a: mf(*a)[1]
    if method == 'cdf':
        return lambda *a: mf(*a)[0]
    if method == 'logsf':
        def g(*a):
            lo, up = mf(*a)
            return mp.log(up) if up < 0.5 else mp.log1p(-lo)
        return g
    if method == 'logcdf':
        def g(*a):
            lo, up = mf(*a)
            return mp.log(lo) if lo < 0.5 else mp.log1p(-up)
        return g
    raise ValueError(method)


def a_dist(dist, method):
    af = DISTS[dist][1]
    if method == 'sf':
        return [lambda *a: af(*a)[1]]
    if method == 'cdf':
        return [lambda *a: af(*a)[0]]
    if method == 'logsf':
        def g(*a):
            lo, up = af(*a)
            return up.log() if lt_half(up) else (-lo).log1p()
        return [g]
    if method == 'logcdf':
        def g(*a):
            lo, up = af(*a)
            return lo.log() if lt_half(lo) else (-up).log1p()
        return [g]
    raise ValueError(method)


# =======================================================================================
# Inverse functions
# =======================================================================================


def _midpoint(lo, hi):
    if lo > 0 and hi / lo > 4:
        return mpmath.sqrt(lo * hi)
    if hi < 0 and lo / hi > 4:
        return -mpmath.sqrt(lo * hi)
    if lo < 0 < hi and (hi - lo) > 4:
        return mpf(0)
    return (lo + hi) / 2


def solve_monotone(F, y, lo, hi, increasing=True, x0=None, fprime=None, max_iter=3000):
    """Solve F(x) = y on (lo, hi) for monotone F at the current precision.

    Safeguarded Newton (when fprime is given) inside a shrinking bracket; bisection is in
    geometric means when the bracket spans decades.  Returns mpf or None."""
    sgn = 1 if increasing else -1
    G = lambda x: sgn * (F(x) - y)
    # a root outside the bracket: the brackets reach beyond the double range, so the
    # correctly rounded answer is then +-inf (|root| > 1e320), 0 (root < 1e-330) or 1
    if G(lo) > 0:
        if lo == -BIG:
            return -mp.inf
        if lo == SMALL:
            return SMALL / 10
        return None
    if G(hi) < 0:
        if hi == BIG:
            return mp.inf
        if hi == 1 - SMALL:
            return 1 - SMALL / 10
        return None
    tol = mpf(2) ** (-mp.prec + 12)
    x = None
    if x0 is not None and math.isfinite(x0):
        x0 = mpf(x0)
        if lo < x0 < hi:
            x = x0
    for _ in range(max_iter):
        if x is None:
            x = _midpoint(lo, hi)
        g = G(x)
        if g == 0:
            return x
        if g < 0:
            lo = x
        else:
            hi = x
        if abs(hi - lo) <= tol * max(abs(lo), abs(hi)):
            return (lo + hi) / 2
        xn = None
        if fprime is not None:
            d = sgn * fprime(x)
            if d and mpmath.isfinite(d) and d > 0:
                xn = x - g / d
                if not (lo < xn < hi):
                    xn = None
                elif abs(xn - x) <= tol * abs(x) / 64:
                    return xn
        x = xn if xn is not None else _midpoint(lo, hi)
    return x


BIG = mpf(10) ** 320
SMALL = mpf(10) ** -330


def inv_ndtri(p, x0=None):
    return solve_monotone(m_ndtr, mpf(p), mpf(-40), mpf(40), True, x0, lambda x: mp.npdf(x))


def inv_erfcinv(y, x0=None):
    return solve_monotone(mp.erfc, mpf(y), mpf(-30), mpf(30), False, x0,
                          lambda x: -2 / mp.sqrt(mp.pi) * mp.exp(-x * x))


def inv_gammaincinv(a, y, x0=None):
    a = mpf(a)
    lg = mp.loggamma(a)
    return solve_monotone(lambda x: m_gammainc(a, x), mpf(y), SMALL, BIG, True, x0,
                          lambda x: mp.exp((a - 1) * mp.log(x) - x - lg))


def inv_betaincinv(a, b, y, x0=None):
    a, b = mpf(a), mpf(b)
    lb = mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
    return solve_monotone(lambda x: m_betainc(a, b, x), mpf(y), SMALL, 1 - SMALL, True, x0,
                          lambda x: mp.exp((a - 1) * mp.log(x) + (b - 1) * mp.log1p(-x) - lb))


def _t_pdf(df):
    lc = mp.loggamma((df + 1) / 2) - mp.loggamma(df / 2) - mp.log(df * mp.pi) / 2
    return lambda t: mp.exp(lc - (df + 1) / 2 * mp.log1p(t * t / df))


def inv_stdtrit(df, p, x0=None):
    df = mpf(df)
    return solve_monotone(lambda t: m_stdtr(df, t), mpf(p), -BIG, BIG, True, x0, _t_pdf(df))


def _pdf_m(dist, x, *s):
    if dist == 'norm':
        return mp.npdf(x)
    if dist == 't':
        return _t_pdf(s[0])(x)
    if dist in ('chi2', 'gamma'):
        a, sc = (s[0] / 2, 2) if dist == 'chi2' else (s[0], 1)
        z = x / sc
        return mp.exp((a - 1) * mp.log(z) - z - mp.loggamma(a)) / sc
    if dist == 'f':
        d1, d2 = s
        lb = mp.loggamma(d1 / 2) + mp.loggamma(d2 / 2) - mp.loggamma((d1 + d2) / 2)
        return mp.exp(d1 / 2 * mp.log(d1 / d2) + (d1 / 2 - 1) * mp.log(x)
                      - (d1 + d2) / 2 * mp.log1p(d1 * x / d2) - lb)
    if dist == 'beta':
        a, b = s
        lb = mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
        return mp.exp((a - 1) * mp.log(x) + (b - 1) * mp.log1p(-x) - lb)
    raise ValueError(dist)


ISF_BRACKET = {'norm': (mpf(-40), mpf(40)), 't': (-BIG, BIG), 'chi2': (SMALL, BIG),
               'gamma': (SMALL, BIG), 'f': (SMALL, BIG), 'beta': (SMALL, 1 - SMALL)}


def inv_isf(dist):
    mf = DISTS[dist][0]

    def solve(q, *shape, x0=None):
        s = [mpf(v) for v in shape]
        lo, hi = ISF_BRACKET[dist]
        return solve_monotone(lambda x: mf(x, *s)[1], mpf(q), lo, hi, False, x0,
                              lambda x: -_pdf_m(dist, x, *s))
    return solve


def binom_isf_ref(q, n, p):
    """Smallest integer k >= -1 with P(X > k) <= q (scipy's discrete isf convention).
    Returns (k, sf(k-1), sf(k)) as mpf, by bisection on k."""
    q, n, p = mpf(q), mpf(n), mpf(p)
    sf = lambda k: _m_binom(mpf(k), n, p)[1]
    lo, hi = -1, int(n)  # sf(-1) = 1 > q (q < 1), sf(n) = 0 <= q
    if sf(lo) <= q:
        return lo, None, sf(lo)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if sf(mid) <= q:
            hi = mid
        else:
            lo = mid
    return hi, sf(hi - 1), sf(hi)
