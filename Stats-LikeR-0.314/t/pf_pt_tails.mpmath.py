#!/usr/bin/env python3
"""Regenerate the frozen tables in t/pf_pt_tails.R.mpmath.t.

    python3 t/pf_pt_tails.mpmath.py

and paste the output into the first half of that file's __DATA__ block (the
part before the %%R%% marker).  The test reads only those frozen literals; it
never runs this script, and the distribution needs neither python nor mpmath
to install or test.

WHY A THIRD OPINION IS NEEDED HERE
----------------------------------
LikeR and R disagreed on both of these, in opposite directions:

  * pf(q, df1, df2, lower.tail = FALSE) for small q.  R had
    0.99999920211563864 at q = 1e-12, df = (1, 1e6); LikeR returned exactly 1.
  * qf(p, df1, df2) in the far lower tail.  R had 1.5765166949677223e-12 at
    p = 1e-6, df = (1, 100); LikeR had 1.5786698446096233e-12 -- a 1.4e-3
    relative gap, so one of them is wrong in the third digit.

mpmath at mp.dps = 80 settles both: LikeR was wrong on pf (fixed by passing
the exact complement into incbeta_xy(); see LikeR.xs), and *R* is the wrong
one on qf, by up to 1.4e-3.  Both verdicts are frozen below so that changing
either is a deliberate act rather than an accident.

METHOD
------
I_x(a, b) comes from the Numerical-Recipes continued fraction evaluated by
Lentz's method at mp.dps = 80, taking the branch on x < (a+1)/(a+b+2) and
receiving x and y = 1-x as separate arguments so no cancellation can occur.
Quantiles come from bisection on that, never from a library inverse.

That is the same continued fraction LikeR uses, so what this pins is
*precision*, not the choice of algorithm -- which is exactly the point: the
bug being regression-tested was a cancellation in forming 1-x, not an error
in the continued fraction.  To keep the algorithm honest as well, every value
is cross-checked here against a genuinely different form,

    I_x(a, b) = x^a (1-x)^b / (a B(a, b)) * 2F1(a+b, 1; a+1; x)

wherever mpmath's hyp2f1 converges (it does not for large a+b with x away
from 0, which is why it cannot be the primary route).  The two agree to
between 75 and 80 digits on every case where both are available; the script
aborts if any pair disagrees by more than 1e-60 relative.

Quadrature of the defining integral was tried first and rejected: for
a = 50, b = 5e5 the beta density is a spike narrower than tanh-sinh can
resolve, and mp.quad returned 0.51864617 where the two forms above and both
implementations agree on 0.51880550.
"""

import sys

from mpmath import mp, mpf, hyp2f1, beta as Beta, exp, log, loggamma

mp.dps = 80

XCHECK_TOL = mpf(10) ** -60      # 2F1-vs-CF agreement demanded, when both work


def _cf(a, b, x):
    """Lentz evaluation of the NR incomplete-beta continued fraction."""
    tiny = mpf(10) ** (-mp.dps - 20)
    eps = mpf(10) ** (-mp.dps - 10)
    qab, qap, qam = a + b, a + 1, a - 1
    c = mpf(1)
    d = 1 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1 / d
    h = d
    for m in range(1, 400000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < eps:
            return h
    raise RuntimeError("incomplete-beta continued fraction did not converge")


def _lbeta(a, b):
    return loggamma(a) + loggamma(b) - loggamma(a + b)


def _ibeta_2f1(a, b, x):
    """The independent cross-check form.  Raises when it will not converge."""
    return x ** a * (1 - x) ** b / (a * Beta(a, b)) * hyp2f1(a + b, 1, a + 1, x,
                                                            maxterms=10 ** 7)


def ibeta_xy(a, b, x, y):
    """I_x(a, b) with the complement y = 1-x supplied, never re-formed."""
    a, b, x, y = mpf(a), mpf(b), mpf(x), mpf(y)
    if x <= 0:
        return mpf(0)
    if y <= 0:
        return mpf(1)
    if x < (a + 1) / (a + b + 2):
        v = exp(a * log(x) + b * log(y) - log(a) - _lbeta(a, b)) * _cf(a, b, x)
    else:
        v = 1 - exp(b * log(y) + a * log(x) - log(b) - _lbeta(a, b)) * _cf(b, a, y)
    try:
        r = _ibeta_2f1(a, b, x)
    except Exception:
        return v                                  # 2F1 unavailable here
    if r != 0 and abs(v - r) / abs(r) > XCHECK_TOL:
        sys.exit(f"CF and 2F1 disagree at a={a} b={b} x={x}: {v} vs {r}")
    return v


def pf_tails(f, d1, d2):
    """(lower, upper) for the F CDF.  x = d1*f/D and y = d2/D are exact
    complements over the one denominator D, which is the whole point."""
    f, d1, d2 = mpf(f), mpf(d1), mpf(d2)
    D = d1 * f + d2
    x, y = d1 * f / D, d2 / D
    lo = ibeta_xy(d1 / 2, d2 / 2, x, y)
    up = ibeta_xy(d2 / 2, d1 / 2, y, x)
    return lo, up


def pt_tails(t, df):
    """(lower, upper) for Student's t.  P(|T| > |t|) = I_z(df/2, 1/2) with
    z = df/(df + t^2) and exact complement t^2/(df + t^2)."""
    t, df = mpf(t), mpf(df)
    tt = t * t
    D = df + tt
    two = ibeta_xy(df / 2, mpf(1) / 2, df / D, tt / D)     # P(|T| > |t|)
    upper = two / 2 if t > 0 else 1 - two / 2
    lower = 1 - two / 2 if t > 0 else two / 2
    return lower, upper


def qf_lower(p, d1, d2):
    """F quantile at lower-tail p, by bisection on pf_tails -- not by any
    library inverse."""
    p = mpf(p)
    lo, hi = mpf("1e-80"), mpf(1)
    while pf_tails(hi, d1, d2)[0] < p:
        hi *= 2
    for _ in range(280):
        mid = (lo + hi) / 2
        if pf_tails(mid, d1, d2)[0] < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def s(v):
    return mp.nstr(v, 21, strip_zeros=False)


PF_DF = [(1, 1), (1, 5), (1, 100), (1, 1e6), (2, 2), (2, 100),
         (5, 5), (5, 1e6), (20, 20), (100, 100)]
PF_Q = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1, 2, 10, 100]
PT_DF = [1, 2, 3, 5, 10, 100, 1e6]
PT_T = [-100, -10, -2, -1, -1e-2, -1e-4, -1e-6, -1e-8, -1e-10, -1e-12,
        1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 2, 10, 100]
QF_DF = [(1, 1), (1, 2), (1, 5), (1, 20), (1, 100), (2, 1), (2, 5),
         (2, 20), (2, 100), (5, 5), (5, 20), (20, 100)]
QF_P = [1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.5, 0.95, 0.999]

print("# Generated by t/pf_pt_tails.mpmath.py -- see that file to re-run.")
print("# mpmath mp.dps = 80; incomplete beta by Lentz continued fraction with")
print("# the complement passed in, cross-checked against the 2F1 form.")
print("#kind\targ\tdf1\tdf2\tlower\tupper")
for (d1, d2) in PF_DF:
    for q in PF_Q:
        lo, up = pf_tails(q, d1, d2)
        print(f"pf\t{q:.17g}\t{d1:g}\t{d2:g}\t{s(lo)}\t{s(up)}")
for df in PT_DF:
    for t in PT_T:
        lo, up = pt_tails(t, df)
        print(f"pt\t{t:.17g}\t{df:g}\t-\t{s(lo)}\t{s(up)}")
for (d1, d2) in QF_DF:
    for p in QF_P:
        print(f"qf\t{p:.17g}\t{d1:g}\t{d2:g}\t{s(qf_lower(p, d1, d2))}\t-")
