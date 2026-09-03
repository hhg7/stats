#!/usr/bin/env python3
"""Third opinion for the rows of t/distributions.R.scipy.t where Stats::LikeR
and R disagree by more than the block tolerance.

CLAUDE.md's rule for a disagreement is to ask mpmath at mp.dps = 60, solving the
*defining* equation by bisection rather than calling a library inverse -- so
that the answer does not come from a third implementation of the same inverse
that might share a weakness. That is what this does: every quantile below is
found by bisecting its own CDF identity, and every CDF is mpmath's
regularized incomplete beta, evaluated at 60 digits.

Re-run with:  python3 t/distributions.mpmath.py
It prints, for each case, the 60-digit truth and the relative error of both R's
value and this module's. The .t file quotes the conclusions; nothing in the test
suite calls this script.
"""

from mpmath import mp, mpf, betainc, findroot

mp.dps = 60


def pf_lower(f, d1, d2):
    """P(F <= f) for F(d1, d2) = I_x(d1/2, d2/2), x = d1 f / (d1 f + d2)."""
    f, d1, d2 = mpf(f), mpf(d1), mpf(d2)
    x = d1 * f / (d1 * f + d2)
    return betainc(d1 / 2, d2 / 2, 0, x, regularized=True)


def pf_upper(f, d1, d2):
    """P(F > f), evaluated as its own integral rather than 1 - pf_lower."""
    f, d1, d2 = mpf(f), mpf(d1), mpf(d2)
    x = d1 * f / (d1 * f + d2)
    return betainc(d1 / 2, d2 / 2, x, 1, regularized=True)


def qf(p, d1, d2, lower):
    """Bisect the defining equation for f.  Bracket by doubling in log space,
    then bisect to 55 significant digits -- far past anything a double can
    hold, so the printed value is the truth to every digit the .t file uses."""
    p = mpf(p)
    cdf = (lambda f: pf_lower(f, d1, d2)) if lower else (lambda f: pf_upper(f, d1, d2))
    # F is increasing in f for the lower tail and decreasing for the upper.
    lo, hi = mpf(2) ** -400, mpf(2) ** 400
    for _ in range(4000):
        mid = mp.sqrt(lo * hi) if lo > 0 else hi / 2
        v = cdf(mid)
        if (v < p) if lower else (v > p):
            lo = mid
        else:
            hi = mid
        if hi - lo <= mpf(10) ** -55 * hi:
            break
    return mp.sqrt(lo * hi)


# The rows where the two disagree.  Each is (p, df1, df2, lower, R, LikeR),
# copied from the failing diagnostics.
CASES = [
    # (fn, p_or_q, df1, df2, lower, R's value, this module's value)
    ("qf", 9.5367431640625e-07,  1, 1,  1, "2.2442048219772914e-12", "2.2440882278497347e-12"),
    ("qf", 1.52587890625e-05,    1, 1,  1, "5.7448668044912665e-10", "5.7448658654869556e-10"),
    ("qf", 0.0009765625,         1, 1,  1, "2.3531007489197009e-06", "2.3531007489843549e-06"),
    ("qf", 0.9990234375,         1, 1,  0, "2.3531007489197009e-06", "2.3531007489843549e-06"),
    ("qf", 0.9999847412109375,   1, 1,  0, "5.7448668044912665e-10", "5.7448658654869556e-10"),
    ("qf", 9.5367431640625e-07,  1, 10, 1, "1.503241975342462e-12",  "1.5017547184536e-12"),
    ("qf", 1.52587890625e-05,    1, 10, 1, "3.8445024941324846e-10", "3.8444920797810381e-10"),
    ("qf", 0.0009765625,         1, 10, 1, "1.5747048642822392e-06", "1.5747048648775778e-06"),
    ("qf", 0.9990234375,         1, 10, 0, "1.5747048642822392e-06", "1.5747048648775778e-06"),
    ("qf", 0.9999847412109375,   1, 10, 0, "3.8445024941324846e-10", "3.8444920797810381e-10"),
    ("qf", 9.5367431640625e-07,  2, 5,  1, "9.5367495289711002e-07", "9.5367495305302613e-07"),
    ("qf", 1.52587890625e-05,    2, 5,  1, "1.5258952045793528e-05", "1.5258952045940082e-05"),
    ("qf", 0.9999847412109375,   2, 5,  0, "1.5258952045793528e-05", "1.5258952045940082e-05"),
    ("pf", 9.5367431640625e-07,  1, 10, 0, "0.999240022820441",      "0.99924002282056479"),
]


def rel(value, truth):
    v = mpf(value)
    return abs(v - truth) / abs(truth) if truth != 0 else abs(v)


print(f"{'fn':>3} {'p':>22} {'df1':>5} {'df2':>5} {'lo':>3}  "
      f"{'rel(R)':>10} {'rel(LikeR)':>10}   winner   truth")
worst_liker = mpf(0)
for fn, p, d1, d2, lower, r_val, l_val in CASES:
    truth = (qf(p, d1, d2, bool(lower)) if fn == "qf"
             else (pf_lower if lower else pf_upper)(p, d1, d2))
    rr, rl = rel(r_val, truth), rel(l_val, truth)
    worst_liker = max(worst_liker, rl)
    winner = "LikeR" if rl < rr else ("R" if rr < rl else "tie")
    print(f"{fn:>3} {p:>22} {d1:>5} {d2:>5} {lower:>3}  "
          f"{mp.nstr(rr, 4):>10} {mp.nstr(rl, 4):>10}   {winner:>6}   "
          f"{mp.nstr(truth, 20)}")

print()
print("worst relative error, Stats::LikeR, over these cases:",
      mp.nstr(worst_liker, 4))

# pbinom(0, 1000, 0.5) is exactly 2^-1000, so it needs no bisection at all.
exact = mpf(2) ** -1000
print()
print("pbinom(0, 1000, 0.5) is exactly 2^-1000 =", mp.nstr(exact, 22))
for who, v in (("R", "9.3326361850325597e-302"),
               ("LikeR", "9.3326361850335731e-302")):
    print(f"  rel({who}) = {mp.nstr(rel(v, exact), 4)}")
