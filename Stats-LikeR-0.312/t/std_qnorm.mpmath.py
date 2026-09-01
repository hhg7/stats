#!/usr/bin/env python3
"""Arbiter and generator for t/qnorm.crit.R.scipy.t.

Two jobs, both about the critical value z with P(-z < Z < z) = conf.level that
every Wald confidence limit in Stats::LikeR is built from:

  1. It prints the 60-digit truth for each conf.level in that test, evaluated
     at the *exact double* the C expression 1 - (1 - conf.level)/2 forms, so
     the frozen table pins the algorithm and not the rounding of conf.level.
  2. It scores R 4.6.1, SciPy 1.18.0 and Stats::LikeR against that truth, and
     against Moro's approximation alone -- which is what LikeR.xs's
     inverse_normal_cdf() returns, and what these nine call sites used through
     0.302.  That is the measurement the tolerance in the .t file is justified
     against.

CLAUDE.md's rule for a disagreement is to ask mpmath at mp.dps = 60 and to
solve the *defining* equation by bisection rather than call a library inverse,
so that the answer cannot come from a third implementation of the same inverse
sharing the same weakness.  The truth below is therefore found by bisecting

    erfc(z / sqrt(2)) / 2 = 1 - p

for z, at 60 digits, exactly as t/distributions.mpmath.py does for qf.

Re-run with:  python3 t/std_qnorm.mpmath.py
Nothing in the test suite calls this script; the .t file quotes its output.
"""

import math

from mpmath import erfc, mp, mpf, nstr, sqrt

mp.dps = 60

CONF_LEVELS = (0.80, 0.90, 0.95, 0.99, 0.999, 0.9999)

# Printed by:
#   Rscript -e 'options(digits=17); for (cl in c(0.80,0.90,0.95,0.99,0.999,0.9999))
#               cat(sprintf("%.17g\n", qnorm(1-(1-cl)/2)))'
R_4_6_1 = ("1.2815515655446008", "1.6448536269514715", "1.9599639845400536",
           "2.5758293035488999", "3.2905267314919255", "3.8905918864131199")

# Printed by scipy.stats.norm.ppf(1-(1-cl)/2), SciPy 1.18.0 / NumPy 2.5.2.
SCIPY_1_18_0 = ("1.2815515655446004", "1.6448536269514722", "1.959963984540054",
                "2.5758293035489004", "3.2905267314919255", "3.8905918864131204")


def truth(p_double):
    """z with P(Z <= p) = p, by bisecting erfc(-z/sqrt(2))/2 = p at 60 digits.

    mpf() of a python float is that double exactly, so what is inverted here is
    the number C actually holds and not the decimal it was written as."""
    p = mpf(p_double)
    lo, hi = mpf(-40), mpf(40)
    for _ in range(400):
        mid = (lo + hi) / 2
        if erfc(-mid / sqrt(2)) / 2 < p:
            lo = mid
        else:
            hi = mid
        if hi - lo <= mpf(10) ** -55 * max(abs(hi), 1):
            break
    return (lo + hi) / 2


def moro(p):
    """inverse_normal_cdf() from LikeR.xs, in double: Moro's approximation,
    the seed normal_quantile_hp() polishes and the whole answer these call
    sites used through 0.302."""
    a = (2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637)
    b = (-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833)
    c = (0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
         0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
         0.0000321767881768, 0.0000002888167364, 0.0000003960315187)
    y = p - 0.5
    if abs(y) < 0.42:
        r = y * y
        return (y * (((a[3] * r + a[2]) * r + a[1]) * r + a[0])
                / ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0))
    r = 1.0 - p if y > 0 else p
    r = math.log(-math.log(r))
    x = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * (
        c[5] + r * (c[6] + r * (c[7] + r * c[8])))))))
    return -x if y < 0 else x


def newton(p):
    """normal_quantile_hp(): Moro seeded, then Newton on approx_pnorm(), which
    is erfc and good to a few ulp at every NV width.  Double throughout, and
    NV_EPSILON is the double one, so this reproduces the default build."""
    z = moro(p)
    for _ in range(4):
        dens = 0.39894228040143267794 * math.exp(-0.5 * z * z)
        if not dens > 0.0:
            break
        step = (0.5 * math.erfc(-z * 0.70710678118654752440) - p) / dens
        z -= step
        if step * step <= 2.0 * 2.220446049250313e-16:
            break
    return z


def std_qnorm(p):
    """The reflection LikeR.xs applies above the median."""
    return -newton(1.0 - p) if p > 0.5 else newton(p)


def dyad(x):
    """x as the exact ratio M / 2^K the .t file freezes, so no atof anywhere
    has to round it back.  perl 5.10.1's atof does not always round correctly,
    which is why the table is written this way and not as a decimal."""
    m, e = math.frexp(x)
    m, e = int(m * 2 ** 53), 53 - e
    while m and m % 2 == 0:
        m //= 2
        e -= 1
    return f"{m}/2^{e}"


ULP = 2.220446049250313e-16


def rel(value, t):
    return abs(mpf(value) - t) / abs(t)


print("The frozen @CRIT table of t/qnorm.crit.R.scipy.t: conf.level, p, this\n"
      "module's z, Moro's z, R's z, SciPy's z -- every one as the exact ratio\n"
      "M / 2^K of the double that project holds.  p is frozen too because\n"
      "1 - (1 - conf.level)/2 rounds differently at each NV width; the .t file\n"
      "explains why that matters.\n")
for cl, r_v, s_v in zip(CONF_LEVELS, R_4_6_1, SCIPY_1_18_0):
    p = 1.0 - (1.0 - cl) / 2.0
    print(f"\t[ '{dyad(cl)}', '{dyad(p)}', '{dyad(std_qnorm(p))}', "
          f"'{dyad(moro(p))}', '{dyad(float(r_v))}', '{dyad(float(s_v))}' ],")

print(f"\n{'conf.level':>10} {'p as a double':>22} "
      f"{'Moro (<=0.302)':>15} {'R 4.6.1':>10} {'SciPy 1.18':>11} "
      f"{'LikeR 0.303':>12}")
print(f"{'':>10} {'relative error ->':>22} "
      f"{'':>15} {'ulp':>10} {'ulp':>11} {'ulp':>12}")
worst = {"moro": mpf(0), "r": mpf(0), "scipy": mpf(0), "liker": mpf(0)}
for cl, r_v, s_v in zip(CONF_LEVELS, R_4_6_1, SCIPY_1_18_0):
    p = 1.0 - (1.0 - cl) / 2.0
    t = truth(p)
    e = {"moro": rel(moro(p), t), "r": rel(r_v, t),
         "scipy": rel(s_v, t), "liker": rel(std_qnorm(p), t)}
    for k, v in e.items():
        worst[k] = max(worst[k], v)
    print(f"{cl:>10} {p:>22.17g} {nstr(e['moro'], 4):>15} "
          f"{float(e['r']) / ULP:>10.1f} {float(e['scipy']) / ULP:>11.1f} "
          f"{float(e['liker']) / ULP:>12.1f}")

print("\nworst relative error over these conf.levels:")
for k, label in (("moro", "Moro alone, i.e. <= 0.302"),
                 ("r", "R 4.6.1 qnorm"),
                 ("scipy", "SciPy 1.18.0 norm.ppf"),
                 ("liker", "Stats::LikeR 0.303 std_qnorm")):
    print(f"  {label:<30} {nstr(worst[k], 4):>10}"
          f"  ({float(worst[k]) / ULP:.1f} ulp)")
print(f"\nimprovement over 0.302: "
      f"{float(worst['moro'] / worst['liker']):.3g}x")

print("\nWhat the reflection buys -- same Newton polish, without inverting "
      "min(p, 1-p):")
for cl in CONF_LEVELS:
    p = 1.0 - (1.0 - cl) / 2.0
    t = truth(p)
    print(f"  conf.level {cl:<8} unreflected {float(rel(newton(p), t)) / ULP:>7.1f} ulp"
          f"   reflected {float(rel(std_qnorm(p), t)) / ULP:>5.1f} ulp")
p = 1.0 - 2.0 ** -16
t = truth(p)
print(f"  p = 1 - 2^-16     unreflected {nstr(rel(newton(p), t), 4):>9}"
      f"   reflected {nstr(rel(std_qnorm(p), t), 4):>9}")
