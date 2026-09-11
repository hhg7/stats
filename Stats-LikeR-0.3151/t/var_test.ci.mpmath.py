# Exact p-values for the table in t/var_test.ci.R.t, at mp.dps = 60.
#   python3 t/var_test.ci.mpmath.py
# The test never runs this; it reads the literals that were pasted from it.
#
# Why not R's p-values.  var.test() in R forms the upper tail as 1 - pf(F),
# which loses everything below the ulp of 1.  Stats::LikeR evaluates each tail
# where it lies (pf_upper()), so in the far tail the two genuinely differ and
# R is the less accurate side -- see the header of the .t file.  The variances
# are computed as exact rationals from the integer data, so the F ratio here
# is the true one and the only error being measured is the tail routine's.
from fractions import Fraction as F
from mpmath import mp, mpf, betainc
mp.dps = 60
# dyadic conf.levels; see the note in t/var_test.ci.R.R
CLS = [0.5, 0.75, 0.9375, 1 - 2**-8, 1 - 2**-12, 1 - 2**-20]
xs = [[1,2,3,4,5,6,7,8,9,10],
      [3,1,4,1,5,9,2,6,5,3,5],
      [2,4,6,8,10,12,14,16],
      [1,1,2,3,5,8,13,21,34]]
ys = [[2,4,1,8,3,9,5,7,6,10],
      [1,2,4,8,16,32,64],
      [5,5,6,6,7,7,8,8,9,9,10,10],
      [1,3,7,15,31,63,127]]
def var(v):
    n = len(v); m = sum(v, F(0))/n
    return sum((a-m)**2 for a in v)/(n-1)
def tails(Fstat, d1, d2):
    d1, d2 = mpf(d1), mpf(d2)
    x = d1*Fstat/(d1*Fstat + d2)
    lo = betainc(d1/2, d2/2, 0, x, regularized=True)
    up = betainc(d2/2, d1/2, 0, d2/(d1*Fstat + d2), regularized=True)
    return lo, up
for xi, x in enumerate(xs, 1):
    vx = var([F(a) for a in x]); d1 = len(x) - 1
    for yi, y0 in enumerate(ys, 1):
        for sc in (F(1), F(1,1024), F(1024)):
            y = [F(a)*sc for a in y0]
            vy = var(y); d2 = len(y) - 1
            r = vx/vy
            Fs = mpf(r.numerator)/mpf(r.denominator)
            lo, up = tails(Fs, d1, d2)
            two = min(2*min(lo, up), mpf(1))
            for ci in range(len(CLS)):
                for alt, p in (("two.sided", two), ("less", lo), ("greater", up)):
                    print("%d\t%d\t%.17g\t%d\t%s\t%s" %
                          (xi, yi, float(sc), ci, alt, mp.nstr(p, 20)))
