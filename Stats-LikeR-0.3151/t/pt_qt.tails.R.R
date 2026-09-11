# Regenerates the frozen table in t/pt_qt.tails.R.t.
#   Rscript t/pt_qt.tails.R.R
# The test never runs this; it reads the literals that were pasted from it.
#
# FLOOR = 1e-290 is not cosmetic, and it is why the far-tail pt() rows stop
# where they do.  Two separate things break below it:
#
#  * perl 5.10.1 cannot parse the literal.  Its numeric parser is accurate to
#    the last ulp or two down to about 1e-290 and then collapses:
#    3.1830988618377598e-295 comes back as 2.9999999999999982e-295, an error
#    of 5.75e-02, and anything past 1e-310 comes back as 0.  A frozen table is
#    decimal literals, so an expected value below the floor is simply not the
#    number the test compares against on that perl.  (5.44.0 parses all of
#    them exactly; the mantissa-1 inputs like 1e-300 survive 5.10.1 to 6.4e-16
#    and so stay in.)
#
#  * below DBL_MIN = 2.2e-308 a double is subnormal and holds only a few
#    significant bits, so R's answer there is not the same number a
#    long-double or __float128 build computes.  pt(-1e160, 2) is
#    4.999944335913415e-321 in R and 5e-321 on the wider builds -- a relative
#    1.1e-05 that is entirely the double's missing mantissa, and no tolerance
#    that admits it would still be testing anything.
#
# The floor keeps both out.  It costs no coverage of the bug this pins: the
# collapse being tested set in at |t| > 1.3407807929942597e154, and rows from
# there down to 1e-290 exercise it several times over.
options(digits = 17, scipen = 500)
FLOOR <- 1e-290
cat("# pt: t, df, expected\n")
for (df in c(0.5, 1, 2, 5, 30, 100)) {
  for (t in c(-1, -10, -1e3, -1e10, -1e50, -1e100, -1e150, -1e154,
              -1.3407807929942597e154, -1e155, -1e160, -1e200, -1e250,
              -1e280, -1e300)) {
    v <- pt(t, df)
    if (v >= FLOOR) cat(sprintf("pt\t%.17g\t%.17g\t%.17g\n", t, df, v))
  }
}
cat("# qt: p, df, expected\n")
for (df in c(0.1, 0.5, 1, 2, 5, 10, 30)) {
  for (p in c(1e-300, 1e-200, 1e-160, 1e-155, 1e-100, 1e-50, 1e-16,
              1e-8, 0.001, 0.025, 0.25)) {
    v <- qt(p, df)
    if (is.finite(v) && abs(v) >= FLOOR) cat(sprintf("qt\t%.17g\t%.17g\t%.17g\n", p, df, v))
  }
}
