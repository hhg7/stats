# Regenerates the frozen table in t/var_test.ci.R.t.
#   Rscript t/var_test.ci.R.R
# The test never runs this; it reads the literals that were pasted from it.
#
# Data are integers scaled by powers of two, so x, y and every variance are
# exactly representable at each NV width and the interval is the only thing
# that can move between builds.
options(digits = 17, scipen = 500)
# conf.level must be DYADIC.  It is read as an NV, so a decimal like 0.999999
# is a different number on a double build than on a long-double one, and
# beta = (1 - conf.level)/2 turns that 5.5e-17 gap into a relative 5.5e-11 in
# a quantity of size 1e-6 -- which then lands in the interval.  The frozen
# table would be R's double reading of the literal, and the wider builds would
# fail against it for a reason that has nothing to do with the code.  These
# six are exact at every NV width, and the test spells them as the same
# arithmetic rather than as decimals.
CLS <- c(0.5, 0.75, 0.9375, 1 - 2^-8, 1 - 2^-12, 1 - 2^-20)
xs <- list(c(1,2,3,4,5,6,7,8,9,10),
           c(3,1,4,1,5,9,2,6,5,3,5),
           c(2,4,6,8,10,12,14,16),
           c(1,1,2,3,5,8,13,21,34))
ys <- list(c(2,4,1,8,3,9,5,7,6,10),
           c(1,2,4,8,16,32,64),
           c(5,5,6,6,7,7,8,8,9,9,10,10),
           c(1,3,7,15,31,63,127))
for (xi in seq_along(xs)) for (yi in seq_along(ys)) {
  x <- xs[[xi]]
  for (sc in c(1, 1/1024, 1024)) {           # dyadic: exact at every NV width
    y <- ys[[yi]] * sc
    for (ci in seq_along(CLS)) {
      cl <- CLS[ci]
      for (alt in c("two.sided", "less", "greater")) {
        r <- var.test(x, y, conf.level = cl, alternative = alt)
        cat(sprintf("%d\t%d\t%.17g\t%d\t%s\t%.17g\t%.17g\t%.17g\t%.17g\n",
                    xi, yi, sc, ci, alt, r$estimate, r$p.value,
                    r$conf.int[1], r$conf.int[2]))
      }
    }
  }
}
