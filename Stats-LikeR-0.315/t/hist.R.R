# Regenerates the frozen table in t/hist.R.t.
#   Rscript t/hist.R.R
# The test never runs this; it reads the literals that were pasted from it.
#
# Written against R 4.6.1 (2026-06-24) -- graphics::hist.default, which is
#   breaks <- pretty(range(x), n = breaks, min.n = 1)
# with `breaks` defaulting to nclass.Sturges(x) = ceiling(log2(n) + 1), then
# counts from C_BinCount with right = TRUE and include.lowest = TRUE, mids the
# midpoints of consecutive breaks, and density = counts / (n * diff(breaks)).
# pretty() itself is R_pretty() in src/appl/pretty.c.
#
# Every dataset is written out in full rather than generated, so the test reads
# the same numbers R saw.  The values are exact in binary at every NV width
# except n100/n1000, which are 6-decimal-place randoms -- those two are here
# for the bin *counts* and are compared at a tolerance for that reason.
options(digits = 17, scipen = 500)

sets <- list(
  int13  = c(1,2,2,3,4,7,9,10,11,15,15,15,20),
  tiny   = c(0.001,0.002,0.0035,0.004),
  neg    = c(-5,-3,-1,0,2,4,6,-7,8,9),
  big    = c(1e6,2e6,3.5e6,4e6,9e6),
  const2 = c(3,3,3,3,3,3,7),
  flat   = rep(7, 5),
  single = 10,
  two    = c(1,2),
  wide   = c(-1000,1000,0,500,-250,750),
  frac   = c(0.1,0.2,0.25,0.7,0.9,1.1,1.4,1.9),
  cross0 = c(-0.4,-0.2,0,0.1,0.35,0.5),
  sturges8 = c(1:8),
  sturges13 = c(1:13)
)

for (nm in names(sets)) {
  x <- sets[[nm]]
  for (nb in c(0, 1, 2, 3, 5, 8, 12, 20)) {
    h <- if (nb == 0) hist(x, plot = FALSE) else hist(x, breaks = nb, plot = FALSE)
    cat(sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\n", nm, nb,
                paste(sprintf("%.17g", x),         collapse = ","),
                paste(sprintf("%.17g", h$breaks),  collapse = ","),
                paste(h$counts,                    collapse = ","),
                paste(sprintf("%.17g", h$mids),    collapse = ","),
                paste(sprintf("%.17g", h$density), collapse = ",")))
  }
}
