# Regenerates the frozen table in t/cor_test.kendall.R.t.
#   Rscript t/cor_test.kendall.R.R
# The test never runs this; it reads the literals that were pasted from it.
#
# Data are small integers throughout, so the tie structure -- which is what
# decides between R's exact branch and its normal approximation -- is identical
# at every NV width.
options(digits = 17, scipen = 500)

cases <- list()
add <- function(name, x, y) cases[[length(cases)+1]] <<- list(name=name, x=x, y=y)

# tie-free: the exact branch, on both sides of R's n < 50 default
set.seed(11)
for (n in c(4, 5, 6, 7, 8, 9, 10, 12, 16, 24, 49, 50, 60)) {
  add(sprintf("perm%d", n), 1:n, sample.int(n))
}
add("perfect10",  1:10, 1:10)
add("reverse10",  1:10, 10:1)

# ties: the normal approximation, and the tie correction in var_S.  Group sizes
# 2 and 3 are both here because an odd group averages to a whole rank.
add("pairs10",   1:10, c(1,1,2,2,3,3,4,4,5,5))
add("triples12", 1:12, c(1,1,1,2,2,2,3,3,3,4,4,4))
add("tiedx12",   c(1,1,1,2,2,2,3,3,3,4,4,4), 1:12)
add("both15",    c(3,2,0,0,1,0,3,1,2,1,2,2,0,2,3), c(0,1,1,3,1,3,1,0,3,1,1,0,3,2,0))
add("binary20",  rep(c(0,1), each=10), c(0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1))
add("heavy30",   rep(1:3, each=10), rep(c(1,2,3,4,5), 6))
add("onetie",    1:8, c(1,2,3,4,5,6,7,7))

for (cs in cases) {
  for (alt in c("two.sided", "greater", "less")) {
    for (ex in c(NA, TRUE, FALSE)) {
      for (cc in c(FALSE, TRUE)) {
        r <- suppressWarnings(
               if (is.na(ex)) cor.test(cs$x, cs$y, method="kendall",
                                       alternative=alt, continuity=cc)
               else cor.test(cs$x, cs$y, method="kendall", alternative=alt,
                             exact=ex, continuity=cc))
        cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%.17g\t%.17g\t%.17g\n",
                    cs$name,
                    paste(cs$x, collapse=","), paste(cs$y, collapse=","),
                    alt, if (is.na(ex)) "d" else if (ex) "1" else "0",
                    if (cc) "1" else "0",
                    r$statistic, r$p.value, r$estimate))
      }
    }
  }
}

# The same datasets through cov() and cor(), which share the counting routine
# with cor_test() since 0.312.  R's cov(method="kendall") sums the sign
# products over the whole n x n space, so it is 2(C - D).
cat("--\n")
for (cs in cases) {
  cat(sprintf("%s\t%s\t%s\t%.17g\t%.17g\n",
              cs$name, paste(cs$x, collapse=","), paste(cs$y, collapse=","),
              cov(cs$x, cs$y, method="kendall"),
              cor(cs$x, cs$y, method="kendall")))
}
