# Regenerates the frozen table in t/cor_test.spearman.R.t.
#   Rscript t/cor_test.spearman.R.R
# The test never runs this; it reads the literals that were pasted from it.
#
# Data are small integers throughout, so the ranks -- and therefore the tie
# structure the exact branch depends on -- are identical at every NV width.
options(digits = 17, scipen = 500)
mk <- function(n, seed) { set.seed(seed); sample.int(n) }
cases <- list()
for (n in c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 20, 24, 32, 50, 100)) {
  for (s in 1:2) cases[[length(cases)+1]] <- list(n = n, y = mk(n, n*100+s))
}
# perfect and near-perfect, both signs, across the n = 9/10 branch boundary
for (n in c(8, 9, 10, 11, 24)) {
  cases[[length(cases)+1]] <- list(n = n, y = 1:n)
  cases[[length(cases)+1]] <- list(n = n, y = rev(1:n))
  yy <- 1:n; yy[1:2] <- yy[2:1]
  cases[[length(cases)+1]] <- list(n = n, y = yy)
}
# ties: forced onto the approximation, in R and here alike
cases[[length(cases)+1]] <- list(n = 10, y = c(1,1,2,2,3,3,4,4,5,5))
cases[[length(cases)+1]] <- list(n = 12, y = c(1,1,1,2,2,2,3,3,3,4,4,4))
for (cs in cases) {
  n <- cs$n; x <- 1:n; y <- cs$y
  for (alt in c("two.sided", "greater", "less")) {
    for (ex in c(NA, TRUE, FALSE)) {
      r <- suppressWarnings(
             if (is.na(ex)) cor.test(x, y, method = "spearman", alternative = alt)
             else cor.test(x, y, method = "spearman", alternative = alt, exact = ex))
      cat(sprintf("%d\t%s\t%s\t%s\t%.17g\t%.17g\n",
                  n, paste(y, collapse = ","), alt,
                  if (is.na(ex)) "d" else if (ex) "1" else "0",
                  r$statistic, r$p.value))
    }
  }
}
