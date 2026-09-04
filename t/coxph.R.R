# Generator for the frozen expected values in t/coxph.R.t.
#
# Re-run with:
#
#     Rscript t/coxph.R.R
#
# and paste the printed numbers into the @expect table in t/coxph.R.t. The
# test itself never invokes this script, or R at all -- it must pass on a
# machine with neither installed.
#
# Produced with R 4.6.1 and survival 3.8.9. Two corpora:
#
#  * `ovarian`, the dataset that package ships and that ?coxph uses. It has no
#    tied event times, so Efron and Breslow agree on it exactly -- itself
#    worth pinning.
#
#  * a 16-observation set built here, with deliberate ties at t = 5 and t = 8
#    so that Efron and Breslow genuinely differ. It is built rather than taken
#    from a package because survival ships nothing that has tied event times
#    *and* values that are all small integers, and a tie structure that
#    depended on rounding would be a different test at each NV width.

suppressMessages(library(survival))
options(digits = 17)

emit <- function(prefix, tm, st, covs) {
  for (tie in c("efron", "breslow")) {
    for (k in seq_along(covs)) {
      x <- covs[seq_len(k)]
      df <- data.frame(tm = tm, st = st, x)
      form <- as.formula(paste("Surv(tm, st) ~", paste(names(x), collapse = " + ")))
      f <- coxph(form, data = df, ties = tie)
      cf <- coef(f); se <- sqrt(diag(vcov(f)))
      cat(sprintf("%s.%s.%d.coef %s\n", prefix, tie, k,
                  paste(sprintf("%.17g", unname(cf)), collapse = " ")))
      cat(sprintf("%s.%s.%d.se %s\n", prefix, tie, k,
                  paste(sprintf("%.17g", unname(se)), collapse = " ")))
      cat(sprintf("%s.%s.%d.loglik %.17g\n",  prefix, tie, k, f$loglik[2]))
      cat(sprintf("%s.%s.%d.loglik0 %.17g\n", prefix, tie, k, f$loglik[1]))
    }
  }
}

# ---- corpus 1: survival's own `ovarian` ---------------------------------
ov <- survival::ovarian
cat("ovarian.n", nrow(ov), "\n")
cat("ovarian.time",   paste(ov$futime, collapse = ","), "\n")
cat("ovarian.status", paste(ov$fustat, collapse = ","), "\n")
cat("ovarian.age",    paste(ov$age,    collapse = ","), "\n")
cat("ovarian.rx",     paste(ov$rx,     collapse = ","), "\n")
emit("ovarian", ov$futime, ov$fustat, list(age = ov$age, rx = ov$rx))

# ---- corpus 2: deliberate ties, all values exactly representable --------
tm <- c(2, 3, 5, 5, 5, 6, 8, 8, 9, 11, 12, 14, 15, 15, 17, 20)
st <- c(1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0)
x1 <- c(1, 4, 2, 6, 3, 8, 5, 7, 2, 9, 4, 6, 8, 1, 3, 5)
x2 <- c(0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1)
cat("tied.tm", paste(tm, collapse = ","), "\n")
cat("tied.st", paste(st, collapse = ","), "\n")
cat("tied.x1", paste(x1, collapse = ","), "\n")
cat("tied.x2", paste(x2, collapse = ","), "\n")
emit("tied", tm, st, list(x1 = x1, x2 = x2))
