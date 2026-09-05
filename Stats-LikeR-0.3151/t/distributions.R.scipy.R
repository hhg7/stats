# Generator for the frozen table in t/distributions.R.scipy.t.
#
# Re-run with:   Rscript t/distributions.R.scipy.R > /tmp/dist.tsv
# then paste the rows into the __DATA__ section of t/distributions.R.scipy.t.
# The test never invokes this script; it reads only the frozen copy.
#
# Every input is chosen dyadic -- exactly representable in binary -- so that the
# double R computes at is bit-identical to the long-double and __float128 the
# Perl side may parse it into. A decimal like 0.975 is three different numbers
# on those three builds, and in a distribution tail that difference shows up in
# the answer; see the "familiar" block at the end, which keeps a few of those on
# purpose and pays for them with a looser tolerance.
#
# Being dyadic is not enough on its own: the value has to survive the round trip
# through decimal text, and perl 5.10.1's string-to-number conversion does not
# parse 0.9999847412109375 exactly -- it lands one ulp away, and 1 - p amplifies
# that ulp to 7.3e-12 relative, which is 400 times the tolerance. So every
# argument is written as an exact dyadic RATIO, "M/2^K", which perl rebuilds
# with one integer division on any build and any version. The expected values
# are still decimal, where a one-ulp parse is 2e-16 against a 1e-13 tolerance.
#
# The smallest dyadic probability is 2^-20; the far tails are covered instead by
# the mpmath references from SciPy, whose own tolerance is loose enough not to
# care about the last bit.

options(digits = 17)

# Render an exactly-representable number as "M/2^K", the form the .t file
# rebuilds by integer division.  Anything that is not dyadic within 60 bits
# falls back to %.17g, and would then be at the mercy of the reader's atof --
# so nothing in the grids below is allowed to reach that branch.
# The 1100-bit cap is for the VALUE column, not the arguments: a binomial tail
# at p = 1/2 is an integer times a power of two, so pbinom(0, 1000, 0.5) is
# exactly 2^-1000 -- and perl 5.10.1's atof turns its decimal, 9.3326...e-302,
# into 8.999999999999996e-302, a 3.7% error. As a ratio it is exact everywhere.
# A value that is not dyadic (nearly all of them) falls back to %.17g, which is
# fine: those sit in the normal range where a one-ulp parse is 2e-16.
dyad <- function(x) {
	if (!is.finite(x)) return(sprintf("%.17g", x))
	k <- 0; y <- x
	while (y != floor(y) && k < 1100) { y <- y * 2; k <- k + 1 }
	if (y != floor(y) || abs(y) >= 2^53) return(sprintf("%.17g", x))
	sprintf("%.0f/2^%d", y, k)
}

emit <- function(fn, args, value) {
	cat(fn, paste(sapply(args, dyad), collapse = ","), dyad(value), sep = "\t")
	cat("\n")
}

dyadic_p  <- c(2^-20, 2^-16, 2^-10, 2^-5, 2^-4, 2^-3, 2^-2, 0.5,
               0.75, 0.875, 0.96875, 1 - 2^-10, 1 - 2^-16)
dyadic_q  <- c(0.0625, 0.5, 1, 1.5, 2.5, 3.5, 5, 10, 20, 100, 1024)
dfs       <- c(1, 2, 3, 5, 10, 30, 100, 4.5, 12.5)
f_dfs     <- list(c(1, 1), c(1, 10), c(2, 5), c(5, 5), c(10, 100), c(4.5, 12.5))

# qnorm: both tails, and one shifted/scaled case per probability
for (p in dyadic_p) {
	emit("qnorm", c(p, 0, 1, 1), qnorm(p))
	emit("qnorm", c(p, 0, 1, 0), qnorm(p, lower.tail = FALSE))
	emit("qnorm", c(p, 10, 0.5, 1), qnorm(p, mean = 10, sd = 0.5))
}

# pt: both tails, symmetric grid of q
for (df in dfs) for (q in c(-rev(dyadic_q), 0, dyadic_q)) {
	emit("pt", c(q, df, 1), pt(q, df))
	emit("pt", c(q, df, 0), pt(q, df, lower.tail = FALSE))
}

# qt: both tails
for (df in dfs) for (p in dyadic_p) {
	emit("qt", c(p, df, 1), qt(p, df))
	emit("qt", c(p, df, 0), qt(p, df, lower.tail = FALSE))
}

# pchisq / qchisq: both tails.  q = 0 and q = 2^-20 pin the lower tail where
# 1 - upper would have nothing left of it.
for (df in dfs) {
	for (q in c(0, 2^-20, dyadic_q)) {
		emit("pchisq", c(q, df, 1), pchisq(q, df))
		emit("pchisq", c(q, df, 0), pchisq(q, df, lower.tail = FALSE))
	}
	for (p in dyadic_p) {
		emit("qchisq", c(p, df, 1), qchisq(p, df))
		emit("qchisq", c(p, df, 0), qchisq(p, df, lower.tail = FALSE))
	}
}

# pf / qf: both tails
for (d in f_dfs) {
	for (q in c(0, 2^-20, dyadic_q)) {
		emit("pf", c(q, d[1], d[2], 1), pf(q, d[1], d[2]))
		emit("pf", c(q, d[1], d[2], 0), pf(q, d[1], d[2], lower.tail = FALSE))
	}
	for (p in dyadic_p) {
		emit("qf", c(p, d[1], d[2], 1), qf(p, d[1], d[2]))
		emit("qf", c(p, d[1], d[2], 0), qf(p, d[1], d[2], lower.tail = FALSE))
	}
}

# pbinom: both tails, over the whole support and past both ends, at dyadic
# probabilities plus the two degenerate ones.
for (sz in c(1, 2, 10, 55, 1000)) for (pr in c(0, 2^-10, 0.25, 0.5, 0.75, 1)) {
	for (k in unique(c(-1, 0, 1, floor(sz/2), sz - 1, sz, sz + 1))) {
		emit("pbinom", c(k, sz, pr, 1), pbinom(k, sz, pr))
		emit("pbinom", c(k, sz, pr, 0), pbinom(k, sz, pr, lower.tail = FALSE))
	}
}

# The familiar non-dyadic probabilities every applied script actually writes.
# 0.975 is a different number on a double, a long double and a __float128, so
# these are checked at a looser tolerance than the dyadic block above.
for (p in c(0.001, 0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99, 0.999)) {
	emit("qnorm_familiar", c(p), qnorm(p))
	for (df in c(1, 5, 10, 30)) {
		emit("qt_familiar",     c(p, df), qt(p, df))
		emit("qchisq_familiar", c(p, df), qchisq(p, df))
	}
	emit("qf_familiar", c(p, 2, 10), qf(p, 2, 10))
}
