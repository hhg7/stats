# Generator for the frozen expected values in t/var_sd_cov.R.t.
#
# Re-run it from the distribution root whenever that corpus changes:
#
#     Rscript t/var_sd_cov.R.R
#
# and paste the two blocks it prints over the __DATA__ section of the .t file.
# The test never invokes this script: t/var_sd_cov.R.t must pass on a machine
# with no R installed, so the numbers live in the .t file as literals.
#
# Written against R 4.6.1 (2026-06-24) -- stats::var, stats::sd, stats::cov,
# stats::cor, base::sum, base::min, base::max.  R's variance and covariance are
# src/library/stats/src/cov.c: the MEAN macro takes a first mean and refines it
# with a second pass (tmp = tmp + sum(x - tmp)/n), and the COV macros then sum
# the centred products about the refined mean.  That is the algorithm the XS
# side reproduces, so these are the numbers it has to match.
#
# The corpus is not written out.  Every vector here is defined by exact integer
# arithmetic over a power-of-two denominator, so t/var_sd_cov.R.t rebuilds the
# identical vectors in perl instead of carrying 5000 literals; the n/sum/min/max
# columns below are what pin the two constructions to each other.  Dyadic values
# also mean the vectors are the same numbers at every NV width -- a corpus of
# tenths would leave R computing on the double nearest 0.1 and a quadmath perl
# on a closer one, and that disagreement would read as a bug in the variance
# rather than in the test data.

# %.40g, not format(digits = 17).
#
# 17 significant digits round-trip a double back to the same double, which is
# all that is needed to reproduce R's value on a double-NV perl.  It is not
# enough for a long-double or __float128 perl: it reads the decimal to its own
# width, and 100000000004.83398 is not the number R had -- the exact sum is
# 100000000004.833984375, and the 17-digit form is short by 4.4e-17 relative.
# That was enough to fail the exact columns on both wide builds.  %.40g prints
# the double's exact decimal expansion, which every NV width then reads back as
# the identical number.
options(digits = 17)
f <- function(v) sprintf("%.40g", v)

cases <- list(
	int10       = 1:10,
	int2        = 1:2,
	eighths     = (0:63) / 8,
	sixteenths  = ((0:99) %% 17) / 16 - 0.5,
	# Mean 1e9 with a spread of 0.05: the case that separates a two-pass
	# variance from the naive sum-of-squares, and the reason R refines its
	# mean at all.  1e9 is itself exactly representable.
	offset_1e9  = 1e9 + (0:99) / 1024,
	offset_1e6  = 1e6 + ((0:499) %% 29) / 4,
	identical   = rep(0.25, 50),
	two         = c(1.5, -2.25),
	negpos      = c(-8.5, 4.25, -2.125, 16, -32.75, 0.5),
	# Long enough that the scan loop runs well past its seeding element.
	long_dyadic = ((0:4999) %% 8191) / 512 - 8
)

# The partner vector for cov()/cor().  x reversed would be a linear function of
# x for every arithmetic sequence here, pinning cor() at exactly -1 and testing
# nothing; adding a short repeating dyadic sawtooth breaks the collinearity
# while keeping every value exact.
partner <- function(x) {
	n <- length(x)
	x * 0.5 + ((seq_len(n) - 1) %% 13) / 4 - 1.5
}

cat("# name\tn\tsum\tmin\tmax\tmean_naive\tvar\tsd\n")
for (nm in names(cases)) {
	x <- as.numeric(cases[[nm]])
	# base::mean() refines its mean with a second pass, exactly as var() does.
	# Stats::LikeR's mean() is the single-pass sum/n, so the column it is
	# checked against is sum(x)/length(x) and not mean(x); see the note in
	# t/var_sd_cov.R.t.
	cat(sprintf("STAT\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", nm, length(x),
		f(sum(x)), f(min(x)), f(max(x)), f(sum(x) / length(x)),
		if (length(x) > 1) f(var(x)) else "NA",
		if (length(x) > 1) f(sd(x))  else "NA"))
}

cat("# name\tn\tcov\tcor\n")
for (nm in names(cases)) {
	x <- as.numeric(cases[[nm]])
	if (length(x) < 2) next
	y <- partner(x)
	# cor() is undefined when either column has zero variance; R warns and
	# returns NA there, and Stats::LikeR croaks, so that case is checked
	# separately in the .t file rather than tabulated here.
	cat(sprintf("PAIR\t%s\t%d\t%s\t%s\n", nm, length(x), f(cov(x, y)),
		if (sd(x) > 0) f(cor(x, y)) else "NA"))
}
