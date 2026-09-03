# Regenerate the frozen table in t/p_adjust.R.t:
#
#     Rscript t/p_adjust.R.R
#
# and paste the output into that file's data section.  Written against
# R 4.6.1 (2026-06-24).  The test reads only the frozen literals; installing or
# testing Stats::LikeR needs no R.
#
# Provenance: stats::p.adjust, driven over R's own p.adjust.methods vector so
# the method list can never drift out of step with R's, plus the cases R's own
# material covers:
#
#   src/library/stats/man/p.adjust.Rd -- the example applies every method in
#       p.adjust.methods to one p-vector, and its comment block records that
#       "The smallest 3 Bonferroni values are smaller than the 'BY' ones,
#       (John Maindonald, PR#17136)", which is reproduced.
#
#       That example's other assertion, p.adj <= p.adj.60, needs R's n=
#       argument to correct a family as if it were larger.  p_adjust() has no
#       n=, so there is nothing to generate for it; the .t file records the gap
#       rather than pretending to test it.
#   tests/reg-tests-1d.R:4339 -- p.adjust(<empty>, n = 0), PR#18002, and
#       "p.adjust(pp, 'holm') worked always but was not strictly tested".
#
# The sorted, reverse-sorted and all-tied vectors are here because holm,
# hochberg and hommel each walk the sorted order differently and only disagree
# when the order is degenerate; the n = 1 and the exact 0/1 vectors are the
# boundaries of every method's cumulative min/max.

options(digits = 17, warn = -1)

f <- function(k, v) cat(sprintf("%s\t%s\n", k,
	paste(sapply(v, function(z) sprintf("%.17g", z)), collapse = ",")))

SETS <- list(
	seq5     = c(0.01, 0.02, 0.03, 0.04, 0.05),          # already sorted
	rev5     = c(0.05, 0.04, 0.03, 0.02, 0.01),          # reverse sorted
	mix10    = c(0.6, 0.02, 0.001, 0.4, 0.31, 0.99, 0.05, 0.0001, 0.5, 0.12),
	tied4    = c(0.05, 0.05, 0.05, 0.05),                # every value tied
	one      = c(0.5),                                   # n = 1
	zeroone  = c(0, 1, 0.5),                             # the exact boundaries
	# p.adjust.Rd's own example vector, reproduced exactly:
	#     set.seed(123); x <- rnorm(50, mean = c(rep(0,25), rep(3,25)))
	#     p <- 2*pnorm(sort(-abs(x)))
	# It is what the page's printed table and both of its documented facts
	# refer to -- the PR#17136 Bonferroni/BY inversion, and the rejection
	# counts 11 11 11 11 20 12 22.
	manpage  = local({
		set.seed(123)
		x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
		2 * pnorm(sort(-abs(x)))
	})
)

cat(sprintf("# p.adjust.Rd example vector: %s\n",
            paste(sprintf("%.17g", SETS$manpage), collapse = ",")))
cat(sprintf("# p.adjust.Rd rejection counts at 0.05 (holm hochberg hommel bonferroni BH BY none): %s\n",
            paste(sapply(c("holm","hochberg","hommel","bonferroni","BH","BY","none"),
                         function(m) sum(p.adjust(SETS$manpage, m) < 0.05)),
                  collapse = ",")))

cat(sprintf("# R's p.adjust.methods: %s\n", paste(p.adjust.methods, collapse = ",")))
for (nm in names(SETS))
	for (m in p.adjust.methods)
		f(sprintf("%s.%s", nm, m), p.adjust(SETS[[nm]], m))
