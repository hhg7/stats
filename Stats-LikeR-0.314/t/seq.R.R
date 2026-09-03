# Generator for the frozen expected values in t/seq.R.t.
#
#   Rscript t/seq.R.R
#
# Prints the two Perl data tables that file carries; the test itself never
# runs this script and never calls R.  Re-run it only when a table has to be
# regenerated, and paste the output back into the .t file.
#
# Every case is base::seq() -- seq.default(), not the seq.int() primitive,
# because seq.default() is what R's seq() dispatches to and the two disagree
# (see the header of seq() in LikeR.xs).  The cases marked "R's own" are
# taken from R 4.6.1's regression suite and reference manual, cited
# individually in t/seq.R.t; the rest walk the branches of seq.default() that
# those do not reach.
#
# The second column of each row is Perl source for the argument list, kept
# textually as close to the R call as Perl allows, so that both languages
# parse the same decimal literals: R's `pi` is spelled 4*atan2(1,1) and an
# absent `by` is spelled undef.  A sequence of more than 200 values is
# emitted to the ENDS table -- length plus its first and last four values --
# rather than written out in full.
#
# Three tables come out, chosen per case by that third column:
#
#   full    every value, as %.17g decimals
#   ends    length, first four values, last four
#   dyadic  every value as an exact integer multiple of a power of two.
#           For a case at either end of the exponent range that is the only
#           portable spelling: perl-5.10.1 reads the decimal literal
#           7.9050503334599447e-323 as 0 and 1e307 five ulp low, so a
#           decimal table there would be testing perl's string-to-double
#           rather than seq().  The generator checks that R's values really
#           are exact multiples before emitting them this way.

options(scipen = 300)

# label (an R call, verbatim) | Perl argument source | how much to emit
cases <- list(
	# --- R's own: tests/reg-tests-1a.R, the "## seq" block
	c("seq(3,3,by=pi)",                 "3, 3, 4*atan2(1,1)",      "full"),
	c("seq(3,3.1,by=pi)",               "3, 3.1, 4*atan2(1,1)",    "full"),
	c("seq(1,6,by=3)",                  "1, 6, 3",                 "full"),
	c("seq(10,4.05,by=-3)",             "10, 4.05, -3",            "full"),
	# --- R's own: tests/reg-tests-1a.R, Don MacQueen 2002-03-26
	c("seq(1024902010,1024902025,by=1)","1024902010, 1024902025, 1","full"),
	# --- R's own: tests/reg-tests-1a.R, the "## Round" block's sample
	c("seq(-2,4,by=0.5)",               "-2, 4, 0.5",              "full"),
	# --- R's own: tests/reg-tests-1b.R, deliberate overshoot from fuzz
	c("seq(0,1,0.00025+5e-16)",         "0, 1, 0.00025+5e-16",     "ends"),
	# --- R's own: src/library/base/man/seq.Rd examples
	c("seq(1,9,by=2)",                  "1, 9, 2",                 "full"),
	c("seq(1,9,by=pi)",                 "1, 9, 4*atan2(1,1)",      "full"),
	c("seq(1.575,5.125,by=0.05)",       "1.575, 5.125, 0.05",      "full"),
	# --- from:to, i.e. 'by' absent.  Looser fuzz, both directions, no clamp.
	c("seq(1,5)",                       "1, 5, undef",             "full"),
	c("seq(5,1)",                       "5, 1, undef",             "full"),
	c("seq(1,1)",                       "1, 1, undef",             "full"),
	c("seq(2,-2)",                      "2, -2, undef",            "full"),
	c("seq(-1,-5)",                     "-1, -5, undef",           "full"),
	c("seq(1.5,4.5)",                   "1.5, 4.5, undef",         "full"),
	c("seq(4.5,1.5)",                   "4.5, 1.5, undef",         "full"),
	c("seq(-2.5,2.5)",                  "-2.5, 2.5, undef",        "full"),
	c("seq(1,4.99)",                    "1, 4.99, undef",          "full"),
	c("seq(1,4.9999999)",               "1, 4.9999999, undef",     "full"),
	# --- the same endpoints with by = 1: the tighter fuzz, two values fewer
	c("seq(1,4.9999999,by=1)",          "1, 4.9999999, 1",         "full"),
	c("seq(4.9999999,1,by=-1)",         "4.9999999, 1, -1",        "full"),
	c("seq(1,5,by=1)",                  "1, 5, 1",                 "full"),
	# --- fractional and negative steps
	c("seq(0,1,0.1)",                   "0, 1, 0.1",               "full"),
	c("seq(0,100,0.1)",                 "0, 100, 0.1",             "ends"),
	c("seq(1,2,0.25)",                  "1, 2, 0.25",              "full"),
	c("seq(10,5,-1)",                   "10, 5, -1",               "full"),
	c("seq(1,1000,0.5)",                "1, 1000, 0.5",            "ends"),
	c("seq(0.1,0.9,by=0.2)",            "0.1, 0.9, 0.2",           "full"),
	c("seq(-3,3,by=1.5)",               "-3, 3, 1.5",              "full"),
	c("seq(-2,4,by=1.5)",               "-2, 4, 1.5",              "full"),
	c("seq(-5,5,by=2.5)",               "-5, 5, 2.5",              "full"),
	c("seq(1,3,by=0.9999999999)",       "1, 3, 0.9999999999",      "full"),
	c("seq(1,3,by=0.99999999999)",      "1, 3, 0.99999999999",     "full"),
	# --- degenerate: one value out, by whichever of R's four routes
	c("seq(0,0,by=0)",                  "0, 0, 0",                 "full"),
	c("seq(0,0,by=5)",                  "0, 0, 5",                 "full"),
	c("seq(3,3,by=0)",                  "3, 3, 0",                 "full"),
	c("seq(3,3,by=1)",                  "3, 3, 1",                 "full"),
	c("seq(0.5,0.5,by=0.5)",            "0.5, 0.5, 0.5",           "full"),
	# --- endpoints indistinguishable at working precision -> just 'from'
	c("seq(1e15,1e15+20,by=2)",         "1e15, 1e15+20, 2",        "full"),
	c("seq(1e17,1e17+1,by=0.5)",        "1e17, 1e17+1, 0.5",       "full"),
	c("seq(2,2.00000000000001,by=1e-16)","2, 2.00000000000001, 1e-16","full"),
	# --- distinguishable by a hair: 200/1e15 clears 100 * DBL_EPSILON
	c("seq(1e15,1e15+200,by=2)",        "1e15, 1e15+200, 2",       "ends"),
	# --- at and past 2**53, where the IV fast path in seq() has to give up:
	# the first of these collapses (8/2**53 is under 100 * DBL_EPSILON), the
	# second does not and steps in the double's own spacing of 2
	c("seq(9007199254740992,9007199254741000,by=2)",
	                             "9007199254740992, 9007199254741000, 2", "full"),
	c("seq(9007199254740992,9007199254741992,by=2)",
	                             "9007199254740992, 9007199254741992, 2", "ends"),
	# --- subnormal step, and the extreme of the exponent range generally.
	# Powers of two, not decimal literals: perl-5.10.1's string-to-double
	# reads 1e-321 as 0 and 1e307 five ulp low, so a decimal here would be
	# testing perl's atof rather than seq().  2^k is exact on every build
	# and on every one of these perls.
	c("seq(0,2^-1066,by=2^-1070)",      "0, 2**-1066, 2**-1070",   "dyadic"),
	# --- (to - from) overflows the NV: R's 4*(from/4 + i*(by/4)) rewrite.
	# Only a double build takes that path -- 2^1024 is finite in a long
	# double or __float128 -- but every value here is a power of two, so
	# both routes produce the same numbers and the case is checkable on all
	# three widths.
	c("seq(-(2^1023),2^1023,by=2^1019)","-(2**1023), 2**1023, 2**1019", "dyadic"),
	c("seq(2^1023,-(2^1023),by=-(2^1019))",
	                                    "2**1023, -(2**1023), -(2**1019)", "dyadic")
)

# One value per comma-separated field, wrapped so no emitted line runs past
# 78 columns; continuation lines are indented to sit under the opening "[ ".
fmt <- function(v, indent = "\t    ") fmtw(sprintf("%.17g", v), indent)

fmtw <- function(tok, indent = "\t    ") {
	out <- character(0)
	line <- ""
	for (i in seq_along(tok)) {
		piece <- paste0(tok[i], if (i < length(tok)) "," else "")
		if (nchar(line) && nchar(line) + 1 + nchar(piece) > 70) {
			out <- c(out, line)
			line <- piece
		} else {
			line <- if (nchar(line)) paste(line, piece) else piece
		}
	}
	paste(c(out, line), collapse = paste0("\n", indent))
}

cat("# label | [ from, to, by (undef = absent) ] | every value, from R\n")
cat("my @FULL = (\n")
for (cs in cases) {
	if (cs[3] != "full") next
	v <- eval(parse(text = cs[1]))
	cat(sprintf("\t[ '%s', [ %s ],\n\t  [ %s ] ],\n", cs[1], cs[2], fmt(v)))
}
cat(");\n")

cat("# label | [ from, to, by ] | length | first four values | last four\n")
cat("my @ENDS = (\n")
for (cs in cases) {
	if (cs[3] != "ends") next
	v <- eval(parse(text = cs[1]))
	n <- length(v)
	cat(sprintf("\t[ '%s', [ %s ], %d,\n\t  [ %s ],\n\t  [ %s ] ],\n",
	            cs[1], cs[2], n, fmt(v[1:4]), fmt(v[(n - 3):n])))
}
cat(");\n")

cat("# label | [ from, to, by ] | exponent k | every value as an integer * 2**k\n")
cat("my @DYADIC = (\n")
for (cs in cases) {
	if (cs[3] != "dyadic") next
	v <- eval(parse(text = cs[1]))
	# the step is a power of two in every case routed here, so it is the
	# scale on which all of from, to and every element are whole numbers
	by <- eval(parse(text = sub("^.*by *= *", "", sub("\\)$", "", cs[1]))))
	k <- round(log2(abs(by)))
	m <- v / 2^k
	stopifnot(abs(by) == 2^k, all(m == round(m)), all(abs(m) < 2^53))
	cat(sprintf("\t[ '%s', [ %s ], %d,\n\t  [ %s ] ],\n",
	            cs[1], cs[2], k, fmtw(sprintf("%.0f", m))))
}
cat(");\n")
