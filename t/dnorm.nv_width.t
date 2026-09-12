#!/usr/bin/env perl
#
# dnorm(): the density is computed at the build's own NV width.
#
# c_dnorm() decides where the density has underflowed from the exponent range
# of the floating-point type, and up to 0.315 it asked <float.h> about a
# *double* -- DBL_MAX, DBL_MIN_EXP, DBL_MANT_DIG -- whatever perl's NV was.  On
# a long-double or __float128 build that cut the tail off at |x| ~ 38.57, where
# a double's subnormals run out, and returned a flat 0 beyond it: dnorm(-100) is
# 1.4e-2174, perfectly representable on a quadmath NV and four thousand orders
# of magnitude inside its range, and came back 0.  The file's own rule is that
# floating point is NV and the constants are the NV_* ones; these three were
# the exception.
#
# Provenance:
#
#   * R 4.6.1 stats::dnorm() at options(digits=17), for @R_DNORM, @R_DNORM_LOG
#     and @R_DNORM_ARGS below.  R computes in double, so those values pin the
#     double build exactly -- including the zeros past its underflow boundary,
#     which are correct there and which this module must keep producing.
#   * The width-dependent section takes no reference value.  It asserts the
#     identity dnorm(x) == exp(-x^2/2) / sqrt(2*pi), evaluated by perl's own
#     NV arithmetic, which is true at every width: on a double build both sides
#     underflow to 0 together, and on a wider one both sides are the same
#     non-zero number.  That is exactly the property the DBL_ constants broke,
#     and it needs no knowledge of which perl is running it.

require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use Config;
use Stats::LikeR qw(dnorm);

# R 4.6.1: [ x, dnorm(x) ]
my @R_DNORM = (
	[ 0,       0.3989422804014327      ],
	[ 1,       0.24197072451914337     ],
	[ 5,       1.4867195147342977e-06  ],
	[ 10,      7.6945986267064199e-23  ],
	[ 20,      5.5209483621597635e-88  ],
	[ 30,      1.4736461348785476e-196 ],
	[ 37,      2.1200065515246056e-298 ],
	[ 38,      1.0972210519949712e-314 ],   # subnormal on a double
	[ 38.5,    5.434722104253712e-323  ],
	[ 38.56,   4.9406564584124654e-324 ],   # the last subnormal
	[ 38.567,  4.9406564584124654e-324 ],
	[ -1,      0.24197072451914337     ],
	[ -38.5,   5.434722104253712e-323  ],
);
# R 4.6.1: [ x, dnorm(x, log = TRUE) ] -- log_p carries all the way out, at
# every width, and so is the same on every build.
my @R_DNORM_LOG = (
	[ -200, -20000.918938533203  ],
	[ -100, -5000.9189385332047  ],
	[ -40,  -800.91893853320471  ],
	[ -5,   -13.418938533204672  ],
	[ 0,    -0.91893853320467278 ],
	[ 5,    -13.418938533204672  ],
	[ 40,   -800.91893853320471  ],
	[ 100,  -5000.9189385332047  ],
	[ 200,  -20000.918938533203  ],
);
# R 4.6.1: [ x, mean, sd, dnorm(x, mean, sd) ]
my @R_DNORM_ARGS = (
	[ 3,    1,    2,     0.12098536225957168     ],
	[ 1e5,  1e5,  0.001, 398.9422804014327       ],
	[ 0,    0,    1,     0.3989422804014327      ],
	[ -7,   2,    0.5,   3.5174990851902079e-71  ],
);

# 8 ulp of a double.  The body of the density is one exp() and one multiply, so
# the error is a couple of ulp; the subnormal values at the far end have fewer
# bits than that and are compared exactly below instead.  Worst relative
# disagreement observed on a double build: 0.
my $TOL = 8 * 2.220446049250313e-16;

sub rel_ok {
	my ($got, $exp, $label) = @_;
	if ($exp == 0) { return is($got, 0, $label) }
	# A subnormal double carries only a handful of significant bits, so
	# anything at or below the smallest normal is compared for equality: there
	# is no relative tolerance to speak of down there.
	if (abs($exp) < 2.2250738585072014e-308) {
		return ok($got == $exp, $label) || diag("got $got, expected $exp");
	}
	return ok(abs($got - $exp) <= $TOL * abs($exp), $label)
		|| diag("got $got, expected $exp");
}

rel_ok(dnorm($_->[0]), $_->[1], "dnorm($_->[0]) matches R") for @R_DNORM;
rel_ok(dnorm($_->[0], 'log' => 1), $_->[1], "dnorm($_->[0], log) matches R")
	for @R_DNORM_LOG;
rel_ok(dnorm($_->[0], mean => $_->[1], sd => $_->[2]), $_->[3],
       "dnorm($_->[0], $_->[1], $_->[2]) matches R") for @R_DNORM_ARGS;

# --------------------------------------------------------------- the tail
#
# The identity, at whatever width this perl carries.  `exp(-x*x/2)` underflows
# to 0 on a double exactly where dnorm must, and does not on a wider NV exactly
# where dnorm must not, so one assertion covers every build in the matrix.
{
	my $inv_sqrt_2pi = 1 / sqrt(8 * atan2(1, 1));   # 1/sqrt(2*pi) at NV width
	for my $x (30, 38, 38.5, 38.6, 39, 45, 60, 100, 150) {
		my $want = exp(-$x * $x / 2) * $inv_sqrt_2pi;
		my $got  = dnorm(-$x);
		if ($want == 0) {
			is($got, 0, "dnorm(-$x) is 0, as exp(-x^2/2) is at this NV width");
		} else {
			ok($got != 0, "dnorm(-$x) is not flattened to 0 (exp(-x^2/2) = $want)")
				or next;
			ok(abs($got - $want) <= 1e-9 * $want,
			   "dnorm(-$x) agrees with exp(-x^2/2)/sqrt(2*pi)")
				or diag("got $got, wanted $want");
		}
	}
	# Say which side of the matrix this run exercised, so a smoker report is
	# readable: on a double NV the loop above checks zeros past x ~ 38.57, and
	# on a wider one it checks live values there instead.
	diag(sprintf 'NV is %s (nvsize %d): dnorm(-100) = %s',
	     $Config{nvtype}, $Config{nvsize}, dnorm(-100));
}

# The symmetric and degenerate cases, unchanged by any of this.
is(dnorm(3), dnorm(-3), 'dnorm is symmetric');
is(dnorm(0, mean => 0, sd => 0), 9**9**9, 'sd = 0 at the mean is Inf');
is(dnorm(1, mean => 0, sd => 0), 0,       'sd = 0 away from the mean is 0');
is(dnorm(1, mean => 0, sd => 9**9**9), 0, 'an infinite sd gives 0');

done_testing();
