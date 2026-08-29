#!/usr/bin/env perl
# min() and max() propagate NaN, as R's min()/max() do.
#
# Provenance: expected values are R 4.6.1 (2026-06-24), min(v) / max(v) over
# the vectors named in each test, run at options(digits = 17).  Every NaN case
# below is NaN in R for both functions; the finite and infinite cases are
# there so the fix cannot be mistaken for "always NaN".  No value here needs
# more precision than an integer, so there is nothing to freeze from a
# generator.
#
# What it pins.  Through 0.31 both functions folded with a bare
#     if (v < mn) mn = v;
# and every comparison against a NaN is false, so a NaN that was not the first
# element was skipped and one that was stayed put.  The answer depended on
# where in the array the NaN happened to sit:
#
#     min([NaN,1,2]) was NaN      min([1,2,NaN]) was 1      min([1,NaN,2]) was 1
#     max([NaN,1,2]) was NaN      max([1,2,NaN]) was 2
#
# R gives NaN for all of them, and sum(), mean(), var(), sd() and median() in
# this module already did.  The position sweep is the point of the file: a
# regression that only reads the first element would still pass a test that
# put the NaN at index 0.
#
# Both code paths are covered.  min()/max() run av_scan_min()/av_scan_max()
# over the leading run of plain numeric SVs and finish the rest through
# av_fetch(); an element that is a string forces the switch, so the arrays
# marked "split path" below carry the NaN on the av_fetch() side of it and the
# others carry it inside the scan.

require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use lib 'blib/lib', 'blib/arch';
use Stats::LikeR qw(min max);

my $INF = 9**9**9;
my $NAN = $INF - $INF;

sub is_nan { my $x = shift; return $x != $x }

# [name, arrayref]  -- NaN anywhere, at every position, in both paths
my @NAN_CASES = (
	['NaN first',              [$NAN, 1, 2, 4]],
	['NaN middle',             [1, 2, $NAN, 4]],
	['NaN last',               [1, 2, 4, $NAN]],
	['NaN alone',              [$NAN]],
	['NaN with +Inf',          [$NAN, $INF]],
	['NaN between infinities', [-$INF, $NAN, $INF]],
	['all NaN',                [($NAN) x 5]],
	['NaN early, long tail',   [5, $NAN, (1) x 50]],
	['NaN late, long head',    [(1) x 50, $NAN, 5]],
	['split path, NaN before the string', [1, 2, $NAN, '3', 4]],
	['split path, NaN after the string',  [1, 2, '3', $NAN, 4]],
);

# [name, arrayref, min, max] -- nothing here may become NaN
my @FINITE_CASES = (
	['plain',        [1, 2, 3],        1,     3],
	['with +Inf',    [1, $INF, 3],     1,     $INF],
	['with -Inf',    [1, -$INF, 3],   -$INF,  3],
	['both Inf',     [-$INF, 0, $INF], -$INF, $INF],
	['single',       [7],              7,     7],
	['descending',   [9, 5, 2, 1],     1,     9],
	['split path',   [4, 1, '2', 9],   1,     9],
);

plan tests => 2 * scalar(@NAN_CASES) + 2 * scalar(@FINITE_CASES) + 6;

for my $c (@NAN_CASES) {
	my ($name, $a) = @$c;
	ok(is_nan(min($a)), "min: $name is NaN")
		or diag("got " . min($a));
	ok(is_nan(max($a)), "max: $name is NaN")
		or diag("got " . max($a));
}

for my $c (@FINITE_CASES) {
	my ($name, $a, $mn, $mx) = @$c;
	is(min($a), $mn, "min: $name");
	is(max($a), $mx, "max: $name");
}

# The flat argument list and the multiple-arrayref forms take their own
# branches through the XSUB, so each needs its own NaN.
ok(is_nan(min(1, $NAN, 2)), 'min: NaN in a flat argument list');
ok(is_nan(max(1, $NAN, 2)), 'max: NaN in a flat argument list');
ok(is_nan(min([1, 2], [3, $NAN])), 'min: NaN in the second of two array refs');
ok(is_nan(max([1, 2], [3, $NAN])), 'max: NaN in the second of two array refs');
is(min([1, 2], [3, 4], 0.5), 0.5, 'min: mixed refs and scalars, no NaN');
is(max([1, 2], [3, 4], 9),   9,   'max: mixed refs and scalars, no NaN');
