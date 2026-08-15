#!/usr/bin/env perl
# Fixture generator for the scaling benchmarks (scale.pl / scale.py / scale.R).
#
#     perl scale.data.pl            # writes into /tmp/likeR.scaling
#     SCALE_DIR=/somewhere perl scale.data.pl
#
# The three benchmark scripts must read *byte-identical* files, otherwise the
# read_table panels compare three different parsing jobs, so exactly one
# program writes them and the other two only read.  Run this first; scale.pl
# runs it for you if the directory is missing, and scale.py and scale.R stop
# with instructions rather than inventing their own copy.
#
# Three shapes are written at every row count in @rows:
#
#   num.<n>.csv   a,b,c,d              four numeric columns, nothing to quote
#   mix.<n>.csv   id,x,y,cat1,cat2     an integer, two numerics, two strings
#   mix.<n>.tsv   the same table, tab separated
#
# Numbers are printed with a fixed "%.6f", and no field ever contains the
# separator or a quote, so every reader sees the same characters and none of
# them is dragged into a quoting slow path the others avoid.
#
# The pseudo-random numbers come from an xorshift64 written out longhand rather
# than from perl's rand(), so the fixtures are reproducible on any perl and any
# platform: drand48 and the various rand() implementations do not agree, and a
# file whose contents depend on the build is not a fixture.
#
# Nothing here uses Stats::LikeR: the generator has to run before the module
# under test is involved at all.
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Path ();
use File::Spec ();

# Row counts.  These are the same list scale.pl, scale.py and scale.R use for
# their I/O panels; if you change one, change all four.
my @rows = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000);

my $dir = $ENV{SCALE_DIR} || '/tmp/likeR.scaling';

# ---------------------------------------------------------------------------
# xorshift64, and Box-Muller on top of it
# ---------------------------------------------------------------------------
# The state is kept as two 32-bit halves so this works unchanged on a perl
# without 64-bit integers, which is the one thing a plain "$s ^= $s << 13"
# would quietly get wrong.
my ($hi, $lo) = (0x12345678, 0x9ABCDEF0);

sub xorshift {
	# s ^= s << 13
	my $h = (($hi << 13) | ($lo >> 19)) & 0xFFFFFFFF;
	my $l = ($lo << 13) & 0xFFFFFFFF;
	($hi, $lo) = ($hi ^ $h, $lo ^ $l);
	# s ^= s >> 7
	$h = $hi >> 7;
	$l = (($lo >> 7) | ($hi << 25)) & 0xFFFFFFFF;
	($hi, $lo) = ($hi ^ $h, $lo ^ $l);
	# s ^= s << 17
	$h = (($hi << 17) | ($lo >> 15)) & 0xFFFFFFFF;
	$l = ($lo << 17) & 0xFFFFFFFF;
	($hi, $lo) = ($hi ^ $h, $lo ^ $l);
	return ($hi, $lo);
}

# a uniform in (0, 1), from the top 32 bits
sub unif {
	my ($h) = xorshift();
	return ($h + 0.5) / 4_294_967_296.0;
}

# one standard normal per call; the second Box-Muller value is cached
my $spare;
sub norm {
	if (defined $spare) {
		my $s = $spare;
		undef $spare;
		return $s;
	}
	my $u = unif();
	my $v = unif();
	my $r = sqrt(-2 * log($u));
	$spare = $r * sin(2 * 3.14159265358979323846 * $v);
	return $r * cos(2 * 3.14159265358979323846 * $v);
}

# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
File::Path::make_path($dir) unless -d $dir;

my @cat1 = qw(A B C);
my @cat2 = qw(X Y);

foreach my $n (@rows) {
	# --- num.<n>.csv: four numeric columns ---------------------------------
	my $num = File::Spec->catfile($dir, "num.$n.csv");
	open my $fh, '>', $num or die "cannot write \"$num\": $!\n";
	print {$fh} "a,b,c,d\n";
	for my $i (1 .. $n) {
		printf {$fh} "%.6f,%.6f,%.6f,%.6f\n",
			norm(), norm() * 2 + 5, unif() * 100, norm() * 0.5 - 3;
	}
	close $fh or die "cannot close \"$num\": $!\n";

	# --- mix.<n>.csv and mix.<n>.tsv: the same table, two separators -------
	# Both files are written from one pass so the two are the same table and
	# not two independent draws.
	my $csv = File::Spec->catfile($dir, "mix.$n.csv");
	my $tsv = File::Spec->catfile($dir, "mix.$n.tsv");
	open my $c, '>', $csv or die "cannot write \"$csv\": $!\n";
	open my $t, '>', $tsv or die "cannot write \"$tsv\": $!\n";
	print {$c} "id,x,y,cat1,cat2\n";
	print {$t} "id\tx\ty\tcat1\tcat2\n";
	for my $i (1 .. $n) {
		my $x  = norm();
		my $y  = norm() * 2 + 5;
		my $c1 = $cat1[ int(unif() * 3) ];
		my $c2 = $cat2[ int(unif() * 2) ];
		printf {$c} "%d,%.6f,%.6f,%s,%s\n",   $i, $x, $y, $c1, $c2;
		printf {$t} "%d\t%.6f\t%.6f\t%s\t%s\n", $i, $x, $y, $c1, $c2;
	}
	close $c or die "cannot close \"$csv\": $!\n";
	close $t or die "cannot close \"$tsv\": $!\n";

	printf "%8d rows: %s, %s, %s\n", $n, $num, $csv, $tsv;
}

printf "Fixtures written to %s\n", $dir;
