#!/usr/bin/env perl

require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use Test::LeakTrace;
use Stats::LikeR;

sub is_approx {
	my ($got, $exp, $msg, $tol) = @_;
	$tol = 1e-9 unless defined $tol;
	ok(abs($got - $exp) < $tol, $msg)
		or diag("got $got, expected $exp (tolerance $tol)");
}

# --- quartile cutpoints match numpy/pandas linear interpolation -------------
{
	my ($codes, $edges) = qcut([1 .. 10], 4);
	my @want_edge = (1, 3.25, 5.5, 7.75, 10);
	is(scalar @$edges, scalar @want_edge, 'four bins -> five edges');
	for my $i (0 .. $#want_edge) {
		is_approx($edges->[$i], $want_edge[$i], "edge $i matches pandas");
	}
	is_deeply($codes, [0, 0, 0, 1, 1, 2, 2, 3, 3, 3], 'quartile codes match pandas');
}

# --- equal-frequency counts -------------------------------------------------
{
	my $codes = qcut([1 .. 100], 4);
	my %n;
	$n{$_}++ for @$codes;
	is($n{$_}, 25, "bin $_ holds 25 of 100") for 0 .. 3;
}

# --- explicit probability vector (top-5% tranche) ---------------------------
{
	my ($codes, $edges) = qcut([1 .. 100], [0, 0.5, 0.95, 1]);
	my %n;
	$n{$_}++ for @$codes;
	is($n{0}, 50, 'lower-half tranche');
	is($n{1}, 45, 'mid tranche');
	is($n{2},  5, 'top-5% tranche');
}

# --- named and interval labels ----------------------------------------------
{
	my $lab = qcut([1 .. 10], 4, labels => [qw/Q1 Q2 Q3 Q4/]);
	is_deeply($lab, [qw/Q1 Q1 Q1 Q2 Q2 Q3 Q3 Q4 Q4 Q4/], 'named labels applied');

	my $iv = qcut([1 .. 10], 4, labels => 'interval');
	is($iv->[0], '[1, 3.25]',    'first interval is closed-closed');
	is($iv->[-1], '(7.75, 10]',   'last interval is open-closed');
}

# --- NA (undef) passes through ----------------------------------------------
{
	my $codes = qcut([1, 2, undef, 4, 5, 6, 7, 8, 9, 10], 4);
	ok(!defined $codes->[2], 'undef stays undef');
	is($codes->[0], 0, 'value before NA binned correctly');
	is($codes->[9], 3, 'value after NA binned correctly');
}

# --- duplicate edges: raise by default, drop on request ---------------------
{
	my @tied = ((0) x 8, 1, 2, 3, 4);
	my $err = !eval { qcut(\@tied, 4); 1 };
	ok($err && $@ =~ /not unique/, 'tied data raises on duplicate edges');

	my ($codes, $edges) = qcut(\@tied, 4, duplicates => 'drop');
	ok(scalar @$edges < 5, 'dropping merges duplicate edges');
}

# --- leak checks (assignments hoisted out for Devel::Cover) ------------------
{
	my @data  = map { $_ / 7 } 1 .. 500;
	my @tied  = ((0) x 50, 1, 2, 3, 4, 5);
	my @probs = (0, 0.25, 0.5, 0.75, 1);
	unless ($INC{'Devel/Cover.pm'}) {
		no_leaks_ok { eval { my $x = qcut(\@data, 10) } } 'qcut: no leaks (integer q)';
		no_leaks_ok { eval { my $x = qcut(\@data, \@probs) } } 'qcut: no leaks (probability vector)';
		no_leaks_ok { eval { my $x = qcut(\@data, 4, labels => [qw/a b c d/]) } } 'qcut: no leaks (labels)';
		no_leaks_ok { eval { my $x = qcut(\@tied, 4, duplicates => 'drop') } } 'qcut: no leaks (drop dups)';
		no_leaks_ok { eval { my $x = qcut(\@tied, 4) } } 'qcut: no leaks (croak path)';
	}
}

done_testing();
