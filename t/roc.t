#!/usr/bin/env perl
require 5.010;
use warnings FATAL => 'all';
use Stats::LikeR;
use Test::Exception;
use Test::More;
use Test::LeakTrace 'no_leaks_ok';

sub is_approx {
	my ($got, $expected, $test_name, $epsilon) = @_;
	$epsilon = 1e-7 if not defined $epsilon;
	my $diff = abs($got - $expected);
	if ($diff <= $epsilon) { pass("$test_name: within $epsilon"); return 1; }
	fail($test_name);
	diag("         got: $got\n    expected: $expected; diff = $diff");
	return 0;
}

# Reference values from an independent R computation:
#   AUC via Mann-Whitney concordance (ties = 0.5), SE via the DeLong midrank
#   formula, CI = AUC +/- qnorm(0.975)*SE clamped to [0,1].
my @s = (0.1,0.4,0.35,0.8,0.7,0.6,0.2,0.9,0.55,0.55,0.3,0.75,0.45,0.85,0.5);
my @l = (0,  0,  1,   1,  1,  0,  0,  1,  0,   1,   0,  1,   1,   1,  0);

my $r = roc(\@s, \@l);
is_approx($r->{auc},       0.8482142857, 'AUC matches Mann-Whitney concordance');
is_approx($r->{auc_se},    0.1016521083, 'DeLong standard error');
is_approx($r->{auc_ci}[0], 0.6489798145, 'AUC CI lower');
is_approx($r->{auc_ci}[1], 1.0000000000, 'AUC CI upper (clamped to 1)');
is($r->{n_pos}, 8, 'positive count');
is($r->{n_neg}, 7, 'negative count');

# auc() convenience function agrees with roc()->{auc}
is_approx(auc(\@s, \@l), $r->{auc}, 'auc() scalar matches roc AUC');

# Youden's J optimal operating point
my $y = $r->{youden};
is_approx($y->{threshold},   0.70,   'Youden threshold');
is_approx($y->{sensitivity}, 0.625,  'Youden sensitivity');
is_approx($y->{specificity}, 1.0,    'Youden specificity');
is_approx($y->{j},           0.625,  "Youden's J");

# Curve is monotone: sensitivity non-decreasing, specificity non-increasing
# as the threshold falls (walking the returned points in order).
my @pts = @{ $r->{curve} };
ok(@pts >= 2, 'curve has operating points');
my $mono = 1;
for my $i (1 .. $#pts) {
	$mono = 0 if $pts[$i]{sensitivity} < $pts[$i-1]{sensitivity} - 1e-12;
	$mono = 0 if $pts[$i]{specificity} > $pts[$i-1]{specificity} + 1e-12;
}
ok($mono, 'curve is monotone (sens up, spec down as threshold falls)');
is_approx($pts[0]{sensitivity}, 0.0, 'curve starts at sensitivity 0');
is_approx($pts[0]{specificity}, 1.0, 'curve starts at specificity 1');

# Perfectly separable -> AUC 1; reversed labels -> AUC 0.
is_approx(auc([1,2,3,4],[0,0,1,1]),           1.0, 'separable data => AUC 1');
is_approx(auc([4,3,2,1],[0,0,1,1]),           0.0, 'anti-correlated => AUC 0');
# direction => '<' (lower marker = positive) flips it back to 1.
is_approx(auc([4,3,2,1],[0,0,1,1],direction=>'<'), 1.0, "direction '<' flips scoring");

# String labels via positive =>
is_approx(auc(['0.9','0.1','0.8','0.2'], ['case','ctrl','case','ctrl'], positive=>'case'),
          1.0, 'string labels with positive=>');

# error paths
dies_ok { roc(\@s) }                       'roc: missing labels dies';
dies_ok { roc(\@s, [1,0,1]) }              'roc: length mismatch dies';
dies_ok { roc([1,2,3],[1,1,1]) }           'roc: single-class dies';
dies_ok { roc(\@s, \@l, bogus=>1) }        'roc: unknown option dies';

unless ($INC{'Devel/Cover.pm'}) {
	no_leaks_ok { eval { roc(\@s, \@l) } } 'roc: no leaks';
	no_leaks_ok { eval { auc(\@s, \@l) } } 'auc: no leaks';
}

done_testing();
