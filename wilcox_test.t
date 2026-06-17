use strict;
use warnings;
use Test::More;

# Stats::LikeR exports wilcox_test.
use Stats::LikeR qw(wilcox_test);

# Prefer Test::Exception; fall back to a tiny shim so this file is
# self-contained and runs even where Test::Exception is not installed.
BEGIN {
	if (eval { require Test::Exception; Test::Exception->import; 1 }) {
		# real Test::Exception in scope
	} else {
		no warnings 'redefine';
		*throws_ok = sub (&$;$) {
			my ($code, $re, $name) = @_;
			my $ok = !eval { $code->(); 1 };
			my $err = $@;
			$ok &&= ($err =~ $re);
			ok($ok, $name) or diag("got: " . (defined $err ? $err : '(no exception)'));
		};
		*dies_ok = sub (&;$) {
			my ($code, $name) = @_;
			ok(!eval { $code->(); 1 }, $name);
		};
		*lives_ok = sub (&;$) {
			my ($code, $name) = @_;
			my $ok = eval { $code->(); 1 };
			ok($ok, $name) or diag("died: $@");
		};
	}
}

# Compare two numbers within an absolute tolerance.
sub is_approx {
	my ($got, $exp, $tol, $name) = @_;
	$tol = 1e-6 unless defined $tol;
	my $ok = defined($got) && abs($got - $exp) <= $tol;
	ok($ok, $name)
		or diag(sprintf "got %.12g, expected %.12g (tol %.3g)",
			(defined $got ? $got : 'nan'), $exp, $tol);
	return $ok;
}

# Run wilcox_test capturing any warnings it emits.
sub wt_warns {
	my @args = @_;
	my @w;
	my $r = do {
		local $SIG{__WARN__} = sub { push @w, $_[0] };
		wilcox_test(@args);
	};
	return ($r, \@w);
}

# R man-page data set (see ?wilcox.test).
my @x = (1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30);
my @y = (0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29);

# ---------------------------------------------------------------------------
# Agreement with R's wilcox.test on known values
# ---------------------------------------------------------------------------

subtest 'two-sample rank sum (ties -> normal approximation)' => sub {
	my ($r, $w) = wt_warns(\@x, \@y);
	is_approx($r->{statistic}, 58,        1e-9,  'W statistic = 58');
	is_approx($r->{p_value},   0.1329189, 1e-6,  'p-value matches R (0.1329)');
	is($r->{method}, 'Wilcoxon rank sum test with continuity correction', 'method string');
	is($r->{alternative}, 'two.sided', 'alternative echoed');
	done_testing();
};

subtest 'two-sample exact (no ties, fully separated)' => sub {
	my $r = wilcox_test([1,2,3,4], [5,6,7,8]);
	is_approx($r->{statistic}, 0,          1e-9, 'W = 0');
	is_approx($r->{p_value},   0.02857143, 1e-7, 'exact two-sided p = 2/70');
	is($r->{method}, 'Wilcoxon rank sum exact test', 'exact method selected');
	done_testing();
};

subtest 'paired signed-rank, one-sided greater (exact)' => sub {
	my $r = wilcox_test(\@x, \@y, paired => 1, alternative => 'greater');
	is_approx($r->{statistic}, 40,         1e-9, 'V statistic = 40');
	is_approx($r->{p_value},   0.01953125, 1e-8, 'exact p = 10/512');
	is($r->{method}, 'Wilcoxon signed rank exact test', 'signed-rank exact method');
	is($r->{alternative}, 'greater', 'alternative echoed');
	done_testing();
};

subtest 'one-sample signed-rank against mu' => sub {
	# Shifting x down by mu must equal testing (x - mu) against 0.
	my $a = wilcox_test([1,2,3,4,5,6,7,8,9,10,11], mu => 6);
	my @shift = map { $_ - 6 } (1..11);
	my $b = wilcox_test(\@shift);
	is_approx($a->{statistic}, $b->{statistic}, 1e-9, 'mu shift == pre-shifted data (V)');
	is_approx($a->{p_value},   $b->{p_value},   1e-9, 'mu shift == pre-shifted data (p)');
	done_testing();
};

# ---------------------------------------------------------------------------
# Option handling
# ---------------------------------------------------------------------------

subtest 'continuity correction toggles the p-value' => sub {
	my $on  = wilcox_test(\@x, \@y, correct => 1);
	my $off = wilcox_test(\@x, \@y, correct => 0);
	ok(abs($on->{p_value} - $off->{p_value}) > 1e-9, 'correct=>0 differs from correct=>1');
	is($off->{method}, 'Wilcoxon rank sum test', 'method drops "continuity correction"');
	done_testing();
};

subtest 'exact can be forced on and off' => sub {
	my $forced_off = wilcox_test([1,2,3,4], [5,6,7,8], exact => 0);
	is($forced_off->{method}, 'Wilcoxon rank sum test with continuity correction',
		'exact=>0 forces the approximation');
	my ($forced_on, $w) = wt_warns([1,2,2,3], [4,5,5,6], exact => 1);
	ok(scalar(@$w) >= 1, 'exact=>1 with ties warns');
	like($w->[0], qr/ties/, 'warning mentions ties');
	is($forced_on->{method}, 'Wilcoxon rank sum test with continuity correction',
		'falls back to approximation despite exact=>1');
	done_testing();
};

subtest 'named x/y override positionals' => sub {
	my $pos   = wilcox_test(\@x, \@y);
	my $named = wilcox_test(x => \@x, y => \@y);
	is_approx($named->{statistic}, $pos->{statistic}, 1e-9, 'named form == positional form');
	done_testing();
};

subtest 'non-numeric and undef cells are dropped' => sub {
	my $clean = wilcox_test([1,2,3,4], [5,6,7,8]);
	my $dirty = wilcox_test([1,2,undef,3,'NA',4], [5,'x',6,7,8,undef]);
	is_approx($dirty->{statistic}, $clean->{statistic}, 1e-9, 'NA/undef removed, same W');
	is_approx($dirty->{p_value},   $clean->{p_value},   1e-9, 'NA/undef removed, same p');
	done_testing();
};

# ---------------------------------------------------------------------------
# Bug regressions
# ---------------------------------------------------------------------------

subtest 'zero variance (all values identical) does not divide by zero' => sub {
	my ($r, $w) = wt_warns([5,5,5], [5,5,5]);
	is_approx($r->{p_value}, 1.0, 1e-12, 'identical samples => p = 1 (not 0 or NaN)');
	ok($r->{p_value} == $r->{p_value}, 'p-value is not NaN');
	ok(scalar(@$w) >= 1, 'a warning is emitted for zero variance');
	done_testing();
};

subtest 'statistic exactly at the mean uses zero continuity correction' => sub {
	# two.sided, variance > 0, z == 0: correction must be 0 (R uses sign(z)*0.5).
	my $r = wilcox_test([1,4], [2,3], exact => 0);
	is_approx($r->{p_value}, 1.0, 1e-12, 'z==0 => p = 1');
	done_testing();
};

subtest 'invalid alternative is rejected' => sub {
	throws_ok { wilcox_test(\@x, \@y, alternative => 'twosided') }
		qr/alternative/, 'typo in alternative dies instead of silently running two-sided';
	for my $alt (qw(two.sided less greater)) {
		lives_ok { wilcox_test(\@x, \@y, alternative => $alt) } "alternative '$alt' accepted";
	}
	done_testing();
};

# ---------------------------------------------------------------------------
# Argument validation / exceptions
# ---------------------------------------------------------------------------

subtest 'argument errors croak' => sub {
	throws_ok { wilcox_test(x => 42) }         qr/ARRAY reference/, 'non-arrayref x dies';
	throws_ok { wilcox_test() }                qr/required/,        'missing x dies';
	throws_ok { wilcox_test([]) }              qr/Not enough/,      'empty x dies';
	throws_ok { wilcox_test(\@x, \@y, 'paired') } qr/Usage/,        'odd trailing args die';
	throws_ok { wilcox_test(\@x, 'bogus' => 1) }  qr/unknown/,      'unknown named arg dies';
	throws_ok { wilcox_test([1,2,3], [1,2], paired => 1) }
		qr/same length/, 'paired length mismatch dies';
	done_testing();
};

subtest 'output shape' => sub {
	my $r = wilcox_test(\@x, \@y);
	is(ref $r, 'HASH', 'returns a hashref');
	ok(exists $r->{$_}, "has '$_'") for qw(statistic p_value method alternative);
	done_testing();
};

done_testing();
