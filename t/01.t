#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use Test::More;
use Test::Exception; # die_ok
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
use Scalar::Util 'looks_like_number';
use JSON qw(decode_json encode_json);

# Gemini helped to write some of the tests
# Custom helper for floating-point comparisons
sub is_approx ($got, $expected, $test_name, $epsilon = 10**-7) {
	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
	foreach my ($i, $arg) (indexed ($got, $expected, $test_name)) {
		next if defined $arg;
		die "\$arg[$i] (see subroutine signature for name) isn't defined in $current_sub";
	}
	my $diff = abs($got - $expected);
	if ($diff < $epsilon) {
		pass($test_name);
		return 1;
	} else {
		fail($test_name);
		diag("         got: $got\n    expected: $expected; diff = $diff");
		return 0;
	}
}

sub json_file_to_ref ($json_filename) {
	die "$json_filename doesn't exist or isn't a file" unless -f $json_filename;
	die "$json_filename has 0 size" if (-s $json_filename == 0);
#	say "Reading $json_filename" if defined $_[0];
	open my $fh, '<:raw', $json_filename; # Read it unmangled
	local $/;                     # Read whole file
	my $json = <$fh>;             # This is UTF-8
#	$json =~ s/NaN/"NaN"/g;
	return decode_json($json); # This produces decoded text
}
#----------------------
#		min
#----------------------
is_approx( min(1,2,2.33,3), 1, 'min of scalars');
my @test_data = (-5..5);
is_approx(min(@test_data), $test_data[0], 'min of array');
my $test_data = \@test_data;
is_approx(min($test_data), $test_data[0], 'min of array reference');
my %h = (A => -3, B => 4, C => 9);
is_approx( min(values %h), -3, 'min takes values from hash');
#----------------------
#		max
#----------------------
is_approx( max(1,2,3), 3, 'max of scalars');
is_approx(max(@test_data), $test_data[-1], 'max of array');
is_approx(max($test_data), $test_data[-1], 'max of array reference');
is_approx( max(values %h), 9, 'max takes values from hash');
#----------------------
#		mean
#----------------------
my $mean_ok = 0;
my $mean = mean(1,2,3);
if ($mean == 2) {
	$mean_ok = 1;
} else {
	die "mean = $mean";
}
ok($mean_ok, 'Verified: Simple mean works');
my @arr = 1..8;
if (mean(@arr, 4, 5) == 4.5) {
	ok(1, 'Arrays can be given to mean and mixed');
} else {
	fail('Arrays given to mean and mixed failed');
}
if (mean([1,1], [2,2]) == 1.5) {
	ok(1, 'Arrays can be given as references');
} else {
	fail('Arrays as references cannot be given');
}
subtest 'Numerical Stability: Catastrophic Cancellation' => sub {
	my @large_data = (1000000000.1, 1000000000.2, 1000000000.3);
	# The variance of (0.1, 0.2, 0.3) is exactly 0.01.
	is_approx( var(@large_data), 0.01, 'var: handles large magnitude data cleanly' );
	is_approx( sd(@large_data), 0.1, 'sd: handles large magnitude data cleanly' );
};
# Exceptional cases for mean
eval { mean() };
like( $@, qr/mean needs >= 1 element/, 'mean: dies when given empty input' );

# -------------------------------
# standard deviation
# -------------------------------
my $stdev = sd(2,4,4,4,5,5,7,9);
my $correct = 2.1380899352994;
if (abs($stdev - $correct) < 10**-14) {
	ok(1, 'stdev works');
} else {
	my $diff = $correct - $stdev;
	fail("stdev does not work, got $stdev with an error of $diff");
}

# Exceptional cases for sd / var
eval { sd(1) };
like( $@, qr/stdev needs >= 2 elements/, 'sd: dies when given < 2 elements' );

dies_ok { var(1) } 'var: dies when given < 2 elements' ;

#------------------
# t.test
#------------------
@test_data = (
[
	[27.5,21.0,19.0,23.6,17.0,17.9,16.9,20.1,21.9,22.6,23.1,19.6,19.0,21.7,21.4],
	[27.1,22.0,20.8,23.4,23.4,23.5,25.8,22.0,24.8,20.2,21.9,22.1,22.9,20.5,24.4],
],
[
	[17.2,20.9,22.6,18.1,21.7,21.4,23.5,24.2,14.7,21.8],
	[21.5, 22.8, 21.0, 23.0, 21.6, 23.6, 22.5, 20.7, 23.4, 21.8, 20.7, 21.7, 21.5, 22.5, 23.6, 21.5, 22.5, 23.5,21.5,21.8]
],
[
	[19.8,20.4,19.6,17.8,18.5,18.9,18.3,18.9,19.5,22.0],
	[28.2,26.6,20.1,23.3,25.2,22.1,17.7,27.6,20.6,13.7,23.2,17.5,20.6,18.0,23.9,21.6,24.3,20.4, 24.0,13.2]
],
[
	[30.02,29.99,30.11,29.97,30.01,29.99],
	[29.89,29.93,29.72,29.98,30.02,29.98]
],
[
	[3.0,4.0,1.0,2.1],
	[490.2,340.0,433.9]
],
[
	[0.010268,0.000167,0.000167],
	[0.159258,0.136278,0.122389]
],
[
	[1.0/15,10.0/62.0],
	[1.0/10,2/50.0]
],
[
	[9/23.0,21/45.0,0/38.0],
	[0/44.0,42/94.0,0/22.0]
]);
foreach my $i (0..$#test_data) { # single sample t-tests
	foreach my $j (0,1) {
		my $t_test = t_test( 'x' => $test_data[$i][$j], mu => mean( $test_data[$i][$j] ));
		is_approx( $t_test->{p_value}, 1,            "t_test: Testing set $i/$j p-value");
		is_approx( $t_test->{df}, scalar @{ $test_data[$i][$j] } - 1, "t_test: df $i/$j");
		is_approx( $t_test->{statistic}, 0, "t_test: t $i/$j");
	}
}
my @correct_t = (
	{ # default
		conf_int     => [
			-3.98409625405368, -0.349237079279662
		],
		df           => 24.9885292902309,
		'estimate_x' => 20.82,
		'estimate_y' => 22.9866666666666,
		p_value      => 0.021378001462867,
		statistic    => -2.45535639828601
	},
	{ # var.equal = True (Student's t-test)
		conf_int     => [
			-3.0124986, -0.0375014
		],
		df           => 28,
		'estimate_x' => 20.610,
		'estimate_y' => 22.135,
		p_value      => 0.04485852,
		statistic    => -2.10004963761047
	},
	{ # paired = true
		conf_int     => [
			-0.06672889, 0.25672889
		],
		df        => 5,
		estimate  => 0.095,
		p_value   => 0.19143688433660,
		statistic => 1.50996688705414
	}
);
my $t_test = t_test(
	'x' => $test_data[0][0],
	'y' => $test_data[0][1]
);
foreach my $key (grep {ref $correct_t[0]{$_} eq ''} keys %{ $correct_t[0] }) {
	if (not defined $t_test->{$key}) {
		die "$key isn't defined in test";
	}
	is_approx( $t_test->{$key}, $correct_t[0]{$key}, "t_test var_equal = true; $key");
}
foreach my $j (0,1) {
	is_approx( $t_test->{'conf_int'}[$j], $correct_t[0]{'conf_int'}[$j], "Conf. interval index $j");
}
$t_test = t_test(
	'x'       => $test_data[1][0],
	'y'       => $test_data[1][1],
	var_equal => true
);
foreach my $key (grep {ref $correct_t[1]{$_} eq ''} keys %{ $correct_t[1] }) {
	if (not defined $t_test->{$key}) {
		die "$key isn't defined in test";
	}
	is_approx( $t_test->{$key}, $correct_t[1]{$key}, "t_test var_equal = true; $key");
}
foreach my $j (0,1) {
	is_approx( $t_test->{'conf_int'}[$j], $correct_t[1]{'conf_int'}[$j], "Conf. interval index $j");
}
# start new test
$t_test = t_test(
	'x'    => $test_data[3][0],
	'y'    => $test_data[3][1],
	paired => true
);
foreach my $key (grep {ref $correct_t[2]{$_} eq ''} keys %{ $correct_t[2] }) {
	if (not defined $t_test->{$key}) {
		die "$key isn't defined in test";
	}
	is_approx( $t_test->{$key}, $correct_t[2]{$key}, "t_test var_equal = true; $key");
}
foreach my $j (0,1) {
	is_approx( $t_test->{'conf_int'}[$j], $correct_t[2]{'conf_int'}[$j], "Conf. interval index $j");
}

# t_test exceptions & alternative hypotheses tests
eval { t_test(y => [1..5]) };
like( $@, qr/must be an ARRAY reference/, 't_test: dies when x is missing' );

eval { t_test('x' => [1..5], paired => 1) };
like( $@, qr/'y' must be provided for paired or two-sample tests/, 't_test: dies paired without y' );

eval { t_test('x' => [1..5], 'y' => [1..4], paired => 1) };
like( $@, qr/Paired arrays must be same length/, 't_test: dies on mismatched paired arrays' );

eval { t_test('x' => [1..5], conf_level => 1.5) };
like( $@, qr/'conf_level' must be between 0 and 1/, 't_test: dies on invalid conf_level' );

$t_test = t_test('x' => [5, 6, 7, 8, 9], mu => 2, alternative => 'greater');
ok( $t_test->{p_value} < 0.05, 't_test alternative greater works (small p_value)' );

$t_test = t_test('x' => [5, 6, 7, 8, 9], mu => 20, alternative => 'less');
ok( $t_test->{p_value} < 0.05, 't_test alternative less works (small p_value)' );

dies_ok {
	t_test( 'x' => [3,3,3,3] )
} '"t_test" dies when data is constant';

#----------------------
#		p ajdust
#----------------------
my @pvalues = (4.533744e-01, 7.296024e-01, 9.936026e-02, 9.079658e-02, 1.801962e-01,
8.752257e-01, 2.922222e-01, 9.115421e-01, 4.355806e-01, 5.324867e-01,
4.926798e-01, 5.802978e-01, 3.485442e-01, 7.883130e-01, 2.729308e-01,
8.502518e-01, 4.268138e-01, 6.442008e-01, 3.030266e-01, 5.001555e-02,
3.194810e-01, 7.892933e-01, 9.991834e-01, 1.745691e-01, 9.037516e-01,
1.198578e-01, 3.966083e-01, 1.403837e-02, 7.328671e-01, 6.793476e-02,
4.040730e-03, 3.033349e-04, 1.125147e-02, 2.375072e-02, 5.818542e-04,
3.075482e-04, 8.251272e-03, 1.356534e-03, 1.360696e-02, 3.764588e-04,
1.801145e-05, 2.504456e-07, 3.310253e-02, 9.427839e-03, 8.791153e-04,
2.177831e-04, 9.693054e-04, 6.610250e-05, 2.900813e-02, 5.735490e-03);

my %correct = (
	'Benjamini-Hochberg' => [6.126681e-01, 8.521710e-01, 1.987205e-01, 1.891595e-01, 3.217789e-01,
9.301450e-01, 4.870370e-01, 9.301450e-01, 6.049731e-01, 6.826753e-01,
6.482629e-01, 7.253722e-01, 5.280973e-01, 8.769926e-01, 4.705703e-01,
9.241867e-01, 6.049731e-01, 7.856107e-01, 4.887526e-01, 1.136717e-01,
4.991891e-01, 8.769926e-01, 9.991834e-01, 3.217789e-01, 9.301450e-01,
2.304958e-01, 5.832475e-01, 3.899547e-02, 8.521710e-01, 1.476843e-01,
1.683638e-02, 2.562902e-03, 3.516084e-02, 6.250189e-02, 3.636589e-03,
2.562902e-03, 2.946883e-02, 6.166064e-03, 3.899547e-02, 2.688991e-03,
4.502862e-04, 1.252228e-05, 7.881555e-02, 3.142613e-02, 4.846527e-03,
2.562902e-03, 4.846527e-03, 1.101708e-03, 7.252032e-02, 2.205958e-02],
	'Benjamini-Yekutieli' => [1.000000e+00, 1.000000e+00, 8.940844e-01, 8.510676e-01, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 5.114323e-01,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.754486e-01, 1.000000e+00, 6.644618e-01,
7.575031e-02, 1.153102e-02, 1.581959e-01, 2.812089e-01, 1.636176e-02,
1.153102e-02, 1.325863e-01, 2.774239e-02, 1.754486e-01, 1.209832e-02,
2.025930e-03, 5.634031e-05, 3.546073e-01, 1.413926e-01, 2.180552e-02,
1.153102e-02, 2.180552e-02, 4.956812e-03, 3.262838e-01, 9.925057e-02],
	'Bonferroni' => [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 7.019185e-01, 1.000000e+00, 1.000000e+00,
2.020365e-01, 1.516674e-02, 5.625735e-01, 1.000000e+00, 2.909271e-02,
1.537741e-02, 4.125636e-01, 6.782670e-02, 6.803480e-01, 1.882294e-02,
9.005725e-04, 1.252228e-05, 1.000000e+00, 4.713920e-01, 4.395577e-02,
1.088915e-02, 4.846527e-02, 3.305125e-03, 1.000000e+00, 2.867745e-01],

	'Hochberg' => [9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 4.632662e-01, 9.991834e-01, 9.991834e-01,
1.575885e-01, 1.383967e-02, 3.938014e-01, 7.600230e-01, 2.501973e-02,
1.383967e-02, 3.052971e-01, 5.426136e-02, 4.626366e-01, 1.656419e-02,
8.825610e-04, 1.252228e-05, 9.930759e-01, 3.394022e-01, 3.692284e-02,
1.023581e-02, 3.974152e-02, 3.172920e-03, 8.992520e-01, 2.179486e-01],
	'Holm' => [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
1.000000e+00, 1.000000e+00, 4.632662e-01, 1.000000e+00, 1.000000e+00,
1.575885e-01, 1.395341e-02, 3.938014e-01, 7.600230e-01, 2.501973e-02,
1.395341e-02, 3.052971e-01, 5.426136e-02, 4.626366e-01, 1.656419e-02,
8.825610e-04, 1.252228e-05, 9.930759e-01, 3.394022e-01, 3.692284e-02,
1.023581e-02, 3.974152e-02, 3.172920e-03, 8.992520e-01, 2.179486e-01],

	'Hommel' => [9.991834e-01, 9.991834e-01, 9.991834e-01, 9.987624e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.595180e-01,
9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01, 9.991834e-01,
9.991834e-01, 9.991834e-01, 4.351895e-01, 9.991834e-01, 9.766522e-01,
1.414256e-01, 1.304340e-02, 3.530937e-01, 6.887709e-01, 2.385602e-02,
1.322457e-02, 2.722920e-01, 5.426136e-02, 4.218158e-01, 1.581127e-02,
8.825610e-04, 1.252228e-05, 8.743649e-01, 3.016908e-01, 3.516461e-02,
9.582456e-03, 3.877222e-02, 3.172920e-03, 8.122276e-01, 1.950067e-01]);

foreach my $method ('Hochberg','Benjamini-Hochberg','Benjamini-Yekutieli', 'Bonferroni', 'Holm', 'Hommel') {
	my @q = p_adjust(\@pvalues, $method);
	my $error = 0.0;
	foreach my $q (0..$#q) {
		$error += abs($q[$q] - $correct{$method}[$q]);
	}
	if ($error < 10**-6) {
		ok(1, "$method works with cumulative error of $error");
	} else {
		fail("$method doesn't work for FDR correction with error = $error");
	}
}

# p_adjust exceptions
eval { p_adjust("not array") };
like( $@, qr/first argument must be an ARRAY reference/, 'p_adjust: dies on non-array' );

eval { p_adjust(\@pvalues, "invalid_method") };
like( $@, qr/Unknown p-value adjustment method/, 'p_adjust: dies on invalid method' );

my @empty_p = p_adjust([]);
is_deeply( \@empty_p, [], 'p_adjust handles empty arrayref gracefully' );

#----------------------
#		var
#----------------------
my @ans = (2.5, 8.3);
foreach my ($i, $arr) (indexed([1..5], [2, 4, 5, 8, 9])) {
	my $var = var( $arr );
	is_approx( $var, $ans[$i], "Test index $i for var");
}
#----------------------
#		median
#----------------------
@ans = (21, 21.55, 19.2);
foreach my ($i, $ans) (indexed @ans) {
	my $median = median( $test_data[$i][0] );
	is_approx( $median, $ans, "Median test $i");
}

eval { median() };
like( $@, qr/median needs >= 1 element/, 'median: dies when given empty input' );

#----------------------
#		cor
#----------------------
$test_data[0] = [1, 2, 3, 4, 5,  5, 6,  7,   8];
$test_data[1] = [2, 4, 6, 8, 10, 9, 12, 14, 16];
%correct = (
	pearson  => 0.9973649,	spearman => 0.9958246,	kendall  => 0.9860133
);
foreach my $method (sort keys %correct) {
	is_approx(
		cor($test_data[0], $test_data[1], $method),
		$correct{$method},
		"cor: method = \"$method\""
	);
}

# cor exceptions and matrix support tests
eval { cor([1,2,3], [1,2], 'pearson') };
like( $@, qr/x and y must have the same length/, 'cor: dies on mismatched flat vector lengths' );

eval { cor([1,2,3], undef, 'unknown_method') };
like( $@, qr/unknown method/, 'cor: dies on unknown method' );

my $mat_x = [[1, 2], [3, 4], [5, 6]];
my $mat_y = [[6, 5], [4, 3], [2, 1]];
my $cor_matrix = cor($mat_x, $mat_y);
is( ref($cor_matrix), 'ARRAY', 'cor with matrices returns an array reference' );

# It flattens the input and returns a standard array
#----------------------
#  SCALE
#----------------------
my @scaled_results = scale(1..5);
my @correct_scaled = (-1.2649111, -0.6324555, 0.0000000, 0.6324555, 1.2649111);
foreach my ($i, $correct) (indexed @correct_scaled) {
	is_approx($scaled_results[$i], $correct, "index $i for scale");
}
@scaled_results = scale(1..5, { center => false });
@correct_scaled = (0.269679944985297, 0.539359889970594, 0.80903983495589, 1.07871977994119, 1.34839972492648);
foreach my ($i, $correct) (indexed @correct_scaled) {
	is_approx($scaled_results[$i], $correct, "index $i for scale");
}
@scaled_results = scale(1..5, {center => true, scale => false});
#p @scaled_results;
@correct_scaled = (-2..2);
foreach my ($i, $correct) (indexed @correct_scaled) {
	is_approx($scaled_results[$i], $correct, "index $i for scale");
}

# scale matrix tests and exceptions
eval { scale(1) };
like( $@, qr/scale needs >= 2 elements to calculate SD/, 'scale: dies with 1 element for default SD scale' );

my $scaled_mat = scale([[1, 2], [3, 4], [5, 6]]);
is( ref($scaled_mat), 'ARRAY', 'scale on matrix returns an array reference' );

#-----------------------
#			MATRIX
#-----------------------
my $matrix_correct = '[[1,3,5],[2,4,6]]';
my $mat1 = matrix(
	data => [1..6],
	nrow => 2
);
if (encode_json($mat1) eq $matrix_correct) {
	pass('simple "matrix" works');
} else {
	fail('simple "matrix" fails');
}

# matrix exceptions
eval { matrix(data => "string", nrow => 2) };
like( $@, qr/must be an array reference/, 'matrix: dies on non-arrayref data' );

eval { matrix(data => [1,2], nrow => 0, ncol => 0) };
like( $@, qr/Dimensions must be greater than 0/, 'matrix: dies on 0 dimensions' );

eval { matrix(data => []) };
like( $@, qr/Data array cannot be empty/, 'matrix: dies on empty data array' );

#$matrix_correct = '[[-0.707106781186547,-0.707106781186547,-0.707106781186547],[0.707106781186547,0.707106781186547,0.707106781186547]]';

#@scaled_results = scale(
#	[1,2,3],
#	[4,5,6]
#);
#p @scaled_results;
# Output: -1.46385, -0.87831, -0.29277, 0.29277, 0.87831, 1.46385
#---------------------------
#       lm
#----------------------------
my $mtcars = json_file_to_ref('mtcars.hoh.json');
my $lm = lm(formula =>  'mpg ~ wt * hp^2', data => $mtcars);
#p $lm;
%correct = (
	coefficients => {
		Intercept => 49.8084234287587,	hp        => -0.120102090978019,
		wt        => -8.21662429724302,	'wt:hp'   => 0.0278481483187383
	},
	'df.residual' => 28,
	'fitted.values' => {
		'Mazda RX4' => 23.09547, 				'Mazda RX4 Wag' => 21.78138,
		'Datsun 710'    => 25.58488, 			'Hornet 4 Drive' => 20.02924,
		'Hornet Sportabout' => 17.28996, 	Valiant             => 18.88542,
		'Duster 360'        => 15.40745, 	'Merc 240D'         => 21.65887,
		'Merc 450SE'        => 15.14994, 	'Merc 450SL'          => 16.23929,
		'Merc 450SLC'         =>16.07909, 	'Cadillac Fleetwood'  => 12.02179,
		'Lincoln Continental' =>11.89490, 	'Chrysler Imperial'   => 12.50221,
		'Fiat 128'            => 27.84866, 	'Honda Civic'   => 32.63195,
		'Toyota Corolla'   => 30.24587, 		'Toyota Corona'   => 24.56317,
		'Dodge Challenger'   => 17.57441, 	'AMC Javelin'   => 17.91776,
		'Camaro Z28'   => 15.03111, 			'Pontiac Firebird'   => 15.93596,
		'Fiat X1-9'   => 29.53900, 			'Porsche 914-2'   => 26.71871,
		'Lotus Europa'   => 28.56630, 		'Ford Pantera L'   =>15.36033,
		'Ferrari Dino'   => 19.52990, 		'Maserati Bora'   => 13.54587,
		'Volvo 142E'      => 22.31363
	},
	rank => 4,
		residuals  => {
		'AMC Javelin'        =>   -2.7177637422554,	'Cadillac Fleetwood' =>   -1.62178684578001,
		'Camaro Z28'         =>   -1.73111177599938, 'Chrysler Imperial'  =>   2.19779322930961,
		'Datsun 710'         =>   -2.78487707944995,	'Dodge Challenger'   =>   -2.07441456805362,
		'Duster 360'         =>   -1.10744532497053,	'Ferrari Dino'       =>   0.17010189824968,
		'Fiat 128'           =>   4.55133689384449,	'Fiat X1-9'           =>  -2.23900443083033,
		'Ford Pantera L'     => 0.439669246713233,	'Honda Civic'         =>  -2.23195395366212,
		'Hornet 4 Drive'     =>   1.37075604153847,	'Hornet Sportabout'  =>   1.41004478703069,
		'Lincoln Continental' =>  -1.4949003236173,	'Lotus Europa'       =>   1.83369534357951,
		'Maserati Bora'      =>   1.45413280824033,	'Mazda RX4'           =>  -2.09547410785999,
		'Mazda RX4 Wag'       =>  -0.781375472403504,'Merc 230'            =>  1.95008336608675,
		'Merc 240D'           =>  2.74113094560431,	'Merc 280'            =>  0.646212827429729,
		'Merc 280C'           =>  -0.75378717257027,	'Merc 450SE'          =>  1.25006037875687,
		'Merc 450SL'          =>  1.06071479480091,	'Merc 450SLC'         =>  -0.879087325205568,
		'Pontiac Firebird'    =>  3.26404011532368,	'Porsche 914-2'       =>  -0.718705557249944,
		'Toyota Corolla'      =>  3.65413017953585,	'Toyota Corona'       =>  -3.06317321493851,
		Valiant               =>  -0.785416091802773,'Volvo 142E'          =>  -0.913625869362747
	}
);
foreach my $key ('Intercept', 'hp', 'wt', 'wt:hp') {
	unless (defined $lm->{coefficients}{$key}) {
		#p $lm;
		die "\"$key\" isn't defined" ;
	}
	is_approx( $lm->{coefficients}{$key}, $correct{coefficients}{$key}, "Checking lm's $key" );
}
foreach my $key ('df.residual', 'rank') {
	unless (defined $lm->{$key}) {
		#p $lm;
		die "\"$key\" isn't defined" ;
	}
	is_approx( $lm->{$key}, $correct{$key}, "Checking \"$key\"");
}
foreach my $key ('fitted.values', 'residuals') {
	foreach my $car (keys %{ $correct{$key} }) {
		unless (defined $lm->{$key}{$car}) {
			#p $lm;
			die "\"$car\" isn't defined in \"fitted.values\"" ;
		}
		is_approx(
			$lm->{$key}{$car},
			$correct{$key}{$car},
			"Checking $key \"$car\"",
			10**-5
		);
	}
}
$lm = lm(formula =>  'mpg ~ wt + hp', data => $mtcars);
%correct = (
	coefficients => {
		Intercept => 37.22727,
		hp        => -0.03177,
		wt        => -3.87783,
	},
	rank => 3,
	'df.residual' => 29,
	summary => {
		hp => {
  			Estimate      => -0.03177295,
      	'Pr(>|t|)'    => 0.00145122851187347,
       	'Std. Error'  => 0.0090297096758557,
         't value'     => -3.518712
      },
      wt => {
      	Estimate      => -3.87783074,
      	'Pr(>|t|)'    => 1.119647e-06,
       	'Std. Error'  => 0.63273349,
         't value'     => -6.128695
      },
      Intercept => {
      	Estimate      => 37.22727012,
      	'Pr(>|t|)'    => 2.565459e-20,
       	'Std. Error'  => 1.59878754,
         't value'     => 23.284689
      }
	}
);
foreach my $key ('Intercept', 'hp', 'wt') {
	unless (defined $lm->{coefficients}{$key}) {
		#p $lm;
		die "\"$key\" isn't defined" ;
	}
	is_approx(
		$lm->{coefficients}{$key},
		$correct{coefficients}{$key},
		"Checking lm's $key",
		0.1
	);
	unless (defined $lm->{summary}{$key}) {
		die "$key is missing summary";
	}
	foreach my $val ('Estimate', 'Pr(>|t|)', 'Std. Error', 't value') {
		my $e = 10**-7;
		if ($correct{summary}{$key}{$val} =~ m/\.(\d+)$/) {
			$e = 10**(2-length $1);
		}
		is_approx( $lm->{summary}{$key}{$val}, $correct{summary}{$key}{$val}, "lm: Summary $key & $val", $e);
	}
}
foreach my $key ('df.residual', 'rank') {
	unless (defined $lm->{$key}) {
		die "\"$key\" isn't defined" ;
	}
	is_approx( $lm->{$key}, $correct{$key}, "Checking \"$key\"");
}

# lm exceptions and additional tests
eval { lm(data => $mtcars) };
like( $@, qr/formula is required/, 'lm: dies without a formula' );

dies_ok {
	lm(formula => 'mpg wt'); # missing ~
} 'lm: dies on bad formula lacking ~';
dies_ok {
	lm(formula => 'mpg ~ wt', data => 'not_a_hash');
} 'lm: dies when given a non-hash';

dies_ok {
	my $lm_no_int = lm(formula => 'mpg ~ wt -1', data => $mtcars);
	ok( !defined($lm_no_int->{coefficients}{Intercept}), 'lm: formula -1 correctly suppresses Intercept' );
} 'lm: formula -1 correctly suppresses Intercept';

#-------------------------------------------------------------------
#  lm: Matching R's Statistical Behaviors
#-------------------------------------------------------------------

subtest 'lm: R-Parity edge cases (NAs, Strings, Collinearity)' => sub {
	my $messy_data = {
		'R1' => { 'y' => 10, 'x1' => 2,  'x2' => 4 },
		'R2' => { 'y' => 12, 'x1' => 3,  'x2' => 6 },
		'R3' => { 'y' => 15, 'x1' => undef, x2 => 8 }, # Missing predictor
		'R4' => { 'y' => 20, 'x1' => 5,  x2 => 10 },
		'R5' => { 'y' => undef, x1 => 6, x2 => 12 },   # Missing response
		'R6' => { 'y' => 22, x1 => 7,  x2 => "blue" }, # String/Factor
		'R7' => { 'y' => 25, x1 => 8,  x2 => 16 },
	};
	my $lm_res;
	lives_ok {
		$lm_res = lm(formula => 'y ~ x1 + x2', data => $messy_data);
	} 'lm: Handles missing data, strings, and collinearity without croaking';
	# 1. Check Listwise Deletion (NAs and Strings should drop the row)
	# Only R1, R2, R4, R7 are strictly valid numeric rows. (N = 4)
	is( $lm_res->{'df.residual'}, 2, 'lm: correctly drops rows with undef or strings (Listwise Deletion)' );
	ok( !exists $lm_res->{'fitted.values'}{'R3'}, 'lm: omitted row R3 has no fitted value' );
	
	# 2. Check Collinearity (x2 is exactly x1 * 2)
	# R sets aliased coefficients to NA/undef.
	ok( exists $lm_res->{coefficients}{x2}, 'lm: collinear coefficient key exists' );
	
	# In Perl, a C 'NaN' translates to a string "NaN" or undef depending on the platform.
	my $x2_coef = $lm_res->{coefficients}{x2};
	ok( !defined($x2_coef) || $x2_coef eq 'NaN' || $x2_coef != $x2_coef, 
	    'lm: aliased/collinear variable x2 is correctly flagged as NaN/undef (matches R)' );
	
	# 3. Check Rank
	is( $lm_res->{rank}, 2, 'lm: rank is reduced to 2 (Intercept + x1)' );
};
#---------------------------
#   rnorm
#----------------------------
my ($rmean, $sd, $n) = (10, 2, 9999);
my $normals = rnorm( n => $n, mean => $rmean, sd => $sd);
is_approx(scalar @{ $normals }, $n, 'rnorm sample size');
is_approx(mean($normals), $rmean, 'rnorm mean', 0.1);
is_approx(sd($normals), $sd, 'rnorm sd', 0.1);
# rnorm exceptions
eval { rnorm(n => 10, sd => -1) };
like( $@, qr/standard deviation must be non-negative/, 'rnorm: dies on negative sd' );

eval { rnorm(n => 10, mean => 0, 'missing_value_key') };
like( $@, qr/must be even key\/value pairs/, 'rnorm: dies on odd argument count' );
#----------------------
#    quantile
#----------------------
my $quantile = quantile('x' => [1..99], probs => [0.05, 0.1, 0.25]);
# R equivalent: quantile(1:99, probs=c(0.01,0.1,0.25))
@ans = (5.9, 10.80, 25.50);
my @quantile_keys = ('5%', '10%', '25%');
foreach my $idx (0..$#ans) {
	is_approx($quantile->{$quantile_keys[$idx]}, $ans[$idx], "quantile: $quantile_keys[$idx]", 10**-14);
}
p $quantile;
#----------------------
#    Fisher's Test
#----------------------
my $array_data = [
	[10, 2],
	[3, 15]
];
my $ft = fisher_test($array_data);
p $ft; # R equivalent: fisher.test( matrix(c(10,2,3,15), nrow = 2)))
%correct = (
	p_value => 0.0005367241
);
is_approx( $ft->{p_value}, $correct{p_value}, 'Fisher\'s test p-value' );
my $conf_int_range = abs $ft->{conf_int}[0] - $ft->{conf_int}[1];
my $correct_conf_int_range = 301.462337971516 - 2.75338278824932;
if (0.99*$correct_conf_int_range < $conf_int_range < 1.01* $correct_conf_int_range) {
	pass('Fisher\'s test is within 1% of correct: ');
} else {
	fail('Fisher\'s test is *NOT* within 1% of correct: ');
}
is_approx( $ft->{estimate}{'odds ratio'}, 21.30533, 'Fisher\'s test odds ratio', 10**-3);
#----------------------
#    hist
#----------------------
subtest 'hist functionality' => sub {
	# 1. Basic properties with a simple dataset
	# Data: 1, 2, 2, 3, 3, 3, 4, 4, 5 (9 elements)
	my $h_data = [1, 2, 2, 3, 3, 3, 4, 4, 5];
	my $breaks = 4;
	my $res = hist($h_data, breaks => $breaks);
	is(ref $res, 'HASH', 'hist: returns a hash reference');
	is(scalar @{$res->{counts}}, $breaks, 'hist: correct number of bins (counts)');
	is(scalar @{$res->{breaks}}, $breaks+1, 'hist: correct number of breaks (n+1)');
	is(scalar @{$res->{mids}},   $breaks, 'hist: correct number of midpoints');

	# 2. Verify counts sum to total elements
	my $total_counts = 0;
	$total_counts += $_ for @{$res->{counts}};
	is($total_counts, scalar @$h_data, 'hist: sum of counts matches input size');

	# 3. Verify midpoints are mathematically correct
	my $mid_ok = 1;
	for my $i (0 .. $#{$res->{mids}}) {
	  my $expected_mid = ($res->{breaks}->[$i] + $res->{breaks}->[$i+1]) / 2;
	  $mid_ok = 0 if abs($res->{mids}->[$i] - $expected_mid) > 1e-12;
	}
	ok($mid_ok, 'hist: midpoints are correctly calculated between breaks');

	# 4. Density check: Total area (sum of density * bin_width) must be 1
	my $area = 0;
	for my $i (0 .. $#{$res->{counts}}) {
	  my $width = $res->{breaks}->[$i+1] - $res->{breaks}->[$i];
	  $area += $res->{density}->[$i] * $width;
	}
	ok(abs($area - 1.0) < 1e-12, 'hist: total area under density curve is 1.0');

	# 5. Test with a single value (Edge case)
	my $single = hist([10], breaks => 1);
	is($single->{counts}->[0], 1, 'hist: handles single-element array');
	is($single->{breaks}->[0], 10, 'hist: single-element break starts at value');

	# Optional: Visual inspection if needed
	p $res;
};
#-------------------------------------------------------------------
#  Performance & Edge Cases for hist()
#-------------------------------------------------------------------
subtest 'hist: O(N) Performance and Edge Cases' => sub {
	# 1. Test standard binning boundaries
	# Data spans 0 to 10 exactly.
	my $data = [0, 2.5, 4.9, 5.0, 7.5, 10];
	my $res = hist($data, breaks => 2);
	
	# With 2 bins spanning 0 to 10, step is 5.
	# Bin 0: [0, 5] -> should capture 0, 2.5, 4.9, 5.0 (4 items)
	# Bin 1: (5, 10] -> should capture 7.5, 10 (2 items)
	is_deeply($res->{counts}, [4, 2], 'hist: correctly assigns boundary values in O(N) logic');
	is_deeply($res->{breaks}, [0, 5, 10], 'hist: correct breaks generation');

	# 2. Test zero-variance edge case (all identical values)
	my $flat_data = [7, 7, 7, 7, 7];
	my $flat_res = hist($flat_data, breaks => 1);
	
	is($flat_res->{counts}->[0], 5, 'hist: properly handles flatline data (step = 0)');
	is($flat_res->{density}->[0], 1, 'hist: density is 1.0 for a single zero-width bin');
	
	# 3. Scale test: Generate a slightly larger array to ensure no memory segfaults
	# and quick processing.
	my @large_uniform = seq(1, 1000, 0.5); 
	my $large_res = hist(\@large_uniform, breaks => 10);
	
	my $total_counted = 0;
	$total_counted += $_ for @{ $large_res->{counts} };
	
	is($total_counted, scalar @large_uniform, 'hist: sum of counts perfectly matches input size on larger datasets');
};
#----------------------
#    hist exceptions
#----------------------
subtest 'hist exceptions' => sub {
	# Should die if not an array ref
	dies_ok { hist("not an array") } 'hist: dies on string input';
	dies_ok { hist({ a => 1 }) }     'hist: dies on hash ref input';
	# Should die on empty array
	dies_ok { hist([]) }             'hist: dies on empty array ref';
	# Should die on non-numeric data (depending on your SVNV strictness)
	dies_ok { hist([qw(a b c)]) }    'hist: dies on non-numeric array content';
};
#----------------------
#   runif
#----------------------
my $unif = runif( n => $n, min => 0, max => 1);
if (scalar @{ $unif } == $n) {
	pass('random uniform distribution has the correct # of elements');
} else {
	fail('random uniform distribution does NOT have the correct # of elements');
}
is_approx( min(@{ $unif }), 0, 'Approximately correct minimum', 10**-3);
is_approx( max(@{ $unif }), 1, 'Approximately correct maximum', 10**-3);
{
	my $unif2 = runif( n => $n, min => 0, max => 1);
	p $unif2;
	my @identical_idx = grep { $unif->[$_] == $unif2->[$_] } 0..$n-1;
	if (scalar @identical_idx == 0) {
		pass('runif does not repeat');
	} else { # > 1 identical value
		fail('runif repeats ' . scalar @identical_idx . "/$n values");
	}
}
#----------------------
#      rbinom
#----------------------
my $binom = rbinom( n => $n, prob => 0.5, size => 9);
if (scalar @{ $binom } == $n) {
	pass('binom has the correct # of elements');
} else {
	fail("binom should have $n elements, but has " . scalar @{ $binom } . ' elements');
}
subtest 'rbinom: Exceptions and Validations' => sub {
	eval { rbinom(n => 10, size => 10) };
	like( $@, qr/'size' and 'prob' are required arguments/, 'rbinom: dies when prob is missing' );

	dies_ok {
		rbinom(n => 10, prob => 0.5)
	} 'qr/"size" and "prob" are required arguments: rbinom: dies when size is missing';

	dies_ok {
		rbinom(n => 10, size => 10, prob => 1.5)
	} 'qr/prob must be between 0 and 1: rbinom: dies when prob > 1';

	dies_ok {
		rbinom(n => 10, size => 10, prob => -0.1) 
	} 'qr/prob must be between 0 and 1: rbinom: dies when prob < 0';
	dies_ok {
		rbinom(n => 10, size => 10, prob => 0.5, 'extra_arg');
	} 'rbinom: dies on odd argument count';
};

subtest 'rbinom: Edge Cases (Deterministic)' => sub {
	my $n_edge = 100;
	# Prob = 0 should strictly return 0
	my $binom_p = rbinom(n => $n_edge, size => 10, prob => 0);
	is_approx( max($binom_p), 0, 'rbinom: prob=0 always returns 0' );
	is_approx( min($binom_p), 0, 'rbinom: prob=0 always returns 0' );
	
	# Prob = 1 should strictly return 'size'
	$binom_p = rbinom(n => $n_edge, size => 15, prob => 1);
	is_approx( min($binom_p), 15, 'rbinom: prob=1 always returns size' );
	is_approx( max($binom_p), 15, 'rbinom: prob=1 always returns size' );

	# Size = 0 should strictly return 0
	$binom_p = rbinom(n => $n_edge, size => 0, prob => 0.5);
	is_approx( max($binom_p), 0, 'rbinom: size=0 always returns 0' );
};

subtest 'rbinom: Statistical Properties' => sub {
	my $n_stat = 10000;
	my $size   = 20;
	my $prob   = 0.3;
	my $binom_stats = rbinom(n => $n_stat, size => $size, prob => $prob);
	
	is_approx( scalar @{ $binom_stats }, $n_stat, 'rbinom: returns correct number of elements' );
	
	my $expected_mean = $size * $prob;
	my $expected_var  = $size * $prob * (1 - $prob);
	
	# Allow a small epsilon (margin of error) for stochasticity
	is_approx( mean($binom_stats), $expected_mean, 'rbinom: empirical mean matches theoretical mean', 0.1 );
	is_approx( var($binom_stats), $expected_var, 'rbinom: empirical variance matches theoretical variance', 0.5 );
};

subtest 'rbinom: Randomness' => sub {
	my $n_rand = 100;
	my $binom_1 = rbinom(n => $n_rand, size => 10, prob => 0.5);
	my $binom_2 = rbinom(n => $n_rand, size => 10, prob => 0.5);
	
	my @identical_idx = grep { $binom_1->[$_] == $binom_2->[$_] } 0..($n_rand-1);
	
	if (scalar @identical_idx < $n_rand) {
		pass('rbinom: consecutive calls generate different non-repeating sequences');
	} else {
		fail('rbinom: consecutive calls generated identical arrays (PRNG state not updating)');
	}
};
#----------------------
#       seq
#----------------------
# Example 1: Standard integer sequence
say 'seq(1, 5):';
my @seq = seq(1, 5);
say join(', ', @seq), "\n";

foreach my ($idx, $item) (indexed @seq) {
	is_approx( $item, $idx + 1, "seq item $idx");
}

# Example 2: Fractional steps
say 'seq(1, 2, 0.25):';
@seq = seq(1, 2, 0.25);
say join(", ", @seq), "\n";
for (my $idx = 2; $idx >= 1; $idx -= 0.25) { # count down to pop
	is_approx(pop @seq, $idx, "seq item $idx with fractional step");
}

# Example 3: Negative steps
say 'seq(10, 5, -1):';
@seq = seq(10, 5, -1);
say join(", ", @seq), "\n";
for (my $idx = 5; $idx <= 10; $idx++) { # count down to pop
	is_approx(pop @seq, $idx, "seq item $idx with negative step");
}
# Example 4: R-style floating point boundary catch
# (In naive C, 2.0 - 0.2 could cause the last element to drop off)
say 'seq(0, 1, 0.1):';
@seq = seq(0, 1, 0.1);
say join(", ", @seq);

subtest 'seq: Floating-point precision drift' => sub {
	my @seq_drift = seq(0, 100, 0.1);
	# If current += by is used, the 500th element (50.0) might evaluate to 49.999999999999
	is_approx( $seq_drift[500], 50.0, 'seq: 500th fractional step maintains exact expected scale without drifting', 10**-12 );
};

#my $wt_result = wilcox_test( 'x' => [1..4], 'y' => [5..8], {});
#p $wt_result;

#----------------------
#       Shapiro Test
#----------------------
my $shapiro = shapiro_test(
	[1..5]
);
is_approx( $shapiro->{p_value}, 0.9671739, 'Shapiro p-value');
is_approx( $shapiro->{W}, 0.9867622, 'Shapiro W');
#--------
$shapiro = shapiro_test(
	[1..19]
);
is_approx( $shapiro->{p_value}, 0.5896506, 'Shapiro p-value: 19 values');
is_approx( $shapiro->{W}, 0.9608707, 'Shapiro W: 19 values');
#------- cor test

my $x = [1, 2, 3, 4, 5];
my $y = [2, 1, 4, 3, 5];

my @correct = (
{ # cor.test(cx, cy, alternative='two.sided', method = 'spearman', continuity=1)
	alternative => 'two.sided',
	estimate    => 0.8,
	'p.value'   => 0.1333333,
	statistic   => 4,
	method      => 'spearman'
},
{ # cor.test(cx, cy, alternative='two.sided', method = 'kendall', continuity=1)
	alternative => 'two.sided',
	estimate    => 0.6, # tau
	method      => 'kendall',
	'p.value'   => 0.2333333,
	statistic   => 8
},
{# cor.test(cx, cy, alternative='two.sided', method = 'pearson', continuity=1)
	alternative => 'two.sided',
	'conf.int' => [
	 -0.279637693499009, 0.986196123450776
	],
	estimate  =>  0.8,
	method    => 'pearson',
	'p.value' => 0.104088,
	parameter => 3,
	statistic => 2.309401
},
{ # cor.test(cx, cy, alternative='less', method = 'pearson', continuity=1
	alternative => 'less',
	'conf.int' => [
	 -1, 0.9785289
	],
	estimate  =>  0.8,
	method    => 'pearson',
	'p.value' => 0.947956,
	parameter => 3,
	statistic => 2.309401	
},
{# cor.test(cx, cy, alternative='less', method = 'kendall', continuity=1)
	alternative => 'less',
	estimate    => 0.6, # tau
	method      => 'kendall',
	'p.value'   => 0.9583333,
	statistic   => 8
}
);
foreach my $meth (@correct) {
	my $result = cor_test(
		'x'         => $x,
		'y'         => $y,
		alternative => $meth->{alternative}, # so that it matches the test
		method      => $meth->{method},      # so that it matches the test
		continuity  => 1
	);
	my @undef_keys = grep {!defined $meth->{$_}} sort keys %{ $result };
	if (scalar @undef_keys > 0) {
		p @undef_keys;
		die "The above keys aren't defined";
	}
	foreach my $key (sort grep {looks_like_number($result->{$_})} keys %{ $result }) {
		is_approx( $result->{$key}, $meth->{$key}, "cor_test: $meth->{method}/$meth->{alternative} & $key");
	}
}
#--------------------
#  aov
#--------------------
#> yield <- c(5.5, 5.4, 5.8, 4.5, 4.8, 4.2) 
#> ctrl <- c(1, 1, 1, 0, 0, 0)
#> aov(yield ~ ctrl)
#Call:
#   aov(formula = yield ~ ctrl)

#Terms:
#                     ctrl Residuals
#Sum of Squares  1.7066667 0.2666667
#Deg. of Freedom         1         4

#Residual standard error: 0.2581989
#Estimated effects may be unbalanced
#> summary(aov(yield ~ ctrl))
#            Df Sum Sq Mean Sq F value  Pr(>F)   
#ctrl         1 1.7067  1.7067    25.6 0.00718 **
#Residuals    4 0.2667  0.0667                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
my $aov_res = 
aov(
	{
		yield => [5.5, 5.4, 5.8, 4.5, 4.8, 4.2],
		ctrl  => [1,     1,   1,   0,   0,   0]
	},
'yield ~ ctrl');
%correct = (
	ctrl => {
		Df        => 1,
		'F value' => 25.6000000000001,
		'Mean Sq' => 1.70666666666667,
		'Pr(>F)'  => 0.007182329,
		'Sum Sq'  => 1.70666666666667
	 },
	 Residuals => {
		Df        => 4,
		'Mean Sq' => 0.0666666666666729,
		'Sum Sq'  => 0.266666666666692
	 }
);
foreach my $k1 ('ctrl', 'Residuals') {
	foreach my $k2 (sort keys %{ $correct{$k1} }) {
		my $e = 10**-7;
		if ($correct{$k1}{$k2} =~ m/\.(\d+)$/) {
			$e = 10**(2-length $1);
		}
		is_approx( $aov_res->{$k1}{$k2}, $correct{$k1}{$k2}, "AOV: $k1/$k2");
	}
}
#-------------------------------------------------------------------
#  glm: Generalized Linear Models
#-------------------------------------------------------------------
subtest 'glm: Gaussian matches lm' => sub {
    # Check that gaussian glm is mathematically identical to OLS lm
    my $lm_res = lm(formula => 'mpg ~ wt + hp', data => $mtcars);
    my $glm_res = glm(formula => 'mpg ~ wt + hp', data => $mtcars, family => 'gaussian');
    
    is_approx($glm_res->{coefficients}{Intercept}, $lm_res->{coefficients}{Intercept}, 'glm gaussian matches lm intercept');
    is_approx($glm_res->{coefficients}{wt}, $lm_res->{coefficients}{wt}, 'glm gaussian matches lm wt');
    is_approx($glm_res->{deviance}, $lm_res->{rss}, 'glm gaussian deviance matches lm RSS');
    is($glm_res->{family}, 'gaussian', 'glm stored family correctly');
};

#-------------------------------------------------------------------
#  glm: Generalized Linear Models
#-------------------------------------------------------------------
subtest 'glm: Gaussian matches lm' => sub {
	# Check that gaussian glm is mathematically identical to OLS lm
	my $lm_res = lm(formula => 'mpg ~ wt + hp', data => $mtcars);
	my $glm_res = glm(formula => 'mpg ~ wt + hp', data => $mtcars, family => 'gaussian');

	is_approx($glm_res->{coefficients}{Intercept}, $lm_res->{coefficients}{Intercept}, 'glm gaussian matches lm intercept');
	is_approx($glm_res->{coefficients}{wt}, $lm_res->{coefficients}{wt}, 'glm gaussian matches lm wt');
	is_approx($glm_res->{deviance}, $lm_res->{rss}, 'glm gaussian deviance matches lm RSS');
	is($glm_res->{family}, 'gaussian', 'glm stored family correctly');
};

subtest 'glm: Binomial (Logistic Regression)' => sub {
	# Test dataset matching exact output from R's glm(am ~ wt + hp, data=mtcars, family=binomial)
	my $glm_bin = glm(formula => 'am ~ wt + hp', data => $mtcars, family => 'binomial');
	# 1. Convergence & Integrity
	ok($glm_bin->{converged}, 'glm binomial converged');
	ok($glm_bin->{iter} < 25, 'glm binomial converged under iteration limit');
	# 2. Coefficients (matching R precisely)
	my %correct_bin_coefs = (
	  Intercept => 18.86630,
	  wt        => -8.08348,
	  hp        =>  0.03626
	);
	foreach my $coef (keys %correct_bin_coefs) {
	  is_approx(
		   $glm_bin->{coefficients}{$coef},
		   $correct_bin_coefs{$coef},
		   "glm binomial matches R coef for $coef",
		   0.001
	  );
	}
	# 3. Summary standard errors & z-values
	is_approx($glm_bin->{summary}{Intercept}{'Std. Error'}, 7.44356, 'glm binomial Std. Error for Intercept', 0.001);
	is_approx($glm_bin->{summary}{wt}{'z value'}, -2.634, 'glm binomial z value for wt', 0.001);

	# 4. Deviance metrics
	is_approx($glm_bin->{deviance}, 10.059, 'glm binomial residual deviance', 0.001);
	is_approx($glm_bin->{'null.deviance'}, 43.229, 'glm binomial null deviance', 0.001);
	is($glm_bin->{'df.residual'}, 29, 'glm binomial residual degrees of freedom');
	is($glm_bin->{'df.null'}, 31, 'glm binomial null degrees of freedom');
};
my %tooth_growth = (
	dose => [qw(0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
1.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
2.0 2.0 2.0)],
	len  => [qw(4.2 11.5  7.3  5.8  6.4 10.0 11.2 11.2  5.2  7.0 16.5 16.5 15.2 17.3 22.5
17.3 13.6 14.5 18.8 15.5 23.6 18.5 33.9 25.5 26.4 32.5 26.7 21.5 23.3 29.5
15.2 21.5 17.6  9.7 14.5 10.0  8.2  9.4 16.5  9.7 19.7 23.3 23.6 26.4 20.0
25.2 25.8 21.2 14.5 27.3 25.5 26.4 22.4 24.5 24.8 30.9 26.4 27.3 29.4 23.0)],
	supp => [qw(VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC
VC VC VC VC VC OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ
OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ)]
);
%correct = (
	aic => 357.3958,
	coefficients => {
		dose      => 9.763571,
		Intercept => 7.4225
	},
	deviance        => 1227.905,
	'df.null'       => 59,
	'df.residual'   => 58,
	iter            => 2,
	'null.deviance' => 3452.209,
	rank            => 2,
	summary => {
		dose  => {
			Estimate     => 9.7636,
			'Pr(>|t|)'   => 1.23*10**-14,
			'Std. Error' => 0.952532934814911,
			't value'    => 10.250114270819
		},
		Intercept => {
			Estimate     => 7.4225,
			'Pr(>|t|)'   => 2.06e-07,
			'Std. Error' => 1.2601,
			't value'    => 5.89
		}
	}
);
foreach my ($idx, $val) (indexed qw(
-8.1042857 -0.8042857 -5.0042857 -6.5042857 -5.9042857 -2.3042857 -1.1042857 
-1.1042857 -7.1042857 -5.3042857 -0.6860714 -0.6860714 -1.9860714  0.1139286 
 5.3139286  0.1139286 -3.5860714 -2.6860714  1.6139286 -1.6860714 -3.3496429 
-8.4496429  6.9503571 -1.4496429 -0.5496429  5.5503571 -0.2496429 -5.4496429 
-3.6496429  2.5503571  2.8957143  9.1957143  5.2957143 -2.6042857  2.1957143 
-2.3042857 -4.1042857 -2.9042857  4.1957143 -2.6042857  2.5139286  6.1139286 
 6.4139286  9.2139286  2.8139286  8.0139286  8.6139286  4.0139286 -2.6860714 
10.1139286 -1.4496429 -0.5496429 -4.5496429 -2.4496429 -2.1496429  3.9503571 
-0.5496429  0.3503571  2.4503571 -3.9496429
)) {
	$correct{'deviance.resid'}{$idx+1} = $val;
}
my $glm_teeth = glm(
	data    => \%tooth_growth,
	formula => 'len ~ dose',
	family  => 'gaussian'
);

foreach my $term (sort keys %{ $correct{coefficients} }) {
	my $e = 10**-7;
	if ($correct{coefficients}{$term} =~ m/\.(\d+)$/) {
		$e = 10**(1-length $1);
	}
	is_approx( $glm_teeth->{coefficients}{$term}, $correct{coefficients}{$term}, "generalized Linear Models (glm) coefficients->$term within $e", $e);
}
foreach my $term (sort keys %{ $correct{summary} }) {
	foreach my $stat ('Estimate', 'Pr(>|t|)', 'Std. Error', 't value') {
		my $e = 10**-7;
		if ($correct{summary}{$term}{$stat} =~ m/\.(\d+)$/) {
			$e = 10**(1-length $1);
		} else {
			my $sp = sprintf '%.3g', $correct{summary}{$term}{$stat};
			if ($sp =~ m/e\-(\d+)$/) {
				$e = 10**(2-$1);
			} else {
				die "$sp failed regex.";
			}
		}
		is_approx( $glm_teeth->{summary}{$term}{$stat}, $correct{summary}{$term}{$stat}, "generalized Linear Models (glm) coefficients->$term/$stat within $e", $e);
	}
}
foreach my $key (sort grep {ref $correct{$_} eq '' } keys %correct) {
	my $e;
	if ($correct{$key} =~ m/\.(\d+)$/) {
		$e = 10**(-length $1);
	} elsif ($correct{$key} =~ m/^\-?\d+$/) {
		$e = 10**-199;
	} else {
		my $sp = sprintf '%.3g', $correct{$key};
		if ($sp =~ m/e\-(\d+)$/) {
			$e = 10**(2-$1);
		} else {
			die "$sp failed regex.";
		}
	}
	is_approx( $glm_teeth->{$key}, $correct{$key}, "$key within $e", $e);
}
foreach my $key (sort keys %{ $correct{'deviance.resid'} } ) {
	my $e;
	if ($correct{'deviance.resid'}{$key} =~ m/\.(\d+)$/) {
		$e = 10**(-length $1);
	} elsif ($correct{'deviance.resid'}{$key} =~ m/^\-?\d+$/) {
		$e = 10**-199;
	} else {
		my $sp = sprintf '%.3g', $correct{'deviance.resid'}{$key};
		if ($sp =~ m/e\-(\d+)$/) {
			$e = 10**(2-$1);
		} else {
			die "$sp failed regex.";
		}
	}
	is_approx( $glm_teeth->{'deviance.resid'}{$key}, $correct{'deviance.resid'}{$key}, "deviance.resid $key within $e");
}
$glm_teeth = glm(
	data    => \%tooth_growth,
	formula => 'len ~ dose + supp',
	family  => 'gaussian'
);
%correct = (
	aic        => 348.41553291891,
	deviance   => 1022.5550357143,
	'df.residual' => 57,
#	dispersion => 17.93956,
	coefficients => {
		dose      => 9.763571,
		Intercept => 9.272500,
		suppVC    => -3.700000
	},
	summary => {
		dose  => {
			Estimate     => 9.763571,
			'Pr(>|t|)'   => 6.313519e-16,
			'Std. Error' => 0.8768343,
			't value'    => 11.135025
		},
		Intercept => {
			Estimate     => 9.272500,
			'Pr(>|t|)'   => 1.312335e-09,
			'Std. Error' => 1.2823649,
			't value'    => 7.230781
		},
		suppVC    => {
			Estimate     => -3.700000,
			'Pr(>|t|)'   => 1.300662e-03,
			'Std. Error' => 1.0936045,
			't value'    => -3.383307
		}
	}
);
foreach my ($idx, $val) (indexed qw(
-6.2542857  1.0457143 -3.1542857 -4.6542857 -4.0542857 -0.4542857  0.7457143 
 0.7457143 -5.2542857 -3.4542857  1.1639286  1.1639286 -0.1360714  1.9639286 
 7.1639286  1.9639286 -1.7360714 -0.8360714  3.4639286  0.1639286 -1.4996429 
-6.5996429  8.8003571  0.4003571  1.3003571  7.4003571  1.6003571 -3.5996429 
-1.7996429  4.4003571  1.0457143  7.3457143  3.4457143 -4.4542857  0.3457143 
-4.1542857 -5.9542857 -4.7542857  2.3457143 -4.4542857  0.6639286  4.2639286 
 4.5639286  7.3639286  0.9639286  6.1639286  6.7639286  2.1639286 -4.5360714 
 8.2639286 -3.2996429 -2.3996429 -6.3996429 -4.2996429 -3.9996429  2.1003571 
-2.3996429 -1.4996429  0.6003571 -5.7996429)) {
	$correct{'deviance.resid'}{$idx+1} = $val;
}
foreach my $term (sort keys %{ $correct{coefficients} }) {
	my $e = 10**-7;
	if ($correct{coefficients}{$term} =~ m/\.(\d+)$/) {
		$e = 10**(-length $1);
	} elsif ($correct{coefficients}{$term} =~ m/^\-?\d+$/) {
		$e = 10**-199;
	} else {
		my $sp = sprintf '%.3g', $correct{coefficients}{$term};
		if ($sp =~ m/e\-(\d+)$/) {
			$e = 10**(2-$1);
		} else {
			die "$sp failed regex.";
		}
	}
	is_approx( $glm_teeth->{coefficients}{$term}, $correct{coefficients}{$term}, "generalized Linear Models (glm) coefficients->$term = $glm_teeth->{coefficients}{$term} within $e of $correct{coefficients}{$term}", $e);
}
foreach my $key (sort grep {ref $correct{$_} eq '' } keys %correct) {
	my $e;
	if ($correct{$key} =~ m/\.(\d+)$/) {
		$e = 10**(-length $1);
	} elsif ($correct{$key} =~ m/^\-?\d+$/) {
		$e = 10**-199;
	} else {
		my $sp = sprintf '%.3g', $correct{$key};
		if ($sp =~ m/e\-(\d+)$/) {
			$e = 10**(2-$1);
		} else {
			die "$sp failed regex.";
		}
	}
	is_approx( $glm_teeth->{$key}, $correct{$key}, "$key within $e of $correct{$key}", $e);
}
foreach my $term (sort keys %{ $correct{summary} }) {
	foreach my $stat ('Estimate', 'Pr(>|t|)', 'Std. Error', 't value') {
		my $e = 10**-7;
		if ($correct{summary}{$term}{$stat} =~ m/\.(\d+)$/) {
			$e = 10**(1-length $1);
		} else {
			my $sp = sprintf '%.3g', $correct{summary}{$term}{$stat};
			if ($sp =~ m/e\-(\d+)$/) {
				$e = 10**(2-$1);
			} else {
				die "$sp failed regex.";
			}
		}
		is_approx( $glm_teeth->{summary}{$term}{$stat}, $correct{summary}{$term}{$stat}, "glm coefficients->$term/$stat: $glm_teeth->{summary}{$term}{$stat}, $correct{summary}{$term}{$stat}", $e);
	}
}
foreach my $key (sort keys %{ $correct{'deviance.resid'} } ) {
	my $e;
	if ($correct{'deviance.resid'}{$key} =~ m/\.(\d+)$/) {
		$e = 10**(-(length $1));
		say __LINE__;
	} elsif ($correct{'deviance.resid'}{$key} =~ m/^\-?\d+$/) {
		$e = 10**-199;
	} else {
		say __LINE__;
		my $sp = sprintf '%.3g', $correct{'deviance.resid'}{$key};
		if ($sp =~ m/e\-(\d+)$/) {
			$e = 10**(2-$1);
		} else {
			die "$sp failed regex.";
		}
	}
	is_approx( $glm_teeth->{'deviance.resid'}{$key}, $correct{'deviance.resid'}{$key}, "deviance.resid $key within $e");
}
done_testing();
