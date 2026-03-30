#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use Test::More;
use Test::Exception; # die_ok
use Stats::LikeR;
use Util;
use JSON qw(decode_json encode_json);
# written with Gemini's help
#my ($simple_task, $log_write, $stopping, $dry_run, $overwrite) = (0,0,0,0,0);
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
# standard deviation
my $stdev = sd(2,4,4,4,5,5,7,9);
my $correct = 2.1380899352994;
if (abs($stdev - $correct) < 10**-14) {
	ok(1, 'stdev works');
} else {
	my $diff = $correct - $stdev;
	fail("stdev does not work, got $stdev with an error of $diff");
}
# -------------------------------
# Custom helper for floating-point comparisons
sub is_approx ($got, $expected, $test_name, $epsilon = 10**-7) {
	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
	foreach my ($i, $arg) (indexed ($got, $expected, $test_name)) {
		next if defined $arg;
		die "\$arg[$i] (see subroutine signature for name) isn't defined in $current_sub";
	}
	if (abs($got - $expected) < $epsilon) {
	  pass($test_name);
	  return 1;
	} else {
	  fail($test_name);
	  diag("         got: $got\n    expected: $expected");
	  return 0;
	}
}
# --- Tests ---
# 1. Valid Mathematical Outcomes
is_approx( pearson_r([1, 2, 3], [2, 4, 6]), 1, 
	'Perfect positive correlation (r = 1)' );

is_approx( pearson_r([1, 2, 3], [6, 4, 2]), -1, 
	'Perfect negative correlation (r = -1)' );

is_approx( pearson_r([-1, 0, 1, 0], [0, 1, 0, -1]), 0, 
	'Zero correlation (r = 0)' );

is_approx( pearson_r([10, 20, 30, 40, 50], [12, 24, 33, 38, 55]), 0.9857122,
	'Standard positive correlation calculation' );

## 2. Mathematical Edge Cases
is( pearson_r([1, 1, 1], [2, 4, 6]), undef, 
	'Zero variance in one array returns undef (division by zero)' );

# 3. Validation / Exception Handling
eval { pearson_r("string", [1, 2]) };
like( $@, qr/Both arguments must be array references/, 
	"Dies correctly when given a non-array reference" );

eval { pearson_r([1, 2, 3], [1, 2]) };
like( $@, qr/Arrays must have the same number of elements/, 
	'Dies correctly on mismatched array lengths' );

eval { pearson_r([1], [2]) };
like( $@, qr/Need at least 2 elements/, 
	'Dies correctly when given less than 2 elements' );

eval { pearson_r([], []) };
like( $@, qr/Need at least 2 elements/, 
	'Dies correctly when given empty arrays' );

my @test_data = (
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
#------------------
# t.test
#------------------
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

my %correct_q = (
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
		$error += abs($q[$q] - $correct_q{$method}[$q]);
	}
	if ($error < 10**-6) {
		ok(1, "$method works with cumulative error of $error");
	} else {
		fail("$method doesn't work for FDR correction with error = $error");
	}
}
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
#----------------------
#		cor
#----------------------
$test_data[0] = [1, 2, 3, 4, 5,  5, 6,  7,   8];
$test_data[1] = [2, 4, 6, 8, 10, 9, 12, 14, 16];
my %correct_cor = (
	pearson  => 0.9973649,	spearman => 0.9958246,	kendall  => 0.9860133
);
foreach my $method (sort keys %correct_cor) {
	is_approx(
		cor($test_data[0], $test_data[1], $method),
		$correct_cor{$method},
		"cor: method = \"$method\""
	);
}
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
p @scaled_results;
@correct_scaled = (-2..2);
foreach my ($i, $correct) (indexed @correct_scaled) {
	is_approx($scaled_results[$i], $correct, "index $i for scale");
}
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
p $lm;
my %correct_lm = (
	coefficients => {
		Intercept => 49.8084234287587,
		hp        => -0.120102090978019,
		wt        => -8.21662429724302,
		'wt:hp'   => 0.0278481483187383
	},
	'df.residual' => 28,
	'fitted.values' => {
		'Mazda RX4' => 23.09547, 		'Mazda RX4 Wag' => 21.78138,
		'Datsun 710'    => 25.58488, 		'Hornet 4 Drive' => 20.02924,
		'Hornet Sportabout' => 17.28996, 		Valiant             => 18.88542,
		'Duster 360'        => 15.40745, 		'Merc 240D'         => 21.65887,
		'Merc 450SE'        => 15.14994, 		'Merc 450SL'          => 16.23929,
		'Merc 450SLC'         =>16.07909, 		'Cadillac Fleetwood'  => 12.02179,
		'Lincoln Continental' =>11.89490, 		'Chrysler Imperial'   => 12.50221,
		'Fiat 128'            => 27.84866, 		'Honda Civic'   => 32.63195,
		'Toyota Corolla'   => 30.24587, 		'Toyota Corona'   => 24.56317,
		'Dodge Challenger'   => 17.57441, 		'AMC Javelin'   => 17.91776,
		'Camaro Z28'   => 15.03111, 		'Pontiac Firebird'   => 15.93596,
		'Fiat X1-9'   => 29.53900, 		'Porsche 914-2'   => 26.71871,
		'Lotus Europa'   => 28.56630, 		'Ford Pantera L'   =>15.36033,
		'Ferrari Dino'   => 19.52990, 		'Maserati Bora'   => 13.54587,
		'Volvo 142E'      => 22.31363
	},
	rank => 4,
		residuals  => {
		'AMC Javelin'        =>   -2.7177637422554,		'Cadillac Fleetwood' =>   -1.62178684578001,
		'Camaro Z28'         =>   -1.73111177599938, 	'Chrysler Imperial'  =>   2.19779322930961,
		'Datsun 710'         =>   -2.78487707944995,		'Dodge Challenger'   =>   -2.07441456805362,
		'Duster 360'         =>   -1.10744532497053,		'Ferrari Dino'       =>   0.17010189824968,
		'Fiat 128'           =>   4.55133689384449,		'Fiat X1-9'           =>  -2.23900443083033,
		'Ford Pantera L'        => 0.439669246713233,	'Honda Civic'         =>  -2.23195395366212,
		'Hornet 4 Drive'     =>   1.37075604153847,		'Hornet Sportabout'  =>   1.41004478703069,
		'Lincoln Continental' =>  -1.4949003236173,		'Lotus Europa'       =>   1.83369534357951,
		'Maserati Bora'      =>   1.45413280824033,		'Mazda RX4'           =>  -2.09547410785999,
		'Mazda RX4 Wag'       =>  -0.781375472403504,	'Merc 230'            =>  1.95008336608675,
		'Merc 240D'           =>  2.74113094560431,		'Merc 280'            =>  0.646212827429729,
		'Merc 280C'           =>  -0.75378717257027,		'Merc 450SE'          =>  1.25006037875687,
		'Merc 450SL'          =>  1.06071479480091,		'Merc 450SLC'         =>  -0.879087325205568,
		'Pontiac Firebird'    =>  3.26404011532368,		'Porsche 914-2'       =>  -0.718705557249944,
		'Toyota Corolla'      =>  3.65413017953585,		'Toyota Corona'       =>  -3.06317321493851,
		Valiant               =>  -0.785416091802773,	'Volvo 142E'          =>  -0.913625869362747
  }
);
foreach my $key ('Intercept', 'hp', 'wt', 'wt:hp') {
	unless (defined $lm->{coefficients}{$key}) {
		p $lm;
		die "\"$key\" isn't defined" ;
	}
	is_approx( $lm->{coefficients}{$key}, $correct_lm{coefficients}{$key}, "Checking lm's $key" );
}
foreach my $key ('df.residual', 'rank') {
	unless (defined $lm->{$key}) {
		p $lm;
		die "\"$key\" isn't defined" ;
	}
	is_approx( $lm->{$key}, $correct_lm{$key}, "Checking \"$key\"");
}
foreach my $key ('fitted.values', 'residuals') {
	foreach my $car (keys %{ $correct_lm{$key} }) {
		unless (defined $lm->{$key}{$car}) {
			p $lm;
			die "\"$car\" isn't defined in \"fitted.values\"" ;
		}
		is_approx(
			$lm->{$key}{$car},
			$correct_lm{$key}{$car},
			"Checking $key \"$car\"",
			10**-5
		);
	}
}
$lm = lm(formula =>  'mpg ~ wt + hp', data => $mtcars);
p $lm;
my ($rmean, $sd, $n) = (10, 2, 9999);
my $normals = rnorm( n => $n, mean => $rmean, sd => $sd);
is_approx(scalar @{ $normals }, $n, 'rnorm sample size');
is_approx(mean($normals), $rmean, 'rnorm mean', 0.1);
is_approx(sd($normals), $sd, 'rnorm sd', 0.1);
done_testing();
