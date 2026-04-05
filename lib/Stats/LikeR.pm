#!/usr/bin/env perl
use 5.042.2;
no source::encoding;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
package Stats::LikeR;
our $VERSION = 0.01;
require XSLoader;
XSLoader::load('Stats::LikeR', $VERSION);
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use List::Util qw(min sum);
use Exporter 'import';
our @EXPORT_OK = qw(cor cor_test fisher_test hist lm matrix mean median min max p_adjust pearson_r quantile rbinom rnorm runif scale sd seq shapiro_test t_test var);
our @EXPORT = @EXPORT_OK;

#sub mean {
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 1 element in the array" if scalar @n < 1;
#	return sum(@n) / scalar @n;
#}

#sub stdev {
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 2 elements in the array" if scalar @n < 2;
#	my $mean = sum(@n) / scalar @n;
#	my $standard_deviation = 0;
#	foreach my $element (@n) {
#		$standard_deviation += ($element-$mean)**2;
#	}
#	return sqrt($standard_deviation/((scalar @n)-1));
#}

#sub population_variance { # sample variance
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 1 elements in the array" if scalar @n < 1;
#	my $mean = sum(@n) / scalar @n;
#	my $var = 0;
#	foreach my $element (@n) {
#		$var += ($element-$mean)**2;
#	}
#	return $var / scalar @n;
#}

#sub sample_variance { # sample variance
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 1 elements in the array" if scalar @n < 1;
#	my $mean = sum(@n) / scalar @n;
#	my $var = 0;
#	foreach my $element (@n) {
#		$var += ($element-$mean)**2;
#	}
#	return $var / (scalar @n - 1);
#}

#sub pvalue ($array1, $array2) {
#	return 1.0 if scalar @$array1 <= 1;
#	return 1.0 if scalar @$array2 <= 1;
#	my $mean1 = sum(@{ $array1 });
#	my $mean2 = sum(@{ $array2 });
#	return 1.0 if ($mean1 == $mean2);
#	$mean2 /= scalar @$array2;
#	$mean1 /= scalar @$array1;
#	my ($variance1, $variance2) = (0, 0);
#	foreach my $x (@$array1) {
#	$variance1 += ($x-$mean1)*($x-$mean1);
#	}
#	foreach my $x (@$array2) {
#	$variance2 += ($x-$mean2)*($x-$mean2);
#	}
#	if (($variance1 == 0.0) && ($variance2 == 0.0)) {
#	return 1.0;
#	}
#	$variance1 = $variance1/(scalar @$array1-1);
#	$variance2 = $variance2/(scalar @$array2-1);
#	my $array1_size = scalar @$array1;
#	my $array2_size = scalar @$array2;
#	my $WELCH_T_STATISTIC = ($mean1-$mean2)/sqrt($variance1/$array1_size+$variance2/$array2_size);
#	my $DEGREES_OF_FREEDOM = (($variance1/$array1_size+$variance2/(scalar @$array2))**2)
#	/
#	(
#	($variance1*$variance1)/($array1_size*$array1_size*($array1_size-1))+
#	($variance2*$variance2)/($array2_size*$array2_size*($array2_size-1))
#	);
#	my $A = $DEGREES_OF_FREEDOM/2;
#	my $value = $DEGREES_OF_FREEDOM/($WELCH_T_STATISTIC*$WELCH_T_STATISTIC+$DEGREES_OF_FREEDOM);
##from here, translation of John Burkhardt's C
#	my $beta = lgamma($A)+0.57236494292470009-lgamma($A+0.5);
#	my $acu = 10**(-15);
#	my($ai,$cx,$indx,$ns,$pp,$psq,$qq,$rx,$temp,$term,$xx);
## Check the input arguments.
#	return $value if $A <= 0.0;# || $q <= 0.0;
#	return $value if $value < 0.0 || 1.0 < $value;
## Special cases
#	return $value if $value == 0.0 || $value == 1.0;
#	$psq = $A + 0.5;
#	$cx = 1.0 - $value;
#	if ($A < $psq * $value) {
#		($xx, $cx, $pp, $qq, $indx) = ($cx, $value, 0.5, $A, 1);
#	} else {
#		($xx, $pp, $qq, $indx) = ($value, $A, 0.5, 0);
#	}
#	$term = 1.0;
#	$ai = 1.0;
#	$value = 1.0;
#	$ns = int($qq + $cx * $psq);
##Soper reduction formula.
#	$rx = $xx / $cx;
#	$temp = $qq - $ai;
#	$rx = $xx if $ns == 0;
#	while (1) {
#		$term = $term * $temp * $rx / ( $pp + $ai );
#		$value = $value + $term;
#		$temp = abs ($term);
#		if ($temp <= $acu && $temp <= $acu * $value) {
#	   	$value = $value * exp ($pp * log($xx)
#	                          + ($qq - 1.0) * log($cx) - $beta) / $pp;
#	   	$value = 1.0 - $value if $indx;
#	   	last;
#		}
#	 	$ai = $ai + 1.0;
#		$ns = $ns - 1;
#		if (0 <= $ns) {
#			$temp = $qq - $ai;
#			$rx = $xx if $ns == 0;
#		} else {
#			$temp = $psq;
#			$psq = $psq + 1.0;
#		}
#	}
#	return $value;
#}

#sub paired_pvalue ($array1, $array2) {
#	return 1.0 if scalar @$array1 <= 1;
#	return 1.0 if scalar @$array2 <= 1;
#	my $array_size = scalar @$array1;
#	if ($array_size != scalar @$array2) {
#		die 'arrays have different sizes, cannot calculate paired p-value.';
#	}
#	my ($sd, $xd) = (0,0);
#	foreach my $x (0..$array_size-1) {
#		$xd += (@$array2[$x] - @$array1[$x]);
#	}
#	return 1.0 if $xd == 0;
#	$xd /= $array_size;
#	foreach my $x (0..$array_size-1) {
#		$sd += ($xd - (@$array2[$x] - @$array1[$x]))**2;
#	}
#	$sd = sqrt($sd / ($array_size - 1));
#	my $t = $xd / ($sd / sqrt($array_size));
#	my $DEGREES_OF_FREEDOM = $array_size - 1;#http://oak.ucc.nau.edu/rh232/courses/EPS525/Handouts/Understanding%20the%20Dependent%20t%20Test.pdf
#	my $A = $DEGREES_OF_FREEDOM/2;
#	my $value = $DEGREES_OF_FREEDOM/($t*$t+$DEGREES_OF_FREEDOM);
##from here, translation of John Burkhardt's C
#	my $beta = lgamma($A)+0.57236494292470009-lgamma($A+0.5);
#	my $acu = 10**(-15);
#	my($ai,$cx,$indx,$ns,$pp,$psq,$qq,$rx,$temp,$term,$xx);
## Check the input arguments.
#	return $value if $A <= 0.0;# || $q <= 0.0;
#	return $value if $value < 0.0 || 1.0 < $value;
## Special cases
#	return $value if $value == 0.0 || $value == 1.0;
#	$psq = $A + 0.5;
#	$cx = 1.0 - $value;
#	if ($A < $psq * $value) {
#		($xx, $cx, $pp, $qq, $indx) = ($cx, $value, 0.5, $A, 1);
#	} else {
#		($xx, $cx, $pp, $qq, $indx) = ($value, $cx, $A, 0.5, 0);
#	}
#	$term = 1.0;
#	$ai = 1.0;
#	$value = 1.0;
#	$ns = int($qq + $cx * $psq);
##Soper reduction formula.
#	$rx = $xx / $cx;
#	$temp = $qq - $ai;
#	$rx = $xx if $ns == 0;
#	while (1) {
#		$term = $term * $temp * $rx / ( $pp + $ai );
#		$value = $value + $term;
#		$temp = abs ($term);
#		if ($temp <= $acu && $temp <= $acu * $value) {
#			$value = $value * exp ($pp * log($xx)
#				                   + ($qq - 1.0) * log($cx) - $beta) / $pp;
#			$value = 1.0 - $value if $indx;
#			last;
#		}
#		$ai = $ai + 1.0;
#		$ns = $ns - 1;
#		if (0 <= $ns) {
#			$temp = $qq - $ai;
#			$rx = $xx if $ns == 0;
#		} else {
#			$temp = $psq;
#			$psq = $psq + 1.0;
#		}
#	}
#	$value;
#}

1;
