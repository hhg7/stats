#!/usr/bin/env perl

use 5.042.1;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':all';
use stats;

#my $x = [1, 2, 3, 4, 5, 6];
#my $y = [2, 4, 6, 8, 10, 12];

# Welch's Two-Sample t-test (Default)
my $res2 = t_test({
	'x' => [27.5,21.0,19.0,23.6,17.0,17.9,16.9,20.1,21.9,22.6,23.1,19.6,19.0,21.7,21.4],
	'y' => [27.1,22.0,20.8,23.4,23.4,23.5,25.8,22.0,24.8,20.2,21.9,22.1,22.9,20.5,24.4],
});
p $res2;
