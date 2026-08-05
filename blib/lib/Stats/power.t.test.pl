#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

h(*power_t_test);
my $power = power_t_test(
	n  => 30,	delta     => 0.5, 
	sd => 1.0,	sig_level => 0.05
);
p $power;
$power = power_t_test(
	power => 0.9,
	delta => 1
);
