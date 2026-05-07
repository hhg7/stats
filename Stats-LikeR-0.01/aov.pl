#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
use Time::HiRes;

my $data_interact = {
	'y' => [1, 2, 3, 4],
	'A' => ['a', 'b', 'a', 'b'],
	'B' => ['x', 'x', 'y', 'y']
};
my $t0 = Time::HiRes::time();
my $ft = aov($data_interact, 'y ~ A:B');
my $t1 = Time::HiRes::time();
printf("aov calculation in %g seconds.\n", $t1-$t0);

