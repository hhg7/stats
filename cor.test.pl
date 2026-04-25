#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

my $x_na = [1, 2,     3, undef,  5, 5, 6,     undef,7];
my $y_na = [2, undef, 6, 8,     10, 9, undef, 14,  16];
my $r = cor_test($x_na, $y_na);
p $r;
