#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
#use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

say mean([1e16, 1, -1e16, 1]);
say sum([1e16, 1, -1e16, 1]);
say mean(1e16, 1, -1e16, 1);
say sum(1e16, 1, -1e16, 1);
