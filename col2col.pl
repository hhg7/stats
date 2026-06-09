#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

my %data = (
	height => [ 170, 165, 180, 175 ],
	weight => [  70,  60,  85,  77 ],
	age    => [  30,  41,  25,  38 ]
);

my $result = col2col(\%data, 'cor');

p $result;
