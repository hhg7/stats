#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

my %pg = (
	ctrl => [4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14],
	trt1 => [4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69],
	trt2 => [6.31,5.12,5.54,5.50,5.37,5.29,4.92,6.15,5.80,5.26]
#    weight => [
#        ,
#        ,
#        ,
#    ],
#    group => [ (('ctrl') x 10), (('trt1') x 10), (('trt2') x 10) ],
);

my $fit = aov(\%pg);#, 'weight ~ group');
my $hsd = TukeyHSD($fit, data => \%pg);#, formula => 'weight ~ group');
p $hsd;
view( $hsd->{group});
