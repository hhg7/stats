#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
#use Util 

my %h = (a => 1, b => 2, c => 3, d => 4);
my @arr = qw(apple banana cherry date elderberry);

foreach my $n (1..4) {
	my $n = sample(\%h, 1);
	p $n;
}
