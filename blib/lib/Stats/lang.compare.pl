#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
#use Util;

my %d;
foreach my $lang ('perl', 'r', 'python') {
	my $data = read_table( "${lang}_benchmarks.tsv" );
	view($data);
	die;
}
