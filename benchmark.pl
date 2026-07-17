#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use feature 'say';
use autodie ':all';
use Devel::Confess 'color';
use Time::HiRes;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Stats::LikeR;

my $r = read_table('xlsx.read.R.tsv');
my $p = read_table('xlsx.read.Python.tsv');

my @stats_likeR;
for (my $n = 0; $n < 7; $n++) {
	my $t0 = Time::HiRes::time();
	my $x = read_table('titanic.more.complete.xlsx');
	my $t1 = Time::HiRes::time();
	push @stats_likeR, $t1-$t0;
}
#-----------
my $aov = oneway_test({
	Perl => \@stats_likeR,
	R    => vals( $r, 'seconds' ),
	Python => vals( $p, 'seconds' )
});
p $aov;
