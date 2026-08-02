#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

my (%ow, %d);
my @lang = ('perl', 'r', 'python');
foreach my $lang (@lang) {
	$d{$lang} = read_table( "${lang}_benchmarks.tsv" );
}
#view($d);
my @func = uniq(vals($d{'r'}, 'Stats::LikeR function'));
p @func;
foreach my $func (@func) {
	my %val;
	foreach my $lang (@lang) {
		my $tbl = filter( $d{$lang}, col('Stats::LikeR function') eq $func );
		$val{$lang} = vals($tbl, 'time');
	}
	$ow{$func} = oneway_test(\%val);
}
my $i = 0;
foreach my $func (sort {
	$ow{$b}{group_stats}{mean}{perl}
#	$ow{$a}{'Group'}{'Pr(>F)'}
	<=>
	$ow{$a}{group_stats}{mean}{perl}
#	$ow{$b}{'Group'}{'Pr(>F)'}
} keys %ow) {
	say '--------------';
	say $func;
	say '--------------';
	p $ow{$func};
	$i++;
	last if $i == 4;
}
p $ow{t_test};
