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
my @stat;
foreach my $func (@func) {
	my %val;
	foreach my $lang (@lang) {
		my $tbl = filter( $d{$lang}, col('Stats::LikeR function') eq $func );
		$val{$lang} = vals($tbl, 'time');
	}
	my $ow = oneway_test(\%val);
	my %result;
	foreach my $lang (@lang) {
		$result{"$lang time"} = $ow->{group_stats}{mean}{$lang};
#		foreach my $lang2 ( grep {$_ ne $lang} @lang) {
#			$result{"$lang $lang2 diff"} =
#				$ow->{group_stats}{mean}{$lang2}
#				-
#				$ow->{group_stats}{mean}{$lang}
#		}
	}
	$result{'Pr(>F)'} = $ow->{Group}{'Pr(>F)'};
	$result{'Func'} = $func;
	push @stat, \%result;
}
view(\@stat);
write_table(\@stat, 'lang.compare.xlsx');
=my $i = 0;
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
