#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR qw(!hist);
use Matplotlib::Simple 'violin';

my $mtcars = read_table(
	'mtcars.tsv',
	auto_row_names => 'model'
);
view($mtcars);
my $gb = group_by($mtcars, 'mpg', 'cyl');
view($gb);
violin(
	df => $gb,
	show => true,
);
my $anova = aov($gb);
p $anova;
