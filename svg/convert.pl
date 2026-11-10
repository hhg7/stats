#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Cn::Util;
use SimpleFlow 'parallel';

my @tasks;
foreach my $svg (dir(regex => '\.svg$')) {
	my $png = $svg;
	$png =~ s/svg$/png/;
	push @tasks, {
		cmd         => "convert $svg $png",
		output_file => $png,
		overwrite   => true
	};
}
parallel(
	tasks => \@tasks,
	jobs  => 4
);
