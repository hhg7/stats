#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Cn::Util;
use SimpleFlow;

foreach my $svg (dir(regex => '\.svg$')) {
	my $png = $svg;
	$png =~ s/svg$/png/;
	task(
		cmd         => "convert $svg $png",
		output_file => $png
	);
}
