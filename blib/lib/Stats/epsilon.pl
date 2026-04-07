#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';

sub get_epsilon ($val) {
	my $scientific_notation = sprintf '%.3e', $val;
	say $scientific_notation;
}
foreach my $val (-0.03177295,0.00145122851187347,0.0090297096758557,-3.518712, -3.87783074, 1.119647e-06, 0.63273349, -6.128695, 37.22727012, 2.565459e-20, 1.59878754, 23.284689) {
	say $val;
	get_epsilon( $val );
}
