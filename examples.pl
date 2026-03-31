#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':all';
use Stats::LikeR;
use JSON 'encode_json';
use Util;

my $m = matrix(data => [1..9], nrow => 3, ncol => 3);
p $m;
my $mtcars = json_file_to_ref('mtcars.json');
my (@mtcars, %mtcars);
foreach my $row (@{ $mtcars }) { # transform AoH to HoA
	my $row_name = delete $row->{_row};
	foreach my $key (keys %{ $row }) {
		push @{ $mtcars{$row_name}{$key} }, $row->{$key};
	}
}
unless (-f 'mtcars.hoh.json') {
	ref_to_json_file(\%mtcars, 'mtcars.hoh.json');
}
#my $lm = lm(formula => 'mpg ~ wt * hp^2', data => \%mtcars);
#say encode_json( \%mtcars );
# matrix
=my $mat1 = matrix(
	data => [1..6],
	nrow => 2
#	byrow => true
);
my $scale = scale($mat1);
p $scale;
say encode_json( $scale );
=say encode_json( $mat1 );
p $mat1;
# scale
my @x = 1..5;
my @scaled_x = scale($mat1);
p @scaled_x;
