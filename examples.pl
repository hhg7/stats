#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':all';
use Stats::LikeR;
#use JSON 'encode_json';
#use Matplotlib::Simple;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';

my $t = t_test(
	'x' => [3,3,3,3]
);
p $t;
=my $mtcars = json_file_to_ref('mtcars.json');
my %mtcars_hoa;
foreach my $row (@{ $mtcars }) {
	my $car = delete $row->{_row};
	push @{ $mtcars_hoa{car} }, $car;
	foreach my $key (keys %{ $row }) {
		push @{ $mtcars_hoa{$key} }, $row->{$key};
	}
}
p %mtcars_hoa;
ref_to_json_file(\%mtcars_hoa, 'mtcars.hoa.json' );
=my $mtcars = json_file_to_ref('mtcars.hoh.json');
my $lm = lm(formula =>  'mpg ~ wt + hp', data => $mtcars);
p $lm;
printf('%.3g' . "\n", $lm->{summary}{hp}{'Pr(>|t|)'});

=my $x = [1, 2, 3, 4, 5];
my $y = [2, 1, 4, 3, 5];

my @correct = (
{ # cor.test(cx, cy, alternative='two.sided', method = 'spearman', continuity=1)
	estimate  => 0.8,
	'p.value' => 0.1333333,
	statistic => 4,
},
{ #cor.test(cx, cy, alternative='two.sided', method = 'kendall', continuity=1)
	estimate  => 0.6, # tau
	'p.value' => 0.2333333,
	statistic => 
#	statistic, 
},
{	 alternative => 'two.sided',
    'conf.int' => [
       -0.279637693499009, 0.986196123450776
    ],
    estimate  =>  0.8,
    method    => 'pearson',
    p.value   => 0.104088,
    parameter => 3,
    statistic => 2.309401

);
foreach my ($idx, $meth) (indexed('spearman', 'kendall', 'pearson')) {
	my $result = cor_test(
		'x'         => $x,
		'y'         => $y,
		alternative => 'two.sided',
		method      => $meth,
		continuity  => 1
	);
	foreach my $key (sort keys %{ $result }) {
	
	}
	p $result;
}
=my $m = matrix(data => [1..9], nrow => 3, ncol => 3);
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
my $binom = rbinom( n => 99, prob => 0.5, size => 1);
p $binom, array_max => scalar @{ $binom };
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
