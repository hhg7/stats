#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Stats::LikeR;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Time::HiRes;

my $mtcars = read_table('mtcars.tsv', 'output.type' => 'hoh', 'auto.row.names' => 'model');
my $lm_no_int = lm(formula => 'mpg ~ wt -1', data => $mtcars);
p $lm_no_int;
#	ok( !defined($lm_no_int->{coefficients}{Intercept}), 'lm: formula -1 correctly suppresses Intercept' );
my $t0 = Time::HiRes::time();
my $lm = lm(formula =>  'mpg ~ .', data => $mtcars);
my $t1 = Time::HiRes::time();
p $lm;
printf("lm ran in %.4g seconds.\n", $t1-$t0);
