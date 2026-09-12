#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
use Time::HiRes;

#my $titanic = read_table('titanic.csv');
#write_table(
#	$titanic,
#	'titanic.more.complete.xlsx',
#	'xlsx.freeze.rows' => 1
#);

my $t0 = Time::HiRes::time();
my $tmp = read_table('Affinity Dataset(main).xlsx');
my $t1 = Time::HiRes::time();
printf("read the affinity xlsx in %lf seconds.\n", $t1 - $t0);
view($tmp);
