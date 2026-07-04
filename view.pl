#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

my $aoh = [
 { row_name => 'p1', age => 41, sex => 'M', 'Testosterone, total (nmol/L)' => 18.2 },
 { row_name => 'p2', age =>  7, sex => 'F', 'Testosterone, total (nmol/L)' => undef },
 { row_name => 'p3', age => 33, sex => 'F', 'Testosterone, total (nmol/L)' => 1.05 },
 { row_name => 'p4', age => 55, sex => 'M', 'Testosterone, total (nmol/L)' => 22.9 },
 { row_name => 'p5', age => 29, sex => 'M', 'Testosterone, total (nmol/L)' => 14.0 },
 { row_name => 'p6', age => 62, sex => 'F', 'Testosterone, total (nmol/L)' => undef },
 { row_name => 'p7', age => 19, sex => 'M', 'Testosterone, total (nmol/L)' => 9.4 },
];

my $hoa = {
 row_name => [ map { "g$_" } 1..7 ],
 gene     => [ qw(BRCA1 TP53 EGFR KRAS MYC PTEN RB1) ],
 logfc    => [ -2.31, 0.04, 3.882901, undef, 1.2, -0.7, undef ],
};

my $hoh = {
 chrX => { start => 1000, end => 2000, strand => '+' },
 chrY => { start =>  500, end => 1500, strand => '-' },
 chr1 => { start =>   10, end =>  9999, strand => '+' },
};
my $aoa = [ [1, 2, 3], [4, 5, 6] ];
view( $aoa );
view( $aoh );
view( $hoa );
view( $hoh );
