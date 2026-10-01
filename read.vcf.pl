#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Devel::Confess 'color';
use Stats::LikeR;

my $tbl = read_table(
	'/home/con/Documents/Genomics/56001801066929_WGZ.snp.vcf.gz',
);
view($tbl);
