#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Cn::Util;
use Stats::LikeR;

my %hoa = (
	'я' => [1..3],
	'ъ' => [4..9],
);
write_table(
	\%hoa, '/tmp/write.cols.xlsx',
	col_names => [qw(я ъ)],
	row_names => false
);
=my @aoh = [
	{
		a => 1,
		b => 2,
	},
	{
		a => 3,
		b => 4,
	}
];
write_table(
	\@aoh, '/tmp/hoa.xlsx',
	'row.names' => false
);
