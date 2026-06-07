#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use Devel::Confess;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Stats::LikeR;
use File::Temp;

my $table = read_table(
	't/HepatitisCdata.csv',
	'output.type' => 'hoh',
	filter => {
		Sex => sub {$_ eq 'f'}
	}
);
p $table;
my %t = (
	A => {
		a => 1,
		b => 2
	},
	B => {
		a => 2,
		c => 3
	}
);
my $fh = File::Temp->new( DIR => '/tmp', SUFFIX => '.tsv');
close $fh;
write_table(\%t, $fh->filename, 'undef.val' => '');
my $t = read_table($fh->filename);
p $t;
