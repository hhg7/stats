#!/usr/bin/env perl
# read_table has two ways of turning a parsed CSV into the shape the caller
# asked for, and this pins them to each other.
#
# The default is the fast path: once the header is fixed, _parse_csv_file is
# handed a plan (csv_plan in LikeR.xs) and assembles every remaining row in C.
# The older path hands each row to a perl closure instead, and is still what
# runs for a 'filter' (and, until 0.319, for 'hoh'). So every case below is read
# twice -- once as-is, and once with a filter that accepts every row, which is
# the shortest way to force the closure -- and the two results must be
# identical. A no-op filter cannot change what comes back: a key of 0 gives the
# callback the whole row, and read_table only writes a mutated $_ back for keys
# above 0.
#
# There is no R or SciPy suite to take these from: output.type, the filter
# callbacks and auto.row.names are this module's own surface, and what is being
# checked is that two of its internal paths agree, not that either matches a
# reference. The values themselves are the fixtures written below.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use File::Temp ();
use File::Spec ();
use Stats::LikeR qw(read_table);

my $HAVE_LEAKTRACE;
BEGIN {
	$HAVE_LEAKTRACE = eval {
		require Test::LeakTrace;
		Test::LeakTrace->import('no_leaks_ok');
		1;
	} ? 1 : 0;
}

my $dir = File::Temp::tempdir(CLEANUP => 1);

sub fixture {
	my ($name, $text) = @_;
	my $path = File::Spec->catfile($dir, $name);
	open my $fh, '>', $path or die "cannot write \"$path\": $!\n";
	print {$fh} $text;
	close $fh or die "cannot close \"$path\": $!\n";
	return $path;
}

# A row count above one field-buffer's worth of rows is not needed: the plan is
# picked up after the first data row, so anything with two or more data rows
# exercises it. Three are used so that "later values win" and the ragged-row
# message both have a row to be wrong about.
my %f = (
	plain   => fixture('plain.csv',   "a,b\n1,2\n3,4\n5,6\n"),
	dup     => fixture('dup.csv',     "a,b,a\n1,2,3\n4,5,6\n7,8,9\n"),
	empty   => fixture('empty.csv',   "a,b\n1,\n,2\n,\nNA,x\n"),
	quoted  => fixture('quoted.csv',
		qq{a,b,c\n"x,1","he said ""hi""",3\n"multi\nline",2,3\n}),
	hdronly => fixture('hdronly.csv', "a,b\n"),
	onerow  => fixture('onerow.csv',  "a,b\n1,2\n"),
	autorn  => fixture('autorn.csv',  "x,y\nA,1,2\nB,3,4\nC,5,6\n"),
	cmt     => fixture('cmt.tsv',     "# a\tb\n1\t2\n3\t4\n5\t6\n"),
	tabs    => fixture('tabs.tsv',    "a\tb\n1\t2\n3\t4\n"),
	crlf    => fixture('crlf.csv',    "a,b\r\n1,2\r\n3,4\r\n"),
	onecol  => fixture('onecol.csv',  "a\n1\n2\n3\n"),
	ragged  => fixture('ragged.csv',  "a,b,c\n1,2,3\n4,5\n6,7,8\n"),
	duprn   => fixture('duprn.csv',   "id,v,w\nA,1,2\nB,3,4\nA,5,\n"),
	undefrn => fixture('undefrn.csv', "id,v\nA,1\n,2\nC,3\n"),
	nark    => fixture('nark.csv',    "id,v\nA,1\nNA,2\n"),
);

# Each case is (label, file, extra read_table options).
my @cases = (
	[ 'plain',                   $f{plain},   [] ],
	[ 'duplicate column names',  $f{dup},     [] ],
	[ 'empty fields',            $f{empty},   [] ],
	[ 'na.strings',              $f{empty},   [ 'na.strings' => [ 'NA', 'x' ] ] ],
	[ 'na.strings, one string',  $f{empty},   [ 'na.strings' => 'NA' ] ],
	[ 'quoted, embedded comma',  $f{quoted},  [] ],
	[ 'header only',             $f{hdronly}, [] ],
	[ 'one data row',            $f{onerow},  [] ],
	[ 'auto.row.names',          $f{autorn},  [ 'auto.row.names' => 1 ] ],
	[ 'auto.row.names, named',   $f{autorn},  [ 'auto.row.names' => 'rn' ] ],
	[ 'commented-out header',    $f{cmt},     [] ],
	[ 'tab separated',           $f{tabs},    [] ],
	[ 'explicit sep',            $f{tabs},    [ sep => "\t" ] ],
	[ 'CRLF line endings',       $f{crlf},    [] ],
	[ 'single column',           $f{onecol},  [] ],
);

# A hoh keys every row by its row.names column, so the cases that only make
# sense as one -- a repeated row name, a missing one, one that na.strings turns
# into a missing one, and an explicit row.names -- are added here.
push @cases,
	[ 'repeated row name',        $f{duprn},   [] ],
	[ 'missing row name',         $f{undefrn}, [] ],
	[ 'na.strings row name',      $f{nark},    [ 'na.strings' => 'NA' ] ],
	[ 'explicit row.names',       $f{duprn},   [ 'row.names' => 'w' ] ],
	[ 'row.names, repeated col',  $f{dup},     [ 'row.names' => 'a' ] ];

# One read, as (result or undef, $@ or '', the warnings it raised). The
# duplicate-column warning is left out: it is raised once per read, from perl,
# whichever path builds the rows, and t/read_table.t already pins it.
sub read_all {
	my @args = @_;
	my @warn;
	local $SIG{__WARN__} = sub {
		push @warn, $_[0] unless $_[0] =~ /duplicate column name/;
	};
	my $r = eval { read_table(@args) };
	return [ $r, $@, \@warn ];
}

for my $case (@cases) {
	my ($label, $file, $opts) = @$case;
	for my $otype (qw(aoh hoa hoh)) {
		my $fast   = read_all($file, @$opts, 'output.type' => $otype);
		my $closed = read_all($file, @$opts, 'output.type' => $otype,
			filter => { 0 => sub { 1 } });
		is_deeply $fast, $closed,
			"$label ($otype): both paths agree, errors and warnings included";
	}
}

# The hoh cases above agree with each other; these pin what they agree ON, so
# that both paths breaking the same way cannot pass.
{
	my $r = read_all($f{duprn}, 'output.type' => 'hoh');
	is_deeply $r->[0], { A => { v => 5, w => undef }, B => { v => 3, w => 4 } },
		'hoh: a repeated row name keeps the later values';
	is_deeply $r->[2],
		[ "read_table: duplicate row name 'A' in $f{duprn} (later values win)\n" ],
		'hoh: and warns once, naming it';
	$r = read_all($f{undefrn}, 'output.type' => 'hoh');
	is $r->[1],
		"read_table: undefined row name (column 'id') in $f{undefrn} data row 2\n",
		'hoh: a missing row name dies, naming the column and the row';
	$r = read_all($f{nark}, 'output.type' => 'hoh', 'na.strings' => 'NA');
	like $r->[1], qr/^read_table: undefined row name \(column 'id'\) .* data row 2$/,
		'hoh: so does one that na.strings makes missing';
	$r = read_all($f{nark}, 'output.type' => 'hoh');
	is_deeply $r->[0], { A => { v => 1 }, NA => { v => 2 } },
		'hoh: without na.strings, "NA" is an ordinary row name';
}

# The alignment message is produced by whichever path is reading, so the two
# have to spell it the same way -- including the data row number, which the
# fast path continues from wherever the closure stopped rather than restarting.
for my $otype (qw(aoh hoa hoh)) {
	my $fast   = eval { read_table($f{ragged}, 'output.type' => $otype); 1 }
		? '' : $@;
	my $closed = eval { read_table($f{ragged}, 'output.type' => $otype,
		filter => { 0 => sub { 1 } }); 1 } ? '' : $@;
	like $fast,
		qr/\AAlignment error on \Q$f{ragged}\E data row 2 \(2 fields vs 3 headers\)\.$/m,
		"ragged row ($otype): fast path names the row and both counts";
	is $fast, $closed, "ragged row ($otype): both paths word it the same";
}

# undef, not "", is what an empty cell becomes -- the one value transformation
# the fast path makes, so check it by identity rather than through is_deeply.
{
	my $aoh = read_table($f{empty});
	ok !defined $aoh->[0]{b}, 'empty cell is undef in an aoh';
	my $hoa = read_table($f{empty}, 'output.type' => 'hoa');
	ok !defined $hoa->{b}[0], 'empty cell is undef in a hoa';
	my $na = read_table($f{empty}, 'na.strings' => 'NA');
	ok !defined $na->[3]{a}, 'an na.strings cell is undef too';
	is $na->[3]{b}, 'x', 'a cell that is not in na.strings is left alone';
}

# The plan is installed once and then never consulted again, so reading the
# same file twice through the same process has to give the same thing.
{
	my $first  = read_table($f{plain});
	my $second = read_table($f{plain});
	is_deeply $second, $first, 'a second read of the same file is unchanged';
}

SKIP: {
	skip 'Test::LeakTrace not installed', 8 unless $HAVE_LEAKTRACE;
	skip 'running under Devel::Cover', 8 if $INC{'Devel/Cover.pm'};

	no_leaks_ok { read_table($f{plain}) } 'no leaks: aoh fast path';
	no_leaks_ok { read_table($f{plain}, 'output.type' => 'hoa') }
		'no leaks: hoa fast path';
	no_leaks_ok { read_table($f{empty}, 'na.strings' => [ 'NA', 'x' ]) }
		'no leaks: na.strings';
	# a duplicate name leaves one parsed field with no column to go to, which
	# is the only place the fast path frees a cell rather than handing it on
	no_leaks_ok {
		local $SIG{__WARN__} = sub { };
		read_table($f{dup});
	} 'no leaks: duplicate column names';
	# the alignment croak unwinds out of the middle of a read
	no_leaks_ok { eval { read_table($f{ragged}) } }
		'no leaks: alignment error';
	no_leaks_ok { read_table($f{plain}, 'output.type' => 'hoh') }
		'no leaks: hoh fast path';
	# the row-name croak unwinds with the row still full of cells
	no_leaks_ok { eval { read_table($f{undefrn}, 'output.type' => 'hoh') } }
		'no leaks: missing row name';
	# a __WARN__ handler that dies on the repeated-row-name warning, which is
	# raised from inside the fast path
	no_leaks_ok {
		local $SIG{__WARN__} = sub { die $_[0] };
		eval { read_table($f{duprn}, 'output.type' => 'hoh') };
	} 'no leaks: dying on the repeated-row-name warning';
}

done_testing;
