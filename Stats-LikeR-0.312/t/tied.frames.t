#!/usr/bin/env perl
#
# filter(), merge() and drop_duplicates() over frames whose columns and rows
# are tied.
#
# All three read a column or a row by walking AvARRAY() directly, which is the
# whole reason they are fast, and which is wrong for a tied array: its cells do
# not exist until FETCH has run and AvARRAY() on one is not the block the
# elements live in -- on an array that has never held a real element it is a
# null pointer.  av_fetch() is the way in, and what it hands back for a tied
# element is a mortal PVLV that only acquires the value once get magic has run
# on it, so SvOK() or SvROK() on one, tested before mg_get(), reports every
# cell of every tied frame as undef and every row as "not a reference".
# t/hot_path.t pins the same pair of mistakes for the reduction functions;
# this file pins them for the three frame functions.
#
# Up to 0.301 all three were wrong, silently and in four different directions:
#
#   * filter() on a tied HoA returned an empty frame whatever the predicate,
#     because flt_row_hoa() read the wrong block and saw no cell as numeric;
#     on a tied AoH it croaked "AoH element 0 is not a HASH reference".
#   * merge() on a tied HoA matched nothing, because mg_key() read every key
#     cell as undef and an undef key matches nothing by design -- so an inner
#     join came back empty and an outer join came back as two disjoint halves,
#     both of them well-formed and neither of them right.
#   * drop_duplicates() read every cell of a tied frame as undef, which made
#     every row identical and collapsed the frame to one row.
#   * drop_duplicates() on a tied AoH segfaulted, in _aoh_key_union, which
#     indexed AvARRAY() with no guard at all.
#
# The shape of every case is the same: run the call over tied input and over an
# identical plain-array copy, and require the two answers to be equal.  A fixed
# expected value would pin the plain answer as well, which t/filter.t, t/merge.t
# and t/drop_duplicates.t already do; what has to be pinned here is that the two
# routes agree.
#
# Not covered, because it is not supported: a tied *frame* hash (the outer hash
# of a HoA or HoH).  All three read its shape and its columns through
# HeVAL(hv_iternext(...)), which for a tied hash is not the value.  A tied
# column, a tied row array and a tied row hash are what this file is about.
#
# Tie::StdArray and Tie::StdHash are core (Tie::Array, Tie::Hash) and have been
# since well before 5.10, so nothing here is skippable.

require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use Stats::LikeR qw(merge filter col drop_duplicates);

# Imported at compile time so the (&;$) prototype is in scope for the block
# call at the end, as t/merge.t does.  Absent module -> that one test is
# skipped; everything above it still runs.
my $HAVE_LEAKTRACE;
BEGIN {
	$HAVE_LEAKTRACE = eval {
		require Test::LeakTrace;
		Test::LeakTrace->import('no_leaks_ok');
		1;
	};
}

{	package TiedArray;
	require Tie::Array;
	our @ISA = ('Tie::StdArray');
}
{	package TiedHash;
	require Tie::Hash;
	our @ISA = ('Tie::StdHash');
}

# A tied array holding the same elements as @$src.
sub ta {
	my ($src) = @_;
	my @a;
	tie @a, 'TiedArray';
	@a = @$src;
	return \@a;
}

# A tied hash holding the same pairs as %$src.
sub th {
	my ($src) = @_;
	my %h;
	tie %h, 'TiedHash';
	%h = %$src;
	return \%h;
}

# Canonical, order-independent text for a result frame of any shape, so two
# frames can be compared with one is().  Rows are sorted, because a join's row
# order is not what this file is testing.
sub sig {
	my ($df) = @_;
	my @rows;
	if (ref $df eq 'ARRAY') {			# AoH, or AoA read positionally
		@rows = map { my $r = $_;
		              ref $r eq 'ARRAY'
		                ? +{ map { ("[$_]" => $r->[$_]) } 0 .. $#$r }
		                : $r } @$df;
	} elsif (ref $df eq 'HASH') {
		my @cols = sort keys %$df;
		if (@cols && ref $df->{ $cols[0] } eq 'ARRAY') {	# HoA
			my $n = scalar @{ $df->{ $cols[0] } };
			@rows = map { my $i = $_; +{ map { $_ => $df->{$_}[$i] } @cols } }
			        0 .. $n - 1;
		} else {					# HoH: the row key is part of the row
			@rows = map { +{ %{ $df->{$_} }, '#row' => $_ } } @cols;
		}
	}
	return join "\n", sort map {
		my $r = $_;
		join '|', map { "$_=" . (defined $r->{$_} ? $r->{$_} : 'UNDEF') }
		          sort keys %$r;
	} @rows;
}

# The one assertion this file makes, over and over.
sub agree {
	my ($tied, $plain, $name) = @_;
	is sig($tied), sig($plain), $name;
}

# The data.  cat carries an undef, so the undef-key and undef-cell paths are
# exercised as well; id is deliberately not sorted and not gap-free.
my @id  = ( 3, 1, 4, 1, 5 );
my @x   = ( -2.5, 0.5, 7, -1, 0 );
my @cat = ( 'a', 'b', undef, 'b', 'a' );

my $plain = { id => [@id], x => [@x], cat => [@cat] };
my $tied  = { id => ta(\@id), x => ta(\@x), cat => ta(\@cat) };

# filter()
for my $out (undef, 'hoa', 'aoh') {
	my $label = defined $out ? $out : 'default';
	my @opt = defined $out ? ('output.type' => $out) : ();

	agree( filter($tied, col('x') > 0,      @opt),
	       filter($plain, col('x') > 0,     @opt), "filter col() numeric, $label" );
	agree( filter($tied, col('cat') eq 'b', @opt),
	       filter($plain, col('cat') eq 'b',@opt), "filter col() string, $label" );
	agree( filter($tied,  (col('x') > 0) | (col('id') == 1), @opt),
	       filter($plain, (col('x') > 0) | (col('id') == 1), @opt),
	       "filter col() boolean, $label" );
	agree( filter($tied,  sub { $_->{x} > 0 }, @opt),
	       filter($plain, sub { $_->{x} > 0 }, @opt), "filter coderef, $label" );
	agree( filter($tied,  sub { $_[1] < 2 },  @opt),
	       filter($plain, sub { $_[1] < 2 },  @opt), "filter on row index, $label" );
}

# The specific 0.301 symptom: not merely "the same as plain", but non-empty.
{
	my $got = filter($tied, col('x') > 0);
	is scalar @{ $got->{x} }, 2, 'filter over a tied HoA keeps the matching rows';
	is_deeply [ @{ $got->{id} } ], [ 1, 4 ], 'and the right ones, in input order';
}

# A tied AoH: the rows are a tied array, and one of them is a tied hash.
{
	my @rows = map { +{ id => $id[$_], x => $x[$_] } } 0 .. $#id;
	my $paoh = [ map { +{ %$_ } } @rows ];
	my $taoh = ta([ map { th($_) } @rows ]);
	agree( filter($taoh, col('x') > 0), filter($paoh, col('x') > 0),
	       'filter over a tied AoH of tied row hashes' );
	agree( filter($taoh, col('x') > 0, 'output.type' => 'hoa'),
	       filter($paoh, col('x') > 0, 'output.type' => 'hoa'),
	       'filter tied AoH -> HoA' );
}

# A HoH whose rows are tied hashes.  The frame hash itself is not tied: a tied
# *frame* hash is a separate matter and is not supported -- filter() reads the
# shape off HeVAL(hv_iternext(...)), which for a tied hash is not the value at
# all, and merge() reads its HoA columns the same way.  Nothing here pretends
# otherwise; the case that is fixed is a tied column or a tied row.
{
	my %prows = map { ("r$_" => { id => $id[$_], x => $x[$_] }) } 0 .. $#id;
	my $phoh = { map { ($_ => { %{ $prows{$_} } }) } keys %prows };
	my $thoh = { map { ($_ => th($prows{$_})) } keys %prows };
	agree( filter($thoh, col('x') > 0), filter($phoh, col('x') > 0),
	       'filter over a HoH of tied row hashes' );
}

# merge()
my @rid = ( 1, 4, 4, 9 );
my @w   = ( 'p', 'q', 'r', 's' );
my $rplain = { id => [@rid], w => [@w] };
my $rtied  = { id => ta(\@rid), w => ta(\@w) };

# id 4 appears twice on the right, so the chain through next[] is walked; id 9
# is right-only and id 3/5 are left-only, so every outer branch is reached.
for my $how (qw(inner left right outer)) {
	for my $out (qw(aoh hoa)) {
		agree( merge($tied,  $rtied,  how => $how, on => 'id', 'output.type' => $out),
		       merge($plain, $rplain, how => $how, on => 'id', 'output.type' => $out),
		       "merge $how, tied both sides, $out out" );
		agree( merge($tied,  $rplain, how => $how, on => 'id', 'output.type' => $out),
		       merge($plain, $rplain, how => $how, on => 'id', 'output.type' => $out),
		       "merge $how, tied left only, $out out" );
		agree( merge($plain, $rtied,  how => $how, on => 'id', 'output.type' => $out),
		       merge($plain, $rplain, how => $how, on => 'id', 'output.type' => $out),
		       "merge $how, tied right only, $out out" );
	}
}

# The specific 0.301 symptom, again as an absolute rather than a comparison.
{
	# left  ids 3 1 4 1 5, right ids 1 4 4 9, walked in left-row order:
	# 3 misses, 1 pairs once, 4 pairs with both right 4s, 1 pairs once, 5 misses.
	my $got = merge($tied, $rtied, how => 'inner', on => 'id');
	is scalar @{ $got->{id} }, 4, 'inner join over tied frames matches rows';
	is_deeply [ @{ $got->{id} } ], [ 1, 4, 4, 1 ],
		'and pairs each left row with every right row carrying its key';
}

# A multi-column key: mg_key() length-prefixes each cell only when there is
# more than one, so the two-key path renders its cells differently.
{
	my @k1 = qw(a a b b);
	my @k2 = ( 1, 2, 1, 2 );
	my @v  = qw(w x y z);
	my $lp = { k1 => [@k1], k2 => [@k2], v => [@v] };
	my $lt = { k1 => ta(\@k1), k2 => ta(\@k2), v => ta(\@v) };
	my @j1 = qw(a b b);
	my @j2 = ( 2, 1, 9 );
	my @u  = qw(P Q R);
	my $rp = { k1 => [@j1], k2 => [@j2], u => [@u] };
	my $rt = { k1 => ta(\@j1), k2 => ta(\@j2), u => ta(\@u) };
	for my $how (qw(inner outer)) {
		agree( merge($lt, $rt, how => $how, on => ['k1','k2']),
		       merge($lp, $rp, how => $how, on => ['k1','k2']),
		       "merge $how on two tied key columns" );
	}
}

# left.on / right.on, and a cross join, which takes no key at all.
{
	agree( merge($tied, $rtied,  how => 'left', 'left.on' => 'id', 'right.on' => 'id'),
	       merge($plain, $rplain, how => 'left', 'left.on' => 'id', 'right.on' => 'id'),
	       'merge left.on/right.on over tied frames' );
	agree( merge($tied,  $rtied,  how => 'cross'),
	       merge($plain, $rplain, how => 'cross'),
	       'cross join over tied frames' );
}

# A tied AoH and a tied HoH on either side of the join, so the row-frame
# branch of mg_cell() is covered too.
{
	my @lrows = map { +{ id => $id[$_], x => $x[$_] } } 0 .. $#id;
	my $laoh_p = [ map { +{ %$_ } } @lrows ];
	my $laoh_t = ta([ map { th($_) } @lrows ]);
	for my $how (qw(inner left outer)) {
		agree( merge($laoh_t, $rtied,  how => $how, on => 'id'),
		       merge($laoh_p, $rplain, how => $how, on => 'id'),
		       "merge $how with a tied AoH on the left" );
	}
	my %hrows = map { ("r$_" => { id => $id[$_], x => $x[$_] }) } 0 .. $#id;
	my $hoh_p = { map { ($_ => { %{ $hrows{$_} } }) } keys %hrows };
	my $hoh_t = { map { ($_ => th($hrows{$_})) } keys %hrows };
	agree( merge($hoh_t, $rtied,  how => 'inner', on => 'id'),
	       merge($hoh_p, $rplain, how => 'inner', on => 'id'),
	       'merge with a HoH of tied rows on the left' );
}

# drop_duplicates()
# Rows 0 and 3 are duplicates of each other; row 4 duplicates nothing.  That
# separates keep => 'first' from 'last' from 0, which is the whole of what the
# survivor bookkeeping does.
{
	my @dk  = ( 'a', 'b', 'a', 'a', 'c' );
	my @dv  = (  1,   2,   9,   1,   3  );
	my $dp  = { k => [@dk], v => [@dv] };
	my $dt  = { k => ta(\@dk), v => ta(\@dv) };
	for my $keep ('first', 'last', 0) {
		agree( drop_duplicates($dt, keep => $keep),
		       drop_duplicates($dp, keep => $keep),
		       "drop_duplicates HoA tied columns, keep => '$keep'" );
	}
	agree( drop_duplicates($dt, subset => 'k'), drop_duplicates($dp, subset => 'k'),
	       'drop_duplicates HoA tied columns, subset' );

	# the specific 0.301 symptom: every cell read as undef, so one row survived
	my $got = drop_duplicates($dt);
	is scalar @{ $got->{k} }, 4, 'drop_duplicates over a tied HoA keeps four rows';
	is_deeply [ @{ $got->{k} } ], [ qw(a b a c) ], 'and the right ones, in input order';

	# AoH: a tied outer array (which used to segfault in _aoh_key_union) and
	# tied row hashes.
	my @drows = map { +{ k => $dk[$_], v => $dv[$_] } } 0 .. $#dk;
	my $paoh  = [ map { +{ %$_ } } @drows ];
	for my $keep ('first', 'last', 0) {
		agree( drop_duplicates(ta([ map { +{ %$_ } } @drows ]), keep => $keep),
		       drop_duplicates($paoh, keep => $keep),
		       "drop_duplicates tied AoH outer array, keep => '$keep'" );
		agree( drop_duplicates([ map { th($_) } @drows ], keep => $keep),
		       drop_duplicates($paoh, keep => $keep),
		       "drop_duplicates AoH of tied row hashes, keep => '$keep'" );
	}

	# AoA: a tied outer array and tied row arrays.
	my @daoa = map { [ $dk[$_], $dv[$_] ] } 0 .. $#dk;
	my $paoa = [ map { [ @$_ ] } @daoa ];
	for my $keep ('first', 'last', 0) {
		agree( drop_duplicates(ta([ map { [ @$_ ] } @daoa ]), keep => $keep),
		       drop_duplicates($paoa, keep => $keep),
		       "drop_duplicates tied AoA outer array, keep => '$keep'" );
		agree( drop_duplicates([ map { ta($_) } @daoa ], keep => $keep),
		       drop_duplicates($paoa, keep => $keep),
		       "drop_duplicates AoA of tied row arrays, keep => '$keep'" );
	}
	agree( drop_duplicates([ map { ta($_) } @daoa ], subset => 0),
	       drop_duplicates($paoa, subset => 0),
	       'drop_duplicates AoA of tied row arrays, subset' );

	# An empty tied frame of each shape: _aoh_key_union indexed AvARRAY() with
	# no bounds check, and on an array that never held an element AvARRAY() is
	# a null pointer, so this is the case that crashed rather than lied.
	{
		my @none; tie @none, 'TiedArray';
		is_deeply drop_duplicates(\@none), [], 'drop_duplicates on an empty tied array';
		my %nocol; $nocol{k} = ta([]);
		is_deeply drop_duplicates(\%nocol), { k => [] },
			'drop_duplicates on a tied but empty column';
	}

	# What survives is shared, not copied -- except from a tied column, which
	# has no cell SV to share, so those are copied.  Both halves are behaviour,
	# not accident, so both are pinned.
	{
		my $plain = { k => [@dk], v => [@dv] };
		my $out   = drop_duplicates($plain);
		$out->{v}[0] = 'touched';
		is $plain->{v}[0], 'touched', 'a plain HoA survivor shares its cell with the input';

		my @tv = @dv;
		my $ti = { k => ta(\@dk), v => ta(\@tv) };
		my $to = drop_duplicates($ti);
		$to->{v}[0] = 'touched';
		is $ti->{v}[0], $dv[0], 'a tied HoA survivor is copied, so the input is untouched';
	}
}

# The input must come back untouched: FETCH is allowed, STORE is not.
{
	is_deeply [ @{ $tied->{id} } ],  [@id],  'tied id column unchanged by the calls above';
	is_deeply [ @{ $tied->{x} } ],   [@x],   'tied x column unchanged';
	is_deeply [ @{ $rtied->{id} } ], [@rid], 'tied right id column unchanged';
}

# No leaks.  av_row_keep() hands back a fresh reference for a tied row and the
# caller's own for a plain one, which is two ownership rules in one place.
SKIP: {
	skip 'Test::LeakTrace not installed', 1 unless $HAVE_LEAKTRACE;
	no_leaks_ok {
		drop_duplicates({ k => ta([qw(a b a)]), v => ta([1,2,1]) });
		drop_duplicates(ta([ { k => 'a' }, { k => 'a' } ]));
		drop_duplicates(ta([ ['a'], ['a'] ]));
		filter($tied, col('x') > 0);
		filter($tied, col('x') > 0, 'output.type' => 'aoh');
		filter($tied, sub { $_->{x} > 0 });
		merge($tied, $rtied, how => $_, on => 'id') for qw(inner left right outer);
		merge($tied, $rtied, how => 'cross');
	} 'no leaks over the tied paths';
}

done_testing();
