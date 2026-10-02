#!/usr/bin/env perl
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Spec;
use File::Temp 'tempdir';
use Stats::LikeR;
use Test::More;
use Tie::Hash;

# A tied hash frame must give the same answer as the plain one it holds.
# (t/tied.frames.t covers tied arrays: columns and AoH rows.)
#
# Up to 0.3212 a tied hash -- the frame itself, or each row of a HoH -- made
# csort, value_counts, kruskal_test, merge, agg and drop_duplicates segfault,
# and made most of the others refuse valid data as empty or malformed: they
# read HeVAL() off hv_iternext(), which a tied iterator never fills in, tested
# hv_fetch() results before their get magic ran, and took key counts that read
# 0 for a tied hash. Each case here calls a function on plain data and on a
# tied copy of the same data and requires the same result. Every call is
# evaluated in this process, so a regression to a segfault fails the file.
#
# Results are compared with a numeric tolerance because a tied hash hands its
# keys over in another order, and a few functions sum across columns or groups
# in key order: 1e-12 relative is far above the last-bit differences that
# reordering a handful of additions can make, and far below any real error.
# Values invented for the test; R and SciPy have no tied hashes.

my $dir = tempdir( CLEANUP => 1 );

sub tied_copy {	# the same hash, tied (rows left as they are)
	my $h = shift;
	tie my %t, 'Tie::StdHash';
	%t = %$h;
	return \%t;
}
sub tied_deep {	# a HoH with the frame and every row tied
	my $h = shift;
	tie my %t, 'Tie::StdHash';
	%t = map { ( $_ => tied_copy( $h->{$_} ) ) } keys %$h;
	return \%t;
}
sub near {
	my ( $g, $e ) = @_;
	return 0 if ref $g ne ref $e;
	if ( ref $g eq 'HASH' ) {
		return 0 unless join( "\0", sort keys %$g ) eq join( "\0", sort keys %$e );
		near( $g->{$_}, $e->{$_} ) || return 0 for keys %$g;
		return 1;
	}
	if ( ref $g eq 'ARRAY' ) {
		return 0 unless @$g == @$e;
		near( $g->[$_], $e->[$_] ) || return 0 for 0 .. $#$g;
		return 1;
	}
	return !defined $e if !defined $g;
	return 0 if !defined $e;
	if ( $g =~ /^-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/ && $e =~ /^-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/ ) {
		return 1 if $g == $e;
		return abs( $g - $e ) <= 1e-12 * ( abs($e) > 1 ? abs($e) : 1 );
	}
	return "$g" eq "$e";
}
# TIED_HASHES_ONLY=regex runs only the cases whose name matches: with a case
# that segfaults still in the file, it is how the others can be run at all.
sub same {
	my ( $name, $call, $plain, @tied ) = @_;
	return if defined $ENV{'TIED_HASHES_ONLY'} && $name !~ /$ENV{'TIED_HASHES_ONLY'}/;
	my $want = eval { $call->($plain) };
	if ( !defined $want ) { fail("$name: the plain call itself died: $@"); return }
	for my $t (@tied) {
		my ( $label, $data ) = @$t;
		my $got = eval { $call->($data) };
		ok( defined $got && near( $got, $want ), "$name ($label)" ) or diag( $@ ? "died: $@" : explain( { got => $got, want => $want } ) );
	}
}

my %hoa = (
	'x'  => [ 1, 2, 3, 5, 4, 6 ],
	'y'  => [ 2, 1, 4, 3, 6, 5 ],
	'g'  => [qw(a a a b b b)],
	'id' => [ 1 .. 6 ],
);
my %hoh = map {
	my $i = $_;
	( "r$i" => { map { ( $_ => $hoa{$_}[ $i - 1 ] ) } keys %hoa } )
} 1 .. 6;
my %groups = ( 'a' => [ 1, 2, 3 ], 'b' => [ 4, 5, 7 ], 'c' => [ 9, 8, 11 ] );
my %counts = ( 'r1' => { 'a' => 10, 'b' => 20 }, 'r2' => { 'a' => 30, 'b' => 15 } );
my %ft     = ( 'x' => { 'a' => 3, 'b' => 1 }, 'y' => { 'a' => 1, 'b' => 3 } );
my %pv     = ( 'a' => [ 0.01, 0.04, 0.03 ], 'b' => [ 0.2, 0.5, 0.001 ] );

my @A = ( [ 'tied HoA', tied_copy( \%hoa ) ] );
my @H = ( [ 'tied HoH', tied_copy( \%hoh ) ], [ 'tied HoH, tied rows', tied_deep( \%hoh ) ] );

# the eight that segfaulted
same( 'csort HoA',             sub { csort( $_[0], 'x' ) }, \%hoa, @A );
same( 'csort HoH',             sub { csort( $_[0], 'x' ) }, \%hoh, @H );
same( 'value_counts HoA',      sub { value_counts( $_[0], 'g' ) }, \%hoa, @A );
same( 'value_counts HoH',      sub { value_counts( $_[0] ) }, \%hoh, @H );
same( 'kruskal_test',          sub { kruskal_test( $_[0] ) }, \%groups, [ 'tied', tied_copy( \%groups ) ] );
same( 'kruskal_test h =>',     sub { kruskal_test( 'h' => $_[0] ) }, \%groups, [ 'tied', tied_copy( \%groups ) ] );
same( 'merge HoA',             sub { merge( $_[0], { 'id' => [ 1, 2 ], 'z' => [ 5, 6 ] }, 'on' => 'id' ) }, \%hoa, @A );
same( 'merge HoH',             sub { merge( $_[0], { 'id' => [ 1, 2 ], 'z' => [ 5, 6 ] }, 'on' => 'id' ) }, \%hoh, @H );
same( 'agg HoA',               sub { agg( $_[0], 'by' => ['g'], 'agg' => { 'x' => 'mean' } ) }, \%hoa, @A );
same( 'drop_duplicates HoA',   sub { drop_duplicates( $_[0] ) }, \%hoa, @A );

# the ones that refused tied data
same( 'hoa2aoh',               sub { hoa2aoh( $_[0] ) }, \%hoa, @A );
same( 'hoa2hoh',               sub { hoa2hoh( $_[0], 'id' ) }, \%hoa, @A );
same( 'hoh2hoa',               sub { hoh2hoa( $_[0] ) }, \%hoh, @H );
same( 'col2col HoA',           sub { col2col( $_[0], 'cor', [ 'x', 'y' ] ) }, \%hoa, @A );
same( 'col2col HoH',           sub { col2col( $_[0], 'cor', [ 'x', 'y' ] ) }, \%hoh, @H );
same( 'lm HoA',                sub { lm( 'formula' => 'x ~ y', 'data' => $_[0] )->{'coefficients'} }, \%hoa, @A );
same( 'lm HoH',                sub { lm( 'formula' => 'x ~ y', 'data' => $_[0] )->{'coefficients'} }, \%hoh, @H );
same( 'glm HoA',               sub { glm( 'formula' => 'x ~ y', 'data' => $_[0] )->{'coefficients'} }, \%hoa, @A );
{
	my $fit = lm( 'formula' => 'x ~ y', 'data' => \%hoa );
	my %new = ( 'y' => [ 1, 2, 3 ] );
	same( 'predict newdata HoA', sub { predict( $fit, $_[0] ) }, \%new, [ 'tied', tied_copy( \%new ) ] );
	my %newh = ( 'p' => { 'y' => 1 }, 'q' => { 'y' => 2 } );
	same( 'predict newdata HoH', sub { predict( $fit, $_[0] ) }, \%newh, [ 'tied', tied_copy( \%newh ) ], [ 'tied, tied rows', tied_deep( \%newh ) ] );
}
same( 'aov',                   sub { aov( $_[0] ) }, \%groups, [ 'tied', tied_copy( \%groups ) ] );
same( 'oneway_test',           sub { oneway_test( $_[0] ) }, \%groups, [ 'tied', tied_copy( \%groups ) ] );
same( 'p_adjust HoA',          sub { p_adjust( $_[0], 'BH' ) }, \%pv, [ 'tied', tied_copy( \%pv ) ] );
same( 'chisq_test',            sub { chisq_test( $_[0] ) }, \%counts, [ 'tied', tied_copy( \%counts ) ], [ 'tied, tied rows', tied_deep( \%counts ) ] );
same( 'fisher_test',           sub { fisher_test( $_[0] ) }, \%ft, [ 'tied', tied_copy( \%ft ) ], [ 'tied, tied rows', tied_deep( \%ft ) ] );
same( 'sample (count)',        sub { scalar keys %{ sample( $_[0], 2 ) } }, \%groups, [ 'tied', tied_copy( \%groups ) ] );
same( 'group_by HoA',          sub { group_by( $_[0], 'x', 'g' ) }, \%hoa, @A );
same( 'group_by HoH',          sub { group_by( $_[0], 'x', 'g' ) }, \%hoh, @H );

# already right; pinned so they stay so
same( 'filter HoA',            sub { filter( $_[0], col('x') > 2 ) }, \%hoa, @A );
same( 'filter HoH',            sub { filter( $_[0], col('x') > 2 ) }, \%hoh, @H );
same( 'cfilter HoA',           sub { cfilter( $_[0], 'keep' => ['x'] ) }, \%hoa, @A );
same( 'transpose',             sub { transpose( $_[0] ) }, \%hoh, @H );
same( 'vals HoA',              sub { vals( $_[0], 'x' ) }, \%hoa, @A );
same( 'vals HoH',              sub { [ sort { $a <=> $b } @{ vals( $_[0], 'x' ) } ] }, \%hoh, @H );
same( 'select_cols HoA',       sub { select_cols( $_[0], 'x' ) }, \%hoa, @A );
{
	my $n = 0;
	my $wt = sub {
		my $f = File::Spec->catfile( $dir, 'wt' . $n++ . '.tsv' );
		write_table( $_[0], $f, 'quiet' => 1, 'row.names' => 0 );
		open my $fh, '<', $f or die "$f: $!";
		local $/;
		my $c = <$fh>;
		return $c;
	};
	same( 'write_table HoA', $wt, \%hoa, @A );
}

done_testing();
