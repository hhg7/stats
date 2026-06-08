#!/usr/bin/env perl
require 5.010;
use warnings FATAL => 'all';
use Stats::LikeR;
use Test::Exception; # dies_ok
use Test::More;
use Test::LeakTrace 'no_leaks_ok';

# cfilter selects columns (the inner/2nd-level keys of a HoH or AoH, or the
# outer keys of a HoA) and returns the data in the same shape. The selector is
# either keep => [names] / remove => [names], or keep/remove => a predicate
# (a CODE ref or function name) called per column as $pred->($col_values, $name)
# where $col_values is the column's defined cells. Exactly one of keep/remove.

# The same three-column table (z is constant -> sd 0) in all three shapes.
my %hoa = ( 'x' => [ 1, 2, 3 ], 'y' => [ 4, 5, 6 ], 'z' => [ 0, 0, 0 ] );
my %hoh = (
	'r1' => { 'x' => 1, 'y' => 4, 'z' => 0 },
	'r2' => { 'x' => 2, 'y' => 5, 'z' => 0 },
	'r3' => { 'x' => 3, 'y' => 6, 'z' => 0 },
);
my @aoh = (
	{ 'x' => 1, 'y' => 4, 'z' => 0 },
	{ 'x' => 2, 'y' => 5, 'z' => 0 },
	{ 'x' => 3, 'y' => 6, 'z' => 0 },
);

# 1. cfilter is defined.
ok( defined &Stats::LikeR::cfilter, 'cfilter is defined in Stats::LikeR' );

# 2. keep / remove by name on a hash of arrays (columns are the outer keys).
is_deeply( cfilter( \%hoa, 'keep' => [ 'x', 'y' ] ), { 'x' => [ 1, 2, 3 ], 'y' => [ 4, 5, 6 ] }, 'HoA: keep selects the named columns' );
is_deeply( cfilter( \%hoa, 'remove' => [ 'z' ] ), { 'x' => [ 1, 2, 3 ], 'y' => [ 4, 5, 6 ] }, 'HoA: remove drops the named column' );
is( ref cfilter( \%hoa, 'keep' => [ 'x' ] ), 'HASH', 'HoA stays a hash ref' );

# 3. keep / remove by name on a hash of hashes (columns are the inner keys).
is_deeply( cfilter( \%hoh, 'keep' => [ 'x' ] ), { 'r1' => { 'x' => 1 }, 'r2' => { 'x' => 2 }, 'r3' => { 'x' => 3 } }, 'HoH: keep trims each row to the named column' );
is_deeply( cfilter( \%hoh, 'remove' => [ 'y', 'z' ] ), { 'r1' => { 'x' => 1 }, 'r2' => { 'x' => 2 }, 'r3' => { 'x' => 3 } }, 'HoH: remove drops the named columns from every row' );

# 4. keep / remove by name on an array of hashes.
is_deeply( cfilter( \@aoh, 'keep' => [ 'x', 'z' ] ), [ { 'x' => 1, 'z' => 0 }, { 'x' => 2, 'z' => 0 }, { 'x' => 3, 'z' => 0 } ], 'AoH: keep trims each row hash' );
is( ref cfilter( \@aoh, 'keep' => [ 'x' ] ), 'ARRAY', 'AoH stays an array ref' );

# 5. A predicate selects columns by their values: keep the constant ones (sd 0).
my $constant = sub { sd( $_[0] ) == 0 };
is_deeply( [ sort keys %{ cfilter( \%hoa, 'keep' => $constant ) } ], [ 'z' ], 'predicate keep: only the sd==0 column (z) survives' );
is_deeply( [ sort keys %{ cfilter( \%hoa, 'remove' => $constant ) } ], [ 'x', 'y' ], 'predicate remove: the sd==0 column (z) is dropped' );
# The predicate works the same on the row-major shapes.
is_deeply( cfilter( \@aoh, 'keep' => $constant ), [ { 'z' => 0 }, { 'z' => 0 }, { 'z' => 0 } ], 'predicate keep on AoH keeps the constant column' );

# 6. The predicate's second argument is the column name.
my %named;
cfilter( \%hoa, 'keep' => sub { $named{ $_[1] } = scalar @{ $_[0] }; 1 } );
is_deeply( \%named, { 'x' => 3, 'y' => 3, 'z' => 3 }, 'predicate receives the column name and its values' );

# 7. A function name (with a package) is resolved and called as a predicate.
sub main::is_const { return sd( $_[0] ) == 0 }
is_deeply( [ sort keys %{ cfilter( \%hoa, 'keep' => 'main::is_const' ) } ], [ 'z' ], 'a function name works as a predicate' );

# 8. Only the column's DEFINED cells reach the predicate (undef/missing dropped).
my %gap = ( 'a' => [ 1, 2, undef, 4 ], 'b' => [ 1, 1, 1, 1 ] );
is_deeply( [ sort keys %{ cfilter( \%gap, 'keep' => sub { scalar( @{ $_[0] } ) == 4 } ) } ], [ 'b' ], 'predicate sees defined cells only (a has 3, b has 4)' );

# 9. Reconstruction preserves undef cells inside a kept column.
is_deeply( cfilter( \%gap, 'keep' => [ 'a' ] ), { 'a' => [ 1, 2, undef, 4 ] }, 'a kept column keeps its undef cells in place' );

# 10. Ragged rows: a kept column simply stays absent from rows that lack it.
is_deeply( cfilter( [ { 'a' => 1, 'b' => 2 }, { 'a' => 3 } ], 'keep' => [ 'b' ] ), [ { 'b' => 2 }, {} ], 'AoH: a row missing the kept column yields an empty row' );
is_deeply( cfilter( { 'r1' => { 'a' => 1, 'b' => 2 }, 'r2' => { 'a' => 3 } }, 'keep' => [ 'b' ] ), { 'r1' => { 'b' => 2 }, 'r2' => {} }, 'HoH: a row missing the kept column yields an empty row' );

# 11. The input is copied, never mutated.
my %orig = ( 'x' => [ 1, 2 ], 'y' => [ 3, 4 ] );
cfilter( \%orig, 'keep' => [ 'x' ] );
is_deeply( \%orig, { 'x' => [ 1, 2 ], 'y' => [ 3, 4 ] }, 'the original structure is left untouched' );

# 12. Bad inputs and bad selectors die with a clear message.
dies_ok { cfilter( \%hoa ) } 'no keep/remove dies';
dies_ok { cfilter( \%hoa, 'keep' => [ 'x' ], 'remove' => [ 'y' ] ) } 'both keep and remove dies';
dies_ok { cfilter( \%hoa, 'keep' => [ 'nope' ] ) } 'keeping an unknown column dies';
dies_ok { cfilter( \%hoa, 'remove' => [ 'nope' ] ) } 'removing an unknown column dies';
dies_ok { cfilter( \%hoa, 'keep' => {} ) } 'a hash-ref selector dies';
dies_ok { cfilter( \%hoa, 'keep' => 'no_such_function' ) } 'an unknown function name dies';
dies_ok { cfilter( \%hoa, 'bogus' => [ 'x' ] ) } 'an unknown option dies';
dies_ok { cfilter( \%hoa, 'keep' ) } 'an option without a value dies (odd args)';
dies_ok { cfilter( 42, 'keep' => [ 'x' ] ) } 'non-reference data dies';
dies_ok { cfilter( [ 1, 2, 3 ], 'keep' => [ 'x' ] ) } 'array of plain scalars dies';

# 13. No memory leaks across the name and predicate paths, all three shapes.
no_leaks_ok { cfilter( \%hoa, 'keep' => [ 'x', 'y' ] ) } 'no leaks: HoA keep by name' unless $INC{'Devel/Cover.pm'};
no_leaks_ok { cfilter( \%hoh, 'remove' => [ 'z' ] ) } 'no leaks: HoH remove by name' unless $INC{'Devel/Cover.pm'};
no_leaks_ok { cfilter( \@aoh, 'keep' => [ 'x' ] ) } 'no leaks: AoH keep by name' unless $INC{'Devel/Cover.pm'};
no_leaks_ok { cfilter( \%hoa, 'keep' => $constant ) } 'no leaks: predicate selector' unless $INC{'Devel/Cover.pm'};
use Devel::Peek;
END { Devel::Peek::Dump(1 == 1) }   # dumps PL_sv_yes
done_testing();
