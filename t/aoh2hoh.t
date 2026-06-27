#!/usr/bin/env perl
require 5.010;
use warnings FATAL => 'all';
use Scalar::Util 'reftype';
use Stats::LikeR;
use Test::Exception;	# dies_ok, throws_ok
use Test::More;
use Test::LeakTrace 'no_leaks_ok';

#--------
# basic: index by key, key column retained in inner hash
#--------
{
	my $in  = [ { id => 'a', x => 1 }, { id => 'b', x => 2 } ];
	my $out = aoh2hoh($in, 'id');
	is(reftype($out), 'HASH', 'aoh2hoh: returns a hashref');
	is_deeply(
		$out,
		{ a => { id => 'a', x => 1 }, b => { id => 'b', x => 2 } },
		'aoh2hoh: basic AoH indexed by id'
	);
}

#--------
# duplicate key value is fatal (primary-key enforcement)
#--------
{
	my $in = [ { id => 'a', x => 1 }, { id => 'a', x => 9 } ];
	throws_ok { aoh2hoh($in, 'id') } qr/duplicate key/,
		'aoh2hoh: duplicate key value dies';
}

#--------
# the die fires on the SECOND occurrence, so a single key never dies
#--------
{
	my $in = [ { id => 'a', x => 1 }, { id => 'b', x => 2 }, { id => 'c', x => 3 } ];
	lives_ok { aoh2hoh($in, 'id') } 'aoh2hoh: all-distinct keys do not die';
}

#--------
# rows skipped before the key check do NOT count toward a collision:
# two rows both missing the key are both skipped, not a duplicate
#--------
{
	my $in  = [ { id => 'a', x => 1 }, { x => 5 }, { x => 6 } ];
	my $out;
	lives_ok { $out = aoh2hoh($in, 'id') }
		'aoh2hoh: two keyless rows are skipped, not a duplicate';
	is_deeply(
		$out,
		{ a => { id => 'a', x => 1 } },
		'aoh2hoh: keyless rows contribute nothing'
	);
}

#--------
# non-hashref rows are skipped, not fatal
#--------
{
	my $in  = [ { id => 'a', x => 1 }, 42, [qw/not a hash/], { id => 'b', x => 2 } ];
	my $out = aoh2hoh($in, 'id');
	is_deeply(
		$out,
		{ a => { id => 'a', x => 1 }, b => { id => 'b', x => 2 } },
		'aoh2hoh: non-hashref rows skipped'
	);
}

#--------
# rows with missing / undef key value are skipped
#--------
{
	my $in  = [ { id => 'a', x => 1 }, { x => 5 }, { id => undef, x => 6 } ];
	my $out = aoh2hoh($in, 'id');
	is_deeply(
		$out,
		{ a => { id => 'a', x => 1 } },
		'aoh2hoh: missing/undef key rows skipped'
	);
}

#--------
# empty input -> empty hash
#--------
{
	my $out = aoh2hoh([], 'id');
	is_deeply($out, {}, 'aoh2hoh: empty arrayref -> empty hashref');
}

#--------
# shallow copy: mutating the output inner hash does not touch input
#--------
{
	my $in  = [ { id => 'a', x => 1 } ];
	my $out = aoh2hoh($in, 'id');
	$out->{a}{x} = 999;
	$out->{a}{new} = 'added';
	is($in->[0]{x}, 1, 'aoh2hoh: input scalar value unchanged after output edit');
	ok(!exists $in->[0]{new}, 'aoh2hoh: added key did not leak back into input');
}

#--------
# shallow copy: a value that is itself a ref is SHARED, like $h{$k} = $row->{$k}
#--------
{
	my $shared = [ 1, 2, 3 ];
	my $in     = [ { id => 'a', data => $shared } ];
	my $out    = aoh2hoh($in, 'id');
	push @{ $out->{a}{data} }, 4;
	is(scalar(@$shared), 4, 'aoh2hoh: nested ref value is shallow-copied (shared referent)');
}

#--------
# numeric and string key values collide (1 vs "1") -> duplicate -> dies
#--------
{
	my $in = [ { id => 1, x => 'int' }, { id => '1', x => 'str' } ];
	throws_ok { aoh2hoh($in, 'id') } qr/duplicate key/,
		'aoh2hoh: numeric/string key collision is a duplicate';
}

#--------
# utf8 / SV key safety on the key column itself
#--------
{
	my $k   = "\x{2603}";	# SNOWMAN
	my $in  = [ { id => $k, x => 1 } ];
	my $out = aoh2hoh($in, 'id');
	ok(exists $out->{$k}, 'aoh2hoh: utf8 key value preserved');
	is($out->{$k}{x}, 1, 'aoh2hoh: row reachable under utf8 key');
}

#--------
# argument validation
#--------
dies_ok { aoh2hoh({ not => 'an array' }, 'id') } 'aoh2hoh: dies on non-arrayref first arg';
dies_ok { aoh2hoh([], undef) }                   'aoh2hoh: dies on undef key';
throws_ok { aoh2hoh(\1, 'id') } qr/arrayref/,    'aoh2hoh: scalarref first arg -> arrayref croak';

#--------
# leak checks (skip under Devel::Cover, matching the suite convention)
#--------
unless ($INC{'Devel/Cover.pm'}) {
	no_leaks_ok {
		eval {
			aoh2hoh([ { id => 'a', x => 1 }, { id => 'b', x => 2 } ], 'id')
		}
	} 'aoh2hoh(): no memory leaks on normal input';

	no_leaks_ok {
		eval {
			aoh2hoh([ { id => 'a', x => 1 }, { id => 'a', x => 9 } ], 'id')
		}
	} 'aoh2hoh(): no leaks when dying on a duplicate key (partial result cleaned up)';

	no_leaks_ok {
		eval {
			aoh2hoh([ { id => 'a' }, 42, { x => 5 } ], 'id')
		}
	} 'aoh2hoh(): no leaks when rows are skipped';
}

done_testing();
