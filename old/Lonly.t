# Retired tests for the original two-argument Lonly (backed by set_difference).
#
# Removed from t/set_ops.t when the `Lonly` name was reassigned to the former
# `get_unique`. The `$$`-prototype guards below (exactly-2-refs, left/right
# argument croaks) no longer apply to the new @-prototype Lonly, whose tests now
# live in t/set_ops.t. Preserved here verbatim.
#
# Relies on the same fixtures/helpers as t/set_ops.t:
#   my @a = (1, 2, 3, 4, 5, 4, 3);
#   my @b = (3, 4, 5, 6, 7);
#   my @c = (5, 6, 7, 8);
#   my $LONLY = \&Stats::LikeR::Lonly;

#--------
# Lonly   (values in left ref, absent from right ref)
#--------
is_list( [Lonly(\@a, \@b)], [1, 2],   'Lonly: left-only values, deduped, in left order');
is_list( [Lonly(\@b, \@a)], [6, 7],   'Lonly: not symmetric with the reversed call');
is_list( [Lonly(\@b, \@b)], [],       'Lonly: identical lists -> empty');
is_list( [Lonly([], \@b)],  [],       'Lonly: empty left -> empty');
is_list( [Lonly(\@a, [])],  [1, 2, 3, 4, 5], 'Lonly: empty right -> left distinct');
is( scalar(Lonly(\@a, \@b)), 2, 'Lonly: scalar context returns count');

throws_ok { my @x = Lonly('x', \@b) }        qr/first \(left\) argument/,  'Lonly: croaks on non-ref left';
throws_ok { my @x = Lonly(\@a, 'x') }        qr/second \(right\) argument/,'Lonly: croaks on non-ref right';
throws_ok { my @x = Lonly([undef], \@b) }    qr/undefined value in left/,  'Lonly: croaks on undef in left';
throws_ok { my @x = Lonly(\@a, [undef]) }    qr/undefined value in right/, 'Lonly: croaks on undef in right';
throws_ok { $LONLY->(\@a) }                  qr/exactly 2 array refs/,     'Lonly: croaks on 1 ref';
throws_ok { $LONLY->(\@a, \@b, \@c) }        qr/exactly 2 array refs/,     'Lonly: croaks on 3 refs';

Lonly(\@a, \@b);                            # hoist
no_leaks_ok {
	eval {
		my @l = Lonly(\@a, \@b);
		my $s = Lonly(\@a, \@b);
	}
} 'Lonly(): no memory leaks' unless $INC{'Devel/Cover.pm'};
