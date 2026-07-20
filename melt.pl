#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
#use Util;

my $df = [ { A => 'a', B => 1, C => 2 },
               { A => 'b', B => 3, C => 4 } ];
view( $df );
$df = melt( $df, id_vars => 'A', value_vars => [ 'B', 'C' ] );
view( $df );
