#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
#use Util;

my $ti = read_table('titanic.more.complete.csv');
view($ti);
$ti = cfilter($ti, remove => ['ticketno']);
view($ti);
