#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
#use Util;
#use latex qw(write_2d_array_to_tex_tabular);
# s/[!@#\$\%^&*\(\)\{\}\[\]\<\>,\/'"\-\h;\+=]+/_/g; # annoying chars
#qw(list_regex_files json_file_to_ref ref_to_json_file);
my @r = rank( -7.765, -9.328, -10.326, -9.038, -9.608, -9.779, -9.975, -6.906);
p @r;
