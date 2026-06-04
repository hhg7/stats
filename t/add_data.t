#!/usr/bin/env perl

require 5.010;
use warnings FATAL => 'all';
use Stats::LikeR; # Adjust this if add_data is exported from a different module
use Test::Exception; # dies_ok
use Test::More;
use Test::LeakTrace 'no_leaks_ok';

#----------------------------------------------------------------------
# 1. Input Validation Tests
#----------------------------------------------------------------------
dies_ok {
	add_data(undef, {});
} 'add_data: dies with undefined first argument';

dies_ok {
	add_data([], {});
} 'add_data: dies with array reference for first argument';

dies_ok {
	add_data({}, "string");
} 'add_data: dies with scalar for second argument';

#----------------------------------------------------------------------
# 2. Hash of Hashes (HoH) Merging
#----------------------------------------------------------------------
my $h_hoh = { A => { 'x' => 1 } };
my $i_hoh = { A => { 'y' => 2 }, B => { z => 3 } };

add_data($h_hoh, $i_hoh);

my $n = scalar keys %{ $h_hoh };
if ($n == 2) {
	pass('add_data (HoH): has the correct # of primary hash keys');
} else {
	fail("add_data (HoH): has $n primary hash keys, when it should have 2");
}

if (defined $h_hoh->{A}{'x'} && defined $h_hoh->{A}->{y}) {
	pass('add_data (HoH): correctly merged keys into existing row');
} else {
	fail('add_data (HoH): failed to merge keys into existing row');
}

if (defined $h_hoh->{B}{z}) {
	pass('add_data (HoH): correctly created new row from secondary hash');
} else {
	fail('add_data (HoH): failed to create new row');
}

no_leaks_ok {
	my $h_test = { A => { x => 1 } };
	my $i_test = { A => { y => 2 }, B => { z => 3 } };
	add_data($h_test, $i_test);
} 'add_data: no leaks when merging Hash of Hashes' unless $INC{'Devel/Cover.pm'};

#----------------------------------------------------------------------
# 3. Hash of Arrays (HoA) Merging (New Functionality)
#----------------------------------------------------------------------
my $h_hoa = { A => [1, 2] };
my $i_hoa = { A => [3, 4], B => [5] };

add_data($h_hoa, $i_hoa);

$n = scalar keys %{ $h_hoa };
if ($n == 2) {
	pass('add_data (HoA): has the correct # of primary hash keys');
} else {
	fail("add_data (HoA): has $n primary hash keys, when it should have 2");
}

my $a_len = scalar @{ $h_hoa->{A} };
if ($a_len == 4) {
	pass('add_data (HoA): correctly appended elements to existing array row');
} else {
	fail("add_data (HoA): existing array has $a_len elements, should have 4");
}

my $b_len = scalar @{ $h_hoa->{B} };
if ($b_len == 1 && $h_hoa->{B}[0] == 5) {
	pass('add_data (HoA): correctly created new array row from secondary hash');
} else {
	fail('add_data (HoA): failed to properly create new array row');
}

no_leaks_ok {
	my $h_test = { A => [1, 2] };
	my $i_test = { A => [3, 4], B => [5] };
	add_data($h_test, $i_test);
} 'add_data: no leaks when merging Hash of Arrays' unless $INC{'Devel/Cover.pm'};

#----------------------------------------------------------------------
# 4. Target Inference (Empty Target Hash)
#----------------------------------------------------------------------
my $h_empty_hoa = {};
my $i_source_hoa = { X => [9, 10] };

add_data($h_empty_hoa, $i_source_hoa);

if (ref($h_empty_hoa->{X}) eq 'ARRAY') {
	pass('add_data (Inference): correctly inferred HoA intent from empty target');
} else {
	fail("add_data (Inference): failed to infer HoA intent, created " . ref($h_empty_hoa->{X}) . " instead");
}

my $h_empty_hoh = {};
my $i_source_hoh = { Y => { k => 'v' } };

add_data($h_empty_hoh, $i_source_hoh);

if (ref($h_empty_hoh->{Y}) eq 'HASH') {
	pass('add_data (Inference): correctly inferred HoH intent from empty target');
} else {
	fail("add_data (Inference): failed to infer HoH intent, created " . ref($h_empty_hoh->{Y}) . " instead");
}

no_leaks_ok {
	my $h_test = {};
	my $i_test = { X => [9, 10] };
	add_data($h_test, $i_test);
} 'add_data: no leaks during empty target inference' unless $INC{'Devel/Cover.pm'};

#----------------------------------------------------------------------
# 5. Legacy Fallback (Target is HoH, Source row is Array)
#----------------------------------------------------------------------
my $h_legacy = { A => { key1 => 'val1' } };
my $i_legacy = { A => [ key2 => 'val2', key3 => 'val3' ] }; # Array acting as kv pairs

add_data($h_legacy, $i_legacy);

if (defined $h_legacy->{A}->{key2} && $h_legacy->{A}->{key3} eq 'val3') {
	pass('add_data (Legacy): correctly processed array as key-value pairs for HoH target');
} else {
	fail('add_data (Legacy): failed to process array as key-value pairs');
}

done_testing();
