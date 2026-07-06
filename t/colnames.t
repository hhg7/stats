#!/usr/bin/env perl

require 5.010;
use warnings FATAL => 'all';
use Stats::LikeR;
use Test::Exception; # throws_ok / dies_ok
use Test::More;
use Scalar::Util 'reftype';
use Test::LeakTrace 'no_leaks_ok';

sub colnames {
	my ($data) = @_;
	unless (defined $data) {
		die 'colnames received undefined data';
	}
	my $type = reftype $data;

	die 'colnames: expected an ARRAY or HASH ref (got '
	  . (defined $data ? (ref($data) || 'non-ref scalar') : 'undef') . ")\n"
	  unless defined $type && ($type eq 'ARRAY' || $type eq 'HASH');

	if ($type eq 'ARRAY') {
	  return wantarray ? () : [] unless @$data;               # empty frame
	  my $r0 = reftype $data->[0];                            # element 0 decides the form

	  # AoH: columns = keys per row
	  if (defined $r0 && $r0 eq 'HASH') {
		   my @cols = keys %{ $data->[0] };
		   my $ncol = scalar @cols;
		   for my $i (1 .. $#$data) {
		       my $row = $data->[$i];
		       die "colnames: AoH row $i is not a hash ref\n"
		           unless defined(reftype $row) && reftype($row) eq 'HASH';
		       my $k = scalar keys %$row;
		       die "colnames: ragged AoH — row $i has $k columns, but row 0 has $ncol\n"
		           if $k != $ncol;
		   }
		   return wantarray ? @cols : \@cols;
	  }

	  # AoA: columns = row length (Mimic R's V1..Vn default naming)
	  if (defined $r0 && $r0 eq 'ARRAY') {
		   my $ncol = scalar @{ $data->[0] };
		   for my $i (1 .. $#$data) {
		       my $row = $data->[$i];
		       die "colnames: AoA row $i is not an array ref\n"
		           unless defined(reftype $row) && reftype($row) eq 'ARRAY';
		       my $len = scalar @$row;
		       die "colnames: ragged AoA — row $i has $len columns, but row 0 has $ncol\n"
		           if $len != $ncol;
		   }
		   my @cols = map { "V$_" } 1 .. $ncol;
		   return wantarray ? @cols : \@cols;
	  }

	  # plain 1-D vector: one column
	  if (not defined $r0) {
		   my @cols = ('V1');
		   return wantarray ? @cols : \@cols; 
	  }

	  die "colnames: array element 0 is a $r0 ref; expected HASH (AoH), ARRAY (AoA), or plain scalars (vector)\n";
	}

	# HASH: HoA (keys are columns) or HoH (keys are rows)
	return wantarray ? () : [] unless %$data;

	my $probe;                                                  # first defined value decides the form
	foreach my $k (keys %$data) {
	  next unless defined $data->{$k};
	  $probe = $data->{$k};
	  last;
	}
	my $vtype = reftype $probe;

	# HoA: keys ARE the columns.
	if (defined $vtype && $vtype eq 'ARRAY') {
	  foreach my $col (keys %$data) {
		   die "colnames: HoA column '$col' is not an array ref\n"
		       unless defined(reftype $data->{$col}) && reftype($data->{$col}) eq 'ARRAY';
	  }
	  my @cols = keys %$data;
	  return wantarray ? @cols : \@cols;
	}

	# HoH: keys are rows; columns = keys of a row hash
	if (defined $vtype && $vtype eq 'HASH') {
	  my ($ncol, $ref_row, @cols);
	  foreach my $row_key (keys %$data) {
		   my $row = $data->{$row_key};
		   die "colnames: HoH row '$row_key' is not a hash ref\n"
		       unless defined(reftype $row) && reftype($row) eq 'HASH';
		   
		   my @current_keys = keys %$row;
		   my $k = scalar @current_keys;
		   
		   if (not defined $ncol) { 
		       $ncol = $k; 
		       $ref_row = $row_key;
		       @cols = @current_keys;
		   }
		   elsif ($k != $ncol) {
		       die "colnames: ragged HoH — row '$row_key' has $k columns, but '$ref_row' has $ncol\n";
		   }
	  }
	  return wantarray ? @cols : \@cols;
	}

	die "colnames: HASH values are neither ARRAY refs (HoA) nor HASH refs (HoH)\n";
}


my %d = (
	A => ['a','b'],
	B => [1,2]
);

is_deeply( colnames(\%d), ('A', 'B'), 'simple HoA works');

done_testing();
