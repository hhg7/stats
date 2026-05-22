#!/usr/bin/env perl

require 5.010;
use strict;
use feature 'say';
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;
use Scalar::Util 'looks_like_number';

sub summary {
	my ($data, %args);
	my $current_sub = (split(/::/,(caller(0))[3]))[-1];

	if (@_ && ref $_[0]) {
	  # Handles: summary(\@arr) or summary(\@arr, nrows => 5) or summary(\%h, nrow => 3)
	  $data = shift;
	  %args = @_; # capture any trailing key/value pairs
	} else {
	  # Handles: summary(@runif) or summary(@runif, nrows => 2)
	  # Extract known trailing named arguments from the flat list
	  while (@_ >= 2 && defined $_[-2] && !ref($_[-2]) && $_[-2] =~ /^(?:nrows|nrow)$/) {
	  	  my $val = pop @_;
		  my $key = pop @_;
		  $args{$key} = $val;
	  }
	  # The remaining items in @_ make up the actual data array
	  my @list = @_;
	  $data = \@list;
	}

	# Normalize nrow -> nrows, default to 10
	$args{nrows} //= delete($args{nrow}) // 10;

	my $ref_type = ref $data;
	if (($ref_type ne 'ARRAY') && ($ref_type ne 'HASH')) {
		die "$current_sub' data must either be a hash or an array, not \"$ref_type\"";
	}
	
	my $single_arr = 0;
	if (($ref_type eq 'ARRAY') && (ref $data->[0] eq '')) {
		$single_arr = 1;
	}
	
	my @header = ('# values', 'Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.');
	my @out;
	
	if ($single_arr == 1) {
		push @out, '-' x 75;
		my $header = sprintf('%9s ' x scalar @header, @header);
		push @out, $header;
		push @out, '-' x 75;
		my @undef = grep {!defined $data->[$_]} 0..scalar @{ $data }-1;
		if (scalar @undef > 0) {
			say STDERR join (',', @undef);
			die "The above indices are not defined in $current_sub";
		}
		my @numeric = grep {looks_like_number($_)} @{ $data };
		my $q = quantile(\@numeric, probs => [0.25, 0.75]);
		my $vals = sprintf('%9.4g ' x scalar @header, scalar @numeric, min(\@numeric), $q->{'25%'}, median(\@numeric), mean(\@numeric), $q->{'75%'}, max(\@numeric));
		push @out, $vals;
	} elsif ($ref_type eq 'ARRAY') {
		push @out, '-' x 75;
		my $header = sprintf('%9s ' x scalar @header, @header);
		unshift @header, 'Index';
		$header = 'Index ' . $header;
		push @out, $header;
		push @out, '-' x 75;
		my $rows_printed = 0;
		foreach my $index (0..$#$data) {
			my @undef = grep {!defined $data->[$index][$_]} 0..scalar @{ $data->[$index] }-1;
			if (scalar @undef > 0) {
				say STDERR join (',', @undef);
				die "The above indices are not defined for index $index in $current_sub";
			}
			my @numeric = grep {looks_like_number($_)} @{ $data->[$index] };
			my $q = quantile(\@numeric, probs => [0.25, 0.75]);
			my $vals = sprintf('%6.4g', $index) . sprintf('%9.4g ' x (scalar @header - 1), scalar @numeric, min(\@numeric), $q->{'25%'}, median(\@numeric), mean(\@numeric), $q->{'75%'}, max(\@numeric));
			push @out, $vals;
			$rows_printed++;
			last if $rows_printed >= $args{nrows}; # Changed to >= just to be safe
		}
	} elsif ($ref_type eq 'HASH') {
		push @out, '-' x 78;
		my $header = sprintf('%9s ' x scalar @header, @header);
		unshift @header, 'Key';
		$header = '  Key    ' . $header;
		push @out, $header;
		push @out, '-' x 78;
		my $rows_printed = 0;
		foreach my $key (sort {lc $a cmp lc $b} keys %{ $data }) {
			say $key;
			my @undef = grep {!defined $data->{$key}[$_]} 0..scalar @{ $data->{$key} }-1;
			if (scalar @undef > 0) {
				say STDERR join (',', @undef);
				die "The above indices are not defined for key $key in $current_sub";
			}
			my @numeric = grep {looks_like_number($_)} @{ $data->{$key} };
			my $q = quantile(\@numeric, probs => [0.25, 0.75]);
			my $print_key;
			if ($key =~ m/^(.{,9})/) { # take the first 9 characters of the key
				$print_key = $1;
			} else {
				die "\"$key\" failed regex";
			}
			if ((length $print_key) < 9) { # make sure that short keys line up correctly
				$print_key .= ' ' x (9 - length $print_key);
			}
			my $vals = $print_key . sprintf('%9.4g ' x (scalar @header - 1), scalar @numeric, min(\@numeric), $q->{'25%'}, median(\@numeric), mean(\@numeric), $q->{'75%'}, max(\@numeric));
			push @out, $vals;
			$rows_printed++;
			last if $rows_printed >= $args{nrows};
		}
	}
	say join ("\n", @out);
	return \@out;
}
my @runif = (0.216281301648454, 0.109465371442155, 0.152664169241813, 0.000945096692653635, 0.297535893169954, 0.139163636355065, 0.433281186173499, 0.408562817186144, 0.407467355710114, 0.544592780352787, 0.487398883576855, 0.643468442596237, 0.68575522492846, 0.151846994960366, 0.0108012662535621, 0.765504103474193, 0.170995624940421, 0.100078161688572, 0.167327253677694, 0.178543828268637, 0.767033648977208, 0.0661950819672228, 0.581462013571265, 0.584800690627297, 0.762539213217881, 0.233645411264945, 0.693534299360277, 0.513290613560038, 0.41433325215603, 0.73243812858739, 0.478323977378576, 0.798072957187451, 0.237619591881074, 0.0780442614619403, 0.0360511965325365, 0.660791977980871, 0.912043981453014, 0.415870135589202, 0.831491877016528, 0.737746524987607, 0.663143394629547, 0.777232190070094, 0.816913688077346, 0.352381995029283, 0.744148647065789, 0.729401956002121, 0.465347760265214, 0.0785176667616199, 0.181269420249411, 0.679185700779414, 0.953224347579702, 0.567208135290578, 0.292655755357845, 0.105132055128408, 0.659550831920821, 0.260928737252719, 0.0114517904292804, 0.351924227264533, 0.539668158788782, 0.923435653386754, 0.679118775225493, 0.541537731065048, 0.235382321740357, 0.443470864148644, 0.49701302243216, 0.124681475319193, 0.403251186205477, 0.587374376354269, 0.0806932538910878, 0.613866141439061, 0.285459073394659, 0.882170197671563, 0.729358588888918, 0.872760579993155, 0.0726024246860497, 0.599972473528148, 0.857066010638153, 0.767531044559306, 0.877534848570345, 0.520403080150906, 0.115952349478963, 0.0624610171846882, 0.869999228452524, 0.294535850510563, 0.735723449504025, 0.727797725687921, 0.232053652861307, 0.486724559407229, 0.497430051763761, 0.65156677164174, 0.456347032400441, 0.785195872302019, 0.120408844445638, 0.45376514163452, 0.198314702590377, 0.144783732275236, 0.064735910938797, 0.30123682582493, 0.437664391094597);
summary(\@runif);
summary(@runif, nrows => 2);
my (%h, @arr);
foreach my $i (0..18) {
	push @arr, runif(22);
	$h{"AAAAAAAA$i"} = runif(22);
}
$h{short} = runif(9);
summary(\@arr, nrows => 19);
summary(\%h);
summary(\%h, nrows => 3);
