#!/usr/bin/env perl
use 5.042.2;
no source::encoding;
package Stats::LikeR;
our $VERSION = 0.01;
require XSLoader;
#use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
#use Devel::Confess 'color';
use 5.042.2;
no source::encoding;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use warnings FATAL => 'all';
use autodie ':default';
use Exporter 'import';
XSLoader::load('Stats::LikeR', $VERSION);
our @EXPORT_OK = qw(aov cor cor_test fisher_test glm hist lm matrix mean median min max p_adjust quantile rbinom rnorm runif scale sd seq shapiro_test t_test var write_table);
our @EXPORT = @EXPORT_OK;
package Stats::LikeR;

require XSLoader;

sub read_table { # mimics R's read.table; returns array of hash
	my ($args) = @_;
	my $args_ref = ref $args;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
	my @req_args = ('filename');#, 'row name');
	if ($args_ref ne 'HASH') {
		say STDERR 'Required args:';
		p @req_args;
		die "\"$args_ref\" was entered, but a \"HASH\" is needed by $current_sub";
	}
	my @undef_args = grep { !defined $args->{$_}} @req_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'the above args are necessary, but were not defined.';
	}
	my @defined_args = (@req_args,
		'comment.char', #character to skip lines after the header
		'sep', # by default ","
		'substitutions',
		'output.type'
	);
	my @bad_args = grep { # which args aren't defined?
								my $key = $_;
								not grep {$_ eq $key} @defined_args
							  } keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args;
		say 'the above arguments are not recognized.';
		p @defined_args;
		die 'The above args are accepted.'
	}
	$args->{sep} = $args->{sep} // ',';
	$args->{'output.type'} = $args->{'output.type'} // 'aoh';
	my (@data, @header, %data);
	open my $txt, '<', $args->{filename};
	while (<$txt>) {
		if (
				($. > 1) # past header
				&& 
				(defined $args->{'comment.char'})
				&&
				($_ =~ m/^$args->{'comment.char'}/)
			) {
			next	
		}
		chomp;
# CSV will often have commas inside the fields, which incorrectly splits some cells
		foreach my $sub (@{ $args->{substitutions} }) {
			$_ =~ s/$sub->[0]/$sub->[1]/g;
		}
		my @line = grep {$_ ne ''} split /$args->{sep}/;
		if ($. == 1) { # header
#			if ((defined $args->{'row name'}) && (!grep { $_ eq $args->{'row name'}} @line)) {
#				p @line;
#				die "the above line is missing $args->{'row name'}"
#			}
			foreach my $cell (@line) {
				$cell =~ s/^#//;
				if ($cell =~ m/^"(.+)"$/) { # remove quotes
					$cell = $1;
				}
			}
			@header = grep {$_ ne ''} @line;
			next;
		}
		if (scalar @line != scalar @header) {
			p @line;
			p @header;
			warn "the number of elements on $args->{filename} line $. != the header.";
			next;
		}
		next if grep {not defined} @line;
		my %line;
		@line{@header} = @line; # hash slice based on header
		if (grep {!defined} values %line) {
			p %line;
			die "got undefined values at line $. in $args->{filename}";
		}
		if ($args->{'output.type'} eq 'aoh') {
			push @data, \%line;
		} elsif ($args->{'output.type'} eq 'hoa') {
			foreach my $col (@header) {
				push @{ $data{$col} }, $line{$col};
			}
		}
#		my $key = $line{$args->{'row name'}};
#		delete $line{$args->{'row name'}};
#		$data{$key} = \%line;
	}
	if ($args->{'output.type'} eq 'aoh') {
		return \@data;
	} elsif ($args->{'output.type'} eq 'hoa') {
		return \%data;
	}
}

sub write_table {
	# 1. Allow passing the data hash ref as the first positional argument
	my $data_ref = (ref($_[0]) eq 'HASH') ? shift : undef;

	my $current_sub = (split(/::/,(caller(0))[3]))[-1];

	my %args = (
	  sep         => ',',
	  'row.names' => 1,      # Default to true, similar to R
	  @_,
	);

	# Map the positional data ref to the 'data' key if provided
	$args{data} //= $data_ref;

	# 2. Argument validation
	my %allowed = map { $_ => 1 } qw(data file row.names sep);
	my @err = grep { !$allowed{$_} } keys %args;
	if (@err > 0) {
	  die "$current_sub: Unknown arguments passed: " . join(", ", @err) . "\n";
	}

	die "$current_sub: 'file' argument is required\n" unless defined $args{file};
	die "$current_sub: 'data' must be a HASH reference\n" 
	  unless defined $args{data} && ref $args{data} eq 'HASH';

	my $data         = $args{data};
	my $sep          = $args{sep};
	my $file         = $args{file};
	my $inc_rownames = $args{'row.names'};

	my @rows = keys %$data;
	return if @rows == 0; # Nothing to write

	# 3. Determine if it's a Hash of Hashes or Hash of Arrays
	my $first_type = ref $data->{$rows[0]};
	die "$current_sub: Data values must be either all HASHes or all ARRAYs\n"
	  unless $first_type eq 'HASH' || $first_type eq 'ARRAY';

	# Validate that all elements match the first element's type
	my @type_err = grep { ref $data->{$_} ne $first_type } @rows;
	if (@type_err > 0) {
	  die "$current_sub: Mixed data types detected. Ensure all values are $first_type references.\n";
	}
	# 4. Open file for writing
	open my $fh, '>', $file or die "$current_sub: Could not open '$file' for writing: $!\n";
	if ($first_type eq 'HASH') {
	  # --- HASH OF HASHES ---
	  # Get all unique inner keys to act as column headers
	  my %col_map;
	  for my $r (@rows) {
		   $col_map{$_} = 1 for keys %{ $data->{$r} };
	  }
	  my @headers = sort keys %col_map;

	  # Print header
	  my @header_row = @headers;
	  unshift @header_row, "" if $inc_rownames; # Leave top-left cell blank for row.names
	  print $fh join($sep, @header_row) . "\n";
	  # Print data
	  for my $r (sort @rows) {
		   # Use "NA" for missing inner keys, mimicking R
		   my @row_data = map { defined $data->{$r}{$_} ? $data->{$r}{$_} : "NA" } @headers;
		   unshift @row_data, $r if $inc_rownames;
		   print $fh join($sep, @row_data) . "\n";
	  }
	} else {
	  # --- HASH OF ARRAYS ---
	  # Find the maximum array length to know how many columns we have
	  my $max_cols = 0;
	  for my $r (@rows) {
		   my $len = scalar @{ $data->{$r} };
		   $max_cols = $len if $len > $max_cols;
	  }

	  # Generate headers (e.g., V1, V2, V3... mimicking R's default behavior)
	  my @headers = map { "V$_" } 1 .. $max_cols;
	  my @header_row = @headers;
	  unshift @header_row, "" if $inc_rownames;
	  print $fh join($sep, @header_row) . "\n";

	  # Print data
	  for my $r (sort @rows) {
		   my @row_data = @{ $data->{$r} };
		   # Pad with "NA" if some arrays are shorter than others
		   push @row_data, "NA" while @row_data < $max_cols;
		   unshift @row_data, $r if $inc_rownames;
		   print $fh join($sep, @row_data) . "\n";
	  }
	}
}
1;
#sub mean {
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 1 element in the array" if scalar @n < 1;
#	return sum(@n) / scalar @n;x
#}

#sub stdev {
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 2 elements in the array" if scalar @n < 2;
#	my $mean = sum(@n) / scalar @n;
#	my $standard_deviation = 0;
#	foreach my $element (@n) {
#		$standard_deviation += ($element-$mean)**2;
#	}
#	return sqrt($standard_deviation/((scalar @n)-1));
#}

#sub population_variance { # sample variance
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 1 elements in the array" if scalar @n < 1;
#	my $mean = sum(@n) / scalar @n;
#	my $var = 0;
#	foreach my $element (@n) {
#		$var += ($element-$mean)**2;
#	}
#	return $var / scalar @n;
#}

#sub sample_variance { # sample variance
#	my @n = map { ref($_) eq 'ARRAY' ? @$_ : $_ } @_;
#	my $current_sub = ( split( /::/, ( caller(0) )[3] ) )[-1];
#	die "$current_sub needs >= 1 elements in the array" if scalar @n < 1;
#	my $mean = sum(@n) / scalar @n;
#	my $var = 0;
#	foreach my $element (@n) {
#		$var += ($element-$mean)**2;
#	}
#	return $var / (scalar @n - 1);
#}

#sub pvalue ($array1, $array2) {
#	return 1.0 if scalar @$array1 <= 1;
#	return 1.0 if scalar @$array2 <= 1;
#	my $mean1 = sum(@{ $array1 });
#	my $mean2 = sum(@{ $array2 });
#	return 1.0 if ($mean1 == $mean2);
#	$mean2 /= scalar @$array2;
#	$mean1 /= scalar @$array1;
#	my ($variance1, $variance2) = (0, 0);
#	foreach my $x (@$array1) {
#	$variance1 += ($x-$mean1)*($x-$mean1);
#	}
#	foreach my $x (@$array2) {
#	$variance2 += ($x-$mean2)*($x-$mean2);
#	}
#	if (($variance1 == 0.0) && ($variance2 == 0.0)) {
#	return 1.0;
#	}
#	$variance1 = $variance1/(scalar @$array1-1);
#	$variance2 = $variance2/(scalar @$array2-1);
#	my $array1_size = scalar @$array1;
#	my $array2_size = scalar @$array2;
#	my $WELCH_T_STATISTIC = ($mean1-$mean2)/sqrt($variance1/$array1_size+$variance2/$array2_size);
#	my $DEGREES_OF_FREEDOM = (($variance1/$array1_size+$variance2/(scalar @$array2))**2)
#	/
#	(
#	($variance1*$variance1)/($array1_size*$array1_size*($array1_size-1))+
#	($variance2*$variance2)/($array2_size*$array2_size*($array2_size-1))
#	);
#	my $A = $DEGREES_OF_FREEDOM/2;
#	my $value = $DEGREES_OF_FREEDOM/($WELCH_T_STATISTIC*$WELCH_T_STATISTIC+$DEGREES_OF_FREEDOM);
##from here, translation of John Burkhardt's C
#	my $beta = lgamma($A)+0.57236494292470009-lgamma($A+0.5);
#	my $acu = 10**(-15);
#	my($ai,$cx,$indx,$ns,$pp,$psq,$qq,$rx,$temp,$term,$xx);
## Check the input arguments.
#	return $value if $A <= 0.0;# || $q <= 0.0;
#	return $value if $value < 0.0 || 1.0 < $value;
## Special cases
#	return $value if $value == 0.0 || $value == 1.0;
#	$psq = $A + 0.5;
#	$cx = 1.0 - $value;
#	if ($A < $psq * $value) {
#		($xx, $cx, $pp, $qq, $indx) = ($cx, $value, 0.5, $A, 1);
#	} else {
#		($xx, $pp, $qq, $indx) = ($value, $A, 0.5, 0);
#	}
#	$term = 1.0;
#	$ai = 1.0;
#	$value = 1.0;
#	$ns = int($qq + $cx * $psq);
##Soper reduction formula.
#	$rx = $xx / $cx;
#	$temp = $qq - $ai;
#	$rx = $xx if $ns == 0;
#	while (1) {
#		$term = $term * $temp * $rx / ( $pp + $ai );
#		$value = $value + $term;
#		$temp = abs ($term);
#		if ($temp <= $acu && $temp <= $acu * $value) {
#	   	$value = $value * exp ($pp * log($xx)
#	                          + ($qq - 1.0) * log($cx) - $beta) / $pp;
#	   	$value = 1.0 - $value if $indx;
#	   	last;
#		}
#	 	$ai = $ai + 1.0;
#		$ns = $ns - 1;
#		if (0 <= $ns) {
#			$temp = $qq - $ai;
#			$rx = $xx if $ns == 0;
#		} else {
#			$temp = $psq;
#			$psq = $psq + 1.0;
#		}
#	}
#	return $value;
#}

#sub paired_pvalue ($array1, $array2) {
#	return 1.0 if scalar @$array1 <= 1;
#	return 1.0 if scalar @$array2 <= 1;
#	my $array_size = scalar @$array1;
#	if ($array_size != scalar @$array2) {
#		die 'arrays have different sizes, cannot calculate paired p-value.';
#	}
#	my ($sd, $xd) = (0,0);
#	foreach my $x (0..$array_size-1) {
#		$xd += (@$array2[$x] - @$array1[$x]);
#	}
#	return 1.0 if $xd == 0;
#	$xd /= $array_size;
#	foreach my $x (0..$array_size-1) {
#		$sd += ($xd - (@$array2[$x] - @$array1[$x]))**2;
#	}
#	$sd = sqrt($sd / ($array_size - 1));
#	my $t = $xd / ($sd / sqrt($array_size));
#	my $DEGREES_OF_FREEDOM = $array_size - 1;#http://oak.ucc.nau.edu/rh232/courses/EPS525/Handouts/Understanding%20the%20Dependent%20t%20Test.pdf
#	my $A = $DEGREES_OF_FREEDOM/2;
#	my $value = $DEGREES_OF_FREEDOM/($t*$t+$DEGREES_OF_FREEDOM);
##from here, translation of John Burkhardt's C
#	my $beta = lgamma($A)+0.57236494292470009-lgamma($A+0.5);
#	my $acu = 10**(-15);
#	my($ai,$cx,$indx,$ns,$pp,$psq,$qq,$rx,$temp,$term,$xx);
## Check the input arguments.
#	return $value if $A <= 0.0;# || $q <= 0.0;
#	return $value if $value < 0.0 || 1.0 < $value;
## Special cases
#	return $value if $value == 0.0 || $value == 1.0;
#	$psq = $A + 0.5;
#	$cx = 1.0 - $value;
#	if ($A < $psq * $value) {
#		($xx, $cx, $pp, $qq, $indx) = ($cx, $value, 0.5, $A, 1);
#	} else {
#		($xx, $cx, $pp, $qq, $indx) = ($value, $cx, $A, 0.5, 0);
#	}
#	$term = 1.0;
#	$ai = 1.0;
#	$value = 1.0;
#	$ns = int($qq + $cx * $psq);
##Soper reduction formula.
#	$rx = $xx / $cx;
#	$temp = $qq - $ai;
#	$rx = $xx if $ns == 0;
#	while (1) {
#		$term = $term * $temp * $rx / ( $pp + $ai );
#		$value = $value + $term;
#		$temp = abs ($term);
#		if ($temp <= $acu && $temp <= $acu * $value) {
#			$value = $value * exp ($pp * log($xx)
#				                   + ($qq - 1.0) * log($cx) - $beta) / $pp;
#			$value = 1.0 - $value if $indx;
#			last;
#		}
#		$ai = $ai + 1.0;
#		$ns = $ns - 1;
#		if (0 <= $ns) {
#			$temp = $qq - $ai;
#			$rx = $xx if $ns == 0;
#		} else {
#			$temp = $psq;
#			$psq = $psq + 1.0;
#		}
#	}
#	$value;
#}

1;
