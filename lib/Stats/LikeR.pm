#!/usr/bin/env perl
use 5.042.2;
no source::encoding;
package Stats::LikeR;
our $VERSION = 0.01;
require XSLoader;
use 5.042.2;
no source::encoding;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use warnings FATAL => 'all';
use autodie ':default';
use Exporter 'import';
XSLoader::load('Stats::LikeR', $VERSION);
our @EXPORT_OK = qw(aov cor cor_test fisher_test glm hist lm matrix mean median min max p_adjust quantile rbinom read_table rnorm runif scale sd seq shapiro_test t_test var write_table);
our @EXPORT = @EXPORT_OK;

require XSLoader;

sub read_table {
	my $file = shift;
	die "\"$file\" is either unreadable or not a file" unless -r -f $file;
	my %args = (
		sep => ',', comment => '#',
		@_,
	);
	my %allowed_args = map {$_ => 1} (
		'comment', #character to skip lines after the header
		'row.name',
		'sep', # by default ","
		'substitutions',
		'output.type' # aoh, hoa, or hoh
	);
	my @undef_args = sort grep {!$allowed_args{$_}} keys %allowed_args;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1];
	if (scalar @undef_args > 0) {
		p @undef_args;
		die "the above args aren't defined for $current_sub";
	}
	# ... (argument validation logic from your original script) ...
	$args{'output.type'} = $args{'output.type'} // 'aoh';
	if ($args{'output.type'} !~ m/^(?:aoh|hoa|hoh)$/) {
		die "\"$args{'output.type'}\" isn't allowed";
	}
	$args{comment} = $args{comment} // '#';
	my (@data, %data, @header);
	open my $txt, '<', $file;
	while (<$txt>) {
		next if $_ =~ m/^$args{comment}/;
		next if /^\h*$/; # Skip empty lines
		$_ =~ s/\r?\n$//; # chomp with annoying and invisible Windows "\r"
		# Apply substitutions if any
		foreach my $sub (@{ $args{substitutions} // [] }) {
			$_ =~ s/$sub->[0]/$sub->[1]/g;
		}
		# Use -1 to keep trailing empty fields
		my @line = split /$args{sep}/, $_;
		if ($. == 1) {
			# --- HEADER PROCESSING ---
			foreach my $cell (@line) {
				$cell =~ s/^#//;      # Remove comment prefix if present
				$cell =~ s/"$//; # Strip surrounding quotes 
				$cell =~ s/^"//;
				$cell =~ s/\"\"/\"/g;  # Un-escape doubled quotes 
			}
			# FIX: Instead of grep, only remove trailing empty fields
			# that might exist due to a trailing separator.
			while (@line && $line[-1] eq '') { pop @line }
			@header = @line; 
			# R-LIKE BEHAVIOR: If the first header is blank (like in HepatitisCdata.csv),
			# give it a name so it can be used as a hash key for the index column.
			if ((scalar @header > 0) && ($header[0] eq '')) {
				$header[0] = 'row_name'; 
			}
			if (($args{'output.type'} eq 'hoh') && (not defined $args{'row.name'})) {
				$args{'row.name'} = $header[0];
			}
			if ((defined $args{'row.name'}) && (!grep {$_ eq $args{'row.name'}} @header)) {
				die "\"$args{'row.name'}\" isn't in the header of $file";
			}
			next;
		}
		# Check for column alignment
		if (scalar @line != scalar @header) {
			warn "Alignment error on $file line $. (" . scalar(@line) . " fields vs " . scalar(@header) . " headers).";
			next;
		}
		# --- DATA PROCESSING ---
		my %line;
		for my $i (0 .. $#header) {
			my $cell = $line[$i];
			# Strip quotes and handle un-escaping for data fields 
			if (defined $cell) {
				$cell =~ s/^\"|\"$//g;
				$cell =~ s/\"\"/\"/g;
				$cell =~ s/"$//;
				$cell =~ s/^"//;
			}
			# R-like behavior: Treat empty strings as 'NA' 
			$line{$header[$i]} = ($cell eq '') ? 'NA' : $cell;
		}
		if ($args{'output.type'} eq 'aoh') {
			push @data, \%line;
		} elsif ($args{'output.type'} eq 'hoa') {
			foreach my $col (@header) {
				push @{ $data{$col} }, $line{$col};
			}
		} elsif ($args{'output.type'} eq 'hoh') {
			my $row_name = $line{$args{'row.name'}};
			foreach my $col (@header) {
				$data{$col}{$row_name} = $line{$col};
			}
		}
	}
	close $txt;
	if ($args{'output.type'} eq 'aoh') {
		return \@data;
	} elsif ($args{'output.type'} =~ m/^(?:hoa|hoh)$/) {
		return \%data;
	}
}

sub write_table {
	my $data_ref = (ref($_[0]) eq 'HASH' || ref($_[0]) eq 'ARRAY') ? shift : undef;
	my $file = shift;
	my %args = (
		sep         => ',',
		'row.names' => 1,      
		@_,
	);
	my $current_sub = (split(/::/,(caller(0))[3]))[-1];
	$args{data} //= $data_ref;
	my %allowed = map { $_ => 1 } qw(data file row.names sep col.names);
	my @err = grep { !$allowed{$_} } keys %args;
	if (@err > 0) {
		die "$current_sub: Unknown arguments passed: " . join(", ", @err) . "\n";
	}
	die "$current_sub: 'data' must be a HASH or ARRAY reference\n" 
		unless defined $args{data} && (ref($args{data}) eq 'HASH' || ref($args{data}) eq 'ARRAY');

	my $col_names = $args{'col.names'};
	if (defined $col_names && ref($col_names) ne 'ARRAY') {
		die "$current_sub: 'col.names' must be an ARRAY reference\n";
	}
	my $data         = $args{data};
	my $sep          = $args{sep};
	my $inc_rownames = $args{'row.names'};
	my $quote_field = sub {
		my ($val, $sep) = @_;
		return '' unless defined $val;
		die "$current_sub: Cannot write nested reference types to table\n" if ref($val);
		
		my $str = "$val";  
		if (index($str, $sep) != -1 || index($str, '"') != -1 || $str =~ /[\r\n]/) {
			$str =~ s/"/""/g;
			$str = qq{"$str"};
		}
		return $str;
	};
	my $data_type = ref $data;
	my ($is_hoh, $is_hoa, $is_aoh) = (0, 0, 0);
	my @rows;
	if ($data_type eq 'HASH') {
		@rows = keys %$data;
		return if @rows == 0; 
		my $first_type = ref $data->{$rows[0]};
		die "$current_sub: Data values must be either all HASHes or all ARRAYs\n"
			unless $first_type eq 'HASH' || $first_type eq 'ARRAY';

		my @type_err = grep { ref $data->{$_} ne $first_type } @rows;
		if (@type_err > 0) {
			die "$current_sub: Mixed data types detected. Ensure all values are $first_type references.\n";
		}
		$is_hoh = ($first_type eq 'HASH');
		$is_hoa = ($first_type eq 'ARRAY');
	} else {
		return if @$data == 0;

		my $first_elem = $data->[0];
		die "$current_sub: For ARRAY data, all elements must be HASH references (Array of Hashes)\n"
			unless defined $first_elem && ref($first_elem) eq 'HASH';

		my @type_err = grep { !defined($_) || ref($_) ne 'HASH' } @$data;
		if (@type_err > 0) {
			die "$current_sub: Mixed data types detected in Array of Hashes. All elements must be HASH references.\n";
		}
		$is_aoh = 1;
	}
	open my $fh, '>', $file or die "$current_sub: Could not open '$file' for writing: $!\n";
	if ($is_hoh) {
		my @headers;
		if ($col_names) {
			@headers = @$col_names;
		} else {
			my %col_map;
			for my $r (@rows) {
				$col_map{$_} = 1 for keys %{ $data->{$r} };
			}
			@headers = sort keys %col_map;
		}

		my @header_row = @headers;
		unshift @header_row, '' if $inc_rownames;
		@header_row = map { $quote_field->($_, $sep) } @header_row;
		print $fh join($sep, @header_row) . "\n";

		for my $r (sort @rows) {
			my @row_data = map { defined $data->{$r}{$_} ? $data->{$r}{$_} : "NA" } @headers;
			unshift @row_data, $r if $inc_rownames;
			my @quoted = map { $quote_field->($_, $sep) } @row_data;
			print $fh join($sep, @quoted) . "\n";
		}
	} elsif ($is_hoa) {
		my $max_cols = 0;
		for my $r (@rows) {
			my $len = scalar @{ $data->{$r} };
			$max_cols = $len if $len > $max_cols;
		}
		my @headers;
		if ($col_names) {
			@headers = @$col_names;
			$max_cols = scalar @headers;
		} else {
			@headers = map { "V$_" } 1 .. $max_cols;
		}
		my @header_row = @headers;
		unshift @header_row, "" if $inc_rownames;
		@header_row = map { $quote_field->($_, $sep) } @header_row;
		say $fh join($sep, @header_row);
		for my $r (sort @rows) {
			my @row_data = map { defined $_ ? $_ : "NA" } @{ $data->{$r} }[0 .. $max_cols - 1];
			# Ensure we pad with NA if the row data is shorter than max_cols
			push @row_data, 'NA' while @row_data < $max_cols;
			unshift @row_data, $r if $inc_rownames;
			my @quoted = map { $quote_field->($_, $sep) } @row_data;
			print $fh join($sep, @quoted) . "\n";
		}
	} elsif ($is_aoh) {
		my @headers;
		if ($col_names) {
			@headers = @$col_names;
		} else {
			my %col_map;
			for my $row_hash (@$data) {
				$col_map{$_} = 1 for keys %$row_hash;
			}
			@headers = sort keys %col_map;
		}
		my @header_row = @headers;
		unshift @header_row, "" if $inc_rownames;
		@header_row = map { $quote_field->($_, $sep) } @header_row;
		say $fh join($sep, @header_row);
		for my $i (0 .. $#$data) {
			my $row_hash = $data->[$i];
			my @row_data = map { defined $row_hash->{$_} ? $row_hash->{$_} : "NA" } @headers;
			unshift @row_data, $i + 1 if $inc_rownames;
			my @quoted = map { $quote_field->($_, $sep) } @row_data;
			print $fh join($sep, @quoted) . "\n";
		}
	}
	close $fh;
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
