#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':all';
use Stats::LikeR;
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Time::HiRes;
use File::Temp;

sub old_read_table {
	my $file = shift;
	die "\"$file\" is either unreadable or not a file" unless -r -f $file;
	my %args = (
		sep => ',', comment => '#',
		@_,
	);
	my %allowed_args = map {$_ => 1} (
		'comment', #character to skip lines after the header
		'row.names',
		'sep', # by default ","
		'substitutions',
		'output.type' # aoh, hoa, or hoh
	);
	my @undef_args = sort grep {!$allowed_args{$_}} keys %args;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1];
	if (scalar @undef_args > 0) {
		p @undef_args;
		die "the above args aren't defined for $current_sub";
	}
	$args{'output.type'} = $args{'output.type'} // 'aoh';
	if ($args{'output.type'} !~ m/^(?:aoh|hoa|hoh)$/) {
		die "\"$args{'output.type'}\" isn't allowed";
	}
	my (@data, %data, @header);
	open my $txt, '<', $file;
	# State machine variables for CSV/TSV parsing
	my $in_quotes = 0;
	my @line;
	my $field = '';
	my $sep_len = length($args{sep});
	while (my $line_str = <$txt>) {
		# Skip comments and empty lines ONLY if we are not inside a multiline quoted field
		if (!$in_quotes) {
			next if $line_str =~ m/^\Q$args{comment}\E/;
			next if $line_str =~ m/^\h*[\r\n]+$/; 
		}
		$line_str =~ s/\r?\n$//; # chomp with annoying invisible Windows "\r"
		# Apply substitutions if any (only on raw unquoted lines to prevent breaking the parser)
		if (!$in_quotes) {
			foreach my $sub (@{ $args{substitutions} // [] }) {
				$line_str =~ s/$sub->[0]/$sub->[1]/g;
			}
		}
		# --- PARSING MACHINE ---
		# Safely extracts fields, handling quotes, doubled-quotes, and inner separators.
		my $len = length($line_str);
		for (my $i = 0; $i < $len; $i++) {
			my $char = substr($line_str, $i, 1);
			
			if ($char eq '"') {
				if ($in_quotes && $i + 1 < $len && substr($line_str, $i + 1, 1) eq '"') {
					$field .= '"';
					$i++; # Skip the escaped second quote
				} elsif ($in_quotes) {
					$in_quotes = 0; # Close quotes
				} else {
					$in_quotes = 1; # Open quotes
				}
			} elsif (!$in_quotes && substr($line_str, $i, $sep_len) eq $args{sep}) {
				push @line, $field;
				$field = '';
				$i += $sep_len - 1; # Advance past multi-char separators
			} else {
				$field .= $char;
			}
		}
		if ($in_quotes) {
			# Line ended but quotes are still open! Append newline and fetch next line
			$field .= "\n";
			next; 
		}
		# Push the final field of the record
		push @line, $field;
		$field = '';
		if (!@header) {
			# --- HEADER PROCESSING ---
			$line[0] =~ s/^\Q$args{comment}\E// if @line && defined $line[0];
			
			while (@line && $line[-1] eq '') { pop @line }
			@header = @line; 
			
			# R-LIKE BEHAVIOR
			if ((scalar @header > 0) && ($header[0] eq '')) {
				$header[0] = 'row_name'; 
			}
			if (($args{'output.type'} eq 'hoh') && (not defined $args{'row.names'})) {
				$args{'row.names'} = $header[0];
			}
			if ((defined $args{'row.names'}) && (!grep {$_ eq $args{'row.names'}} @header)) {
				die "\"$args{'row.names'}\" isn't in the header of $file";
			}
			@line = (); # Reset for the first data line
			next;
		}
		# Check for column alignment
		if (scalar @line != scalar @header) {
			p @line;
			p @header;
			die "Alignment error on $file (" . scalar(@line) . " fields vs " . scalar(@header) . " headers).";
		}
		# --- DATA PROCESSING ---
		my %line_hash;
		for my $i (0 .. $#header) {
			my $cell = $line[$i];
			# R-like behavior: Treat completely empty fields as 'NA' 
			$line_hash{$header[$i]} = (defined($cell) && $cell eq '') ? 'NA' : $cell;
		}
		if ($args{'output.type'} eq 'aoh') {
			push @data, \%line_hash;
		} elsif ($args{'output.type'} eq 'hoa') {
			foreach my $col (@header) {
				push @{ $data{$col} }, $line_hash{$col};
			}
		} elsif ($args{'output.type'} eq 'hoh') {
			my $row_name = $line_hash{$args{'row.names'}};
			foreach my $col (@header) {
				next if $col eq $args{'row.names'}; # Bug 2: don't store the key column as a regular column
				$data{$col}{$row_name} = $line_hash{$col};
			}
		}
		@line = (); # Reset array for next record
	}
	close $txt;
	
	# Sanity check: Ensure the file didn't abruptly end with an unclosed quote
	if ($in_quotes) {
		die "Error parsing $file: Reached EOF while inside quotes.";
	}

	if ($args{'output.type'} eq 'aoh') {
		return \@data;
	} elsif ($args{'output.type'} =~ m/^(?:hoa|hoh)$/) {
		return \%data;
	}
}
my (@old, @new);
foreach my $x (0..9) {
	my $t0 = Time::HiRes::time();
	my $x = old_read_table('lib/Stats/HepatitisCdata.csv', 'output.type' => 'hoa');
	my $t1 = Time::HiRes::time();
	push @old, $t1-$t0;
}
foreach my $x (0..9) {
	my $t0 = Time::HiRes::time();
	my $x = read_table('lib/Stats/HepatitisCdata.csv', 'output.type' => 'hoa');
	my $t1 = Time::HiRes::time();
	push @new, $t1-$t0;
}
p @old;
p @new;
my $t = t_test('x' => \@old, 'y' => \@new);
p $t;
=my $fh = File::Temp->new(DIR => '/tmp', SUFFIX => '.csv', UNLINK => 1);
close $fh;
my $d = read_table('lib/Stats/HepatitisCdata.csv', 'output.type' => 'hoa');
my (@gemini, @compiled);
foreach my $x (0..9) {
	my $t0 = Time::HiRes::time();
	my $x = write_table($d, $fh->filename);
	my $t1 = Time::HiRes::time();
	push @compiled, $t1-$t0;
}
say 'Compiled mean: ' . mean(\@compiled);
foreach my $x (0..9) {
	my $t0 = Time::HiRes::time();
	my $x = gemini_write_table($d, $fh->filename);
	my $t1 = Time::HiRes::time();
	push @gemini, $t1-$t0;
}
say 'Gemini mean: ' . mean(\@gemini);
my $t = t_test('x' => \@compiled, 'y' => \@gemini);
p $t;
=my %hoa = (
	'r1' => [42, 'hello,world', undef, undef],
	'r2' => [99, undef, 'quote"here', undef],
	'r3' => [undef, "tab\tin", undef, undef],
);

write_table(
	\%hoa,
	'/tmp/hoa.tsv',
	sep => "\t"
);

my $table = read_table( '/tmp/hoa.tsv', sep => "\t", 'output.type' => 'hoa');
p $table;
foreach my $key (sort keys %hoa) {	
	my $max_i_hoa = scalar @{ $hoa{$key} } - 1;
	my $max_i_table = scalar @{ $table->{$key} } - 1;
	if ($max_i_hoa == $max_i_table) {
		say "$key has the same number of elements";
	}
	foreach my $i (0..$max_i_hoa) {
		if (
				(defined $hoa{$key}[$i])
				&&
				(defined $table->{$key}[$i])
				&&
				($hoa{$key}[$i] eq $table->{$key}[$i])
			) {
			say 'pass';
		}
	}
}
=my @s = seq(10,1,1);
p @s;
die;
my %tooth_growth = (
	dose => [qw(0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
1.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
2.0 2.0 2.0)],
	len  => [qw(4.2 11.5  7.3  5.8  6.4 10.0 11.2 11.2  5.2  7.0 16.5 16.5 15.2 17.3 22.5
17.3 13.6 14.5 18.8 15.5 23.6 18.5 33.9 25.5 26.4 32.5 26.7 21.5 23.3 29.5
15.2 21.5 17.6  9.7 14.5 10.0  8.2  9.4 16.5  9.7 19.7 23.3 23.6 26.4 20.0
25.2 25.8 21.2 14.5 27.3 25.5 26.4 22.4 24.5 24.8 30.9 26.4 27.3 29.4 23.0)],
	supp => [qw(VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC
VC VC VC VC VC OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ
OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ)]
);

my $glm_teeth = glm(
	data    => \%tooth_growth,
	formula => 'len ~ dose + supp',
	family  => 'gaussian'
);
p $glm_teeth;
=my $mtcars = json_file_to_ref('mtcars.json');
my %mtcars_hoa;
foreach my $row (@{ $mtcars }) {
	my $car = delete $row->{_row};
	push @{ $mtcars_hoa{car} }, $car;
	foreach my $key (keys %{ $row }) {
		push @{ $mtcars_hoa{$key} }, $row->{$key};
	}
}
p %mtcars_hoa;
ref_to_json_file(\%mtcars_hoa, 'mtcars.hoa.json' );
=my $mtcars = json_file_to_ref('mtcars.hoh.json');
my $lm = lm(formula =>  'mpg ~ wt + hp', data => $mtcars);
p $lm;
printf('%.3g' . "\n", $lm->{summary}{hp}{'Pr(>|t|)'});

=my $x = [1, 2, 3, 4, 5];
my $y = [2, 1, 4, 3, 5];

my @correct = (
{ # cor.test(cx, cy, alternative='two.sided', method = 'spearman', continuity=1)
	estimate  => 0.8,
	'p.value' => 0.1333333,
	statistic => 4,
},
{ #cor.test(cx, cy, alternative='two.sided', method = 'kendall', continuity=1)
	estimate  => 0.6, # tau
	'p.value' => 0.2333333,
	statistic => 
#	statistic, 
},
{	 alternative => 'two.sided',
    'conf.int' => [
       -0.279637693499009, 0.986196123450776
    ],
    estimate  =>  0.8,
    method    => 'pearson',
    p.value   => 0.104088,
    parameter => 3,
    statistic => 2.309401

);
foreach my ($idx, $meth) (indexed('spearman', 'kendall', 'pearson')) {
	my $result = cor_test(
		'x'         => $x,
		'y'         => $y,
		alternative => 'two.sided',
		method      => $meth,
		continuity  => 1
	);
	foreach my $key (sort keys %{ $result }) {
	
	}
	p $result;
}
=my $m = matrix(data => [1..9], nrow => 3, ncol => 3);
p $m;
my $mtcars = json_file_to_ref('mtcars.json');
my (@mtcars, %mtcars);
foreach my $row (@{ $mtcars }) { # transform AoH to HoA
	my $row_name = delete $row->{_row};
	foreach my $key (keys %{ $row }) {
		push @{ $mtcars{$row_name}{$key} }, $row->{$key};
	}
}
unless (-f 'mtcars.hoh.json') {
	ref_to_json_file(\%mtcars, 'mtcars.hoh.json');
}
my $binom = rbinom( n => 99, prob => 0.5, size => 1);
p $binom, array_max => scalar @{ $binom };
#my $lm = lm(formula => 'mpg ~ wt * hp^2', data => \%mtcars);
#say encode_json( \%mtcars );
# matrix
=my $mat1 = matrix(
	data => [1..6],
	nrow => 2
#	byrow => true
);
my $scale = scale($mat1);
p $scale;
say encode_json( $scale );
=say encode_json( $mat1 );
p $mat1;
# scale
my @x = 1..5;
my @scaled_x = scale($mat1);
p @scaled_x;
