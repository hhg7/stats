#!/usr/bin/env perl
# How Stats::LikeR's XS functions scale with the size of their input: write the
# fixtures, take the measurements, draw the curves.
#
# Where benchmark.pl asks "how long does each function take on one 10,000-row
# frame", this asks "what shape is that curve".  The same call is timed at a
# ladder of sizes an order of magnitude apart, seven times at each size, and the
# result is drawn against the same measurements from scale.py and scale.R.  Read
# on a log-log axis, the slope of the line is the exponent: flat is O(1), 45
# degrees is O(n), steeper is worse.
#
# Three stages.  With no switch all three run, in this order; name one or more
# to run only those.
#
#     --data      write the fixtures the I/O panels read
#     --measure   time Stats::LikeR                 -> perl_scaling.tsv
#     --plot      draw the three .tsv files         -> scaling.*.svg
#
# Python and R are measured by their own programs, so a full run is:
#
#     perl plot.scaling.pl --data                  # once, writes the fixtures
#     perl -Iblib/arch -Iblib/lib \
#          plot.scaling.pl --measure               # -> perl_scaling.tsv
#     python3 scale.py                             # -> python_scaling.tsv
#     Rscript scale.R                              # -> r_scaling.tsv
#     perl plot.scaling.pl --plot                  # -> scaling.*.svg
#
# --measure writes any fixture that is missing before it starts, so --data is
# only worth naming on its own when scale.py or scale.R will run first: those
# two stop with instructions rather than inventing their own copy of the data.
#
# Environment:
#
#     SCALE_DIR    where the fixtures live (/tmp/likeR.scaling)
#     SCALE_RUNS   runs per (function, size); default 7
#     SCALE_CAP    seconds; once one run of a function takes longer than this,
#                  that function is not tried at any larger size.  Default 4.
#     SCALE_MAX_N  hard ceiling on the row count, for a quick partial run
#     SCALE_TARGET seconds a single measurement should span; a call faster than
#                  this is repeated until it does.  Default 0.002.
#     SCALE_CPU    the CPU --measure pins itself to; default 0, "none" to let
#                  the scheduler place it.
#
# On a hybrid CPU -- Intel's P-core/E-core parts, and the big.LITTLE ARM
# designs -- a forked child landing on an E-core reads 1.5 to 2 times slower
# than the identical loop on a P-core, which is a wider band than most of the
# differences these plots are drawn to show.  It is worst in --measure, because
# the fork per run gives the scheduler a fresh chance to place the work
# somewhere else every time, so --measure now pins itself to one CPU rather
# than asking to be started under "taskset"; see pin_to_one_cpu() below.
#
# scale.py and scale.R pin themselves the same way, and all three default to
# CPU 0, so the three languages land on the same core without anyone having to
# remember.  SCALE_CPU moves it, and has to be given the same value in all
# three.
#
# Beyond Stats::LikeR itself, --data and --measure need only core modules;
# --plot also needs Matplotlib::Simple.
#
# What is being compared, and what is not:
#
#   * Every panel asks the three languages for the *same result* by each one's
#     idiomatic route, which is the rule benchmark.pl already follows.  Where
#     that is impossible the panel is simply missing that language rather than
#     being filled with something else: R has no skew(), no kurtosis() and no
#     row-record frame to build, so it does not appear in those four panels.
#   * The read_table and write_table panels all read the byte-identical files
#     --data wrote, so no reader is handed a quoting or type-guessing job the
#     others were spared.
#   * Absolute heights are worth less than slopes here.  A panel where one line
#     sits below another by a constant factor over four decades is telling you
#     about per-element cost; a panel where the lines converge or cross is
#     telling you about fixed overhead, which is the more interesting finding
#     and the one a single-size benchmark cannot show.
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Path ();
use File::Spec ();
use Getopt::Long ();
use IO::Compress::Zip ();
use POSIX ();
use Time::HiRes ();
use Stats::LikeR;

# Matplotlib::Simple is loaded by the --plot stage rather than here; see the
# note beside the require in plot_all().

# ---------------------------------------------------------------------------
# 0. Stages and settings
# ---------------------------------------------------------------------------
my $usage = <<"USAGE";
usage: perl [-Iblib/arch -Iblib/lib] $0 [--data] [--measure] [--plot]

  --data      write the read_table/write_table fixtures into \$SCALE_DIR
  --measure   time Stats::LikeR at every size  -> perl_scaling.tsv
  --plot      draw perl/python/r_scaling.tsv   -> scaling.*.svg

With no switch, all three run in that order.  See the comment at the top of
this file for the environment variables and for how scale.py and scale.R fit
in.
USAGE

my ($do_data, $do_measure, $do_plot, $want_help) = (0, 0, 0, 0);
Getopt::Long::GetOptions(
	'data!'    => \$do_data,
	'measure!' => \$do_measure,
	'plot!'    => \$do_plot,
	'help|h'   => \$want_help,
) or die $usage;
if ($want_help) {
	print $usage;
	exit 0;
}
($do_data, $do_measure, $do_plot) = (1, 1, 1)
	unless $do_data || $do_measure || $do_plot;

my $dir      = $ENV{SCALE_DIR}  || '/tmp/likeR.scaling';
my $runs     = $ENV{SCALE_RUNS} || 7;
my $cap      = defined $ENV{SCALE_CAP} ? $ENV{SCALE_CAP} : 4;
my $max_n    = $ENV{SCALE_MAX_N} || 0;
my $target   = defined $ENV{SCALE_TARGET} ? $ENV{SCALE_TARGET} : 0.002;
my $max_reps = 10_000;

# ---------------------------------------------------------------------------
# 1. Size ladders
# ---------------------------------------------------------------------------
# Half-decade steps, so seven points span three decades and the log-log slope
# is estimated from evenly spaced x.  scale.py and scale.R carry the same three
# lists; change one, change all three.
my @vec_n   = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000);
my @io_n    = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000);
my @frame_n = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000);

# The fixture row counts are @io_n before SCALE_MAX_N is applied, and stay that
# way: the files are shared with scale.py and scale.R, which have their own
# ladders and their own idea of a partial run, so a truncated --measure must not
# leave them a truncated corpus.
my @fixture_n = @io_n;

if ($max_n) {
	@vec_n   = grep { $_ <= $max_n } @vec_n;
	@io_n    = grep { $_ <= $max_n } @io_n;
	@frame_n = grep { $_ <= $max_n } @frame_n;
}

# ---------------------------------------------------------------------------
# 2. The fixtures
# ---------------------------------------------------------------------------
# The three languages must read *byte-identical* files, otherwise the
# read_table panels compare three different parsing jobs, so exactly one
# program writes them and the other two only read.
#
# Four shapes are written at every row count in @fixture_n:
#
#   num.<n>.csv   a,b,c,d              four numeric columns, nothing to quote
#   mix.<n>.csv   id,x,y,cat1,cat2     an integer, two numerics, two strings
#   mix.<n>.tsv   the same table, tab separated
#   mix.<n>.xlsx  the same table again, as a spreadsheet
#
# Numbers are printed with a fixed "%.6f", and no field ever contains the
# separator or a quote, so every reader sees the same characters and none of
# them is dragged into a quoting slow path the others avoid.
#
# Nothing in this stage calls Stats::LikeR: the fixtures describe the job, not
# the module under test, and writing them through the code being timed would
# make the corpus move whenever that code did.  That rule is what decides how
# the .xlsx is written.  write_table() can write one, but write_table's .xlsx
# output is itself one of the panels below, so a fixture written that way would
# move the read_table (xlsx) curve every time the writer changed -- and the two
# would stop being independent measurements.  So the .xlsx is assembled here
# instead, from core IO::Compress::Zip, as the five-member package
# xlsx_writer() describes.  Verified read by all three: Stats::LikeR's
# read_table, pandas' read_excel (openpyxl) and R's readxl::read_excel.
#
# The pseudo-random numbers come from an xorshift64 written out longhand rather
# than from perl's rand(), so the fixtures are reproducible on any perl and any
# platform: drand48 and the various rand() implementations do not agree, and a
# file whose contents depend on the build is not a fixture.
#
# The state is kept as two 32-bit halves so this works unchanged on a perl
# without 64-bit integers, which is the one thing a plain "$s ^= $s << 13"
# would quietly get wrong.
my ($hi, $lo) = (0x12345678, 0x9ABCDEF0);

sub xorshift {
	# s ^= s << 13
	my $h = (($hi << 13) | ($lo >> 19)) & 0xFFFFFFFF;
	my $l = ($lo << 13) & 0xFFFFFFFF;
	($hi, $lo) = ($hi ^ $h, $lo ^ $l);
	# s ^= s >> 7
	$h = $hi >> 7;
	$l = (($lo >> 7) | ($hi << 25)) & 0xFFFFFFFF;
	($hi, $lo) = ($hi ^ $h, $lo ^ $l);
	# s ^= s << 17
	$h = (($hi << 17) | ($lo >> 15)) & 0xFFFFFFFF;
	$l = ($lo << 17) & 0xFFFFFFFF;
	($hi, $lo) = ($hi ^ $h, $lo ^ $l);
	return ($hi, $lo);
}

# a uniform in (0, 1), from the top 32 bits
sub unif {
	my ($h) = xorshift();
	return ($h + 0.5) / 4_294_967_296.0;
}

# one standard normal per call; the second Box-Muller value is cached
my $spare;
sub norm {
	if (defined $spare) {
		my $s = $spare;
		undef $spare;
		return $s;
	}
	my $u = unif();
	my $v = unif();
	my $r = sqrt(-2 * log($u));
	$spare = $r * sin(2 * 3.14159265358979323846 * $v);
	return $r * cos(2 * 3.14159265358979323846 * $v);
}

# The four files one row count needs, in the order the panels use them.
sub fixture_files {
	my ($n) = @_;
	return (
		File::Spec->catfile($dir, "num.$n.csv"),
		File::Spec->catfile($dir, "mix.$n.csv"),
		File::Spec->catfile($dir, "mix.$n.tsv"),
		File::Spec->catfile($dir, "mix.$n.xlsx"),
	);
}

# --- the .xlsx fixture ---------------------------------------------------
# An .xlsx is a zip of XML parts, and the five below are the whole of a
# single-sheet workbook: the content-type map, the package relationships, the
# workbook and its relationships, and the sheet.  It is the same set
# Stats::LikeR's own writer emits, minus its docProps/core.xml -- that part
# carries a timestamp and the writing program's name, neither of which belongs
# in a fixture.
#
# Everything here is written by hand rather than through a spreadsheet module
# for the reason given at the top of this section: the fixture may not come out
# of the code being timed.  Nothing about it is clever; the point is that all
# three readers see one file.
#
# Strings are inline (t="inlineStr"), not shared-string indices.  Both are
# ordinary xlsx and every reader takes either, and inline strings keep this a
# single streaming pass: a sharedStrings.xml has to know every distinct string
# before the first row can be written, which for 300,000 rows means holding the
# table in memory to write a file that is only being written so it can be read.
my $XLSX_EPOCH = 946_684_800;	# 2000-01-01 UTC

# The zip's own per-member timestamp is the one byte-level difference between
# two machines: DOS times are local, so a fixed epoch still records the
# writer's timezone.  Nothing reads it, and the three languages always read the
# copy written on the machine they are running on, so it costs nothing -- but
# it does mean the .xlsx is only reproducible byte for byte within a timezone,
# where the .csv and .tsv are reproducible everywhere.
my @xlsx_parts = (
	'[Content_Types].xml' =>
		'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
	  . '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
	  . '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
	  . '<Default Extension="xml" ContentType="application/xml"/>'
	  . '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
	  . '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
	  . '</Types>',
	'_rels/.rels' =>
		'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
	  . '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
	  . '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
	  . '</Relationships>',
	'xl/workbook.xml' =>
		'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
	  . '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
	  . ' xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
	  . '<sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>',
	'xl/_rels/workbook.xml.rels' =>
		'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
	  . '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
	  . '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
	  . '</Relationships>',
);

# A1, B1, ... -- the fixture is five columns wide, so a handful is plenty.
my @xlsx_col = ('A' .. 'Z');

# Open $path, write the four fixed parts and the sheet's preamble, and hand
# back the zip stream.  The caller writes rows into it with xlsx_row() and ends
# with xlsx_close(); the sheet is streamed rather than assembled, so writing
# 300,000 rows costs no more memory than writing one.
sub xlsx_writer {
	my ($path) = @_;
	my $zip;
	for (my $i = 0; $i < @xlsx_parts; $i += 2) {
		my ($name, $body) = @xlsx_parts[ $i, $i + 1 ];
		if ($zip) {
			$zip->newStream(Name => $name, Time => $XLSX_EPOCH, ExtAttr => 0)
				or die "cannot write \"$name\" in \"$path\": "
				     . "$IO::Compress::Zip::ZipError\n";
		} else {
			$zip = IO::Compress::Zip->new($path, Name => $name,
				Time => $XLSX_EPOCH, ExtAttr => 0)
				or die "cannot write \"$path\": "
				     . "$IO::Compress::Zip::ZipError\n";
		}
		print {$zip} $body or die "cannot write \"$name\": $!\n";
	}
	$zip->newStream(Name => 'xl/worksheets/sheet1.xml',
		Time => $XLSX_EPOCH, ExtAttr => 0)
		or die "cannot write the worksheet in \"$path\": "
		     . "$IO::Compress::Zip::ZipError\n";
	print {$zip} '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
	  . '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
	  . '<sheetData>';
	return $zip;
}

# One row, 1-based like the spreadsheet.  $str says, per column, whether the
# cell is text: a number goes in bare and a string goes in as an inline string,
# which is the difference every reader keys the column's type off.  Nothing in
# this fixture needs XML escaping -- the values are digits, "." and the
# single letters A/B/C and X/Y -- and the generator below is the only caller.
sub xlsx_row {
	my ($zip, $r, $str, $vals) = @_;
	my $out = qq{<row r="$r">};
	for my $c (0 .. $#$vals) {
		my $ref = $xlsx_col[$c] . $r;
		$out .= $str->[$c]
			? qq{<c r="$ref" t="inlineStr"><is><t>$vals->[$c]</t></is></c>}
			: qq{<c r="$ref"><v>$vals->[$c]</v></c>};
	}
	print {$zip} $out, '</row>' or die "cannot write row $r: $!\n";
	return;
}

sub xlsx_close {
	my ($zip, $path) = @_;
	print {$zip} '</sheetData></worksheet>' or die "cannot write \"$path\": $!\n";
	$zip->close or die "cannot close \"$path\": $!\n";
	return;
}

sub write_fixtures {
	File::Path::make_path($dir) unless -d $dir;

	my @cat1 = qw(A B C);
	my @cat2 = qw(X Y);

	foreach my $n (@fixture_n) {
		my ($num, $csv, $tsv, $xlsx) = fixture_files($n);

		# --- num.<n>.csv: four numeric columns -----------------------------
		open my $fh, '>', $num or die "cannot write \"$num\": $!\n";
		print {$fh} "a,b,c,d\n";
		for my $i (1 .. $n) {
			printf {$fh} "%.6f,%.6f,%.6f,%.6f\n",
				norm(), norm() * 2 + 5, unif() * 100, norm() * 0.5 - 3;
		}
		close $fh or die "cannot close \"$num\": $!\n";

		# --- mix.<n>.{csv,tsv,xlsx}: one table, three containers ----------
		# All three are written from one pass so they are the same table and
		# not three independent draws.  The .xlsx carries the values as the
		# same text the .csv does -- "%.6f" for the two numerics, the integer
		# bare -- so a reader that parses both is parsing identical numbers and
		# only the container differs.
		open my $c, '>', $csv or die "cannot write \"$csv\": $!\n";
		open my $t, '>', $tsv or die "cannot write \"$tsv\": $!\n";
		my $z = xlsx_writer($xlsx);
		# id, x, y are numbers; cat1, cat2 are text.  The header is all text.
		my @is_str = (0, 0, 0, 1, 1);
		print {$c} "id,x,y,cat1,cat2\n";
		print {$t} "id\tx\ty\tcat1\tcat2\n";
		xlsx_row($z, 1, [ (1) x 5 ], [qw(id x y cat1 cat2)]);
		for my $i (1 .. $n) {
			my $x  = norm();
			my $y  = norm() * 2 + 5;
			my $c1 = $cat1[ int(unif() * 3) ];
			my $c2 = $cat2[ int(unif() * 2) ];
			my ($xs, $ys) = (sprintf('%.6f', $x), sprintf('%.6f', $y));
			printf {$c} "%d,%s,%s,%s,%s\n",   $i, $xs, $ys, $c1, $c2;
			printf {$t} "%d\t%s\t%s\t%s\t%s\n", $i, $xs, $ys, $c1, $c2;
			xlsx_row($z, $i + 1, \@is_str, [ $i, $xs, $ys, $c1, $c2 ]);
		}
		close $c or die "cannot close \"$csv\": $!\n";
		close $t or die "cannot close \"$tsv\": $!\n";
		xlsx_close($z, $xlsx);

		printf "%8d rows: %s, %s, %s, %s\n", $n, $num, $csv, $tsv, $xlsx;
	}

	printf "Fixtures written to %s\n", $dir;
}

# ---------------------------------------------------------------------------
# 3. Input builders
# ---------------------------------------------------------------------------
# One builder per figure, called once per size.  Whatever it returns is handed
# to every benchmark in that figure, and is built in the parent so that the
# forked children below all measure the same object rather than each paying to
# construct it.
my %build = (
	vector => sub {
		my ($n) = @_;
		srand(42);
		return {
			x     => rnorm(n => $n, mean => 0, sd => 1),
			y     => rnorm(n => $n, mean => 5, sd => 2),
			label => rbinom(n => $n, prob => 0.5, size => 1),
			n     => $n,
		};
	},
	io => sub {
		my ($n) = @_;
		my ($num_csv, $mix_csv, $mix_tsv, $mix_xlsx) = fixture_files($n);
		my %f = (
			num_csv  => $num_csv,
			mix_csv  => $mix_csv,
			mix_tsv  => $mix_tsv,
			mix_xlsx => $mix_xlsx,
			out      => File::Spec->catfile($dir, "out.perl.$$.tmp"),
			# a separate target: write_table picks the writer off the
			# extension, so the .xlsx panel needs a name ending in .xlsx
			out_xlsx => File::Spec->catfile($dir, "out.perl.$$.xlsx"),
		);
		# The frames handed to write_table are read in, not synthesized, so the
		# three languages all write out the same table.  hoa is the columnar
		# shape pandas and R hold natively; aoa is the shape write_table docs
		# lead with.
		$f{hoa} = read_table($f{mix_csv}, 'output.type' => 'hoa');
		my @cols = qw(id x y cat1 cat2);
		$f{aoa} = [ \@cols,
			map { my $i = $_; [ map { $f{hoa}{$_}[$i] } @cols ] }
			0 .. $#{ $f{hoa}{id} } ];
		return \%f;
	},
	frame => sub {
		my ($n) = @_;
		srand(42);
		my $df = {
			x      => rnorm(n => $n, mean => 0, sd => 1),
			y      => rnorm(n => $n, mean => 5, sd => 2),
			cat1   => [ map { (qw(A B C))[ int rand 3 ] } 1 .. $n ],
			cat2   => [ map { (qw(X Y))[ int rand 2 ] } 1 .. $n ],
			binary => rbinom(n => $n, prob => 0.5, size => 1),
		};
		# merge needs a key column: a Stats::LikeR frame has neither row.names
		# nor an index, so all three scripts join on an explicit id, as
		# benchmark.pl does.
		my $df_id = { %$df, id => [ 1 .. $n ] };
		my @cols  = qw(x y cat1 cat2 binary);
		# The transpose input is the three numeric columns only, because
		# numpy's counterpart is a numeric matrix and R's is as.matrix(); a
		# five-column AoA with strings in it would be a wider job than either.
		my @num = qw(x y binary);
		return {
			df    => $df,
			df_id => $df_id,
			aoh   => [ map { my $i = $_; { map { $_ => $df->{$_}[$i] } @cols } }
			           0 .. $n - 1 ],
			aoa   => [ map { my $i = $_; [ map { $df->{$_}[$i] } @num ] }
			           0 .. $n - 1 ],
			n     => $n,
		};
	},
);

# ---------------------------------------------------------------------------
# 4. The benchmarks
# ---------------------------------------------------------------------------
# figure  which output image the panel belongs to, and which size ladder and
#         builder it uses
# name    the panel title, and the key scale.py and scale.R must agree with
# call    what was actually called, recorded in the file for the reader
# code    the timed body, handed the builder's return value
my @benchmarks = (
	# --- reductions over one numeric vector --------------------------------
	{ figure => 'vector', name => 'sum',      call => 'sum($x)',
	  code => sub { sum($_[0]{x}) } },
	{ figure => 'vector', name => 'min',      call => 'min($x)',
	  code => sub { min($_[0]{x}) } },
	{ figure => 'vector', name => 'max',      call => 'max($x)',
	  code => sub { max($_[0]{x}) } },
	{ figure => 'vector', name => 'mean',     call => 'mean($x)',
	  code => sub { mean($_[0]{x}) } },
	{ figure => 'vector', name => 'median',   call => 'median($x)',
	  code => sub { median($_[0]{x}) } },
	{ figure => 'vector', name => 'sd',       call => 'sd($x)',
	  code => sub { sd($_[0]{x}) } },
	{ figure => 'vector', name => 'var',      call => 'var($x)',
	  code => sub { var($_[0]{x}) } },
	{ figure => 'vector', name => 'quantile', call => 'quantile(x => $x, probs => [.25,.5,.75])',
	  code => sub { quantile(x => $_[0]{x}, probs => [ 0.25, 0.5, 0.75 ]) } },
	{ figure => 'vector', name => 'cor',      call => 'cor($x, $y)',
	  code => sub { cor($_[0]{x}, $_[0]{y}) } },
	{ figure => 'vector', name => 'cov',      call => 'cov($x, $y)',
	  code => sub { cov($_[0]{x}, $_[0]{y}) } },
	# skew and kurtosis have no base-R equivalent (only e1071/moments, whose
	# type= conventions differ), so those two panels have two lines, not three.
	{ figure => 'vector', name => 'skew',     call => 'skew($x)',
	  code => sub { skew($_[0]{x}) } },
	{ figure => 'vector', name => 'kurtosis', call => 'kurtosis($x)',
	  code => sub { kurtosis($_[0]{x}) } },

	# --- transforms that return something the size of their input ----------
	{ figure => 'transform', name => 'rank',   call => 'rank($x)',
	  code => sub { my @r = rank($_[0]{x}); return \@r } },
	{ figure => 'transform', name => 'uniq',   call => 'uniq($x)',
	  code => sub { my @u = uniq($_[0]{x}); return \@u } },
	{ figure => 'transform', name => 'scale',  call => 'scale($x)',
	  code => sub { my @s = scale($_[0]{x}); return \@s } },
	{ figure => 'transform', name => 'sample', call => 'sample($x, n/10)',
	  code => sub { sample($_[0]{x}, int($_[0]{n} / 10) + 1) } },
	{ figure => 'transform', name => 'seq',    call => 'seq(1, n)',
	  code => sub { my @s = seq(1, $_[0]{n}); return \@s } },
	{ figure => 'transform', name => 'auc',    call => 'auc($y, $label)',
	  code => sub { auc($_[0]{y}, $_[0]{label}) } },

	# --- read_table and write_table, over four inputs each -----------------
	# read_table's default output.type is 'aoh', one hash per row; 'hoa' is the
	# columnar shape a DataFrame and a data.frame already are.  Both are timed
	# on the same file so the difference between the two panels is the cost of
	# the shape and nothing else.
	{ figure => 'io', name => 'read_table (csv, numeric)',
	  call => "read_table('num.csv')",
	  code => sub { read_table($_[0]{num_csv}) } },
	{ figure => 'io', name => 'read_table (csv, mixed)',
	  call => "read_table('mix.csv')",
	  code => sub { read_table($_[0]{mix_csv}) } },
	{ figure => 'io', name => 'read_table (tsv, mixed)',
	  call => "read_table('mix.tsv')",
	  code => sub { read_table($_[0]{mix_tsv}) } },
	{ figure => 'io', name => 'read_table (csv, hoa)',
	  call => "read_table('mix.csv', 'output.type' => 'hoa')",
	  code => sub { read_table($_[0]{mix_csv}, 'output.type' => 'hoa') } },
	# The same table as the three panels above, out of a spreadsheet instead of
	# a text file, so what this panel shows against them is what the container
	# costs.  All three languages reach for a different implementation here --
	# read_table unzips and scans the sheet XML itself, pandas goes through
	# openpyxl, R through readxl's C++ parser -- which is the point: there is
	# no "the xlsx reader", and a module that ships one is competing with
	# whichever the other two ecosystems settled on.
	{ figure => 'io', name => 'read_table (xlsx)',
	  call => "read_table('mix.xlsx')",
	  code => sub { read_table($_[0]{mix_xlsx}) } },
	{ figure => 'io', name => 'write_table (csv, hoa)',
	  call => "write_table(\$hoa, file, 'row.names' => 0)",
	  code => sub { write_table($_[0]{hoa}, $_[0]{out}, 'row.names' => 0) } },
	{ figure => 'io', name => 'write_table (tsv, hoa)',
	  call => "write_table(\$hoa, file, sep => \"\\t\", 'row.names' => 0)",
	  code => sub { write_table($_[0]{hoa}, $_[0]{out}, sep => "\t",
	                            'row.names' => 0) } },
	{ figure => 'io', name => 'write_table (csv, aoa)',
	  call => "write_table(\$aoa, file, 'row.names' => 0)",
	  code => sub { write_table($_[0]{aoa}, $_[0]{out}, 'row.names' => 0) } },
	{ figure => 'io', name => 'write_table (csv, row.names)',
	  call => "write_table(\$hoa, file, 'row.names' => 1)",
	  code => sub { write_table($_[0]{hoa}, $_[0]{out}, 'row.names' => 1) } },
	# write_table picks the writer off the extension, so this is the same call
	# as the first write panel with a different file name.  R is absent from
	# this one: base R cannot write an .xlsx and neither can readxl, which is a
	# reader by design, so scale.R has nothing idiomatic to offer -- see the
	# note beside the read panel in that file.
	{ figure => 'io', name => 'write_table (xlsx, hoa)',
	  call => "write_table(\$hoa, file.xlsx, 'row.names' => 0)",
	  code => sub { write_table($_[0]{hoa}, $_[0]{out_xlsx}, 'row.names' => 0) } },

	# --- whole-frame operations --------------------------------------------
	{ figure => 'frame', name => 'filter', call => 'filter($df, col(\'x\') > 0)',
	  code => sub { filter($_[0]{df}, col('x') > 0) } },
	# select_cols on a HoA hands back the same array references rather than
	# copying them, so this panel is flat: it measures picking two keys out of
	# a hash and nothing else.  scale.py's counterpart is therefore the pandas
	# view, not a .copy(); R has no view to offer and its subset materializes,
	# which is the difference the panel is there to show.
	{ figure => 'frame', name => 'select_cols', call => "select_cols(\$df, ['x','cat1'])",
	  code => sub { select_cols($_[0]{df}, [ 'x', 'cat1' ]) } },
	# group_by only splits, so the mean-per-group half is timed with it, which
	# is what dplyr's group_by %>% summarise and pandas' groupby().mean() do.
	# All three average one column, not every numeric one.
	{ figure => 'frame', name => 'group_by + mean', call => 'group_by + mean',
	  code => sub {
		my $groups = group_by($_[0]{df}, 'x', 'cat1');
		return { map { $_ => mean($groups->{$_}) } keys %$groups };
	  } },
	{ figure => 'frame', name => 'merge', call => "merge(\$df_id, \$df_id, how => 'inner', on => 'id')",
	  code => sub { merge($_[0]{df_id}, $_[0]{df_id}, how => 'inner', on => 'id') } },
	{ figure => 'frame', name => 'value_counts', call => "value_counts(\$df, 'cat1')",
	  code => sub { value_counts($_[0]{df}, 'cat1') } },
	{ figure => 'frame', name => 'drop_duplicates', call => 'drop_duplicates($df)',
	  code => sub { drop_duplicates($_[0]{df}) } },
	{ figure => 'frame', name => 'transpose', call => 'transpose($aoa)',
	  code => sub { transpose($_[0]{aoa}) } },
	# aoh2hoa turns one hash per row into one array per column.  numpy and
	# pandas do this all the time (pd.DataFrame(records)); base R has no
	# row-record frame to convert from, so R is absent from this panel.
	{ figure => 'frame', name => 'aoh2hoa', call => 'aoh2hoa($aoh)',
	  code => sub { aoh2hoa($_[0]{aoh}) } },
);

my %ladder = (vector => \@vec_n, transform => \@vec_n,
              io => \@io_n, frame => \@frame_n);
my %builder_for = (vector => 'vector', transform => 'vector',
                   io => 'io', frame => 'frame');

# The order the figures are measured in, and the order their images come out in.
my @figure_order = qw(vector transform io frame);

my %figure_title = (
	vector    => 'Reductions over one numeric vector',
	transform => 'Transforms that return as much as they are given',
	io        => 'read_table and write_table',
	frame     => 'Whole-frame operations',
);

# ---------------------------------------------------------------------------
# 5. Execution engine
# ---------------------------------------------------------------------------
# The CPU the whole stage runs on, fixed before the first measurement.
#
# Saying "run this under taskset" in a comment is not the same as running it
# under taskset.  An unpinned rerun of an identical build read min() at 2.75 ms
# against 1.43 pinned, and mean() at 2.61 against 1.35 -- a factor of 1.9 on
# this 13th-generation Core i7, which is wider than nearly every difference
# these curves are drawn to show, and it lands on whichever functions happened
# to draw an E-core rather than on all of them evenly.  A reading that moves by
# 1.9x depending on where the scheduler put the fork is not a measurement.
#
# sched_setaffinity() is in neither POSIX nor any core module, and its syscall
# number differs by architecture, so this goes through taskset(1) -- the tool
# the header used to ask for by hand -- applied to this process's own pid.  Not
# a re-exec under taskset: that would lose the "-Iblib/arch -Iblib/lib" the
# caller put on the command line, which is how this script is always run
# against a build that is not installed.  The forked children inherit the
# affinity, which is the point; taskset is called once, not once per run.
#
# Every failure here is survivable and none of them stop the run:
#
#   * No /proc/self/status -- macOS, the BSDs, Windows -- means the affinity
#     cannot be read, let alone set, so the stage says so once and carries on.
#     There is no taskset on those platforms either.
#   * A taskset that is absent or refuses leaves the run on the CPUs it had,
#     with its complaint repeated verbatim.
#   * A process that is already pinned to a single CPU is left on it, so
#     "taskset -c 5 perl plot.scaling.pl --measure" keeps CPU 5 and this does
#     not overrule the operator.
#
# CPU 0 is the default because the P-cores are enumerated first on the hybrid
# parts, so it is a fast core on the machines where the choice matters at all,
# and because it is what the other two languages get told to use.
sub cpus_allowed {
	open my $fh, '<', '/proc/self/status' or return undef;
	while (my $line = <$fh>) {
		return $1 if $line =~ /^Cpus_allowed_list:\s*(\S+)/;
	}
	return undef;	#no such field: a kernel too old to report it
}

sub pin_to_one_cpu {
	my $want = defined $ENV{SCALE_CPU} ? $ENV{SCALE_CPU} : 0;
	if ($want eq 'none' || $want eq '') {
		print "SCALE_CPU=none: leaving the scheduler to place the measurements\n";
		return;
	}
	$want =~ /\A[0-9]+\z/
		or die "SCALE_CPU must be a CPU number or 'none', not '$want'\n";

	my $before = cpus_allowed();
	unless (defined $before) {
		print "cannot read this process's CPU affinity, so the measurements "
		    . "are not pinned; on a hybrid CPU, start this under the "
		    . "platform's affinity tool\n";
		return;
	}
	#a list with no comma and no dash is one CPU, so someone has already pinned us
	if ($before !~ /[,\-]/) {
		printf "already pinned to CPU %s; keeping it\n", $before;
		return;
	}

	my $said = qx{taskset -pc $want $$ 2>&1};
	$said = '' unless defined $said;
	my $after = cpus_allowed();
	if (defined $after && $after eq $want) {
		printf "pinned to CPU %s (was %s); run scale.py and scale.R under "
		     . "taskset -c %s to match\n", $after, $before, $want;
	} else {
		$said =~ s/\s+/ /g;
		$said =~ s/\s+\z//;
		printf STDERR "could not pin to CPU %s, staying on %s%s\n",
			$want, $before, length($said) ? ": $said" : '';
	}
}

# Each run happens in a forked child that does nothing else.  benchmark.pl
# forks so that its memory figures mean something; here it is for a different
# reason: a timed loop in one process lets an earlier run leave the later ones
# a grown arena, a warmed allocator and a partly numified input, and those
# effects grow with n, which is exactly the axis being measured.  A child that
# starts from the parent's untouched heap every time cannot do that.
#
# Inside the child the call is made once untimed before the clock starts.  That
# is not politeness towards the cache: a freshly forked process has a
# copy-on-write address space, and the first dozen writes perl makes -- its
# stack, the eval context, the first temporaries -- each take a page fault that
# copies 4 KiB.  Measured, that is a fixed ~50 microseconds sitting on top of
# every reading, which at n = 1,000 is ten times the call itself and would draw
# a Perl line that looks flat up to n = 10,000 for reasons that have nothing to
# do with Stats::LikeR.  The untimed call pays it, and what is then measured is
# the steady-state cost -- the same quantity scale.py and scale.R measure.
#
# A call faster than $target is repeated until the pair of clock readings spans
# $target and the total divided by the count, because Time::HiRes resolves
# microseconds and seq() at n = 1,000 does not take one.
#
# The report comes back down a pipe of its own rather than down the child's
# STDOUT, and the child's STDOUT is sent to the bit bucket, because some of the
# functions being timed are chatty -- write_table announces the file it wrote.
# Sharing one channel would splice that announcement into the timing, and
# leaving it connected to the terminal would put a few hundred lines of noise
# in the middle of the progress report and charge the write for it.
sub measure {
	my ($code, $input) = @_;

	pipe(my $from_child, my $to_parent) or die "pipe failed: $!\n";

	my $pid = fork();
	die "fork failed: $!\n" unless defined $pid;

	if (!$pid) {                                  # the child
		close $from_child;
		select((select($to_parent), $| = 1)[0]);  # _exit does not flush
		open my $saved_out, '>&', \*STDOUT or POSIX::_exit(1);
		open STDOUT, '>', File::Spec->devnull() or POSIX::_exit(1);

		# untimed: page faults, allocator, page cache, and a broken call
		my $ok  = eval { $code->($input); 1 };
		my $err = $ok ? '' : $@;

		my $secs = 0;
		my $reps = 1;
		if ($ok) {
			# A second untimed call, this one clocked, sizes the repeat count.
			# The first cannot: it carries the fork's page faults, so it reads
			# five to ten times high and would choose far too few repeats.
			my $cal = Time::HiRes::time();
			$code->($input);
			my $one = Time::HiRes::time() - $cal;

			$reps = $one > 0 ? int($target / $one) + 1 : $max_reps;
			$reps = $max_reps if $reps > $max_reps;
			my $start = Time::HiRes::time();
			$code->($input) for 1 .. $reps;
			$secs = (Time::HiRes::time() - $start) / $reps;
		}
		$err =~ s/\s+/ /g;

		open STDOUT, '>&', $saved_out or POSIX::_exit(1);
		print {$to_parent} join("\t", $secs, $reps, $err), "\n";
		POSIX::_exit(0);
	}

	close $to_parent;
	my $line = <$from_child>;
	close $from_child;
	waitpid $pid, 0;
	return (undef, undef, 'the child died without reporting') unless defined $line;

	chomp $line;
	my ($secs, $reps, $err) = split /\t/, $line, 3;
	return ($secs, $reps, (defined $err && length $err) ? $err : undef);
}

# ---------------------------------------------------------------------------
# 6. The --measure stage
# ---------------------------------------------------------------------------
my $perl_tsv = 'perl_scaling.tsv';

sub measure_all {
	my @results;
	my %too_slow; # name => 1 once it exceeds $cap, or once it fails

	# Grouped by figure, then by size, so each input is built once and every
	# benchmark that wants it is measured before it is thrown away.  The
	# 1,000,000 element frames are large enough that holding two ladders' worth
	# at once is worth avoiding.
	my %by_figure;
	push @{ $by_figure{ $_->{figure} } }, $_ for @benchmarks;

	foreach my $figure (@figure_order) {
		my $list = $by_figure{$figure} or next;
		foreach my $n (@{ $ladder{$figure} }) {
			my @todo = grep { !$too_slow{ $_->{name} } } @$list;
			next unless @todo;

			my $input = $build{ $builder_for{$figure} }->($n);
			foreach my $b (@todo) {
				my ($slowest, $reps_used, $failed) = (0, 0, 0);
				for my $run (0 .. $runs - 1) {
					my ($secs, $reps, $err) = measure($b->{code}, $input);
					if (defined $err) {
						# Report it once and stop trying this function at every
						# larger size too: a call that dies at n = 1,000 is not
						# going to start working at n = 1,000,000.
						printf STDERR "%s at n=%d: %s\n", $b->{name}, $n, $err;
						$too_slow{ $b->{name} } = 1;
						$failed = 1;
						last;
					}
					$slowest   = $secs if $secs > $slowest;
					$reps_used = $reps;
					push @results, [ $figure, $b->{name}, $b->{call}, $n, $run, $secs ];
				}
				next if $failed;
				printf "%-9s %-30s n=%-8d %.6f s%s\n", $figure, $b->{name}, $n,
					$slowest, $reps_used > 1 ? " (x$reps_used)" : '';
				$too_slow{ $b->{name} } = 1 if $slowest > $cap;
			}
			undef $input;
		}
	}

	# Anything the cap stopped is worth saying out loud: a curve that ends early
	# is not the same as a curve that was never measured, and the plot cannot
	# tell you which one you are looking at.
	if (%too_slow) {
		printf "Stopped early (a run exceeded %g s, or it failed): %s\n",
			$cap, join(', ', sort keys %too_slow);
	}

	# clean up the write_table targets
	for my $out (File::Spec->catfile($dir, "out.perl.$$.tmp"),
	             File::Spec->catfile($dir, "out.perl.$$.xlsx")) {
		unlink $out if -f $out;
	}

	write_table(
		[ [ 'figure', 'function', 'call', 'n', 'run', 'seconds' ], @results ],
		$perl_tsv,
		sep         => "\t",
		'row.names' => 0,
	);
	printf "Done. %d measurements written to %s\n", scalar @results, $perl_tsv;
}

# ---------------------------------------------------------------------------
# 7. Axis labelling and curve fitting
# ---------------------------------------------------------------------------
# Both axes carry log10 values and are relabelled with the units they stand
# for, so a reader never has to do the arithmetic.
sub log10 { return log($_[0]) / log(10) }

# 1000 -> 1k, 1e6 -> 1M: the row counts are round numbers, and "1M" is easier
# to land on than "1e+06".
sub count_label {
	my ($e) = @_;
	return sprintf('%gM', 10 ** ($e - 6)) if $e >= 6;
	return sprintf('%gk', 10 ** ($e - 3)) if $e >= 3;
	return sprintf('%g', 10 ** $e);
}

# -6 -> 1us, -3 -> 1ms, 0 -> 1s.  ASCII only: these strings are pasted into a
# generated python file whose encoding is not ours to assume, and "us" costs
# nothing that a mu would buy.
sub time_label {
	my ($e) = @_;
	return sprintf('%gs',  10 ** $e)        if $e >= 0;
	return sprintf('%gms', 10 ** ($e + 3))  if $e >= -3;
	return sprintf('%gus', 10 ** ($e + 6))  if $e >= -6;
	return sprintf('%gns', 10 ** ($e + 9));
}

# Matplotlib::Simple quotes a title for you only when it holds no comma and no
# apostrophe, and half of these titles hold both -- "read_table (csv, mixed)".
# Quote them here rather than relying on that.
sub pyq {
	my ($s) = @_;
	$s =~ s/'/\\'/g;
	return "'$s'";
}

# Ordinary least squares of log10(seconds) on log10(n), over the mean of the
# runs at each size.  Two of them are reported: over every size, and over the
# largest three, because the small-n end of most of these curves is flat -- it
# is measuring per-call overhead, not the algorithm -- and averaging that in
# understates the exponent.  The tail is the one to quote.
sub fit_slope {
	my ($points) = @_;                      # [ [log10 n, log10 seconds], ... ]
	return undef if @$points < 2;
	my ($sx, $sy, $sxx, $sxy) = (0, 0, 0, 0);
	foreach my $p (@$points) {
		$sx  += $p->[0];
		$sy  += $p->[1];
		$sxx += $p->[0] * $p->[0];
		$sxy += $p->[0] * $p->[1];
	}
	my $k    = scalar @$points;
	my $den  = $k * $sxx - $sx * $sx;
	return undef if $den == 0;
	return ($k * $sxy - $sx * $sy) / $den;
}

# ---------------------------------------------------------------------------
# 8. The --plot stage
# ---------------------------------------------------------------------------
# One image per figure -- scaling.vector.svg, scaling.transform.svg,
# scaling.io.svg, scaling.frame.svg -- and one panel per function inside it,
# with Stats::LikeR, Python and R drawn together.
#
# The panels are Matplotlib::Simple's "wide" plot type, which is what these
# measurements want: every one of the seven runs is drawn as a faint line, the
# mean of the runs as a solid one, and one standard deviation either side as a
# translucent ribbon.  Seven lines per language per panel would be a thicket,
# and the mean alone would hide the run-to-run spread -- which on a hybrid
# P-core/E-core CPU is a factor of two and is the first thing you need to see
# before believing any small difference between two languages.
#
# Both axes are logged by hand rather than with set_xscale('log'), for two
# reasons that both come from how "wide" builds its ribbon.  It interpolates
# each run onto 101 points spaced evenly across the group's x range: on a raw
# n axis running to a million, that leaves about one point in the whole first
# decade, so the left-hand third of a log-scaled panel would be a straight
# guess between two samples.  And it takes the mean and standard deviation of
# the raw y, so a ribbon around a 4-microsecond mean reaches below zero and
# falls off the bottom of a log axis entirely.  Feeding it log10 values makes
# the interpolation grid uniform in the space the picture is actually drawn in,
# and makes the ribbon a multiplicative band, which is the right shape for a
# spread that is proportional rather than absolute.  The ticks are then
# relabelled with real units, so nothing about the axis reads as a logarithm.
#
# What to look for: the slope is the exponent.  A flat line is O(1), 45 degrees
# is O(n), steeper than that is worse than linear, and a line that bends upward
# at the right-hand end is an algorithm changing behaviour once the data stops
# fitting in cache.  A constant vertical gap between two languages is per-call
# overhead if it closes as n grows, and per-element cost if it does not.  The
# fitted exponents are printed at the end and written to scaling.slopes.tsv.

# The three files, in the order their lines are labelled.  A file that is not
# there is skipped with a warning rather than being fatal: one language's
# curves are still worth looking at.
my @sources = (
	{ file => $perl_tsv,            group => 'Stats::LikeR', color => 'blue'  },
	{ file => 'python_scaling.tsv', group => 'Python',       color => 'green' },
	{ file => 'r_scaling.tsv',      group => 'R',            color => 'red'   },
);

sub plot_all {
	# require does not run import(), so plt is never installed into main and is
	# called below by its full name.  That is the point: a "use" at the top of
	# the file would make --measure -- the stage that takes hours -- fail at
	# startup on a machine with no plotting stack.
	require Matplotlib::Simple;

	# $t{$figure}{$function}{$group}{$run} = [ [n, seconds], ... ]
	# $panel_order{$figure} keeps first-appearance order, so the panels come out
	# in the order --measure took them rather than in hash order.
	my (%t, %panel_order, %seen_panel, %color);

	foreach my $src (@sources) {
		unless (-f $src->{file}) {
			warn "$src->{file} is not there; skipping $src->{group}\n";
			next;
		}
		$color{ $src->{group} } = $src->{color};

		my $tab = read_table($src->{file}, sep => "\t", 'output.type' => 'hoa');
		my $rows = @{ $tab->{function} };
		for my $i (0 .. $rows - 1) {
			my $figure = $tab->{figure}[$i];
			my $fn     = $tab->{function}[$i];
			push @{ $panel_order{$figure} }, $fn unless $seen_panel{$figure}{$fn}++;
			push @{ $t{$figure}{$fn}{ $src->{group} }{ $tab->{run}[$i] } },
				[ $tab->{n}[$i], $tab->{seconds}[$i] ];
		}
	}

	die "no *_scaling.tsv found; run --measure, scale.py and scale.R first\n"
		unless %t;

	my @slope_rows;

	foreach my $figure (@figure_order) {
		next unless $t{$figure};

		my (@plots, $lo_x, $hi_x, $lo_y, $hi_y);

		foreach my $fn (@{ $panel_order{$figure} }) {
			my %data;

			foreach my $group (map { $_->{group} } @sources) {
				my $group_runs = $t{$figure}{$fn}{$group} or next;

				# One faint line per run: the x of that line is the size ladder
				# and the y is that run's time at each size.
				foreach my $run (sort { $a <=> $b } keys %$group_runs) {
					my @pts = sort { $a->[0] <=> $b->[0] } @{ $group_runs->{$run} };
					next unless @pts;
					my (@x, @y);
					foreach my $p (@pts) {
						# A measurement of exactly zero cannot be logged.  None
						# has ever been seen -- the repeat loop in the three
						# scripts exists so that no reading is at the clock's
						# resolution -- but a zero would take the whole panel
						# with it.
						next unless $p->[1] > 0;
						push @x, log10($p->[0]);
						push @y, log10($p->[1]);
					}
					next unless @x;
					push @{ $data{$group} }, [ \@x, \@y ];

					for my $v (@x) {
						$lo_x = $v if !defined $lo_x || $v < $lo_x;
						$hi_x = $v if !defined $hi_x || $v > $hi_x;
					}
					for my $v (@y) {
						$lo_y = $v if !defined $lo_y || $v < $lo_y;
						$hi_y = $v if !defined $hi_y || $v > $hi_y;
					}
				}

				# the exponent, from the mean over runs at each size
				my (%sum, %count);
				foreach my $run (keys %$group_runs) {
					foreach my $p (@{ $group_runs->{$run} }) {
						$sum{ $p->[0] }   += $p->[1];
						$count{ $p->[0] } += 1;
					}
				}
				my @means = map { [ log10($_), log10($sum{$_} / $count{$_}) ] }
				            grep { $sum{$_} > 0 }
				            sort { $a <=> $b } keys %sum;
				next unless @means;
				my @tail = @means > 3 ? @means[ -3 .. -1 ] : @means;
				my $all  = fit_slope(\@means);
				my $tail = fit_slope(\@tail);
				push @slope_rows, [
					$figure, $fn, $group, scalar @means,
					defined $all  ? sprintf('%.3f', $all)  : 'NA',
					defined $tail ? sprintf('%.3f', $tail) : 'NA',
				];
			}

			next unless %data;
			push @plots, {
				'plot.type' => 'wide',
				data        => \%data,
				color       => { map { $_ => $color{$_} } keys %data },
				title       => pyq($fn),
			};
		}

		next unless @plots;

		# Ticks: whole decades only, on both axes, bracketing the data.
		# Positions and labels go in one set_xticks(ticks, labels) call rather
		# than in set_xticks followed by set_xticklabels, because
		# Matplotlib::Simple emits its axes methods in the order they appear in
		# its own @ax_methods list -- where set_xticklabels comes first -- and
		# labelling before the locator is set is how you get a tick labelled 100
		# on an axis that starts at 1000.  The nudge is not cosmetic:
		# log(1000)/log(10) is 2.9999999999999996, so a bare floor() puts a tick
		# labelled "100" on an axis whose smallest sample is a thousand rows.
		my $eps = 1e-9;
		my @xt = (POSIX::floor($lo_x + $eps) .. POSIX::ceil($hi_x - $eps));
		my @yt = (POSIX::floor($lo_y + $eps) .. POSIX::ceil($hi_y - $eps));
		my $xticks = '[' . join(',', @xt) . '], ['
			. join(',', map { pyq(count_label($_)) } @xt) . ']';
		my $yticks = '[' . join(',', @yt) . '], ['
			. join(',', map { pyq(time_label($_)) } @yt) . ']';

		# Only the outermost panels get axis labels; four columns of "rows" and
		# "seconds per call" is noise, and the ticks already say what the units
		# are.
		my $ncols = @plots >= 8 ? 4 : (@plots >= 3 ? 3 : scalar @plots);
		for my $i (0 .. $#plots) {
			$plots[$i]{set_xticks} = $xticks;
			$plots[$i]{set_yticks} = $yticks;
			$plots[$i]{xlabel} = pyq('rows / elements')
				if $i >= @plots - $ncols;
			$plots[$i]{ylabel} = pyq('seconds per call')
				if $i % $ncols == 0;
			# One legend is enough for a whole figure, and "wide" draws its own
			# in every panel it is left on for.
			$plots[$i]{'show.legend'} = 0 if $i > 0;
		}
		$plots[0]{suptitle} = pyq(
			$figure_title{$figure} . ' -- Stats::LikeR vs Python vs R'
		);

		my $out = "scaling.$figure.svg";
		Matplotlib::Simple::plt(
			'output.file' => $out,
			ncols         => $ncols,
			scalex        => 2.4,
			scaley        => 1.6,
			p             => \@plots,
		);
		printf "%s: %d panels\n", $out, scalar @plots;
	}

	# "tail" is the slope over the largest three sizes and is the number to
	# quote; "all" includes the flat, overhead-dominated left-hand end of every
	# curve and is only there to show how much of the panel that end takes up.
	write_table(
		[ [ 'figure', 'function', 'language', 'sizes', 'slope.all', 'slope.tail' ],
		  @slope_rows ],
		'scaling.slopes.tsv',
		sep         => "\t",
		'row.names' => 0,
	);

	printf "%-9s %-30s %-13s %6s %6s\n",
		'figure', 'function', 'language', 'all', 'tail';
	printf "%-9s %-30s %-13s %6s %6s\n", @$_[0 .. 2], @$_[4, 5] for @slope_rows;
}

# ---------------------------------------------------------------------------
# 9. Run the requested stages
# ---------------------------------------------------------------------------
write_fixtures() if $do_data;

if ($do_measure) {
	# scale.py and scale.R stop when a fixture is missing, because a benchmark
	# script inventing its own corpus is the one thing that would make the three
	# languages incomparable.  Here the generator is in the same file, so the
	# honest thing is to run it -- but only when something is actually absent,
	# so a rerun does not spend a minute rewriting 300,000 rows it already has.
	# Only the sizes this run will actually read are checked, so a truncated
	# SCALE_MAX_N run is not made to sit through the whole corpus; the generator
	# then writes all of it, because that is the corpus scale.py and scale.R
	# expect to find and there is no partial version of it.
	unless ($do_data) {
		my @missing = grep { !-f $_ } map { fixture_files($_) } @io_n;
		if (@missing) {
			printf "%s is missing; writing the fixtures first\n", $missing[0];
			write_fixtures();
		}
	}
	#after the fixtures: writing 40 MB of them is not a measurement, and there
	#is no reason to make it wait for one core
	pin_to_one_cpu();
	measure_all();
}

plot_all() if $do_plot;
