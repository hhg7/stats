#!/usr/bin/env perl
# How Stats::LikeR's XS functions scale with the size of their input.
#
# Where benchmark.pl asks "how long does each function take on one 10,000-row
# frame", this asks "what shape is that curve".  The same call is timed at a
# ladder of sizes an order of magnitude apart, seven times at each size, and
# the result is written out for plot.scaling.pl to draw against the same
# measurements from scale.py and scale.R.  Read on a log-log axis, the slope of
# the line is the exponent: flat is O(1), 45 degrees is O(n), steeper is worse.
#
#     perl scale.data.pl                          # once, writes the fixtures
#     perl -Iblib/arch -Iblib/lib scale.pl        # -> perl_scaling.tsv
#     python3 scale.py                            # -> python_scaling.tsv
#     Rscript scale.R                             # -> r_scaling.tsv
#     perl plot.scaling.pl                        # -> scaling.*.png
#
# Environment:
#
#     SCALE_DIR    where scale.data.pl put the fixtures (/tmp/likeR.scaling)
#     SCALE_RUNS   runs per (function, size); default 7
#     SCALE_CAP    seconds; once one run of a function takes longer than this,
#                  that function is not tried at any larger size.  Default 4.
#     SCALE_MAX_N  hard ceiling on the row count, for a quick partial run
#     SCALE_TARGET seconds a single measurement should span; a call faster than
#                  this is repeated until it does.  Default 0.002.
#
# On a hybrid CPU -- Intel's P-core/E-core parts, and the big.LITTLE ARM
# designs -- run all three scripts under "taskset -c 0" (or the platform's
# equivalent).  Measured here, a forked child landing on an E-core reads 1.5 to
# 2 times slower than the identical loop on a P-core, which is a wider band
# than most of the differences these plots are drawn to show.  It is worst in
# this script, because the fork per run gives the scheduler a fresh chance to
# place the work somewhere else every time.
#
# Only Time::HiRes and POSIX are needed beyond Stats::LikeR itself, both core.
#
# What is being compared, and what is not:
#
#   * Every panel asks the three languages for the *same result* by each one's
#     idiomatic route, which is the rule benchmark.pl already follows.  Where
#     that is impossible the panel is simply missing that language rather than
#     being filled with something else: R has no skew(), no kurtosis() and no
#     row-record frame to build, so it does not appear in those four panels.
#   * The read_table and write_table panels all read the byte-identical files
#     scale.data.pl wrote, so no reader is handed a quoting or type-guessing
#     job the others were spared.
#   * Absolute heights are worth less than slopes here.  A panel where one line
#     sits below another by a constant factor over four decades is telling you
#     about per-element cost; a panel where the lines converge or cross is
#     telling you about fixed overhead, which is the more interesting finding
#     and the one a single-size benchmark cannot show.
require 5.010;
use strict;
use warnings FATAL => 'all';
use Stats::LikeR;
use File::Spec ();
use POSIX ();
use Time::HiRes ();

my $dir  = $ENV{SCALE_DIR}  || '/tmp/likeR.scaling';
my $runs = $ENV{SCALE_RUNS} || 7;
my $cap  = defined $ENV{SCALE_CAP} ? $ENV{SCALE_CAP} : 4;
my $max_n = $ENV{SCALE_MAX_N} || 0;
my $target = defined $ENV{SCALE_TARGET} ? $ENV{SCALE_TARGET} : 0.002;
my $max_reps = 10_000;

# ---------------------------------------------------------------------------
# 1. Size ladders
# ---------------------------------------------------------------------------
# Half-decade steps, so seven points span three decades and the log-log slope
# is estimated from evenly spaced x.  scale.py and scale.R carry the same three
# lists; scale.data.pl's row counts are @io_n.  Change one, change all four.
my @vec_n   = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000);
my @io_n    = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000);
my @frame_n = (1_000, 3_000, 10_000, 30_000, 100_000, 300_000);

if ($max_n) {
	@vec_n   = grep { $_ <= $max_n } @vec_n;
	@io_n    = grep { $_ <= $max_n } @io_n;
	@frame_n = grep { $_ <= $max_n } @frame_n;
}

# ---------------------------------------------------------------------------
# 2. Input builders
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
		my %f = (
			num_csv => File::Spec->catfile($dir, "num.$n.csv"),
			mix_csv => File::Spec->catfile($dir, "mix.$n.csv"),
			mix_tsv => File::Spec->catfile($dir, "mix.$n.tsv"),
			out     => File::Spec->catfile($dir, "out.perl.$$.tmp"),
		);
		for my $k (qw(num_csv mix_csv mix_tsv)) {
			die "missing fixture \"$f{$k}\"; run \"perl scale.data.pl\" first\n"
				unless -f $f{$k};
		}
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
# 3. The benchmarks
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

# ---------------------------------------------------------------------------
# 4. Execution engine
# ---------------------------------------------------------------------------
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

my @results;
my %too_slow; # name => 1 once it exceeds $cap

# Grouped by figure, then by size, so each input is built once and every
# benchmark that wants it is measured before it is thrown away.  The 1,000,000
# element frames are large enough that holding two ladders' worth at once is
# worth avoiding.
my %by_figure;
push @{ $by_figure{ $_->{figure} } }, $_ for @benchmarks;

foreach my $figure (qw(vector transform io frame)) {
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

# Anything the cap stopped is worth saying out loud: a curve that ends early is
# not the same as a curve that was never measured, and the plot cannot tell you
# which one you are looking at.
if (%too_slow) {
	printf "Stopped early (a run exceeded %g s, or it failed): %s\n",
		$cap, join(', ', sort keys %too_slow);
}

# clean up the write_table target
my $out = File::Spec->catfile($dir, "out.perl.$$.tmp");
unlink $out if -f $out;

# ---------------------------------------------------------------------------
# 5. Output
# ---------------------------------------------------------------------------
my $output_file = 'perl_scaling.tsv';
write_table(
	[ [ 'figure', 'function', 'call', 'n', 'run', 'seconds' ], @results ],
	$output_file,
	sep         => "\t",
	'row.names' => 0,
);
printf "Done. %d measurements written to %s\n", scalar @results, $output_file;
