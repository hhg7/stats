#!/usr/bin/env perl
# The Stats::LikeR side of benchmark.R and benchmark.py: the same 32 operations,
# on the same 10,000-row mock frame, seven runs each, written out as
# perl_benchmarks.tsv with the same four columns the other two scripts write.
#
# Run it from the distribution root against the built module:
#
#     perl -Iblib/arch -Iblib/lib benchmark.pl
#
# Only Time::HiRes and POSIX are needed beyond Stats::LikeR itself, both core.
#
# Where the mapping is not one-to-one, the "Perl sub" column says what was
# actually called, and a comment on the entry says why:
#
#   * cor/cov take two vectors here rather than a two-column frame.  Same
#     computation, one number out instead of a 2x2 matrix.
#   * chisq_test is handed a contingency table, so building that table is inside
#     the timed call -- exactly as R times table() inside chisq.test(table(..)).
#   * group_by only splits, so the summarise half (mean per group) is timed with
#     it, which is what dplyr's group_by %>% summarise does.
#   * merge joins on an id column; R merges by row.names and pandas by index,
#     neither of which Stats::LikeR has an equivalent of.
use strict;
use warnings;
use Stats::LikeR;
use Time::HiRes ();
use POSIX ();

# ---------------------------------------------------------------------------
# 1. Setup Mock Data
# ---------------------------------------------------------------------------
srand(42);
my $n = 10_000;

my $df = {
	x      => rnorm(n => $n, mean => 0, sd => 1),
	y      => rnorm(n => $n, mean => 5, sd => 2),
	cat1   => [ map { (qw(A B C))[ int rand 3 ] } 1 .. $n ],
	cat2   => [ map { (qw(X Y))[ int rand 2 ] } 1 .. $n ],
	binary => rbinom(n => $n, prob => 0.5, size => 1),
};

# the same frame with a hole in x, for dropna and fillna
my $df_missing = { %$df, x => [ @{ $df->{x} } ] };
$df_missing->{x}[$_] = undef for 10 .. 20;

# merge needs a key column: R joins these two frames by row.names and pandas by
# index, and a Stats::LikeR frame has neither.
my $df_id = { %$df, id => [ 1 .. $n ] };

my $y_true   = rbinom(n => $n, prob => 0.5, size => 1);
my $y_scores = runif(n => $n, min => 0, max => 1);

# ---------------------------------------------------------------------------
# 2. Define Function Mappings (Stats::LikeR function -> the call being timed)
# ---------------------------------------------------------------------------
# An array rather than a hash, so the rows come out in this order.
my @benchmarks = (
	# Basic Stats
	[ 'mean',     'Stats::LikeR::mean',     sub { mean($df->{x}) } ],
	[ 'median',   'Stats::LikeR::median',   sub { median($df->{x}) } ],
	[ 'sd',       'Stats::LikeR::sd',       sub { sd($df->{x}) } ],
	[ 'var',      'Stats::LikeR::var',      sub { var($df->{x}) } ],
	[ 'min',      'Stats::LikeR::min',      sub { min($df->{x}) } ],
	[ 'max',      'Stats::LikeR::max',      sub { max($df->{x}) } ],
	[ 'quantile', 'Stats::LikeR::quantile',
		sub { quantile(x => $df->{x}, probs => [ 0.25, 0.5, 0.75 ]) } ],
	[ 'cor',      'Stats::LikeR::cor',      sub { cor($df->{x}, $df->{y}) } ],
	[ 'cov',      'Stats::LikeR::cov',      sub { cov($df->{x}, $df->{y}) } ],

	# Distributions & Tests
	[ 'rnorm', 'Stats::LikeR::rnorm', sub { rnorm(n => 1000, mean => 0, sd => 1) } ],
	[ 'runif', 'Stats::LikeR::runif', sub { runif(n => 1000, min => 0, max => 1) } ],
	[ 't_test', 'filter + vals + t_test', sub {
		t_test(
			vals(filter($df, col('cat2') eq 'X'), 'x'),
			vals(filter($df, col('cat2') eq 'Y'), 'x'),
		);
	} ],
	[ 'wilcox_test', 'filter + vals + wilcox_test', sub {
		wilcox_test(
			vals(filter($df, col('cat2') eq 'X'), 'x'),
			vals(filter($df, col('cat2') eq 'Y'), 'x'),
		);
	} ],
	[ 'chisq_test', 'pivot_table + chisq_test', sub {
		# the cross-tabulation R gets from table(): counts of cat1 x cat2
		my $xt   = pivot_table($df, index => 'cat1', columns => 'cat2',
		                       values => 'x', aggfunc => 'count');
		my @cols = grep { $_ ne 'cat1' } sort keys %$xt;
		chisq_test([ map { my $r = $_; [ map { $xt->{$_}[$r] || 0 } @cols ] }
		             0 .. $#{ $xt->{cat1} } ]);
	} ],
	[ 'shapiro_test', 'Stats::LikeR::shapiro_test',
		sub { shapiro_test([ @{ $df->{x} }[ 0 .. 4999 ] ]) } ],   # as in R: n <= 5000
	[ 'binom_test', 'Stats::LikeR::binom_test',
		sub { binom_test(500, 1000, p => 0.5) } ],

	# Data Manipulation
	[ 'filter',      'filter + col',            sub { filter($df, col('x') > 0) } ],
	[ 'select_cols', 'Stats::LikeR::select_cols', sub { select_cols($df, [ 'x', 'cat1' ]) } ],
	[ 'drop_cols',   'Stats::LikeR::drop_cols',   sub { drop_cols($df, 'y') } ],
	[ 'rename_cols', 'Stats::LikeR::rename_cols', sub { rename_cols($df, cat1 => 'Category_1') } ],
	[ 'dropna',      'Stats::LikeR::dropna',
		sub { dropna($df_missing, cols => [ sort keys %$df_missing ]) } ],
	[ 'fillna',      'Stats::LikeR::fillna',  sub { fillna($df_missing, value => 0) } ],
	[ 'drop_duplicates', 'Stats::LikeR::drop_duplicates', sub { drop_duplicates($df) } ],
	[ 'group_by', 'group_by + mean', sub {
		my $groups = group_by($df, 'x', 'cat1');
		return { map { $_ => mean($groups->{$_}) } keys %$groups };
	} ],
	[ 'concat', 'Stats::LikeR::concat', sub { concat($df, $df) } ],
	[ 'merge',  'Stats::LikeR::merge',
		sub { merge($df_id, $df_id, how => 'inner', on => 'id') } ],
	[ 'value_counts', 'Stats::LikeR::value_counts', sub { value_counts($df, 'cat1') } ],
	[ 'pivot_table',  'Stats::LikeR::pivot_table', sub {
		pivot_table($df, index => 'cat1', columns => 'cat2',
		            values => 'x', aggfunc => 'mean');
	} ],

	# Modeling & Metrics
	[ 'lm',  'Stats::LikeR::lm',  sub { lm(formula => 'y ~ x + cat1', data => $df) } ],
	[ 'glm', 'Stats::LikeR::glm',
		sub { glm(data => $df, formula => 'binary ~ x + y', family => 'binomial') } ],
	[ 'auc', 'Stats::LikeR::auc', sub { auc($y_scores, $y_true) } ],
	[ 'prcomp', 'Stats::LikeR::prcomp', sub {
		# prcomp wants a frame of its own, and building it is part of the call,
		# the way df[, c("x","y")] is part of R's.
		prcomp([ map { [ $df->{x}[$_], $df->{y}[$_] ] } 0 .. $n - 1 ]);
	} ],
);

# ---------------------------------------------------------------------------
# 3. Execution Engine
# ---------------------------------------------------------------------------
# Each run happens in a forked child that does nothing else, and reports back
# the seconds it took and the resident memory it grew by.
#
# The fork is what makes the memory figure mean anything.  Perl hands freed
# memory back to its own arenas rather than to the OS, so measuring in one
# process gives the first run the whole bill and every later run zero.  A child
# starts from the parent's untouched heap every time, so all seven runs of an
# operation are measured against the same baseline and agree with each other.
#
# What the number contains: the memory the result holds, plus whatever the
# operation allocated and freed along the way (the arena keeps it), plus the
# copy-on-write copies of any input the call writes to -- reading a column as a
# number caches the conversion in the SV, which dirties its page.  It is
# page-granular (4 KiB), so small operations round up.  Like the gc() figure in
# benchmark.R, it is the right order of magnitude rather than an allocation
# ledger of the kind Python's tracemalloc keeps.
#
# Forking, timing and reporting cost some memory of their own, the same amount
# every time.  That floor is measured below by running an empty subroutine and
# subtracted from every reading, which is what benchmark.R's base_mem does.
sub rss_bytes {
	open my $fh, '<', '/proc/self/status' or return undef;
	while (my $line = <$fh>) {
		return $1 * 1024 if $line =~ /^VmRSS:\s+(\d+)\s+kB/;
	}
	return undef;
}

# One measured run.  Returns (seconds, bytes, error), any of which may be undef
# if the child could not report.
sub measure {
	my ($code) = @_;

	my $pid = open my $from_child, '-|';
	die "fork failed: $!\n" unless defined $pid;

	if (!$pid) {                                  # the child
		$| = 1;                                   # _exit does not flush
		rss_bytes();                              # first read costs ~300 KiB of
		my $before = rss_bytes();                 # its own; keep it out of the
		                                          # baseline
		my $start = Time::HiRes::time();
		my $ok    = eval { $code->(); 1 };
		my $end   = Time::HiRes::time();
		my $err   = $ok ? '' : $@;

		my $after = rss_bytes();
		$err =~ s/\s+/ /g;
		print join("\t", $end - $start,
		                 (defined $after && defined $before) ? $after - $before : '',
		                 $err), "\n";
		POSIX::_exit(0);
	}

	my $line = <$from_child>;
	close $from_child;
	return (undef, undef, "the child died without reporting") unless defined $line;

	chomp $line;
	my ($secs, $bytes, $err) = split /\t/, $line, 3;
	# the memory field is empty when /proc/self/status is not there to read
	$bytes = undef unless defined $bytes && length $bytes;
	return ($secs, $bytes, (defined $err && length $err) ? $err : undef);
}

my $runs = 7;
my @results;

# what a measured run of nothing at all costs, in bytes
my $floor = 0;
for (1 .. 5) {
	my (undef, $bytes) = measure(sub { });
	$floor = $bytes if defined $bytes && (!$floor || $bytes < $floor);
}
printf "Measurement floor: %d bytes, subtracted from every RAM figure below.\n",
	$floor;

print "Running Perl benchmarks...\n";
for my $b (@benchmarks) {
	my ($liker, $perl_sub, $code) = @$b;
	for my $i (1 .. $runs) {
		my ($secs, $bytes, $err) = measure($code);
		print "Error in $liker : $err\n" if defined $err;
		if (defined $bytes) {
			$bytes -= $floor;
			$bytes = 0 if $bytes < 0;
		}
		push @results, [ $liker, $perl_sub,
		                 defined $secs  ? $secs  : 'NA',
		                 defined $bytes ? $bytes : 'NA' ];
	}
}

# ---------------------------------------------------------------------------
# 4. Output
# ---------------------------------------------------------------------------
my $output_file = 'perl_benchmarks.tsv';
write_table(
	[ [ 'Stats::LikeR function', 'Perl sub', 'time', 'RAM' ], @results ],
	$output_file,
	'row.names' => 0,
);

printf "Done. Perl results written to %s\n", $output_file;
