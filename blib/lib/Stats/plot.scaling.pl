#!/usr/bin/env perl
# Draw the scaling curves that scale.pl, scale.py and scale.R measured.
#
#     perl scale.data.pl                          # once, writes the fixtures
#     perl -Iblib/arch -Iblib/lib scale.pl        # -> perl_scaling.tsv
#     python3 scale.py                            # -> python_scaling.tsv
#     Rscript scale.R                             # -> r_scaling.tsv
#     perl plot.scaling.pl                        # -> scaling.*.svg
#
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
require 5.010;
use strict;
use warnings FATAL => 'all';
use POSIX ();
use Stats::LikeR;
use Matplotlib::Simple 'plt';

# The three files, in the order their panels are laid out and their lines are
# labelled.  A file that is not there is skipped with a warning rather than
# being fatal: one language's curves are still worth looking at.
my @sources = (
	{ file => 'perl_scaling.tsv',   group => 'Stats::LikeR', color => 'blue'  },
	{ file => 'python_scaling.tsv', group => 'Python',       color => 'green' },
	{ file => 'r_scaling.tsv',      group => 'R',            color => 'red'   },
);

my %figure_title = (
	vector    => 'Reductions over one numeric vector',
	transform => 'Transforms that return as much as they are given',
	io        => 'read_table and write_table',
	frame     => 'Whole-frame operations',
);
my @figure_order = qw(vector transform io frame);

# ---------------------------------------------------------------------------
# 1. Read
# ---------------------------------------------------------------------------
# $t{$figure}{$function}{$group}{$run} = [ [n, seconds], ... ]
# @panel_order{$figure} keeps first-appearance order, so the panels come out in
# the order scale.pl measured them rather than in hash order.
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

die "no *_scaling.tsv found; run scale.pl, scale.py and scale.R first\n"
	unless %t;

# ---------------------------------------------------------------------------
# 2. Axis labelling
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

# ---------------------------------------------------------------------------
# 3. Fitted exponents
# ---------------------------------------------------------------------------
# Ordinary least squares of log10(seconds) on log10(n), over the mean of the
# runs at each size.  Two of them: over every size, and over the largest three,
# because the small-n end of most of these curves is flat -- it is measuring
# per-call overhead, not the algorithm -- and averaging that in understates the
# exponent.  The tail is the one to quote.
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

my @slope_rows;

# ---------------------------------------------------------------------------
# 4. Draw
# ---------------------------------------------------------------------------
foreach my $figure (@figure_order) {
	next unless $t{$figure};

	my (@plots, $lo_x, $hi_x, $lo_y, $hi_y);

	foreach my $fn (@{ $panel_order{$figure} }) {
		my %data;

		foreach my $group (map { $_->{group} } @sources) {
			my $runs = $t{$figure}{$fn}{$group} or next;

			# One faint line per run: the x of that line is the size ladder and
			# the y is that run's time at each size.
			foreach my $run (sort { $a <=> $b } keys %$runs) {
				my @pts = sort { $a->[0] <=> $b->[0] } @{ $runs->{$run} };
				next unless @pts;
				my (@x, @y);
				foreach my $p (@pts) {
					# A measurement of exactly zero cannot be logged.  None has
					# ever been seen -- the repeat loop in the three scripts
					# exists so that no reading is at the clock's resolution --
					# but a zero would take the whole panel with it.
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
			foreach my $run (keys %$runs) {
				foreach my $p (@{ $runs->{$run} }) {
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

	# Ticks: whole decades only, on both axes, bracketing the data.  Positions
	# and labels go in one set_xticks(ticks, labels) call rather than in
	# set_xticks followed by set_xticklabels, because Matplotlib::Simple emits
	# its axes methods in the order they appear in its own @ax_methods list --
	# where set_xticklabels comes first -- and labelling before the locator is
	# set is how you get a tick labelled 100 on an axis that starts at 1000.
	# The nudge is not cosmetic: log(1000)/log(10) is 2.9999999999999996, so a
	# bare floor() puts a tick labelled "100" on an axis whose smallest sample
	# is a thousand rows.
	my $eps = 1e-9;
	my @xt = (POSIX::floor($lo_x + $eps) .. POSIX::ceil($hi_x - $eps));
	my @yt = (POSIX::floor($lo_y + $eps) .. POSIX::ceil($hi_y - $eps));
	my $xticks = '[' . join(',', @xt) . '], ['
		. join(',', map { pyq(count_label($_)) } @xt) . ']';
	my $yticks = '[' . join(',', @yt) . '], ['
		. join(',', map { pyq(time_label($_)) } @yt) . ']';

	# Only the outermost panels get axis labels; four columns of "rows" and
	# "seconds per call" is noise, and the ticks already say what the units are.
	my $ncols = @plots >= 8 ? 4 : (@plots >= 3 ? 3 : scalar @plots);
	for my $i (0 .. $#plots) {
		$plots[$i]{set_xticks} = $xticks;
		$plots[$i]{set_yticks} = $yticks;
		$plots[$i]{xlabel} = pyq('rows / elements')
			if $i >= @plots - $ncols;
		$plots[$i]{ylabel} = pyq('seconds per call')
			if $i % $ncols == 0;
		# One legend is enough for a whole figure, and "wide" draws its own in
		# every panel it is left on for.
		$plots[$i]{'show.legend'} = 0 if $i > 0;
	}
	$plots[0]{suptitle} = pyq(
		$figure_title{$figure} . ' -- Stats::LikeR vs Python vs R'
	);

	my $out = "scaling.$figure.svg";
	plt(
		'output.file' => $out,
		ncols         => $ncols,
		scalex        => 2.4,
		scaley        => 1.6,
		p             => \@plots,
	);
	printf "%s: %d panels\n", $out, scalar @plots;
}

# ---------------------------------------------------------------------------
# 5. The exponents, as a table
# ---------------------------------------------------------------------------
# "tail" is the slope over the largest three sizes and is the number to quote;
# "all" includes the flat, overhead-dominated left-hand end of every curve and
# is only there to show how much of the panel that end takes up.
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
