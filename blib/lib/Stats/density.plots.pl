#!/usr/bin/env perl
# Draws the figures that illustrate density() in README.md (and therefore in
# read.me.pod and lib/Stats/LikeR.pm).  Author-only: it is not shipped, and it
# needs Matplotlib::Simple, python3 and matplotlib.
#
#   perl -Iblib/lib -Iblib/arch density.plots.pl
#
# Every figure is written to img/, at the natural size the POD and the README
# embed them at.  The data are R's own: faithful$eruptions (n = 272, the
# example in ?density) and precip (n = 70, the example in ?bw.nrd), copied from
# t/density.R.scipy.t so the pictures and the tests describe the same numbers.
#
# Colours are the validated categorical/sequential palette: identity (kernel,
# weighting) gets categorical slots in their documented order, and anything
# ordered -- the adjust ladder, the cut ladder -- gets the one-hue sequential
# ramp instead, dark meaning more.
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Path 'make_path';
use Stats::LikeR qw(density hist bw_nrd0 bw_nrd bw_ucv bw_bcv bw_sj);
use Matplotlib::Simple 'plt';

my $DIR = 'img';
make_path($DIR) unless -d $DIR;

# --- palette ----------------------------------------------------------------
my $BLUE   = '#2a78d6';    # categorical slot 1
my $ORANGE = '#eb6834';    # categorical slot 2
my $GREY   = '#8c8c88';    # recessive reference marks
my $INK    = '#3d3d3a';    # rug ticks
# the blue sequential ramp, light -> dark, for ordered quantities
my @RAMP = ('#86b6ef', '#3987e5', '#1c5cab', '#0d366b');

# --- the data ---------------------------------------------------------------
my @eruptions = (
	3.6, 1.8, 3.333, 2.283, 4.533, 2.883, 4.7, 3.6, 1.95, 4.35, 1.833, 3.917,
	4.2, 1.75, 4.7, 2.167, 1.75, 4.8, 1.6, 4.25, 1.8, 1.75, 3.45, 3.067,
	4.533, 3.6, 1.967, 4.083, 3.85, 4.433, 4.3, 4.467, 3.367, 4.033, 3.833, 2.017,
	1.867, 4.833, 1.833, 4.783, 4.35, 1.883, 4.567, 1.75, 4.533, 3.317, 3.833, 2.1,
	4.633, 2, 4.8, 4.716, 1.833, 4.833, 1.733, 4.883, 3.717, 1.667, 4.567, 4.317,
	2.233, 4.5, 1.75, 4.8, 1.817, 4.4, 4.167, 4.7, 2.067, 4.7, 4.033, 1.967,
	4.5, 4, 1.983, 5.067, 2.017, 4.567, 3.883, 3.6, 4.133, 4.333, 4.1, 2.633,
	4.067, 4.933, 3.95, 4.517, 2.167, 4, 2.2, 4.333, 1.867, 4.817, 1.833, 4.3,
	4.667, 3.75, 1.867, 4.9, 2.483, 4.367, 2.1, 4.5, 4.05, 1.867, 4.7, 1.783,
	4.85, 3.683, 4.733, 2.3, 4.9, 4.417, 1.7, 4.633, 2.317, 4.6, 1.817, 4.417,
	2.617, 4.067, 4.25, 1.967, 4.6, 3.767, 1.917, 4.5, 2.267, 4.65, 1.867, 4.167,
	2.8, 4.333, 1.833, 4.383, 1.883, 4.933, 2.033, 3.733, 4.233, 2.233, 4.533, 4.817,
	4.333, 1.983, 4.633, 2.017, 5.1, 1.8, 5.033, 4, 2.4, 4.6, 3.567, 4,
	4.5, 4.083, 1.8, 3.967, 2.2, 4.15, 2, 3.833, 3.5, 4.583, 2.367, 5,
	1.933, 4.617, 1.917, 2.083, 4.583, 3.333, 4.167, 4.333, 4.5, 2.417, 4, 4.167,
	1.883, 4.583, 4.25, 3.767, 2.033, 4.433, 4.083, 1.833, 4.417, 2.183, 4.8, 1.833,
	4.8, 4.1, 3.966, 4.233, 3.5, 4.366, 2.25, 4.667, 2.1, 4.35, 4.133, 1.867,
	4.6, 1.783, 4.367, 3.85, 1.933, 4.5, 2.383, 4.7, 1.867, 3.833, 3.417, 4.233,
	2.4, 4.8, 2, 4.15, 1.867, 4.267, 1.75, 4.483, 4, 4.117, 4.083, 4.267,
	3.917, 4.55, 4.083, 2.417, 4.183, 2.217, 4.45, 1.883, 1.85, 4.283, 3.95, 2.333,
	4.15, 2.35, 4.933, 2.9, 4.583, 3.833, 2.083, 4.367, 2.133, 4.35, 2.2, 4.45,
	3.567, 4.5, 4.15, 3.817, 3.917, 4.45, 2, 4.283, 4.767, 4.533, 1.85, 4.25,
	1.983, 2.25, 4.75, 4.117, 2.15, 4.417, 1.817, 4.467
);

my @precip = (
	67, 54.7, 7, 48.5, 14, 17.2, 20.7, 13, 43.4, 40.2, 38.9, 54.5,
	59.8, 48.3, 22.9, 11.5, 34.4, 35.1, 38.7, 30.8, 30.6, 43.1, 56.8, 40.8,
	41.8, 42.5, 31, 31.7, 30.2, 25.9, 49.2, 37, 35.9, 15, 30.2, 7.2,
	36.2, 45.5, 7.8, 33.4, 36.1, 40.2, 42.7, 42.5, 16.2, 39, 35, 37,
	31.4, 37.6, 39.9, 36.2, 42.8, 46.4, 24.7, 49.1, 46, 35.9, 7.8, 48.2,
	15.2, 32.5, 44.7, 42.6, 38.8, 17.4, 40.8, 29.1, 14.6, 59.2
);

# --- helpers ----------------------------------------------------------------

# The staircase outline of a histogram on the density scale, from this module's
# own hist(), so the bars and the curve are directly comparable.  The outline is
# closed down to zero at both ends, the way the bars themselves would be.
sub staircase {
	my ($x, $breaks) = @_;
	my $h = hist($x, breaks => $breaks);
	my @sx = ( $h->{breaks}[0] );
	my @sy = ( 0 );
	for my $i (0 .. $#{ $h->{counts} }) {
		push @sx, $h->{breaks}[$i], $h->{breaks}[ $i + 1 ];
		push @sy, $h->{density}[$i], $h->{density}[$i];
	}
	push @sx, $h->{breaks}[-1];
	push @sy, 0;
	return (\@sx, \@sy);
}

# One vertical tick per observation, drawn as its own two-point line.
sub rug {
	my ($x, $height, $width) = @_;
	$width = 0.6 unless defined $width;
	return {
		'plot.type'   => 'plot',
		data          => [ map { [ [ $_, $_ ], [ 0, $height ] ] } @$x ],
		'show.legend' => 0,
		'set.options' => 'color = "' . $INK . '", linewidth = ' . $width . ', alpha = 0.55',
	};
}

# --- 1. what density() does -------------------------------------------------
# Left: the definition -- one kernel per observation, averaged.  Right: the same
# estimate over a real sample, against the histogram it smooths.
sub fig_what {
	my @small = (1.4, 2, 2.2, 3.6, 4.4, 4.6, 5.4);
	my %grid  = (bw => 0.4, from => 0, to => 7, n => 512);

	my $sum = density(\@small, %grid);
	my @bumps;
	for my $pt (@small) {
		my $k = density([$pt], %grid);
		push @bumps, [ $k->{x}, [ map { $_ / scalar @small } @{ $k->{y} } ] ];
	}
	my $first = shift @bumps;

	my $d = density(\@eruptions);
	my ($sx, $sy) = staircase(\@eruptions, 12);

	plt(
		'output.file' => "$DIR/density.what.png",
		ncol          => 2,
		p             => [
			[
				{
					'plot.type' => 'plot',
					data        => {
						'sum of the kernels = density()' => [ $sum->{x}, $sum->{y} ],
						'one kernel / n'                 => $first,
					},
					'key.order'   => [ 'one kernel / n', 'sum of the kernels = density()' ],
					'set.options' => {
						'sum of the kernels = density()' => 'color = "' . $BLUE . '", linewidth = 2.2',
						'one kernel / n'                 => 'color = "' . $GREY . '", linewidth = 1',
					},
					title         => 'A kernel at every observation',
					xlabel        => 'x',
					ylabel        => 'density',
					legend        => 'loc = "upper left", frameon = False, fontsize = 9',
					set_figwidth  => 11,
					set_figheight => 4.2,
					set_dpi       => 100,
				},
				{
					'plot.type'   => 'plot',
					data          => \@bumps,
					'show.legend' => 0,
					'set.options' => 'color = "' . $GREY . '", linewidth = 1',
				},
				rug(\@small, 0.02, 1.4),
			],
			[
				{
					'plot.type' => 'plot',
					data        => {
						'density(eruptions)' => [ $d->{x}, $d->{y} ],
						'hist(eruptions)'    => [ $sx, $sy ],
					},
					'key.order'   => [ 'hist(eruptions)', 'density(eruptions)' ],
					'set.options' => {
						'density(eruptions)' => 'color = "' . $BLUE . '", linewidth = 2.2',
						'hist(eruptions)'    => 'color = "' . $GREY . '", linewidth = 1',
					},
					title  => '"The same estimate, over R\'s faithful$eruptions"',
					xlabel => 'eruption length (minutes)',
					ylabel => 'density',
					legend => 'loc = "upper left", frameon = False, fontsize = 9',
				},
				rug(\@eruptions, 0.02),
			],
		],
	);
	return;
}

# --- 2. bandwidth ------------------------------------------------------------
# adjust is an ordered quantity, so it gets the sequential ramp: darker = more
# smoothing.
sub fig_bandwidth {
	my @adjust = (0.15, 0.4, 1, 2.5);
	my (%data, %opts, @order);
	for my $i (0 .. $#adjust) {
		my $d     = density(\@eruptions, adjust => $adjust[$i]);
		my $label = sprintf 'adjust = %-4s h = %.3f', $adjust[$i], $d->{bw};
		push @order, $label;
		$data{$label} = [ $d->{x}, $d->{y} ];
		$opts{$label} = 'color = "' . $RAMP[$i] . '", linewidth = 2';
	}

	plt(
		'output.file' => "$DIR/density.bandwidth.png",
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => \%data,
					'key.order'   => \@order,
					'set.options' => \%opts,
					title         => '"bw sets the width of the kernel; adjust scales it"',
					xlabel        => 'eruption length (minutes)',
					ylabel        => 'density',
					legend        => 'loc = "upper left", frameon = False, fontsize = 9',
					set_figwidth  => 8,
					set_figheight => 4.5,
					set_dpi       => 100,
				},
				rug(\@eruptions, 0.025),
			],
		],
	);
	return;
}

# --- 3. the five bandwidth rules --------------------------------------------
# Six panels rather than six lines on one axes: the rules are told apart by the
# panel they are in, not by a colour the reader has to match to a legend.
sub fig_rules {
	my @rule = (
		[ 'nrd0'    => bw_nrd0(\@eruptions),                  q{bw => 'nrd0' (the default)} ],
		[ 'nrd'     => bw_nrd(\@eruptions),                   q{bw => 'nrd'} ],
		[ 'ucv'     => bw_ucv(\@eruptions),                   q{bw => 'ucv'} ],
		[ 'bcv'     => bw_bcv(\@eruptions),                   q{bw => 'bcv'} ],
		[ 'SJ-ste'  => bw_sj(\@eruptions),                    q{bw => 'SJ'} ],
		[ 'SJ-dpi'  => bw_sj(x => \@eruptions, method => 'dpi'), q{bw => 'SJ-dpi'} ],
	);
	my ($sx, $sy) = staircase(\@eruptions, 12);

	my @panels;
	for my $i (0 .. $#rule) {
		my ($name, $h, $call) = @{ $rule[$i] };
		my $d = density(\@eruptions, bw => $h);
		my %first = (
			title  => sprintf('"%s   h = %.3f"', $call, $h),
			xlabel => $i >= 3 ? 'eruption length (minutes)' : undef,
			ylabel => $i % 3 == 0 ? 'density' : undef,
		);
		delete $first{$_} for grep { !defined $first{$_} } keys %first;
		if ($i == 0) {
			$first{set_figwidth}  = 12;
			$first{set_figheight} = 6;
			$first{set_dpi}       = 100;
		}
		push @panels, [
			{
				'plot.type'   => 'plot',
				data          => [ [ $d->{x}, $d->{y} ] ],
				'show.legend' => 0,
				'set.options' => 'color = "' . $BLUE . '", linewidth = 2',
				%first,
			},
			{
				'plot.type'   => 'plot',
				data          => [ [ $sx, $sy ] ],
				'show.legend' => 0,
				'set.options' => 'color = "' . $GREY . '", linewidth = 0.9',
			},
		];
	}

	plt(
		'output.file' => "$DIR/density.bw.rules.png",
		p             => \@panels,
		ncol          => 3,
		sharex        => 'True',
		sharey        => 'True',
	);
	return;
}

# --- 4. the seven kernels ----------------------------------------------------
# Seven panels of kernel shape -- density() of a single observation at 0 with
# bw = 1 is the kernel itself -- and an eighth panel showing what changing the
# kernel does to a real estimate, which is almost nothing.
sub fig_kernels {
	my @kernel = qw(gaussian epanechnikov rectangular triangular biweight cosine optcosine);

	my @panels;
	for my $i (0 .. $#kernel) {
		my $k = $kernel[$i];
		my $d = density([0], bw => 1, kernel => $k, from => -3.5, to => 3.5, n => 512);
		my %first = (
			title  => sprintf('"%s   R(K) = %.4f"', $k, density(kernel => $k, give_rkern => 1)),
			ylabel => $i % 4 == 0 ? 'K(t)' : undef,
			xlabel => $i >= 4 ? 't (standard deviations)' : undef,
		);
		delete $first{$_} for grep { !defined $first{$_} } keys %first;
		if ($i == 0) {
			$first{set_figwidth}  = 13;
			$first{set_figheight} = 6;
			$first{set_dpi}       = 100;
		}
		push @panels, {
			'plot.type'   => 'plot',
			data          => [ [ $d->{x}, $d->{y} ] ],
			'show.legend' => 0,
			'set.options' => 'color = "' . $BLUE . '", linewidth = 2',
			%first,
		};
	}

	# the eighth panel: all seven kernels over the same sample, at one bandwidth
	my $h = bw_nrd0(\@eruptions);
	my $g = density(\@eruptions, bw => $h, kernel => 'gaussian');
	my @others = map { my $d = density(\@eruptions, bw => $h, kernel => $_); [ $d->{x}, $d->{y} ] }
	             grep { $_ ne 'gaussian' } @kernel;
	my $first_other = shift @others;

	push @panels, [
		{
			'plot.type' => 'plot',
			data        => {
				'gaussian'       => [ $g->{x}, $g->{y} ],
				'the other six'  => $first_other,
			},
			'key.order'   => [ 'the other six', 'gaussian' ],
			'set.options' => {
				'gaussian'      => 'color = "' . $BLUE . '", linewidth = 2',
				'the other six' => 'color = "' . $GREY . '", linewidth = 1',
			},
			title  => q{"all seven on one sample", fontsize = 10},
			xlabel => 'eruption length (minutes)',
			legend => 'loc = "upper left", frameon = False, fontsize = 8',
		},
		{
			'plot.type'   => 'plot',
			data          => \@others,
			'show.legend' => 0,
			'set.options' => 'color = "' . $GREY . '", linewidth = 1',
		},
	];

	plt(
		'output.file' => "$DIR/density.kernels.png",
		p             => \@panels,
		ncol          => 4,
		sharey        => 'True',    # equal standard deviation, unequal peak
	);
	return;
}

# --- 5. the grid, and weights ------------------------------------------------
# Left: cut (and from/to) decide where the grid stops, not what the curve is.
# Right: weights change the estimate itself.
sub fig_grid_weights {
	my @cut = (0, 3, 6);
	my (%cdata, %copts, @corder);
	for my $i (0 .. $#cut) {
		my $d     = density(\@precip, cut => $cut[$i]);
		my $label = sprintf 'cut = %d  [%.1f, %.1f]', $cut[$i], $d->{x}[0], $d->{x}[-1];
		push @corder, $label;
		$cdata{$label} = [ $d->{x}, $d->{y} ];
		# the three curves are the same curve, so they are drawn widest first
		# and thinnest last: all three stay visible, nested, and only their ends
		# (the "|" markers) differ
		$copts{$label} = sprintf 'color = "%s", linewidth = %s, marker = "|", markevery = [0, -1], markersize = 13',
			$RAMP[ $i + 1 ], (4.5, 2.6, 1.2)[$i];
	}

	# weights that sum to one, tripling the short eruptions
	my @w = map { $_ < 3 ? 3 : 1 } @eruptions;
	my $tot = 0;
	$tot += $_ for @w;
	@w = map { $_ / $tot } @w;
	my $h  = bw_nrd0(\@eruptions);
	my $u  = density(\@eruptions, bw => $h);
	my $wd = density(\@eruptions, bw => $h, weights => \@w);

	plt(
		'output.file' => "$DIR/density.grid.weights.png",
		ncol          => 2,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => \%cdata,
					'key.order'   => \@corder,
					'set.options' => \%copts,
					title         => '"cut: only the ends of the grid move"',
					xlabel        => 'precipitation (inches)',
					ylabel        => 'density',
					legend        => 'loc = "upper left", frameon = False, fontsize = 9',
					set_figwidth  => 11,
					set_figheight => 4.2,
					set_dpi       => 100,
				},
				rug(\@precip, 0.0015),
			],
			[
				{
					'plot.type' => 'plot',
					data        => {
						'default weights = 1/n'    => [ $u->{x},  $u->{y} ],
						'short eruptions weighted' => [ $wd->{x}, $wd->{y} ],
					},
					'key.order'   => [ 'default weights = 1/n', 'short eruptions weighted' ],
					'set.options' => {
						'default weights = 1/n'    => 'color = "' . $BLUE . '", linewidth = 2',
						'short eruptions weighted' => 'color = "' . $ORANGE . '", linewidth = 2',
					},
					title  => '"weights: one bandwidth, two estimates"',
					xlabel => 'eruption length (minutes)',
					ylabel => 'density',
					legend => 'loc = "upper left", frameon = False, fontsize = 9',
				},
				rug(\@eruptions, 0.025),
			],
		],
	);
	return;
}

fig_what();
fig_bandwidth();
fig_rules();
fig_kernels();
fig_grid_weights();
print "wrote $DIR/density.*.png\n";
