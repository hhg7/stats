#!/usr/bin/env perl
# Draws the figures that illustrate skew() and kurtosis() in README.md (and
# therefore in read.me.pod and lib/Stats/LikeR.pm).  Author-only: it is not
# shipped, and it needs Matplotlib::Simple, python3 and matplotlib.
#
#   perl -Iblib/lib -Iblib/arch skew.kurtosis.plots.pl
#
# Every figure is written to img/, at the natural size the POD and the README
# embed them at.  Nothing outside Stats::LikeR and Matplotlib::Simple is used:
# the samples come from rnorm() and runif(), the curves from density() and
# dnorm(), and every number printed on a panel comes back out of skew() and
# kurtosis() themselves, so the pictures and the module cannot drift apart.
#
# srand(42) is R's set.seed(42): the samples, and therefore the numbers in the
# titles, are the same on every run of this script on this perl.  They are not
# the same on a long-double or quadmath perl -- a different NV width draws a
# different stream -- so the committed images are the ones this script writes
# on the default (double) perl.
#
# Every sample is standardized to mean 0, sd 1 before it is drawn, because skew
# and kurtosis are both scale-free: putting each shape on the same footing as
# N(0, 1) is what makes the comparison to the grey normal curve honest.  The
# panels of a row also share their axis limits, so a taller peak or a fatter
# tail is a real difference and not a difference in scale.
#
# Colours follow density.plots.pl and t.test.plots.pl: identity (a shape) gets
# the categorical slots in their documented order, and the normal reference --
# the thing every panel is measured against, never the subject -- is the
# recessive grey.
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Path 'make_path';
use Stats::LikeR qw(density dnorm kurtosis mean median rnorm runif sd skew);
use Matplotlib::Simple 'plt';

my $DIR = 'img';
make_path($DIR) unless -d $DIR;

# --- palette ----------------------------------------------------------------
my $BLUE   = '#2a78d6';    # categorical slot 1: the positive departure
my $ORANGE = '#eb6834';    # categorical slot 2: the negative departure
my $GREEN  = '#2f9e6b';    # categorical slot 3: the normal sample itself
my $GREY   = '#8c8c88';    # the N(0, 1) reference every panel is measured against
my $INK    = '#3d3d3a';    # annotation

my $N = 50_000;            # large enough that the curves are shape, not noise

srand(42);

# --- the samples -------------------------------------------------------------

# z-scores, so every shape is drawn on the N(0, 1) footing skew() and kurtosis()
# already report it on
sub standardize {
	my ($x) = @_;
	my $m = mean($x);
	my $s = sd($x);
	return [ map { ($_ - $m) / $s } @$x ];
}

# the baseline: a normal sample, skew ~ 0 and excess kurtosis ~ 0
my $normal = standardize(rnorm(n => $N, mean => 0, sd => 1));

# a long right tail: the log-normal, exp() of a normal.  This is R's
# rlnorm(n, meanlog = 0, sdlog = 0.6), whose population skewness is
# (e^s^2 + 2) * sqrt(e^s^2 - 1) = 2.26 at s = 0.6.
my $right = standardize([ map { exp(0.6 * $_) } @{ rnorm(n => $N, mean => 0, sd => 1) } ]);

# a long left tail: the same shape, mirrored, so the two panels differ in the
# sign of the third moment and in nothing else
my $left = standardize([ map { -$_ } @$right ]);

# light tails: the uniform, whose excess kurtosis is -6/5 exactly
my $light = standardize(runif(n => $N, min => -1, max => 1));

# heavy tails: a scale mixture of normals -- one observation in ten drawn with
# three times the spread.  The mixture is still symmetric (skew 0), and its
# population excess kurtosis is 3 * (0.9 + 0.1 * 3^4) / (0.9 + 0.1 * 3^2)^2 - 3
# = 5.33, so the panel isolates the fourth moment from the third.
my $wide   = rnorm(n => $N, mean => 0, sd => 3);
my $narrow = rnorm(n => $N, mean => 0, sd => 1);
my $pick   = runif(n => $N, min => 0, max => 1);
my $heavy  = standardize([ map { $pick->[$_] < 0.1 ? $wide->[$_] : $narrow->[$_] } 0 .. $N - 1 ]);

# --- helpers -----------------------------------------------------------------

# The N(0, 1) curve every panel is measured against, on its own grid.
sub normal_curve {
	my ($lo, $hi, $n) = @_;
	$n = 512 unless defined $n;
	my @x = map { $lo + ($hi - $lo) * $_ / $n } 0 .. $n;
	return (\@x, dnorm(\@x));
}

# One panel: the sample's density() in colour, the standard normal dashed in
# grey behind it.  %extra is everything that is particular to the panel --
# labels, limits, annotation, and (on the first panel of a figure) the size of
# the whole figure.
sub shape_panel {
	my ($x, $label, $colour, $title, $lo, $hi, $extra) = @_;
	my $d = density($x, from => $lo, to => $hi, n => 512);
	my ($nx, $ny) = normal_curve($lo, $hi);
	return {
		'plot.type' => 'plot',
		data        => {
			'normal, for comparison' => [ $nx, $ny ],
			$label                   => [ $d->{x}, $d->{y} ],
		},
		'key.order'   => [ 'normal, for comparison', $label ],
		'set.options' => {
			'normal, for comparison' => 'color = "' . $GREY . '", linewidth = 1.4, linestyle = "--"',
			$label                   => 'color = "' . $colour . '", linewidth = 2.2',
		},
		title    => $title,
		legend   => 'loc = "upper right", frameon = False, fontsize = 8',
		set_xlim => sprintf('%g, %g', $lo, $hi),
		%{ $extra || {} },
	};
}

# --- 1. skew -----------------------------------------------------------------
# Three shapes on one row: a long left tail, no tail either way, a long right
# tail.  Each panel marks the mean and the median, because which side of the
# median the mean falls on is the sign of skew(), read straight off the picture.
sub fig_skew {
	my @shape = (
		[ $left,   'negative skew',   $ORANGE ],
		[ $normal, 'a normal sample', $GREEN ],
		[ $right,  'positive skew',   $BLUE ],
	);
	my ($lo, $hi) = (-5, 5);
	my $TOP = 0.85;    # one y axis for all three panels: the peaks are comparable

	my @panels;
	for my $i (0 .. $#shape) {
		my ($x, $label, $colour) = @{ $shape[$i] };
		my $g  = skew($x);
		my $mn = mean($x);
		my $md = median($x);

		my %first;
		if ($i == 0) {
			%first = (set_figwidth => 14, set_figheight => 4.6, set_dpi => 100);
		}

		push @panels, [
			shape_panel($x, $label, $colour,
				sprintf('"%s: skew = %+.3f"', $label, $g),
				$lo, $hi,
				{
					xlabel   => 'standard deviations from the mean',
					($i == 0 ? (ylabel => 'density') : ()),
					set_ylim => sprintf('%.6g, %.6g', -0.02 * $TOP, $TOP),
					vlines   => [
						sprintf('%.10g, 0, %.6g, color = "%s", linewidth = 1.8', $mn, 0.70 * $TOP, $INK),
						sprintf('%.10g, 0, %.6g, color = "%s", linewidth = 1.8, linestyle = ":"', $md, 0.70 * $TOP, $colour),
					],
					text     => [
						# the readout sits in the corner rather than on the lines,
						# so it cannot collide with them when they nearly coincide
						sprintf('%.10g, %.6g, "mean      %+.3f", fontsize = 8.5, ha = "left", va = "top", color = "%s"',
							$lo + 0.15, 0.97 * $TOP, $mn, $INK),
						sprintf('%.10g, %.6g, "median  %+.3f", fontsize = 8.5, ha = "left", va = "top", color = "%s"',
							$lo + 0.15, 0.90 * $TOP, $md, $colour),
						sprintf('%.10g, %.6g, "%s", fontsize = 8.5, ha = "left", va = "top", color = "%s"',
							$lo + 0.15, 0.34 * $TOP,
							($g > 0.05    ? 'the long tail is on the right,\nand it pulls the mean\nabove the median'
							 : $g < -0.05 ? 'the long tail is on the left,\nand it pulls the mean\nbelow the median'
							 :              'symmetric: the mean and the\nmedian sit on top of each other'),
							$INK),
					],
					%first,
				},
			),
		];
	}

	plt(
		'output.file' => "$DIR/skew.what.png",
		ncol          => 3,
		p             => \@panels,
	);
	return;
}

# --- 2. kurtosis -------------------------------------------------------------
# The same three-shape row, but the departure is in the fourth moment, so it is
# the tails that move and not the side the mass is on.  On a linear axis (top
# row) a heavy tail is nearly invisible -- the whole difference is a number too
# small to see -- so the bottom row repeats each panel with a logarithmic
# density, which is the only way to look at a tail.
sub fig_kurtosis {
	my @shape = (
		[ $light,  'negative excess kurtosis', $ORANGE ],
		[ $normal, 'a normal sample',          $GREEN ],
		[ $heavy,  'positive excess kurtosis', $BLUE ],
	);
	my ($lo, $hi) = (-6, 6);
	my $TOP = 0.62;    # one y axis across the top row: the peaks are comparable

	my (@top, @bottom);
	for my $i (0 .. $#shape) {
		my ($x, $label, $colour) = @{ $shape[$i] };
		my $g2 = kurtosis($x);

		my %first;
		if ($i == 0) {
			%first = (set_figwidth => 14, set_figheight => 8.2, set_dpi => 100);
		}

		push @top, [
			shape_panel($x, $label, $colour,
				sprintf('"%s: kurtosis = %+.3f"', $label, $g2),
				$lo, $hi,
				{
					($i == 0 ? (ylabel => 'density') : ()),
					set_ylim => sprintf('%.6g, %.6g', -0.02 * $TOP, $TOP),
					text     => sprintf('%.10g, %.6g, "%s", fontsize = 8.5, ha = "left", va = "top", color = "%s"',
						$lo + 0.2, 0.97 * $TOP,
						($g2 > 0.05    ? 'a sharper peak, bought with\nmass moved out into the tails'
						 : $g2 < -0.05 ? 'a flat shoulder and short tails:\nno mass left for the extremes'
						 :               'the shape the other two\nare measured against'),
						$INK),
					%first,
				},
			),
		];

		# the same two curves on a log density: three decades of tail that the
		# linear panel above cannot show
		push @bottom, [
			shape_panel($x, $label, $colour,
				'"the same tails, log density"',
				$lo, $hi,
				{
					xlabel => 'standard deviations from the mean',
					($i == 0 ? (ylabel => 'density (log scale)') : ()),
					set_ylim   => '2e-5, 1.6',
					set_yscale => '"log"',
				},
			),
		];
	}

	plt(
		'output.file' => "$DIR/kurtosis.what.png",
		ncol          => 3,
		p             => [ @top, @bottom ],
	);
	return;
}

fig_skew();
fig_kurtosis();
print "wrote $DIR/skew.what.png and $DIR/kurtosis.what.png\n";
