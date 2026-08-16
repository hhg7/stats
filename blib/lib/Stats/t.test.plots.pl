#!/usr/bin/env perl
# Draws the figures that illustrate t_test() in README.md (and therefore in
# read.me.pod and lib/Stats/LikeR.pm).  Author-only: it is not shipped, and it
# needs Matplotlib::Simple, python3 and matplotlib.
#
#   perl -Iblib/lib -Iblib/arch t.test.plots.pl
#
# Every figure is written to img/, at the natural size the POD and the README
# embed them at.  The data are R's own: the two examples in ?t.test -- the
# sleep data (Cushny & Peebles' extra hours of sleep under two drugs, n = 10 per
# group, the example behind "t = -1.8608, df = 17.776, p-value = 0.07939" and
# the paired "t = -4.0621, df = 9, p-value = 0.002833") and t.test(1:10, 7:20)
# ("P = .00001855" in ?t.test).  Every number drawn comes back out of t_test()
# itself, so the pictures and the module cannot drift apart.
#
# Colours follow density.plots.pl: identity (a test, a hypothesis) gets the
# categorical slots in their documented order, and anything ordered -- the
# conf_level ladder -- gets the one-hue sequential ramp instead, dark meaning
# more.
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Path 'make_path';
use POSIX 'lgamma';
use Stats::LikeR qw(t_test density mean sd);
use Matplotlib::Simple 'plt';

my $DIR = 'img';
make_path($DIR) unless -d $DIR;

# --- palette ----------------------------------------------------------------
my $BLUE   = '#2a78d6';    # categorical slot 1: the estimate, the test in hand
my $ORANGE = '#eb6834';    # categorical slot 2: the null, the rejected region
my $GREEN  = '#2f9e6b';    # categorical slot 3: the contrasting test
my $GREY   = '#8c8c88';    # recessive reference marks
my $INK    = '#3d3d3a';    # rug ticks, annotation
# the blue sequential ramp, light -> dark, for ordered quantities
my @RAMP = ('#9dc3f2', '#63a2ec', '#2a78d6', '#1c5cab', '#0d366b', '#04203f');

my $PI  = 4 * atan2(1, 1);
my $INF = 9**9**9;

# --- the data ---------------------------------------------------------------
# R's sleep: extra hours of sleep, drug 1 then drug 2, on the same ten patients
my @drug1 = (0.7, -1.6, -0.2, -1.2, -0.1,  3.4,  3.7, 0.8, 0.0, 2.0);
my @drug2 = (1.9,  0.8,  1.1,  0.1, -0.1,  4.4,  5.5, 1.6, 4.6, 3.4);
my @diff  = map { $drug1[$_] - $drug2[$_] } 0 .. $#drug1;

# ?t.test's other example, t.test(1:10, y = c(7:20)): unequal n, unequal spread
my @small = (1 .. 10);
my @big   = (7 .. 20);

# --- helpers ----------------------------------------------------------------

# Student's t density.  Only used to draw the null distribution a p-value is an
# area under; every number annotated on the figures comes from t_test().
sub t_pdf {
	my ($t, $df) = @_;
	my $ln = lgamma(($df + 1) / 2) - lgamma($df / 2)
		- 0.5 * log($df * $PI)
		- (($df + 1) / 2) * log(1 + $t * $t / $df);
	return exp($ln);
}

sub t_curve {
	my ($df, $lo, $hi, $n) = @_;
	$n = 400 unless defined $n;
	my (@x, @y);
	for my $i (0 .. $n) {
		my $t = $lo + ($hi - $lo) * $i / $n;
		push @x, $t;
		push @y, t_pdf($t, $df);
	}
	return (\@x, \@y);
}

# A filled region under the t density, as an ax.add_patch() argument.  Matplotlib
# ::Simple has no fill_between that lands on the right subplot, so the area is
# drawn as an explicit polygon closed down to y = 0 at both ends.
sub shade_t {
	my ($df, $lo, $hi, $colour, $alpha) = @_;
	$alpha = 0.45 unless defined $alpha;
	my @pt = ( sprintf '[%.6g,0]', $lo );
	for my $i (0 .. 200) {
		my $t = $lo + ($hi - $lo) * $i / 200;
		push @pt, sprintf '[%.6g,%.6g]', $t, t_pdf($t, $df);
	}
	push @pt, sprintf '[%.6g,0]', $hi;
	return sprintf 'plt.Polygon([%s], closed = True, facecolor = "%s", edgecolor = "none", alpha = %s)',
		join(',', @pt), $colour, $alpha;
}

# One vertical tick per observation, drawn as its own two-point line.
sub rug {
	my ($x, $y0, $y1, $colour, $width) = @_;
	$colour = $INK unless defined $colour;
	$width  = 1.1  unless defined $width;
	return {
		'plot.type'   => 'plot',
		data          => [ map { [ [ $_, $_ ], [ $y0, $y1 ] ] } @$x ],
		'show.legend' => 0,
		'set.options' => 'color = "' . $colour . '", linewidth = ' . $width . ', alpha = 0.7',
	};
}

# The standard error t_test() worked in: statistic = (estimate - mu) / SE, so the
# SE is recoverable from the result without recomputing it a second way.
sub se_of {
	my ($r, $mu) = @_;
	$mu = 0 unless defined $mu;
	my $est = defined $r->{estimate} ? $r->{estimate}
		: $r->{estimate_x} - $r->{estimate_y};
	return ($est - $mu) / $r->{statistic};
}

sub est_of {
	my ($r) = @_;
	return defined $r->{estimate} ? $r->{estimate} : $r->{estimate_x} - $r->{estimate_y};
}

# --- 1. what a t-test asks ---------------------------------------------------
# Left: the estimate, mu, and the standard error between them -- the three
# numbers the statistic is built from.  Middle: the null distribution the
# p-value is an area under.  Right: that area, magnified until it is visible,
# which is the only way to see what p = 0.003 looks like.
sub fig_what {
	my $r   = t_test(\@diff, mu => 0);
	my $se  = se_of($r, 0);
	my $est = $r->{estimate};
	my $t   = $r->{statistic};
	my $df  = $r->{df};

	my ($cx, $cy) = t_curve($df, -5.6, 5.6);
	my ($zx, $zy) = t_curve($df, -7, -3.4);
	my $peak = t_pdf(0, $df);
	my $ztop = t_pdf(-3.4, $df);

	plt(
		'output.file' => "$DIR/t.test.what.png",
		ncol          => 3,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ $est - $se, $est + $se ], [ 0.42, 0.42 ] ] ],
					'show.legend' => 0,
					'set.options' => 'color = "' . $ORANGE . '", linewidth = 3, marker = "|", markersize = 11',
					title         => '"estimate, mu, and the SE between them", fontsize = 11',
					xlabel        => '"extra sleep, drug 1 - drug 2 (hours)"',
					set_figwidth  => 15.5,
					set_figheight => 4.4,
					set_dpi       => 100,
					set_xlim      => '-5.4, 1.3',
					set_ylim      => '-0.08, 1.08',
					set_yticks    => '[]',
					vlines        => [
						sprintf('0, 0, 0.9, color = "%s", linewidth = 1.6, linestyle = "--"', $GREY),
						sprintf('%.10g, 0, 0.72, color = "%s", linewidth = 2.2', $est, $BLUE),
					],
					text          => [
						sprintf('0.06, 0.86, "mu = 0\n(the null)", fontsize = 9, color = "%s"', $INK),
						sprintf('%.10g, 0.75, "estimate = %.3g", fontsize = 9, ha = "center", color = "%s"', $est, $est, $BLUE),
						sprintf('%.10g, 0.42, "+/- 1 SE   (SE = %.5f)  ", fontsize = 9, ha = "right", va = "center", color = "%s"', $est - $se, $se, $ORANGE),
						sprintf('-5.2, 0.16, "the 10 paired differences", fontsize = 9, color = "%s"', $INK),
						sprintf('-5.2, 1.02, "statistic = (estimate - mu) / SE = %.4f", fontsize = 10, va = "top", color = "%s"', $t, $INK),
					],
				},
				rug(\@diff, 0, 0.1, $INK, 1.4),
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ $cx, $cy ] ],
					'show.legend' => 0,
					'set.options' => 'color = "' . $BLUE . '", linewidth = 2.2',
					title         => sprintf('"the null distribution: t with df = %g", fontsize = 11', $df),
					xlabel        => 't',
					ylabel        => 'density',
					set_ylim      => sprintf('%.6g, %.6g', -0.02 * $peak, 1.18 * $peak),
					add_patch     => [
						shade_t($df, -5.6, $t,  $ORANGE, 0.9),
						shade_t($df, -$t,  5.6, $ORANGE, 0.9),
					],
					vlines        => [
						sprintf('%.10g, 0, %.6g, color = "%s", linewidth = 2', $t, 1.02 * $peak, $ORANGE),
						sprintf('%.10g, 0, %.6g, color = "%s", linewidth = 1.2, linestyle = "--"', -$t, 1.02 * $peak, $ORANGE),
					],
					text          => [
						sprintf('%.10g, %.6g, "statistic = %.4f", fontsize = 9, ha = "left", color = "%s"', $t + 0.2, 1.04 * $peak, $t, $ORANGE),
						sprintf('-5.4, %.6g, "p_value = %.5f\n= the two shaded tails,\nwhich is why they are invisible", fontsize = 9, va = "top", color = "%s"', 0.95 * $peak, $r->{p_value}, $INK),
					],
				},
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ $zx, $zy ] ],
					'show.legend' => 0,
					'set.options' => 'color = "' . $BLUE . '", linewidth = 2.2',
					title         => sprintf('"the left tail, magnified %.0fx", fontsize = 11', $peak / $ztop),
					xlabel        => 't',
					ylabel        => 'density',
					set_xlim      => '-7, -3.4',
					set_ylim      => sprintf('%.6g, %.6g', -0.04 * $ztop, 1.5 * $ztop),
					add_patch     => [ shade_t($df, -7, $t, $ORANGE, 0.6) ],
					vlines        => sprintf('%.10g, 0, %.6g, color = "%s", linewidth = 2', $t, 1.1 * $ztop, $ORANGE),
					text          => [
						sprintf('-6.9, %.6g, "area = p_value / 2 = %.6f", fontsize = 9, va = "top", color = "%s"', 1.45 * $ztop, $r->{p_value} / 2, $INK),
						sprintf('%.10g, %.6g, "statistic = %.4f  ", fontsize = 9, ha = "right", color = "%s"', $t, 1.13 * $ztop, $t, $ORANGE),
					],
				},
			],
		],
	);
	return;
}

# --- 2. the confidence interval ---------------------------------------------
# Left: what conf_int is made of.  Right: conf_level widens it, and the level at
# which it first touches mu is exactly 1 - p_value.
sub fig_conf_int {
	my $r    = t_test(\@diff, mu => 0);
	my $se   = se_of($r, 0);
	my $est  = $r->{estimate};
	my $half = ($r->{conf_int}[1] - $r->{conf_int}[0]) / 2;
	my $tcrit = $half / $se;

	my @level = (0.80, 0.90, 0.95, 0.99, 1 - $r->{p_value}, 0.999);
	my (@lo, @hi, @ypos, @vline, @ltext);
	for my $i (0 .. $#level) {
		my $c = t_test(\@diff, conf_level => $level[$i])->{conf_int};
		my $y = $#level - $i;
		push @lo, $c->[0];
		push @hi, $c->[1];
		push @ypos, $y;
		push @ltext, sprintf('%.10g, %.10g, "conf_level = %s", fontsize = 9, va = "bottom", ha = "center", color = "%s"',
			$est, $y + 0.13,
			($i == 4 ? sprintf('1 - p_value = %.5f', $level[$i]) : sprintf('%g', $level[$i])),
			$RAMP[$i]);
	}

	my @bars = map {
		{
			'plot.type'   => 'plot',
			data          => [ [ [ $lo[$_], $hi[$_] ], [ $ypos[$_], $ypos[$_] ] ] ],
			'show.legend' => 0,
			'set.options' => sprintf('color = "%s", linewidth = 3.2, marker = "|", markersize = 12', $RAMP[$_]),
		}
	} 0 .. $#level;

	plt(
		'output.file' => "$DIR/t.test.conf.int.png",
		ncol          => 2,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ $r->{conf_int}[0], $r->{conf_int}[1] ], [ 0.5, 0.5 ] ] ],
					'show.legend' => 0,
					'set.options' => 'color = "' . $BLUE . '", linewidth = 3.2, marker = "|", markersize = 14',
					title         => '"conf_int = estimate +/- t(conf_level, df) * SE"',
					xlabel        => '"extra sleep, drug 1 - drug 2 (hours)"',
					set_figwidth  => 12,
					set_figheight => 4.6,
					set_dpi       => 100,
					set_xlim      => '-3.0, 0.65',
					set_ylim      => '0, 1.25',
					set_yticks    => '[]',
					vlines        => [
						sprintf('0, 0, 1.1, color = "%s", linewidth = 1.6, linestyle = "--"', $GREY),
						sprintf('%.10g, 0.34, 0.66, color = "%s", linewidth = 2', $est, $ORANGE),
					],
					text          => [
						sprintf('0.03, 1.0, "mu = 0", fontsize = 9, color = "%s"', $INK),
						sprintf('%.10g, 0.7, "estimate\n%.3g", fontsize = 9, ha = "center", color = "%s"', $est, $est, $ORANGE),
						sprintf('%.10g, 0.38, "%.5f", fontsize = 9, ha = "center", color = "%s"', $r->{conf_int}[0], $r->{conf_int}[0], $BLUE),
						sprintf('%.10g, 0.38, "%.5f", fontsize = 9, ha = "center", color = "%s"', $r->{conf_int}[1], $r->{conf_int}[1], $BLUE),
						sprintf('-2.95, 1.15, "t(0.95, df = %g) = %.5f      SE = %.5f      half-width = %.5f", fontsize = 9, va = "top", color = "%s"',
							$r->{df}, $tcrit, $se, $half, $INK),
						sprintf('-2.95, 0.15, "mu lies outside the interval, and p_value = %.5f < 0.05:\nthe same statement, twice", fontsize = 9, va = "bottom", color = "%s"',
							$r->{p_value}, $INK),
					],
				},
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ $est, $est ], [ -0.5, scalar @level ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1.2, linestyle = ":"', $GREY),
					title         => '"a wider conf_level is a wider conf_int"',
					xlabel        => '"extra sleep, drug 1 - drug 2 (hours)"',
					set_xlim      => '-3.8, 0.75',
					set_ylim      => sprintf('-0.75, %.10g', $#level + 0.8),
					set_yticks    => '[]',
					vlines        => sprintf('0, -0.75, %.10g, color = "%s", linewidth = 1.6, linestyle = "--"', $#level + 0.75, $GREY),
					text          => [
						@ltext,
						sprintf('0.06, -0.6, "mu = 0", fontsize = 9, color = "%s"', $INK),
						sprintf('-3.75, -0.6, "the interval first reaches mu at conf_level = 1 - p_value", fontsize = 9, color = "%s"', $INK),
					],
				},
				@bars,
			],
		],
	);
	return;
}

# --- 3. alternative ----------------------------------------------------------
# Top row: which area of the null distribution the p-value is.  Bottom row: the
# interval that goes with it -- a one-sided alternative gives a one-sided
# interval, with an infinite bound.
sub fig_alternative {
	my @alt = ('two.sided', 'less', 'greater');
	my $lo_x = -3.6;
	my $hi_x = 1.4;

	my (@top, @bottom);
	for my $i (0 .. $#alt) {
		my $r  = t_test(\@drug1, \@drug2, alternative => $alt[$i]);
		my $t  = $r->{statistic};
		my $df = $r->{df};
		my ($cx, $cy) = t_curve($df, -4.4, 4.4);
		my $peak = t_pdf(0, $df);

		my @patch;
		if ($alt[$i] eq 'two.sided') {
			@patch = (shade_t($df, -4.4, $t, $ORANGE), shade_t($df, -$t, 4.4, $ORANGE));
		} elsif ($alt[$i] eq 'less') {
			@patch = (shade_t($df, -4.4, $t, $ORANGE));
		} else {
			@patch = (shade_t($df, $t, 4.4, $ORANGE));
		}

		my %first;
		if ($i == 0) {
			%first = (set_figwidth => 12.5, set_figheight => 5.3, set_dpi => 100);
		}
		push @top, {
			'plot.type'   => 'plot',
			data          => [ [ $cx, $cy ] ],
			'show.legend' => 0,
			'set.options' => 'color = "' . $BLUE . '", linewidth = 2',
			title         => sprintf('"alternative => \'%s\'   p_value = %.5f"', $alt[$i], $r->{p_value}),
			xlabel        => sprintf('t   (df = %.4f)', $df),
			ylabel        => $i == 0 ? 'density' : undef,
			set_ylim      => sprintf('%.6g, %.6g', -0.02 * $peak, 1.2 * $peak),
			add_patch     => \@patch,
			vlines        => sprintf('%.10g, 0, %.6g, color = "%s", linewidth = 2', $t, 1.05 * $peak, $ORANGE),
			text          => sprintf('%.10g, %.6g, "statistic = %.4f", fontsize = 9, ha = "left", color = "%s"',
				$t + 0.15, 1.08 * $peak, $t, $ORANGE),
			%first,
		};
		delete $top[-1]{ylabel} unless defined $top[-1]{ylabel};

		# the interval, on the data scale, with the infinite end drawn as an arrow
		my ($l, $u) = @{ $r->{conf_int} };
		my $dl = $l == -$INF ? $lo_x + 0.12 : $l;
		my $du = $u ==  $INF ? $hi_x - 0.12 : $u;
		my @arrow;
		push @arrow, sprintf('%.10g, 0.5, "<", fontsize = 12, ha = "right", va = "center", color = "%s"', $dl, $BLUE)
			if $l == -$INF;
		push @arrow, sprintf('%.10g, 0.5, ">", fontsize = 12, ha = "left", va = "center", color = "%s"', $du, $BLUE)
			if $u == $INF;

		push @bottom, [
			{
				'plot.type'   => 'plot',
				data          => [ [ [ $dl, $du ], [ 0.5, 0.5 ] ] ],
				'show.legend' => 0,
				'set.options' => 'color = "' . $BLUE . '", linewidth = 3.2, marker = "|", markersize = 13',
				title         => sprintf('"conf_int = [%s, %s]", fontsize = 10',
					($l == -$INF ? '-Inf' : sprintf('%.5f', $l)),
					($u ==  $INF ? 'Inf'  : sprintf('%.5f', $u))),
				xlabel        => '"mean(drug 1) - mean(drug 2), hours of extra sleep"',
				set_xlim      => sprintf('%g, %g', $lo_x, $hi_x),
				set_ylim      => '0, 1',
				set_yticks    => '[]',
				vlines        => [
					sprintf('0, 0, 1, color = "%s", linewidth = 1.6, linestyle = "--"', $GREY),
					sprintf('%.10g, 0.36, 0.64, color = "%s", linewidth = 2', est_of($r), $ORANGE),
				],
				text          => [
					@arrow,
					sprintf('0.05, 0.86, "mu = 0", fontsize = 9, color = "%s"', $INK),
					sprintf('%.10g, 0.68, "estimate_x - estimate_y = %.3g", fontsize = 8.5, ha = "center", color = "%s"',
						est_of($r), est_of($r), $ORANGE),
				],
			},
		];
	}

	plt(
		'output.file' => "$DIR/t.test.alternative.png",
		ncol          => 3,
		p             => [ @top, @bottom ],
	);
	return;
}

# --- 4. mu, p_value and conf_int are one statement ---------------------------
# Left: p_value as a function of mu.  It peaks at 1 over the estimate and falls
# through alpha exactly at the two conf_int bounds -- the interval is the set of
# mu the test does not reject.  Right: moving mu moves the statistic and the
# p-value, and leaves conf_int where it was.
sub fig_duality {
	my $r    = t_test(\@diff, mu => 0);
	my $est  = $r->{estimate};
	my ($lo, $hi) = @{ $r->{conf_int} };

	my ($a, $b) = ($est - 2.6, $est + 2.6);
	my (@mu, @p, @stat);
	for my $i (0 .. 400) {
		my $m = $a + ($b - $a) * $i / 400;
		my $s = t_test(\@diff, mu => $m);
		push @mu,   $m;
		push @p,    $s->{p_value};
		push @stat, $s->{statistic};
	}

	# the conf_int at three values of mu: identical, three times over
	my @mus = (-3, 0, 1);
	my (@ci_bar, @ci_text);
	for my $i (0 .. $#mus) {
		my $s = t_test(\@diff, mu => $mus[$i]);
		my $y = $#mus - $i;
		push @ci_bar, {
			'plot.type'   => 'plot',
			data          => [ [ [ $s->{conf_int}[0], $s->{conf_int}[1] ], [ $y, $y ] ] ],
			'show.legend' => 0,
			'set.options' => sprintf('color = "%s", linewidth = 3, marker = "|", markersize = 12', $BLUE),
		};
		push @ci_text, sprintf('%.10g, %.10g, "mu = %g:  statistic = %.4f,  p_value = %.5g", fontsize = 9, va = "bottom", ha = "center", color = "%s"',
			$est, $y + 0.12, $mus[$i], $s->{statistic}, $s->{p_value}, $INK);
	}

	plt(
		'output.file' => "$DIR/t.test.duality.png",
		ncol          => 2,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ \@mu, \@p ] ],
					'show.legend' => 0,
					'set.options' => 'color = "' . $BLUE . '", linewidth = 2.2',
					title         => '"conf_int is the set of mu the test does not reject"',
					xlabel        => 'mu',
					ylabel        => 'p_value',
					set_figwidth  => 12,
					set_figheight => 4.6,
					set_dpi       => 100,
					set_ylim      => '-0.03, 1.12',
					hlines        => sprintf('0.05, %.10g, %.10g, color = "%s", linewidth = 1.4, linestyle = "--"', $a, $b, $ORANGE),
					vlines        => [
						sprintf('%.10g, 0, 0.05, color = "%s", linewidth = 1.8', $lo, $ORANGE),
						sprintf('%.10g, 0, 0.05, color = "%s", linewidth = 1.8', $hi, $ORANGE),
						sprintf('%.10g, 0, 1, color = "%s", linewidth = 1.4, linestyle = ":"', $est, $GREY),
					],
					text          => [
						sprintf('%.10g, 0.09, "conf_int[0]\n%.5f", fontsize = 9, ha = "right", color = "%s"', $lo - 0.05, $lo, $ORANGE),
						sprintf('%.10g, 0.09, "conf_int[1]\n%.5f", fontsize = 9, ha = "left", color = "%s"', $hi + 0.05, $hi, $ORANGE),
						sprintf('%.10g, 1.02, "p_value = 1 at mu = estimate = %.3g", fontsize = 9, ha = "center", color = "%s"', $est, $est, $INK),
						sprintf('%.10g, 0.06, "1 - conf_level = 0.05", fontsize = 8.5, color = "%s"', $a + 0.05, $ORANGE),
					],
				},
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ $est, $est ], [ -0.6, 3 ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1.2, linestyle = ":"', $GREY),
					title         => '"mu moves the statistic and the p-value, never conf_int"',
					xlabel        => '"extra sleep, drug 1 - drug 2 (hours)"',
					set_xlim      => '-3.6, 1.7',
					set_ylim      => '-0.8, 2.8',
					set_yticks    => '[]',
					vlines        => [ map {
						sprintf('%.10g, -0.8, 2.65, color = "%s", linewidth = 1.4, linestyle = "--"', $_, $GREY)
					} @mus ],
					text          => [
						@ci_text,
						( map { sprintf('%.10g, -0.72, "mu = %g", fontsize = 8.5, ha = "center", color = "%s"', $_, $_, $GREY) } @mus ),
						sprintf('%.10g, 2.72, "one conf_int, drawn three times", fontsize = 9, ha = "center", va = "top", color = "%s"', $est, $INK),
					],
				},
				@ci_bar,
			],
		],
	);
	return;
}

# --- 5. paired, var_equal, and the Welch df ----------------------------------
# Left: the same twenty numbers, three ways.  Middle and right: what var_equal
# costs -- the Welch df as the two spreads are pulled apart, against the pooled
# df that never moves, and the p-value each of them reaches.
sub fig_designs {
	my @test = (
		[ 'paired => 1'                => t_test(\@drug1, \@drug2, paired => 1) ],
		[ 'var_equal => 1'             => t_test(\@drug1, \@drug2, var_equal => 1) ],
		[ 'the default (Welch)'        => t_test(\@drug1, \@drug2) ],
	);

	my (@bar, @label);
	for my $i (0 .. $#test) {
		my ($name, $r) = @{ $test[$i] };
		my $y = $#test - $i;
		push @bar, {
			'plot.type'   => 'plot',
			data          => [ [ [ $r->{conf_int}[0], $r->{conf_int}[1] ], [ $y, $y ] ] ],
			'show.legend' => 0,
			'set.options' => sprintf('color = "%s", linewidth = 3.2, marker = "|", markersize = 13', ($BLUE, $GREEN, $ORANGE)[$i]),
		};
		push @label, sprintf('-3.5, %.10g, "%s\ndf = %.4f,  statistic = %.4f,  p_value = %.5f", fontsize = 9, va = "bottom", color = "%s"',
			$y + 0.1, $name, $r->{df}, $r->{statistic}, $r->{p_value}, ($BLUE, $GREEN, $ORANGE)[$i]);
	}

	# the Welch df as the spread of y is scaled about its own mean
	my $my = mean(\@big);
	my (@ratio, @welch_df, @pool_df, @welch_p, @pool_p);
	for my $i (0 .. 160) {
		my $s = exp( log(0.15) + (log(12) - log(0.15)) * $i / 160 );
		my @y = map { $my + $s * ($_ - $my) } @big;
		my $w = t_test(\@small, \@y);
		my $v = t_test(\@small, \@y, var_equal => 1);
		push @ratio,    sd(\@y) / sd(\@small);
		push @welch_df, $w->{df};
		push @pool_df,  $v->{df};
		push @welch_p,  $w->{p_value};
		push @pool_p,   $v->{p_value};
	}

	my $equal = t_test(\@small, \@big);    # the ratio-1 point of the sweep

	plt(
		'output.file' => "$DIR/t.test.designs.png",
		ncol          => 3,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ 0, 0 ], [ -0.75, 2.65 ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1.6, linestyle = "--"', $GREY),
					title         => '"paired, var_equal, Welch: the same twenty numbers", fontsize = 11',
					xlabel        => '"mean difference, drug 1 - drug 2 (hours)"',
					set_figwidth  => 15.5,
					set_figheight => 4.4,
					set_dpi       => 100,
					set_xlim      => '-3.75, 0.9',
					set_ylim      => '-0.8, 2.95',
					set_yticks    => '[]',
					text          => [
						@label,
						sprintf('0.05, -0.72, "mu = 0", fontsize = 9, color = "%s"', $INK),
					],
				},
				@bar,
			],
			[
				{
					'plot.type' => 'plot',
					data        => {
						'the default (Welch)' => [ \@ratio, \@welch_df ],
						'var_equal => 1'      => [ \@ratio, \@pool_df ],
					},
					'key.order'   => [ 'var_equal => 1', 'the default (Welch)' ],
					'set.options' => {
						'the default (Welch)' => 'color = "' . $ORANGE . '", linewidth = 2.2',
						'var_equal => 1'      => 'color = "' . $GREEN . '", linewidth = 2, linestyle = "--"',
					},
					title      => '"df, as the spreads of t.test(1:10, 7:20) separate", fontsize = 11',
					xlabel     => '"sd(y) / sd(x)"',
					ylabel     => 'df',
					legend     => 'loc = "upper right", frameon = False, fontsize = 9',
					set_xscale => '"log"',
					set_ylim   => '7.5, 25.5',
					hlines     => [
						sprintf('9, %.6g, %.6g, color = "%s", linewidth = 1.1, linestyle = ":"',  $ratio[0], $ratio[-1], $GREY),
						sprintf('13, %.6g, %.6g, color = "%s", linewidth = 1.1, linestyle = ":"', $ratio[0], $ratio[-1], $GREY),
					],
					text       => [
						sprintf('%.6g, 8.3, "floor: n(x) - 1 = 9", fontsize = 8.5, color = "%s"', $ratio[0] * 1.05, $INK),
						sprintf('%.6g, 11.4, "floor: n(y) - 1 = 13", fontsize = 8.5, ha = "right", color = "%s"', $ratio[-1] * 0.95, $INK),
						sprintf('%.6g, 25.1, "equal spreads: Welch df = %.4f", fontsize = 8.5, ha = "left", va = "top", color = "%s"',
							$ratio[0] * 1.05, $equal->{df}, $INK),
					],
				},
			],
			[
				{
					'plot.type' => 'plot',
					data        => {
						'the default (Welch)' => [ \@ratio, \@welch_p ],
						'var_equal => 1'      => [ \@ratio, \@pool_p ],
					},
					'key.order'   => [ 'var_equal => 1', 'the default (Welch)' ],
					'set.options' => {
						'the default (Welch)' => 'color = "' . $ORANGE . '", linewidth = 2.2',
						'var_equal => 1'      => 'color = "' . $GREEN . '", linewidth = 2, linestyle = "--"',
					},
					title      => '"and the p_value each of them reaches", fontsize = 11',
					xlabel     => '"sd(y) / sd(x)"',
					ylabel     => 'p_value',
					legend     => 'loc = "upper left", frameon = False, fontsize = 9',
					set_xscale => '"log"',
					set_yscale => '"log"',
					hlines     => sprintf('0.05, %.6g, %.6g, color = "%s", linewidth = 1.2, linestyle = ":"',
						$ratio[0], $ratio[-1], $GREY),
					text       => sprintf('%.6g, 0.062, "0.05", fontsize = 8.5, color = "%s"', $ratio[0] * 1.05, $INK),
				},
			],
		],
	);
	return;
}

# --- 6. what a falling p-value looks like ------------------------------------
# Two samples from two different distributions, pulled apart until t_test()
# reports each of four p-values five orders of magnitude apart.  Top row: the
# two distributions, by this module's own density().  Bottom row: the conf_int
# that goes with each.  Neither df nor the width of the interval changes -- only
# how far the interval sits from mu.
sub fig_p_and_ci {
	my @target = (0.5, 0.05, 0.001, 1e-5);

	# the shift of x that makes t_test() report exactly that p.  Bisection on
	# t_test's own p_value, so the labels are the module's numbers and not a
	# quantile recomputed here; p rises with delta up to the null, so the
	# bracket stops at the shift that would put the estimate on mu.
	my $null_shift = mean(\@drug2) - mean(\@drug1);
	my @shift;
	for my $target (@target) {
		my ($lo, $hi) = (-8, $null_shift);
		for (1 .. 200) {
			my $mid = ($lo + $hi) / 2;
			my $p   = t_test([ map { $_ + $mid } @drug1 ], \@drug2)->{p_value};
			if ($p < $target) { $lo = $mid } else { $hi = $mid }
		}
		push @shift, ($lo + $hi) / 2;
	}

	# a common grid, so the four panels are directly comparable
	my %grid = (from => -9, to => 9, n => 512);
	my $dy   = density(\@drug2, %grid);

	my (@top, @bottom);
	for my $i (0 .. $#target) {
		my @x = map { $_ + $shift[$i] } @drug1;
		my $r = t_test(\@x, \@drug2);
		my $d = density(\@x, %grid);

		my %first;
		if ($i == 0) {
			%first = (set_figwidth => 16, set_figheight => 5.8, set_dpi => 100);
		}
		# only the first panel carries a legend; the rest are the same two colours
		my $series = $i == 0
			? {
				'plot.type' => 'plot',
				data        => {
					'x = drug 1, shifted' => [ $d->{x},  $d->{y} ],
					'y = drug 2'          => [ $dy->{x}, $dy->{y} ],
				},
				'key.order'   => [ 'x = drug 1, shifted', 'y = drug 2' ],
				'set.options' => {
					'x = drug 1, shifted' => 'color = "' . $BLUE . '", linewidth = 2',
					'y = drug 2'          => 'color = "' . $ORANGE . '", linewidth = 2',
				},
				legend => 'loc = "upper left", frameon = False, fontsize = 8.5',
			}
			: {
				'plot.type'   => 'plot',
				data          => [ [ $d->{x}, $d->{y} ] ],
				'show.legend' => 0,
				'set.options' => 'color = "' . $BLUE . '", linewidth = 2',
			};

		push @top, [
			{
				%$series,
				title      => sprintf('"p_value = %s", fontsize = 11', p_label($r->{p_value})),
				xlabel     => '"extra sleep (hours)"',
				ylabel     => $i == 0 ? 'density' : undef,
				set_xlim   => '-8, 8.5',
				set_ylim   => '-0.008, 0.30',
				%first,
			},
			( $i == 0 ? () : {
				'plot.type'   => 'plot',
				data          => [ [ $dy->{x}, $dy->{y} ] ],
				'show.legend' => 0,
				'set.options' => 'color = "' . $ORANGE . '", linewidth = 2',
			} ),
			rug(\@x,      0, 0.012, $BLUE,   1.3),
			rug(\@drug2,  0, 0.012, $ORANGE, 1.3),
		];
		delete $top[-1][0]{ylabel} unless defined $top[-1][0]{ylabel};

		push @bottom, [
			{
				'plot.type'   => 'plot',
				data          => [ [ [ 0, 0 ], [ 0, 1 ] ] ],
				'show.legend' => 0,
				'set.options' => sprintf('color = "%s", linewidth = 1.6, linestyle = "--"', $GREY),
				title         => sprintf('"conf_int = [%.4f, %.4f]", fontsize = 11', @{ $r->{conf_int} }),
				xlabel        => '"mean(x) - mean(y), hours of extra sleep"',
				set_xlim      => '-7.6, 2.4',
				set_ylim      => '0, 1',
				set_yticks    => '[]',
				vlines        => sprintf('%.10g, 0.36, 0.64, color = "%s", linewidth = 2', est_of($r), $ORANGE),
				text          => [
					sprintf('0.15, 0.88, "mu = 0", fontsize = 8.5, color = "%s"', $INK),
					sprintf('%.10g, 0.66, "estimate = %.4f", fontsize = 8.5, ha = "center", color = "%s"',
						est_of($r), est_of($r), $ORANGE),
					sprintf('%.10g, 0.3, "statistic = %.4f", fontsize = 8.5, ha = "center", va = "top", color = "%s"',
						est_of($r), $r->{statistic}, $INK),
					sprintf('-7.4, 0.08, "df = %.4f,  width = %.4f", fontsize = 8.5, color = "%s"',
						$r->{df}, $r->{conf_int}[1] - $r->{conf_int}[0], $GREY),
				],
			},
			{
				'plot.type'   => 'plot',
				data          => [ [ [ $r->{conf_int}[0], $r->{conf_int}[1] ], [ 0.5, 0.5 ] ] ],
				'show.legend' => 0,
				'set.options' => sprintf('color = "%s", linewidth = 3.2, marker = "|", markersize = 13', $RAMP[ $i + 1 ]),
			},
		];
	}

	plt(
		'output.file' => "$DIR/t.test.p.and.ci.png",
		ncol          => scalar @target,
		suptitle      => 'the sleep data with drug 1 shifted: as the two distributions separate the p-value falls '
			. 'by five orders of magnitude while df and the width of conf_int never move',
		p             => [ @top, @bottom ],
	);
	return;
}

# p-values spanning five orders of magnitude do not all want the same format
sub p_label {
	my ($p) = @_;
	return $p >= 0.001 ? sprintf('%g', $p) : sprintf('%.0e', $p);
}

fig_what();
fig_conf_int();
fig_alternative();
fig_duality();
fig_designs();
fig_p_and_ci();
print "wrote $DIR/t.test.*.png\n";
