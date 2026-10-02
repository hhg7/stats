#!/usr/bin/env perl
# Draws the figures that illustrate aov() and anova() in README.md (and
# therefore in read.me.pod and lib/Stats/LikeR.pm).  Author-only: it is not
# shipped, and it needs Matplotlib::Simple, python3 and matplotlib.
#
#   perl -Iblib/lib -Iblib/arch aov.anova.plots.pl
#
# Every figure is written to img/.  The data are R's own, the same corpora
# t/aov.R.t and t/anova.R.t freeze:
#
#   PlantGrowth       src/library/datasets/man/PlantGrowth.Rd, weight by group
#                     (ctrl, trt1, trt2; n = 10 each), given as a named list so
#                     that aov() stacks it.
#   warpbreaks        src/library/datasets/man/warpbreaks.Rd, breaks ~ wool *
#                     tension: balanced, 9 looms per cell.
#   LifeCycleSavings  src/library/stats/man/anova.lm.Rd's
#                     sr ~ pop15 + pop75 + dpi + ddpi: 50 countries, numeric
#                     regressors, so not orthogonal and the term order matters.
#
# Every number drawn or annotated comes back out of aov() or anova() itself,
# except the grand mean (a weighted mean of group.stats) and the F density
# drawn under Pr(>F), so the pictures and the module cannot drift apart.
#
# Colours follow t.test.plots.pl: categorical slot 1 (blue) is the term in
# hand, grey is the residual -- what no term explains -- and orange is the tail
# a p-value is the area of.  Where terms are compared with one another (the
# term-order figure) each term takes the next categorical slot in order.
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Path 'make_path';
use POSIX 'lgamma';
use Stats::LikeR qw(aov anova);
use Matplotlib::Simple 'plt';

my $DIR = 'img';
make_path($DIR) unless -d $DIR;

# --- palette ----------------------------------------------------------------
my $BLUE   = '#2a78d6';    # categorical slot 1: the term in hand
my $ORANGE = '#eb6834';    # categorical slot 2: the Pr(>F) tail
my $GREEN  = '#2f9e6b';    # categorical slot 3
my $YELLOW = '#eda100';    # categorical slot 4
my $GREY   = '#8c8c88';    # the residual, and recessive reference marks
my $INK    = '#3d3d3a';    # observations, annotation

# --- the data ---------------------------------------------------------------
my %pg = (
	ctrl => [4.17, 5.58, 5.18, 6.11, 4.5,  4.61, 5.17, 4.53, 5.33, 5.14],
	trt1 => [4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69],
	trt2 => [6.31, 5.12, 5.54, 5.5,  5.37, 5.29, 4.92, 6.15, 5.8,  5.26],
);

my %wb = (
	breaks  => [26, 30, 54, 25, 70, 52, 51, 26, 67, 18, 21, 29, 17, 12, 18, 35, 30, 36,
		36, 21, 24, 18, 10, 43, 28, 15, 26, 27, 14, 29, 19, 29, 31, 41, 20, 44,
		42, 26, 19, 16, 39, 28, 21, 39, 29, 20, 21, 24, 17, 13, 15, 15, 16, 28],
	wool    => [ ('A') x 27, ('B') x 27 ],
	tension => [ map { ( ('L') x 9, ('M') x 9, ('H') x 9 ) } 1 .. 2 ],
);

my %lcs = (
	'sr' => [11.43, 12.07, 13.17, 5.75, 12.88, 8.79, 0.6, 11.9, 4.98, 10.78, 16.85, 3.59, 11.24, 12.64, 12.55, 10.67, 3.01, 7.7, 1.27, 9, 11.34, 14.28, 21.1, 3.98, 10.35, 15.48, 10.25, 14.65, 10.67, 7.3, 4.44, 2.02, 12.7, 12.78, 12.49, 11.14, 13.3, 11.77, 6.86, 14.13, 5.13, 2.81, 7.81, 7.56, 9.22, 18.56, 7.72, 9.24, 8.89, 4.71],
	'pop15' => [29.35, 23.32, 23.8, 41.89, 42.19, 31.72, 39.74, 44.75, 46.64, 47.64, 24.42, 46.31, 27.84, 25.06, 23.31, 25.62, 46.05, 47.32, 34.03, 41.31, 31.16, 24.52, 27.01, 41.74, 21.8, 32.54, 25.95, 24.71, 32.61, 45.04, 43.56, 41.18, 44.19, 46.26, 28.96, 31.94, 31.92, 27.74, 21.44, 23.49, 43.42, 46.12, 23.27, 29.81, 46.4, 45.25, 41.12, 28.13, 43.69, 47.2],
	'pop75' => [2.87, 4.41, 4.43, 1.67, 0.83, 2.85, 1.34, 0.67, 1.06, 1.14, 3.93, 1.19, 2.37, 4.7, 3.35, 3.1, 0.87, 0.58, 3.08, 0.96, 4.19, 3.48, 1.91, 0.91, 3.73, 2.47, 3.67, 3.25, 3.17, 1.21, 1.2, 1.05, 1.28, 1.12, 2.85, 2.28, 1.52, 2.87, 4.54, 3.73, 1.08, 1.21, 4.46, 3.43, 0.9, 0.56, 1.73, 2.72, 2.07, 0.66],
	'dpi' => [2329.68, 1507.99, 2108.47, 189.13, 728.47, 2982.88, 662.86, 289.52, 276.65, 471.24, 2496.53, 287.77, 1681.25, 2213.82, 2457.12, 870.85, 289.71, 232.44, 1900.1, 88.94, 1139.95, 1390, 1257.28, 207.68, 2449.39, 601.05, 2231.03, 1740.7, 1487.52, 325.54, 568.56, 220.56, 400.06, 152.01, 579.51, 651.11, 250.96, 768.79, 3299.49, 2630.96, 389.66, 249.87, 1813.93, 4001.89, 813.39, 138.33, 380.47, 766.54, 123.58, 242.69],
	'ddpi' => [2.87, 3.93, 3.82, 0.22, 4.56, 2.43, 2.67, 6.51, 3.08, 2.8, 3.99, 2.19, 4.32, 4.52, 3.44, 6.28, 1.48, 3.19, 1.12, 1.54, 2.99, 3.54, 8.21, 5.81, 1.57, 8.12, 3.62, 7.66, 1.76, 2.48, 3.61, 1.03, 0.67, 2, 7.48, 2.19, 2, 4.35, 3.01, 2.7, 2.96, 1.13, 2.01, 2.45, 0.53, 5.14, 10.23, 1.88, 16.71, 5.08],
);

# --- helpers ----------------------------------------------------------------

# The F density on (d1, d2) degrees of freedom.  Only used to draw the null
# distribution Pr(>F) is an area under; the p-values annotated come from aov().
sub f_pdf {
	my ($x, $d1, $d2) = @_;
	my $ln = lgamma(($d1 + $d2) / 2) - lgamma($d1 / 2) - lgamma($d2 / 2)
		+ ($d1 / 2) * log($d1 / $d2) + ($d1 / 2 - 1) * log($x)
		- (($d1 + $d2) / 2) * log(1 + $d1 * $x / $d2);
	return exp($ln);
}

# The density is evaluated from just above 0, since log(0) appears in it.
sub f_curve {
	my ($d1, $d2, $hi, $lo) = @_;
	$lo = 0 unless defined $lo;
	my (@x, @y);
	for my $i (1 .. 400) {
		my $x = $lo + ($hi - $lo) * $i / 400;
		push @x, $x;
		push @y, f_pdf($x, $d1, $d2);
	}
	return (\@x, \@y);
}

# The area under the F density from $lo to $hi, as an ax.add_patch() argument,
# closed down to y = 0 at both ends.
sub shade_f {
	my ($d1, $d2, $lo, $hi, $colour) = @_;
	my @pt = ( sprintf '[%.6g,0]', $lo );
	for my $i (0 .. 200) {
		my $x = $lo + ($hi - $lo) * $i / 200;
		push @pt, sprintf '[%.6g,%.6g]', $x, f_pdf($x, $d1, $d2);
	}
	push @pt, sprintf '[%.6g,0]', $hi;
	return sprintf 'plt.Polygon([%s], closed = True, facecolor = "%s", edgecolor = "none", alpha = 0.9)',
		join(',', @pt), $colour;
}

# A bar as an ax.add_patch() argument.  The white edge is the 2px surface gap
# that keeps adjacent fills apart.
sub rect {
	my ($x, $y, $w, $h, $colour) = @_;
	return sprintf 'plt.Rectangle((%.10g, %.10g), %.10g, %.10g, facecolor = "%s", edgecolor = "white", linewidth = 2)',
		$x, $y, $w, $h, $colour;
}

# Points, as one series with no line between them.
sub dots {
	my ($x, $y, $colour, $size, $alpha) = @_;
	return {
		'plot.type'   => 'plot',
		data          => [ [ $x, $y ] ],
		'show.legend' => 0,
		'set.options' => sprintf('color = "%s", linestyle = "none", marker = "o", markersize = %g, alpha = %g',
			$colour, $size, $alpha),
	};
}

# Several two-point segments sharing one style.
sub segments {
	my ($seg, $colour, $width, $style) = @_;
	$style = '-' unless defined $style;
	return {
		'plot.type'   => 'plot',
		data          => [ map { [ [ $_->[0], $_->[2] ], [ $_->[1], $_->[3] ] ] } @$seg ],
		'show.legend' => 0,
		'set.options' => sprintf('color = "%s", linewidth = %g, linestyle = "%s"', $colour, $width, $style),
	};
}

# --- 1. what aov() computes ---------------------------------------------------
# Left: the input -- a named list, stacked into Value ~ Group -- with the two
# kinds of deviation the table is built from.  Middle: the Sum Sq split and the
# division by Df that turns it into Mean Sq and then F.  Right: Pr(>F), the
# upper tail of F(Df, Residuals Df) beyond the F value, and that tail magnified,
# since at p = 0.016 it is too thin to see at full scale.
sub fig_aov_what {
	my $r   = aov(\%pg);
	my @lev = @{ $r->{xlevels}{Group} };
	my $gs  = $r->{'group.stats'};
	my ($n, $sum) = (0, 0);
	for my $l (@lev) {
		$n   += $gs->{size}{$l};
		$sum += $gs->{size}{$l} * $gs->{mean}{$l};
	}
	my $grand = $sum / $n;
	my ($g, $e) = @{$r}{qw(Group Residuals)};

	my (@px, @py, @within, @between, @mean_bar, @gtext);
	for my $i (0 .. $#lev) {
		my $l = $lev[$i];
		my $m = $gs->{mean}{$l};
		my @y = @{ $pg{$l} };
		for my $j (0 .. $#y) {
			my $x = $i - 0.27 + 0.5 * $j / $#y;
			push @px, $x;
			push @py, $y[$j];
			push @within, [ $x, $y[$j], $x, $m ];
		}
		push @mean_bar, [ $i - 0.33, $m, $i + 0.33, $m ];
		push @between,  [ $i + 0.4, $grand, $i + 0.4, $m ];
		push @gtext, sprintf('%d, 3.25, "mean %.4g\nsize %d", fontsize = 9, ha = "center", va = "bottom", color = "%s"',
			$i, $m, $gs->{size}{$l}, $INK);
	}

	my $ss_tot = $g->{'Sum Sq'} + $e->{'Sum Sq'};
	my $fv     = $g->{'F value'};
	my ($cx, $cy) = f_curve($g->{Df}, $e->{Df}, 8);
	my $zlo  = 0.8 * $fv;    # the magnified panel starts a little left of F
	my ($zx, $zy) = f_curve($g->{Df}, $e->{Df}, 8, $zlo);
	my $ztop = f_pdf($zlo, $g->{Df}, $e->{Df});
	my $fpk  = f_pdf($fv,  $g->{Df}, $e->{Df});

	plt(
		'output.file' => "$DIR/aov.what.png",
		ncol          => 4,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ -0.5, $#lev + 0.5 ], [ $grand, $grand ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1.4, linestyle = "--"', $GREY),
					title         => '"aov(\\%list): stacked to Value ~ Group", fontsize = 11',
					ylabel        => '"Value (PlantGrowth weight)"',
					set_figwidth  => 21,
					set_figheight => 4.8,
					set_dpi       => 100,
					set_xlim      => sprintf('-0.55, %g', $#lev + 0.6),
					set_ylim      => '2.9, 7.6',
					set_xticks    => sprintf('[%s], [%s]', join(',', 0 .. $#lev), join(',', map {"\"Group = $_\""} @lev)),
					text          => [
						@gtext,
						sprintf('-0.5, 7.5, "Sum Sq(Group) = sum of size * (mean - grand)^2 = %.5g", fontsize = 9, va = "top", color = "%s"',
							$g->{'Sum Sq'}, $BLUE),
						sprintf('-0.5, 7.2, "Sum Sq(Residuals) = sum of (Value - its group mean)^2 = %.5g", fontsize = 9, va = "top", color = "%s"',
							$e->{'Sum Sq'}, $INK),
						sprintf('-0.5, 6.9, "group.stats holds each mean and size", fontsize = 9, va = "top", color = "%s"', $INK),
						sprintf('-0.5, 6.6, "dashed: the grand mean, %.4g", fontsize = 9, va = "top", color = "%s"', $grand, $GREY),
					],
				},
				segments(\@within,  $GREY, 1),
				segments(\@mean_bar, $BLUE, 2.6),
				segments(\@between, $BLUE, 4),
				dots(\@px, \@py, $INK, 4, 0.85),
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ 0, 0 ], [ -0.4, 3.4 ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1', $INK),
					title         => '"Sum Sq / Df = Mean Sq, and their ratio is F", fontsize = 11',
					xlabel        => '"sum of squares"',
					set_xlim      => sprintf('-0.3, %.6g', 1.12 * $ss_tot),
					set_ylim      => '-0.4, 3.6',
					set_yticks    => '[2.75, 1.25, 0.25], ["Sum Sq", "Mean Sq\nGroup", "Mean Sq\nResiduals"]',
					add_patch     => [
						rect(0, 2.5, $g->{'Sum Sq'}, 0.5, $BLUE),
						rect($g->{'Sum Sq'}, 2.5, $e->{'Sum Sq'}, 0.5, $GREY),
						rect(0, 1.0, $g->{'Mean Sq'}, 0.5, $BLUE),
						rect(0, 0.0, $e->{'Mean Sq'}, 0.5, $GREY),
					],
					text          => [
						sprintf('0, 3.1, "Group: %.5g on Df %d", fontsize = 9, va = "bottom", color = "%s"',
							$g->{'Sum Sq'}, $g->{Df}, $INK),
						sprintf('%.6g, 3.1, "Residuals: %.5g on Df %d", fontsize = 9, ha = "right", va = "bottom", color = "%s"',
							$ss_tot, $e->{'Sum Sq'}, $e->{Df}, $INK),
						sprintf('%.6g, 1.25, "  %.5g / %d = %.5g", fontsize = 9, va = "center", color = "%s"',
							$g->{'Mean Sq'}, $g->{'Sum Sq'}, $g->{Df}, $g->{'Mean Sq'}, $INK),
						sprintf('%.6g, 0.25, "  %.5g / %d = %.5g", fontsize = 9, va = "center", color = "%s"',
							$e->{'Mean Sq'}, $e->{'Sum Sq'}, $e->{Df}, $e->{'Mean Sq'}, $INK),
						sprintf('%.6g, 0.75, "F value = %.5g / %.5g = %.5g", fontsize = 10, va = "center", color = "%s"',
							0.42 * $ss_tot, $g->{'Mean Sq'}, $e->{'Mean Sq'}, $fv, $INK),
						sprintf('%.6g, 0.45, "Residuals has no F of its own", fontsize = 9, va = "center", color = "%s"',
							0.42 * $ss_tot, $GREY),
					],
				},
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ $cx, $cy ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 2.2', $BLUE),
					title         => sprintf('"Pr(>F): the F(%d, %d) tail beyond F value", fontsize = 11', $g->{Df}, $e->{Df}),
					xlabel        => 'F',
					ylabel        => 'density',
					set_xlim      => '0, 8',
					set_ylim      => '0, 1.08',
					add_patch     => [ shade_f($g->{Df}, $e->{Df}, $fv, 8, $ORANGE) ],
					vlines        => sprintf('%.10g, 0, 0.55, color = "%s", linewidth = 2', $fv, $ORANGE),
					text          => [
						sprintf('%.10g, 0.57, "F value = %.4f", fontsize = 9, ha = "center", va = "bottom", color = "%s"', $fv, $fv, $INK),
						sprintf('%.10g, 0.2, "Pr(>F) = %.5f\n= the shaded area,\ntoo thin to see here", fontsize = 9, va = "bottom", color = "%s"', $fv + 0.3, $g->{'Pr(>F)'}, $INK),
						sprintf('7.8, 1.0, "the curve is the F density on\nDf = %d (Group) and Df = %d (Residuals):\nwhat F values look like when\nthe groups do not differ", fontsize = 9, ha = "right", va = "top", color = "%s"',
							$g->{Df}, $e->{Df}, $INK),
					],
				},
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ $zx, $zy ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 2.2', $BLUE),
					title         => sprintf('"the same tail, magnified %.0fx", fontsize = 11', 1.08 / (1.3 * $ztop)),
					xlabel        => 'F',
					ylabel        => 'density',
					set_xlim      => sprintf('%.6g, 8', $zlo),
					set_ylim      => sprintf('0, %.6g', 1.3 * $ztop),
					add_patch     => [ shade_f($g->{Df}, $e->{Df}, $fv, 8, $ORANGE) ],
					vlines        => sprintf('%.10g, 0, %.6g, color = "%s", linewidth = 2', $fv, 1.1 * $fpk, $ORANGE),
					text          => [
						sprintf('%.10g, %.6g, "F value = %.4f", fontsize = 9, ha = "left", va = "bottom", color = "%s"', $fv + 0.05, 1.1 * $fpk, $fv, $INK),
						sprintf('%.10g, %.6g, "area = Pr(>F) = %.5f", fontsize = 9, va = "bottom", color = "%s"', $fv + 0.6, 0.3 * $fpk, $g->{'Pr(>F)'}, $INK),
					],
				},
			],
		],
	);
	return;
}

# --- 2. coefficients and fitted.values -----------------------------------------
# breaks ~ wool * tension.  Under treatment contrasts the Intercept is the mean
# of the reference cell -- the first of each factor's sorted xlevels, which
# makes tension "H" the reference, not "L" -- and every other cell's fitted
# value is the Intercept plus the steps its levels add: a main effect for each
# non-reference level, and an interaction for each non-reference pair.  The
# staircase in each cell is that sum, and it lands on the fitted value, which in
# a full factorial is the cell mean.
sub fig_aov_coefficients {
	my $r  = aov(\%wb, 'breaks ~ wool * tension');
	my $b  = $r->{coefficients};
	my $fv = $r->{'fitted.values'};
	my @wl = @{ $r->{xlevels}{wool} };
	my @tl = @{ $r->{xlevels}{tension} };

	# each cell's fitted value, from the row of its first observation in %wb
	my %cell;
	for my $i (0 .. $#{ $wb{breaks} }) {
		my $k = "$wb{wool}[$i]:$wb{tension}[$i]";
		$cell{$k} = $fv->{ $i + 1 } unless exists $cell{$k};
	}

	my $INTC = $INK;    # the Intercept is a baseline, not one of the effects
	my $FITY = 50;    # the row the "fitted" labels sit in, above every staircase
	my (@patch, @link, @fit, @text, @tick, @tlab);
	for my $wi (0 .. $#wl) {
		for my $ti (0 .. $#tl) {
			my ($w, $t) = ($wl[$wi], $tl[$ti]);
			my $c = $wi * (@tl + 0.6) + $ti;
			my @step = ([ 'Intercept', $INTC ]);
			push @step, [ "tension$t", $BLUE ] if $ti;
			push @step, [ "wool$w", $GREEN ] if $wi;
			push @step, [ "wool$w:tension$t", $YELLOW ] if $wi && $ti;
			my $cum = 0;
			for my $k (0 .. $#step) {
				my ($name, $colour) = @{ $step[$k] };
				my $v  = $b->{$name};
				my $x0 = $c - 0.42 + 0.22 * $k;
				my ($lo, $hi) = $v < 0 ? ($cum + $v, $cum) : ($cum, $cum + $v);
				push @patch, rect($x0, $lo, 0.18, $hi - $lo, $colour);
				push @text, sprintf('%.6g, %.6g, "%s", fontsize = 8, ha = "center", va = "%s", color = "%s"',
					$x0 + 0.09, $v < 0 ? $lo - 0.5 : $hi + 0.5, sprintf(abs($v) < 1 ? '%+.2f' : '%+.3g', $v),
					$v < 0 ? 'top' : 'bottom', $INK) if $k;
				$cum += $v;
				push @link, [ $x0 + 0.18, $cum, $x0 + 0.22, $cum ] if $k < $#step;
			}
			push @fit, [ $c - 0.46, $cell{"$w:$t"}, $c + 0.46, $cell{"$w:$t"} ];
			push @text, sprintf('%.6g, %g, "fitted %.4g", fontsize = 9, ha = "center", va = "bottom", color = "%s"',
				$c, $FITY, $cell{"$w:$t"}, $INK);
			push @tick, $c;
			push @tlab, sprintf '"wool %s, tension %s\n%s"', $w, $t, join('\n+ ', map { $_->[0] } @step);
		}
	}
	my @key = (
		[ 'Intercept: the reference cell', $INTC ],
		[ 'tension main effect', $BLUE ],
		[ 'wool main effect', $GREEN ],
		[ 'wool:tension interaction', $YELLOW ],
	);
	my (@swatch, @ktext);
	for my $i (0 .. $#key) {
		my $y = 71 - 3.4 * $i;
		push @swatch, segments([ [ -0.5, $y, -0.38, $y ] ], $key[$i][1], 8);
		push @ktext, sprintf('-0.32, %.6g, "%s", fontsize = 9, va = "center", color = "%s"', $y, $key[$i][0], $INK);
	}

	plt(
		'output.file' => "$DIR/aov.coefficients.png",
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ -0.6, $tick[-1] + 0.6 ], [ 0, 0 ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1', $INK),
					title         => '"aov(\\%wb, \'breaks ~ wool * tension\'): coefficients add up to fitted.values", fontsize = 11',
					ylabel        => 'breaks',
					set_figwidth  => 16,
					set_figheight => 6,
					set_dpi       => 100,
					set_xlim      => sprintf('-0.6, %g', $tick[-1] + 0.6),
					set_ylim      => '0, 75',
					set_xticks    => sprintf('[%s], [%s], fontsize = 8', join(',', @tick), join(',', @tlab)),
					add_patch     => \@patch,
					text          => [ @text, @ktext,
						sprintf('%g, 73, "xlevels are sorted, so the reference levels are\nwool %s and tension %s; the dashes are fitted.values,\nwhich in a full factorial are the cell means", fontsize = 9, ha = "right", va = "top", color = "%s"',
							$tick[-1] + 0.55, $wl[0], $tl[0], $INK),
					],
				},
				segments(\@link, $INK, 1, ':'),
				segments(\@fit,  $INK, 1.6, '--'),
				@swatch,
			],
		],
	);
	return;
}

# --- 3. group.stats and the table ---------------------------------------------
# Left: group.stats, each factor's own (marginal) means, which are not the cell
# means of the figure before.  Right: the table, one row per term in R's order.
sub fig_aov_outputs {
	my $r  = aov(\%wb, 'breaks ~ wool * tension');
	my @wl = @{ $r->{xlevels}{wool} };
	my @tl = @{ $r->{xlevels}{tension} };

	my $grand = 0;
	$grand += $_ for @{ $wb{breaks} };
	$grand /= @{ $wb{breaks} };

	# group.stats: wool's levels, then tension's, along one axis, each labelled
	# on the side away from the grand-mean line
	my $gs = $r->{'group.stats'};
	my (@gx, @gy, @glab, @gtext);
	my $x = 0;
	for my $f ([ wool => \@wl ], [ tension => \@tl ]) {
		for my $l (@{ $f->[1] }) {
			my $m  = $gs->{mean}{ $f->[0] }{$l};
			my $up = $m > $grand;
			push @gx, $x;
			push @gy, $m;
			push @glab, "\"$f->[0]\\n$l\"";
			push @gtext, sprintf('%g, %.6g, "%.4g\nsize %d", fontsize = 9, ha = "center", va = "%s", color = "%s"',
				$x, $up ? $m + 1.2 : $m - 1.2, $m, $gs->{size}{ $f->[0] }{$l}, $up ? 'bottom' : 'top', $INK);
			$x++;
		}
		$x += 0.6;
	}

	# the table, top to bottom in R's term order
	my @terms = ('wool', 'tension', 'wool:tension', 'Residuals');
	my (@bars, @btext);
	for my $i (0 .. $#terms) {
		my $row = $r->{ $terms[$i] };
		my $y   = $#terms - $i;
		my $res = $terms[$i] eq 'Residuals';
		push @bars, rect(0, $y - 0.3, $row->{'Sum Sq'}, 0.6, $res ? $GREY : $BLUE);
		my $lab = $res
			? sprintf('Df %d   Mean Sq %.5g', $row->{Df}, $row->{'Mean Sq'})
			: sprintf('Df %d   Mean Sq %.5g   F value %.4g   Pr(>F) %.3g',
				$row->{Df}, $row->{'Mean Sq'}, $row->{'F value'}, $row->{'Pr(>F)'});
		push @btext, sprintf('%.6g, %g, "  Sum Sq %.5g\n  %s", fontsize = 9, va = "center", color = "%s"',
			$res ? 0 : $row->{'Sum Sq'}, $res ? $y - 0.55 : $y, $row->{'Sum Sq'}, $lab, $INK);
	}

	plt(
		'output.file' => "$DIR/aov.outputs.png",
		ncol          => 2,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ -0.5, $gx[-1] + 0.5 ], [ $grand, $grand ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1.4, linestyle = "--"', $GREY),
					title         => '"group.stats: each factor\'s own means", fontsize = 11',
					ylabel        => 'breaks',
					set_figwidth  => 16,
					set_figheight => 4.8,
					set_dpi       => 100,
					set_xlim      => sprintf('-0.6, %g', $gx[-1] + 0.6),
					set_ylim      => '15, 45',
					set_xticks    => sprintf('[%s], [%s]', join(',', @gx), join(',', @glab)),
					text          => [
						@gtext,
						sprintf('1.8, %.6g, "grand mean %.4g", fontsize = 9, ha = "center", va = "top", color = "%s"',
							$grand - 0.4, $grand, $GREY),
						sprintf('-0.5, 44, "mean and size, keyed by factor and then level,\neach averaged over the other factor:\nnot the cell means", fontsize = 9, va = "top", color = "%s"', $INK),
					],
				},
				dots(\@gx, \@gy, $BLUE, 9, 1),
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ 0, 0 ], [ -0.9, $#terms + 0.5 ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1', $INK),
					title         => '"one row per term, in R\'s order, then Residuals", fontsize = 11',
					xlabel        => '"Sum Sq"',
					set_xlim      => sprintf('0, %.6g', 1.6 * $r->{Residuals}{'Sum Sq'}),
					set_ylim      => sprintf('-1.1, %g', $#terms + 0.6),
					set_yticks    => sprintf('[%s], [%s]', join(',', reverse 0 .. $#terms), join(',', map {"\"$_\""} @terms)),
					add_patch     => \@bars,
					text          => \@btext,
				},
			],
		],
	);
	return;
}

# --- 4. anova(): the table is sequential ---------------------------------------
# Each term's Sum Sq is what it adds to the terms before it.  With numeric,
# correlated regressors (LifeCycleSavings) reversing the formula moves the
# model sum of squares between terms, though its total and Residuals do not
# move; in a balanced design (warpbreaks) the terms are orthogonal and the
# order changes nothing.
sub fig_anova_order {
	my %col = (pop15 => $BLUE, pop75 => $ORANGE, dpi => $GREEN, ddpi => $YELLOW,
		wool => $BLUE, tension => $ORANGE);

	my $panel = sub {
		my ($data, $resp, $order, $title, $first) = @_;
		my (@patch, @text, @yt, @yl);
		my $max = 0;
		for my $i (0 .. $#$order) {
			my $f = "$resp ~ " . join(' + ', @{ $order->[$i] });
			my $r = anova($data, $f);
			my $y = 1.5 * ($#$order - $i);
			my $x = 0;
			for my $k (0 .. $#{ $order->[$i] }) {
				my $t   = $order->[$i][$k];
				my $row = $r->{$t};
				push @patch, rect($x, $y - 0.3, $row->{'Sum Sq'}, 0.6, $col{$t});
				push @text, sprintf('%.6g, %g, "%s\n%.4g\np %.2g", fontsize = 8, ha = "center", va = "%s", color = "%s"',
					$x + $row->{'Sum Sq'} / 2, $k % 2 ? $y - 0.36 : $y + 0.36, $t, $row->{'Sum Sq'}, $row->{'Pr(>F)'},
					$k % 2 ? 'top' : 'bottom', $INK);
				$x += $row->{'Sum Sq'};
			}
			$max = $x if $x > $max;
			push @text, sprintf('%.6g, %g, "  Residuals %.5g\n  on Df %d", fontsize = 8, va = "center", color = "%s"',
				$x, $y, $r->{Residuals}{'Sum Sq'}, $r->{Residuals}{Df}, $GREY);
			push @yt, $y;
			push @yl, "\"$f\"";
		}
		my %p = (
			'plot.type'   => 'plot',
			data          => [ [ [ 0, 0 ], [ -0.9, 1.5 * $#$order + 0.9 ] ] ],
			'show.legend' => 0,
			'set.options' => sprintf('color = "%s", linewidth = 1', $INK),
			title         => "\"$title\", fontsize = 11",
			xlabel        => '"Sum Sq, each term after the terms to its left"',
			set_xlim      => sprintf('0, %.6g', 1.4 * $max),
			set_ylim      => sprintf('-1, %g', 1.5 * $#$order + 1),
			set_yticks    => sprintf('[%s], [%s]', join(',', @yt), join(',', @yl)),
			add_patch     => \@patch,
			text          => \@text,
		);
		%p = (%p, set_figwidth => 16, set_figheight => 4.4, set_dpi => 100) if $first;
		return [ \%p ];
	};

	plt(
		'output.file' => "$DIR/anova.order.png",
		ncol          => 2,
		p             => [
			$panel->(\%lcs, 'sr', [ [qw(pop15 pop75 dpi ddpi)], [qw(ddpi dpi pop75 pop15)] ],
				'correlated regressors: the order moves Sum Sq between terms', 1),
			$panel->(\%wb, 'breaks', [ [qw(wool tension)], [qw(tension wool)] ],
				'a balanced design: the order changes nothing', 0),
		],
	);
	return;
}

# --- 5. anova() comparing nested formulas --------------------------------------
# Left: each model's RSS, and the drop to the next, which is that row's
# Sum of Sq.  Right: F is that drop per Df over the residual mean square of
# the model with the fewest residual Df -- here the largest -- so the chain of
# comparisons reproduces the single-model table of the largest formula.
sub fig_anova_compare {
	my @f = ('breaks ~ 1', 'breaks ~ wool', 'breaks ~ wool + tension', 'breaks ~ wool * tension');
	my @short = ('breaks ~ 1', '+ wool', '+ tension', '+ wool:tension');
	my @terms = ('', 'wool', 'tension', 'wool:tension');
	my $t   = anova(\%wb, @f);
	my $one = anova(\%wb, $f[-1]);
	my $den = $t->[-1]{RSS} / $t->[-1]{'Res.Df'};

	my (@patch, @text);
	for my $i (0 .. $#$t) {
		my $row = $t->[$i];
		push @patch, rect($i - 0.35, 0, 0.7, $row->{RSS}, $GREY);
		push @text, sprintf('%d, %.6g, "RSS %.5g\nRes.Df %d", fontsize = 9, ha = "center", va = "top", color = "white"',
			$i, $row->{RSS} - 150, $row->{RSS}, $row->{'Res.Df'});
		next unless $i;
		push @patch, rect($i - 0.35, $row->{RSS}, 0.7, $row->{'Sum of Sq'}, $BLUE);
		push @text, sprintf('%g, %.6g, "Sum of Sq %.5g\nDf %d", fontsize = 9, va = "center", color = "%s"',
			$i + 0.38, $row->{RSS} + $row->{'Sum of Sq'} / 2, $row->{'Sum of Sq'}, $row->{Df}, $BLUE);
	}

	my (@mpatch, @mtext, @yt, @yl);
	for my $i (1 .. $#$t) {
		my $row = $t->[$i];
		my $y   = $#$t - $i;
		my $ms  = $row->{'Sum of Sq'} / $row->{Df};
		push @mpatch, rect(0, $y - 0.3, $ms, 0.6, $BLUE);
		push @mtext, sprintf('%.6g, %g, "  %.5g / %d = %.5g\n  F = %.5g / %.5g = %.4g,  Pr(>F) %.3g\n  anova(\\\\%%wb, \'%s\') gives %s F = %.4g", fontsize = 8.5, va = "center", color = "%s"',
			$ms, $y, $row->{'Sum of Sq'}, $row->{Df}, $ms, $ms, $den, $row->{F}, $row->{'Pr(>F)'},
			$f[-1], $terms[$i], $one->{ $terms[$i] }{'F value'}, $INK);
		push @yt, $y;
		push @yl, "\"row $i: $short[$i]\"";
	}

	plt(
		'output.file' => "$DIR/anova.compare.png",
		ncol          => 2,
		p             => [
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ -0.5, $#$t + 0.5 ], [ 0, 0 ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1', $INK),
					title         => '"anova(\\%wb, four formulas): one row per model", fontsize = 11',
					ylabel        => 'RSS',
					set_figwidth  => 16,
					set_figheight => 4.8,
					set_dpi       => 100,
					set_xlim      => sprintf('-0.5, %g', $#$t + 0.95),
					set_ylim      => sprintf('0, %.6g', 1.08 * $t->[0]{RSS}),
					set_xticks    => sprintf('[%s], [%s]', join(',', 0 .. $#$t), join(',', map {"\"row $_\\n$short[$_]\""} 0 .. $#$t)),
					add_patch     => \@patch,
					text          => [ @text,
						sprintf('%g, %.6g, "each row after the first adds\nSum of Sq = the drop in RSS\nDf = the drop in Res.Df", fontsize = 9, ha = "right", va = "top", color = "%s"',
							$#$t + 0.9, 1.06 * $t->[0]{RSS}, $INK),
					],
				},
			],
			[
				{
					'plot.type'   => 'plot',
					data          => [ [ [ $den, $den ], [ -0.6, $#$t - 0.4 ] ] ],
					'show.legend' => 0,
					'set.options' => sprintf('color = "%s", linewidth = 1.6, linestyle = "--"', $GREY),
					title         => '"F: each drop per Df, over the largest model\'s RSS / Res.Df", fontsize = 11',
					xlabel        => '"Sum of Sq / Df"',
					set_xlim      => sprintf('0, %.6g', 2.6 * $t->[2]{'Sum of Sq'} / $t->[2]{Df}),
					set_ylim      => sprintf('-0.6, %g', $#$t - 0.3),
					set_yticks    => sprintf('[%s], [%s]', join(',', @yt), join(',', @yl)),
					add_patch     => \@mpatch,
					text          => [ @mtext,
						sprintf('%.6g, %g, " %.5g / %d = %.5g", fontsize = 9, va = "bottom", color = "%s"',
							$den, $#$t - 0.55, $t->[-1]{RSS}, $t->[-1]{'Res.Df'}, $den, $GREY),
					],
				},
			],
		],
	);
	return;
}

fig_aov_what();
fig_aov_coefficients();
fig_aov_outputs();
fig_anova_order();
fig_anova_compare();
