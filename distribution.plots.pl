#!/usr/bin/env perl
# Draws the figures that illustrate the distribution family in README.md (and
# therefore in read.me.pod and lib/Stats/LikeR.pm).  Author-only: it is not
# shipped, and it needs Matplotlib::Simple, python3 and matplotlib.
#
#   perl -Iblib/lib -Iblib/arch distribution.plots.pl
#
# One figure per function, ten in all: dnorm and pnorm, and the eight added in
# 0.303.  Each is written to img/<fn>.what.png at the size the README and the
# POD embed it at.
#
# The point every figure has to make is which end of the same picture the
# function is answering from, because that is the whole difference between a
# p-function and a q-function:
#
#   dnorm   a HEIGHT at a point -- no area anywhere, which is why it is the
#           only figure here with no shading
#   p*      given the BOUNDARY, return the shaded AREA        (the integral)
#   q*      given the AREA, return the dashed BOUNDARY        (its inverse)
#
# So a p* figure puts its answer on the shading and a q* figure puts its answer
# on the vertical line, and both carry the same integral in mathtext so the pair
# reads as one identity seen from two sides. pbinom is the exception: a binomial
# is discrete, so its figure is bars and its annotation is a sum, not an
# integral -- writing the integral there would be a lie about what pbinom does.
#
# Every number drawn comes back out of the module itself -- the densities below
# are only for the curve, never for a quoted value -- so the pictures and the
# code cannot drift apart.
#
# Colours and sizes follow t.test.plots.pl and density.plots.pl: identity gets
# the categorical slots in their documented order, the answer gets the second
# slot so the eye lands on it, and recessive marks get the grey.
require 5.010;
use strict;
use warnings FATAL => 'all';
use File::Path 'make_path';
use POSIX 'lgamma';
use Stats::LikeR qw(dnorm pnorm qnorm pt qt pchisq qchisq pf qf pbinom);
use Matplotlib::Simple 'plt';

my $DIR = 'img';
make_path($DIR) unless -d $DIR;

# --- palette ----------------------------------------------------------------
my $BLUE   = '#2a78d6';    # categorical slot 1: the distribution itself
my $ORANGE = '#eb6834';    # categorical slot 2: the answer
my $GREY   = '#8c8c88';    # recessive reference marks
my $INK    = '#3d3d3a';    # annotation

my $PI = 4 * atan2(1, 1);

# --- densities, for the curves only -----------------------------------------
# None of these is used for a number the figure quotes; every quoted value
# comes from the module.  They exist because the module exposes no d* function
# but dnorm, and a CDF picture needs the density to shade.
sub d_t {
	my ($x, $df) = @_;
	return exp( lgamma( ( $df + 1 ) / 2 ) - lgamma( $df / 2 )
	          - 0.5 * log( $df * $PI )
	          - ( ( $df + 1 ) / 2 ) * log( 1 + $x * $x / $df ) );
}

sub d_chisq {
	my ($x, $df) = @_;
	return 0 if $x <= 0;
	return exp( ( $df / 2 - 1 ) * log($x) - $x / 2
	          - ( $df / 2 ) * log(2) - lgamma( $df / 2 ) );
}

sub d_f {
	my ($x, $d1, $d2) = @_;
	return 0 if $x <= 0;
	my $lbeta = lgamma( $d1 / 2 ) + lgamma( $d2 / 2 ) - lgamma( ( $d1 + $d2 ) / 2 );
	return exp( ( $d1 / 2 ) * log( $d1 / $d2 )
	          + ( $d1 / 2 - 1 ) * log($x)
	          - ( ( $d1 + $d2 ) / 2 ) * log( 1 + $d1 * $x / $d2 )
	          - $lbeta );
}

sub d_binom {
	my ($k, $n, $p) = @_;
	return exp( lgamma( $n + 1 ) - lgamma( $k + 1 ) - lgamma( $n - $k + 1 )
	          + $k * log($p) + ( $n - $k ) * log( 1 - $p ) );
}

# --- plumbing ---------------------------------------------------------------
# curve($pdf, $lo, $hi) -- [x], [y] over 400 points, for 'plot'.
sub curve {
	my ($pdf, $lo, $hi, $n) = @_;
	$n = 400 unless defined $n;
	my (@x, @y);
	for my $i (0 .. $n) {
		my $x = $lo + ( $hi - $lo ) * $i / $n;
		push @x, sprintf '%.8g', $x;
		push @y, sprintf '%.8g', $pdf->($x);
	}
	return ( \@x, \@y );
}

# shade($pdf, $lo, $hi) -- a filled polygon under $pdf between $lo and $hi,
# as the plt.Polygon string that add_patch takes.  Closed down to y = 0 at both
# ends so the fill is the area under the curve and not the curve itself.
sub shade {
	my ($pdf, $lo, $hi, $alpha) = @_;
	$alpha = 0.42 unless defined $alpha;
	my @pt = ( sprintf '[%.8g,0]', $lo );
	for my $i (0 .. 300) {
		my $x = $lo + ( $hi - $lo ) * $i / 300;
		push @pt, sprintf '[%.8g,%.8g]', $x, $pdf->($x);
	}
	push @pt, sprintf '[%.8g,0]', $hi;
	return sprintf 'plt.Polygon([%s], closed = True, facecolor = "%s", '
	             . 'edgecolor = "none", alpha = %s)',
		join( ',', @pt ), $ORANGE, $alpha;
}

# one_panel(%spec) -- write a single-panel figure. Keeps the size, dpi and
# spine choices identical across all ten so they sit together in the README.
sub one_panel {
	my (%s) = @_;
	plt(
		'output.file' => "$DIR/$s{file}",
		'plot.type'   => 'plot',
		data          => $s{data},
		'show.legend' => 0,
		'set.options' => 'color = "' . $BLUE . '", linewidth = 2.2',
		title         => $s{title},
		xlabel        => $s{xlabel},
		ylabel        => '"density"',
		set_figwidth  => 7.6,
		set_figheight => 3.5,
		set_dpi       => 110,
		set_xlim      => $s{xlim},
		set_ylim      => $s{ylim},
		( $s{patch} ? ( add_patch => $s{patch} ) : () ),
		( $s{vlines} ? ( vlines => $s{vlines} ) : () ),
		( $s{hlines} ? ( hlines => $s{hlines} ) : () ),
		( $s{bar}    ? ( bar    => $s{bar} )    : () ),
		( $s{annotate} ? ( annotate => $s{annotate} ) : () ),
		text          => $s{text},
	);
	return;
}


# --- annotation discipline --------------------------------------------------
# Every figure places its text the same way, because the alternative is ten
# separate collisions with the curve:
#
#   top-left, line 1     the identity, in mathtext
#   top-left, line 2     the concrete call, in the colour of whichever side is
#                        the answer
#   below the axis       the boundary, at its own x -- grey when it is given,
#                        orange when it is the answer. Text is not clipped to
#                        the axes, so hanging it under y = 0 can never overlap
#                        anything drawn.
#   one aside            in ink, in whichever corner the density leaves empty
#
# $peak scales the vertical positions so the same numbers work for a density
# whose maximum is 0.4 and one whose maximum is 0.24.
sub ann {
	my (%a) = @_;
	my @t;
	push @t, sprintf '%.8g, %.8g, %s, fontsize = 13, va = "top", ha = "%s", color = "%s"',
		$a{ix}, ( $a{identy} || 1.42 ) * $a{peak}, $a{identity}, ( $a{iha} || 'left' ), $INK;
	# cally lets a panel whose identity is TALL -- pbinom's summation sign,
	# whose limits sit above and below the sigma -- push the second line down
	# without moving the other nine.
	push @t, sprintf '%.8g, %.8g, "%s", fontsize = 10, va = "top", ha = "%s", color = "%s"',
		$a{ix}, ( $a{cally} || 1.22 ) * $a{peak}, $a{call}, ( $a{iha} || 'left' ),
		$a{callcolour};
	push @t, sprintf '%.8g, %.8g, "%s", fontsize = 9.5, ha = "center", va = "top", color = "%s"',
		$a{bx}, -0.055 * $a{peak}, $a{blabel}, $a{bcolour};
	push @t, sprintf '%.8g, %.8g, "%s", fontsize = 9, va = "top", ha = "%s", color = "%s"',
		$a{ax}, $a{ay}, $a{aside}, ( $a{aha} || 'left' ), $INK
		if defined $a{aside};
	return \@t;
}

# The integral, in mathtext, with the limits and the answer filled in. One
# helper so the ten figures cannot disagree about how the identity is written.
sub integral {
	my ($lo, $hi, $var, $value) = @_;
	return sprintf 'r"$\int_{%s}^{%s}\!f(%s)\,d%s = %s$"', $lo, $hi, $var, $var, $value;
}

# ---------------------------------------------------------------------------
# dnorm -- a height, and the only figure here with nothing shaded.
# ---------------------------------------------------------------------------
{
	my $x0   = 1;
	my $h    = dnorm($x0);
	my $peak = dnorm(0);
	my ($cx, $cy) = curve( sub { dnorm( $_[0] ) }, -4, 4 );
	one_panel(
		file   => 'dnorm.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => '"dnorm: the HEIGHT of the normal density at a point", fontsize = 11',
		xlabel => 'x',
		xlim   => '-4, 4',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 2.6', $x0, $h, $ORANGE ),
		],
		hlines => [
			sprintf( '%.8g, -4, %.8g, color = "%s", linewidth = 1.2, linestyle = ":"',
				$h, $x0, $ORANGE ),
		],
		text => ann(
			peak       => $peak,
			ix         => -3.9,
			identity   => 'r"$f(x) = \frac{1}{\sqrt{2\pi\sigma^2}}\,e^{-(x-\mu)^2/2\sigma^2}$"',
			call       => sprintf( 'dnorm(%g) = %.6f', $x0, $h ),
			callcolour => $ORANGE,
			bx         => $x0,
			blabel     => sprintf( 'x = %g', $x0 ),
			bcolour    => $ORANGE,
			ax         => 3.9,
			ay         => 0.72 * $peak,
			aha        => 'right',
			aside      => 'a height, not an area:\nnothing here is integrated.\nEvery other figure below shades one.',
		),
	);
}

# ---------------------------------------------------------------------------
# pnorm / qnorm -- the same picture, answered from opposite ends.
# ---------------------------------------------------------------------------
{
	my $pdf  = sub { dnorm( $_[0] ) };
	my $peak = dnorm(0);
	my ($cx, $cy) = curve( $pdf, -4, 4 );

	my $q = 1.28;
	my $p = pnorm($q);
	one_panel(
		file   => 'pnorm.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => '"pnorm: given the boundary, return the AREA to its left", fontsize = 11',
		xlabel => 'x',
		xlim   => '-4, 4',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, -4, $q ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 1.8, linestyle = "--"',
				$q, $pdf->($q), $GREY ),
		],
		text => ann(
			peak       => $peak,
			ix         => -3.9,
			identity   => integral( '-\infty', sprintf('%g', $q), 'x', sprintf('%.5f', $p) ),
			call       => sprintf( 'pnorm(%g) = %.5f      <- the answer is the shaded area', $q, $p ),
			callcolour => $ORANGE,
			bx         => $q,
			blabel     => sprintf( 'q = %g (given)', $q ),
			bcolour    => $GREY,
			ax         => 3.9,
			ay         => 0.66 * $peak,
			aha        => 'right',
			aside      => 'lower => 0 shades the other\nside instead, and returns\n1 minus this -- computed as\nits own integral, not by\nsubtracting.',
		),
	);

	my $pp = 0.9;
	my $qq = qnorm($pp);
	one_panel(
		file   => 'qnorm.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => '"qnorm: given the area, return the BOUNDARY -- pnorm run backwards", fontsize = 11',
		xlabel => 'x',
		xlim   => '-4, 4',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, -4, $qq ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 3', $qq, 0.62 * $peak, $ORANGE ),
		],
		text => ann(
			peak       => $peak,
			ix         => -3.9,
			identity   => sprintf( 'r"$\int_{-\infty}^{q}\!f(x)\,dx = %.2f \Rightarrow q = ?$"', $pp ),
			call       => sprintf( 'the area %.2f is given; qnorm(%.2f) = %.6f is the answer', $pp, $pp, $qq ),
			callcolour => $ORANGE,
			bx         => $qq,
			blabel     => sprintf( 'q = %.6f\n(the answer)', $qq ),
			bcolour    => $ORANGE,
			ax         => -1.3,
			ay         => 0.40 * $peak,
			aha        => 'center',
			aside      => sprintf( 'shaded area = %.2f\n(given)', $pp ),
		),
	);
}

# ---------------------------------------------------------------------------
# pt / qt
# ---------------------------------------------------------------------------
{
	my $df   = 10;
	my $pdf  = sub { d_t( $_[0], $df ) };
	my $peak = d_t( 0, $df );
	my ($cx, $cy) = curve( $pdf, -5, 5 );

	my $q = 2.5;
	my $p = pt( $q, $df );
	one_panel(
		file   => 'pt.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => sprintf( '"pt: the AREA of the t density left of t, here df = %g", fontsize = 11', $df ),
		xlabel => 't',
		xlim   => '-5, 5',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, -5, $q ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 1.8, linestyle = "--"',
				$q, $pdf->($q), $GREY ),
		],
		text => ann(
			peak       => $peak,
			ix         => -4.9,
			identity   => integral( '-\infty', sprintf('%g', $q), 't', sprintf('%.5f', $p) ),
			call       => sprintf( 'pt(%g, %g) = %.5f      <- the answer is the shaded area', $q, $df, $p ),
			callcolour => $ORANGE,
			bx         => $q,
			blabel     => sprintf( 't = %g (given)', $q ),
			bcolour    => $GREY,
			ax         => 4.9,
			ay         => 0.66 * $peak,
			aha        => 'right',
			aside      => sprintf( 'the upper tail is the\nunshaded part:\npt(%g, %g, lower => 0)\n= %.5f', $q, $df, pt( $q, $df, lower => 0 ) ),
		),
	);

	my $pp = 0.975;
	my $qq = qt( $pp, $df );
	one_panel(
		file   => 'qt.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => sprintf( '"qt: the t that leaves a given area to its left, df = %g", fontsize = 11', $df ),
		xlabel => 't',
		xlim   => '-5, 5',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, -5, $qq ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 3', $qq, 0.62 * $peak, $ORANGE ),
		],
		text => ann(
			peak       => $peak,
			ix         => -4.9,
			identity   => sprintf( 'r"$\int_{-\infty}^{t}\!f(t)\,dt = %.3f \Rightarrow t = ?$"', $pp ),
			call       => sprintf( 'qt(%.3f, %g) = %.6f is the answer', $pp, $df, $qq ),
			callcolour => $ORANGE,
			bx         => $qq,
			blabel     => sprintf( 't = %.6f\n(the answer)', $qq ),
			bcolour    => $ORANGE,
			ax         => -1.6,
			ay         => 0.40 * $peak,
			aha        => 'center',
			aside      => sprintf( 'shaded area = %.3f\n(given)\n\nthis is the 1.96 that a\nt interval is built from', $pp ),
		),
	);
}

# ---------------------------------------------------------------------------
# pchisq / qchisq -- the right tail, because that is the one a test uses.
# ---------------------------------------------------------------------------
{
	my $df   = 3;
	my $pdf  = sub { d_chisq( $_[0], $df ) };
	my $peak = d_chisq( 1.0, $df );
	my ($cx, $cy) = curve( $pdf, 1e-6, 14 );

	my $q  = 7.81;
	my $up = pchisq( $q, $df, lower => 0 );
	one_panel(
		file   => 'pchisq.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => sprintf( '"pchisq: the chi-square tail AREA, here the upper one, df = %g", fontsize = 11', $df ),
		xlabel => '"chi-square statistic"',
		xlim   => '0, 14',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, $q, 14 ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 1.8, linestyle = "--"',
				$q, $pdf->($q), $GREY ),
		],
		text => ann(
			peak       => $peak,
			ix         => 3.6,
			identity   => integral( sprintf('%g', $q), '\infty', 'x', sprintf('%.5f', $up) ),
			call       => sprintf( 'pchisq(%g, %g, lower => 0) = %.5f      <- the answer', $q, $df, $up ),
			callcolour => $ORANGE,
			bx         => $q,
			blabel     => sprintf( 'statistic = %g (given)', $q ),
			bcolour    => $GREY,
			ax         => 13.8,
			ay         => 0.62 * $peak,
			aha        => 'right',
			aside      => 'this shaded tail is the\np-value that every\nchi-square test reports:\nthe mass beyond the\nstatistic it observed.',
		),
	);

	my $pp = 0.95;
	my $qq = qchisq( $pp, $df );
	one_panel(
		file   => 'qchisq.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => sprintf( '"qchisq: the critical value that cuts off a given tail, df = %g", fontsize = 11', $df ),
		xlabel => '"chi-square statistic"',
		xlim   => '0, 14',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, $qq, 14 ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 3', $qq, 0.62 * $peak, $ORANGE ),
		],
		text => ann(
			peak       => $peak,
			ix         => 3.6,
			identity   => sprintf( 'r"$\int_{q}^{\infty}\!f(x)\,dx = %.2f \Rightarrow q = ?$"', 1 - $pp ),
			call       => sprintf( 'qchisq(%.2f, %g) = %.6f is the answer', $pp, $df, $qq ),
			callcolour => $ORANGE,
			bx         => $qq,
			blabel     => sprintf( 'q = %.6f\n(the answer)', $qq ),
			bcolour    => $ORANGE,
			ax         => 13.8,
			ay         => 0.62 * $peak,
			aha        => 'right',
			aside      => sprintf( 'shaded area = %.2f (given)\n\nthe 3.84 of a 1-df test\nis qchisq(0.95, 1)', 1 - $pp ),
		),
	);
}

# ---------------------------------------------------------------------------
# pf / qf
# ---------------------------------------------------------------------------
{
	my ($d1, $d2) = ( 2, 10 );
	my $pdf  = sub { d_f( $_[0], $d1, $d2 ) };
	my $peak = d_f( 0.2, $d1, $d2 );
	my ($cx, $cy) = curve( $pdf, 1e-6, 8 );

	my $q  = 4.1;
	my $up = pf( $q, $d1, $d2, lower => 0 );
	one_panel(
		file   => 'pf.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => sprintf( '"pf: the F tail AREA beyond an observed F, df = %g and %g", fontsize = 11', $d1, $d2 ),
		xlabel => '"F statistic"',
		xlim   => '0, 8',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, $q, 8 ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 1.8, linestyle = "--"',
				$q, $pdf->($q), $GREY ),
		],
		text => ann(
			peak       => $peak,
			ix         => 2.0,
			identity   => integral( sprintf('%g', $q), '\infty', 'F', sprintf('%.5f', $up) ),
			call       => sprintf( 'pf(%g, %g, %g, lower => 0) = %.5f      <- the answer', $q, $d1, $d2, $up ),
			callcolour => $ORANGE,
			bx         => $q,
			blabel     => sprintf( 'F = %g (given)', $q ),
			bcolour    => $GREY,
			ax         => 7.9,
			ay         => 0.62 * $peak,
			aha        => 'right',
			aside      => 'this shaded tail is the\nPr(>F) column of an\nANOVA table.',
		),
	);

	my $pp = 0.95;
	my $qq = qf( $pp, $d1, $d2 );
	one_panel(
		file   => 'qf.what.png',
		data   => [ [ $cx, $cy ] ],
		title  => sprintf( '"qf: the F that cuts off a given upper tail, df = %g and %g", fontsize = 11', $d1, $d2 ),
		xlabel => '"F statistic"',
		xlim   => '0, 8',
		ylim   => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.50 * $peak ),
		patch  => [ shade( $pdf, $qq, 8 ) ],
		vlines => [
			sprintf( '%.8g, 0, %.8g, color = "%s", linewidth = 3', $qq, 0.62 * $peak, $ORANGE ),
		],
		text => ann(
			peak       => $peak,
			ix         => 2.0,
			identity   => sprintf( 'r"$\int_{q}^{\infty}\!f(F)\,dF = %.2f \Rightarrow q = ?$"', 1 - $pp ),
			call       => sprintf( 'qf(%.2f, %g, %g) = %.6f is the answer', $pp, $d1, $d2, $qq ),
			callcolour => $ORANGE,
			bx         => $qq,
			blabel     => sprintf( 'q = %.6f\n(the answer)', $qq ),
			bcolour    => $ORANGE,
			ax         => 7.9,
			ay         => 0.62 * $peak,
			aha        => 'right',
			aside      => sprintf( 'shaded area = %.2f (given)\n\nthe number an F table\nused to be printed for.', 1 - $pp ),
		),
	);
}

# ---------------------------------------------------------------------------
# pbinom -- discrete, so bars and a SUM. The one figure here that must not
# carry an integral sign, because there is nothing to integrate: the answer is
# a finite sum of the bar heights that are shaded.
#
# The bars are Rectangle patches rather than a bar plot, for the same reason the
# shading elsewhere is a Polygon: add_patch takes whatever geometry the figure
# needs, and two colours of bar in one panel is geometry.
# ---------------------------------------------------------------------------
{
	my ($n, $prob, $k) = ( 10, 0.5, 3 );
	my $lower = pbinom( $k, $n, $prob );
	my @kk    = ( 0 .. $n );
	my $peak  = d_binom( $n / 2, $n, $prob );
	my $W     = 0.72;

	my @bars = map {
		sprintf 'plt.Rectangle((%.8g, 0), %.8g, %.8g, facecolor = "%s", '
		      . 'edgecolor = "none", alpha = %s)',
			$_ - $W / 2, $W, d_binom( $_, $n, $prob ),
			( $_ <= $k ? $ORANGE : $GREY ), ( $_ <= $k ? 0.85 : 0.55 )
	} @kk;

	plt(
		'output.file'  => "$DIR/pbinom.what.png",
		'plot.type'    => 'plot',
		# An invisible series: the panel needs a 'plot' to exist at all, and
		# every mark that matters is a patch.
		data           => [ [ [ 0, $n ], [ 0, 0 ] ] ],
		'show.legend'  => 0,
		'set.options'  => 'color = "none", linewidth = 0',
		title          => sprintf( '"pbinom: a SUM of probabilities, not an integral (size = %g, prob = %g)", fontsize = 11', $n, $prob ),
		xlabel         => '"number of successes"',
		ylabel         => '"probability"',
		set_figwidth   => 7.6,
		set_figheight  => 3.5,
		set_dpi        => 110,
		set_xlim       => '-0.8, 10.8',
		set_ylim       => sprintf( '%.6g, %.6g', -0.15 * $peak, 1.88 * $peak ),
		set_xticks     => '[' . join( ',', @kk ) . ']',
		add_patch      => \@bars,
		text           => ann(
			peak       => $peak,
			ix         => 5.4,
			identity   => sprintf( 'r"$\sum_{i=0}^{%g}\binom{%g}{i}p^{i}(1-p)^{%g-i} = %.6f$"', $k, $n, $n, $lower ),
			call       => sprintf( 'pbinom(%g, %g, %g) = %.6f      <- the answer', $k, $n, $prob, $lower ),
			callcolour => $ORANGE,
			iha        => 'center',
			# The bars reach the full peak, so pbinom needs more headroom than
			# the nine density panels: a summation sign is two lines tall on its
			# own and the call has to clear the tallest bar underneath it.
			identy     => 1.72,
			cally      => 1.34,
			bx         => $k,
			blabel     => sprintf( 'k = %g (given)', $k ),
			bcolour    => $ORANGE,
			ax         => 10.4,
			ay         => 0.70 * $peak,
			aha        => 'right',
			aside      => 'no integral sign here:\na binomial is discrete, so\nthe lower tail is a finite\nsum of the shaded bars.\nlower => 0 sums the grey ones.',
		),
	);
}

print "\nall ten figures written to $DIR/\n";
