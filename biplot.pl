#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use utf8;
use Stats::LikeR;
use Matplotlib::Simple 'scatter';

=head2 biplot

An R-style C<biplot.prcomp> for the HashRef returned by C<Stats::LikeR::prcomp>,
drawn through C<Matplotlib::Simple>.

 biplot($prcomp, 'output.file' => 'svg/biplot.png', top => 20);

Options: C<choices> ([1,2] PCs to plot), C<top> (only arrow the N variables with
the longest loading in that plane; 0 = all), C<biplot.scale> (R's C<scale>, 1 by
default), C<trim> (fraction of the most extreme scores that don't get to set the
frame, 0.01 by default; C<< trim =E<gt> 0 >> frames every point), C<labels>
(ArrayRef of row names to print instead of points), plus any Matplotlib::Simple
option, which is passed through to C<scatter> -- so C<scale> still means what it
means everywhere else in Matplotlib::Simple, the figure size.

=cut

sub extent {    # symmetric extent of @v, ignoring the most extreme $trim of it
	my ($trim, @v) = @_;
	my @a = sort { $a <=> $b } map { abs } @v;
	return ($a[ int((1 - $trim) * $#a + 0.5) ] || $a[-1] || 1);
}

sub biplot {
	my ($pca, %opt) = @_;
	foreach my $k ('x', 'rotation', 'sdev') {
		die "biplot: \$prcomp has no \"$k\"; call prcomp with retx => true" unless defined $pca->{$k};
	}
	my $ch    = delete $opt{choices} // [1, 2];
	my $top   = delete $opt{top}     // 0;
	my $sc    = delete $opt{'biplot.scale'} // 1;
	my $lab   = delete $opt{labels};
	my $trim  = delete $opt{trim}    // 0.01;
	my $file  = delete $opt{'output.file'} // 'biplot.png';
	my ($i, $j) = ($ch->[0] - 1, $ch->[1] - 1);
	my @sdev  = @{ $pca->{sdev} };
	my $names = $pca->{varnames} // [map { "V$_" } 1 .. scalar @{ $pca->{rotation} }];

	# --- R's biplot.prcomp: lam = sdev * sqrt(n), scores /= lam^scale, loadings *= lam^scale
	my $n   = scalar @{ $pca->{x} };
	my @lam = map { ($sdev[$_] * sqrt($n)) ** $sc } ($i, $j);
	$_ ||= 1 foreach @lam;    # a zero-variance component must not divide
	my @px  = map { $_->[$i] / $lam[0] } @{ $pca->{x} };    # scores, PC$ch->[0]
	my @py  = map { $_->[$j] / $lam[1] } @{ $pca->{x} };    # scores, PC$ch->[1]
	my @lx  = map { $_->[$i] * $lam[0] } @{ $pca->{rotation} };    # loadings
	my @ly  = map { $_->[$j] * $lam[1] } @{ $pca->{rotation} };

	# --- keep the N longest loadings only, so 100+ arrows aren't a hairball
	my @v = 0 .. $#lx;
	if (($top > 0) && ($top < scalar @v)) {
		my @by_len = sort { ($lx[$b]**2 + $ly[$b]**2) <=> ($lx[$a]**2 + $ly[$a]**2) } @v;
		@v = sort { $a <=> $b } @by_len[0 .. $top - 1];
	}

	# --- the extent of the scores on each axis, taken past the most extreme $trim of
	#     them: a single outlying complex must not decide the scale of everything else
	my ($sx, $sy) = (extent($trim, @px), extent($trim, @py));

	# --- R draws the arrows on their own axes: ratio maps loading units to score units
	my $ratio = max(
		(max map { abs $lx[$_] } @v) / $sx,
		(max map { abs $ly[$_] } @v) / $sy
	) || 1;
	my $var = sum map { $_ ** 2 } @sdev;

	# --- every value below is emitted verbatim into the python, so it must be a
	#     python expression: quote the strings, and let matplotlib see the commas
	my (@arrow, @lbl);
	foreach my $k (@v) {
		my ($dx, $dy) = ($lx[$k] / $ratio, $ly[$k] / $ratio);
		# FancyArrowPatch, not FancyArrow: its head is sized in points, so the heads
		# stay square once the two axes no longer share a scale
		push @arrow, sprintf
		 'plt.matplotlib.patches.FancyArrowPatch((0, 0), (%.6g, %.6g), arrowstyle = "-|>", '
		. 'mutation_scale = 8, lw = 0.8, shrinkA = 0, shrinkB = 0, color = "firebrick", '
		. 'alpha = 0.7, zorder = 3)',
		 $dx, $dy;
		my $name = $names->[$k];
		$name = Encode::decode_utf8($name) unless Encode::is_utf8($name); # ΔG, not ÎG
		push @lbl, {
			x    => 1.04 * $dx,
			y    => 1.04 * $dy,
			ha   => $dx < 0 ? 'right' : 'left',
			name => $name =~ s/(["\\])/\\$1/gr,
		};
	}
	# --- near-collinear loadings would print their names on top of each other, so
	#     nudge each side's labels apart vertically, lowest first
	my $gap = 0.03 * $sy;
	foreach my $side ('left', 'right') {
		my $prev;
		foreach my $l (sort { $a->{y} <=> $b->{y} } grep { $_->{ha} eq $side } @lbl) {
			$l->{y} = $prev + $gap if ((defined $prev) && ($l->{y} - $prev < $gap));
			$prev = $l->{y};
		}
	}

	# --- each axis gets its own limits: one shared square frame (R's default) wastes
	#     most of the plot as soon as the two PCs differ in spread.  Each axis has to
	#     hold its scores, its arrow tips, and the names hanging off those tips.
	my $limx = 1.06 * max($sx, map { abs $_->{x} } @lbl);
	my $limy = 1.08 * max($sy, map { abs $_->{y} } @lbl);
	# a 6pt character is roughly 0.6% of the half-width, so allow that much per
	# character -- but a long name may not push the frame out by more than half
	$limx = min(1.5 * $limx,
	 max($limx, map { abs($_->{x}) + 0.006 * $limx * length $_->{name} } @lbl));
	my $out = grep { (abs($px[$_]) > $limx) || (abs($py[$_]) > $limy) } 0 .. $n - 1;
	say "biplot: $out of $n points lie outside the frame; trim => 0 fits them in" if $out;
	my @text = map { sprintf
	 '%.6g, %.6g, "%s", color = "firebrick", fontsize = 6, zorder = 4, ha = "%s", va = "center"',
	 $_->{x}, $_->{y}, $_->{name}, $_->{ha} } @lbl;
	if (defined $lab) {
		die 'biplot: labels must hold one name per row of $prcomp->{x}' if scalar @{$lab} != $n;
		push @text, sprintf('%.6g, %.6g, "%s", fontsize = 5, ha = "center", va = "center"',
		 $px[$_], $py[$_], $lab->[$_] =~ s/(["\\])/\\$1/gr) foreach 0 .. $n - 1;
	}

	scatter(
		data => {
			X => \@px,
			Y => \@py,
		},
		keys => ['X', 'Y'],           # X is the horizontal axis, Y the vertical one
		'set.options' => defined $lab
		 ? 's = 0, color = "none"'    # labels replace the points, as in R
		 : 's = 9, color = "black", alpha = 0.55, linewidths = 0',
		add_patch    => \@arrow,
		text         => \@text,
		axhline      => '0, color = "gray", lw = 0.5, ls = "--", zorder = 0',
		axvline      => '0, color = "gray", lw = 0.5, ls = "--", zorder = 0',
		set_xlim     => sprintf('%.6g, %.6g', -$limx, $limx),
		set_ylim     => sprintf('%.6g, %.6g', -$limy, $limy),
		# xlabel/ylabel, not set_xlabel/set_ylabel: scatter defaults those two from
		# the data keys, and its default would be emitted after (and so overwrite) mine
		xlabel       => sprintf('PC%d (%.1f%% of variance)', $ch->[0], 100 * $sdev[$i]**2 / $var),
		ylabel       => sprintf('PC%d (%.1f%% of variance)', $ch->[1], 100 * $sdev[$j]**2 / $var),
		# the loading axes, in loading units, on top/right -- R's second frame
		secondary_xaxis => sprintf('"top", functions = (lambda v: v * %.9g, lambda v: v / %.9g)', $ratio, $ratio),
		secondary_yaxis => sprintf('"right", functions = (lambda v: v * %.9g, lambda v: v / %.9g)', $ratio, $ratio),
		tick_params  => 'labelsize = 7',
		'output.file' => $file,
		%opt
	);
}

my $feat220 = read_table('all.feat.tsv');
my @col     = colnames($feat220);
$feat220 = drop_cols($feat220, grep {/_rank$/} @col);
$feat220 = drop_cols($feat220, grep {/rank err$/} colnames($feat220));
$feat220 = drop_cols($feat220, grep {/(?:25|50|75)/} colnames($feat220));
$feat220 = drop_cols($feat220, grep {/(?:bias|mmgbsa|fram)/} colnames($feat220));
my $prcomp  = prcomp($feat220, scale => true);
my $ncol = ncol($feat220);
my $nrow = nrow($feat220);
biplot($prcomp,
	top           => 20,
	scale         => 2.0,
	set_title     => "'PCA biplot: 220 complexes, $ncol features'", # quote it: a comma
	                                                              # in a title isn't auto-quoted
	'output.file' => 'svg/biplot.svg',
);
