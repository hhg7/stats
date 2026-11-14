#!/usr/bin/env perl
# cpantesters.matrix.pl - draw a CPAN Testers perl x OS matrix as a clickable SVG.
#
# matrix.cpantesters.org and the www.cpantesters.org distribution pages answer
# 403 (or 404) as of 2026-11-09, but the JSON API behind them still works, and
# so does each report's own page. This script rebuilds the old matrix from the
# API with Matplotlib::Simple's colored_table: one row per perl, one column per
# OS, each cell colored by the share of its graded reports that pass and
# linked to one of its reports.
#
#   ./cpantesters.matrix.pl                       # newest release of dist.ini's dist
#   ./cpantesters.matrix.pl --version 0.3417
#   ./cpantesters.matrix.pl --dist Foo-Bar --columns platform -o foo.svg
#   ./cpantesters.matrix.pl --json saved.summary.json --version 0.3418
#
# The dist is the top-level "name" in the dist.ini beside this script, unless
# --dist names another.
#
# Open the SVG in a browser directly: links in an SVG do nothing when the file
# is shown through an <img> tag or in most image viewers.
#
# Sources, all checked by hand on 2026-11-09:
# - https://api.cpantesters.org/v3/summary/<dist>[/<version>] returns a JSON
#   array of reports, each with guid, grade (pass, fail, na, unknown), perl
#   (such as "5.42.3"), osname (such as "MSWin32"), platform (such as
#   "MSWin32-x64-multi-thread"), version and date.
# - https://www.cpantesters.org/cpan/report/<guid> is the report page each
#   cell links to; it answered 200 where the matrix answered 403.
# - The links rely on matplotlib's SVG backend writing an artist's set_url()
#   as an <a xlink:href>. colored_table draws its cells with Axes.table(),
#   whose cells and their texts are artists, so the Python prelude below wraps
#   Axes.table() to set each cell's URL and text before the figure is saved.
#   Matplotlib::Simple has no option for per-cell URLs, but it appends its
#   script to a File::Temp handle passed as "fh", so a prelude written there
#   first runs first.
#
# Author tool: not part of the distribution ([PruneFiles] in dist.ini drops
# every root-level script, and Makefile.PL's libscan keeps it out of
# `make install`), so it may use modern perl.

use 5.036;
use warnings FATAL => 'all';
use Getopt::Long 'GetOptions';
use HTTP::Tiny;
use JSON::PP ();
use MIME::Base64 'encode_base64';
use File::Temp ();
use File::Spec ();
use FindBin ();
use version ();
use Matplotlib::Simple 'colored_table';

my %opt = (
	columns => 'osname',	# osname = one column per OS; platform = one per archname
);
GetOptions(\%opt, 'dist=s', 'version=s', 'output|o=s', 'json=s', 'columns=s', 'help|h')
	or die "usage: $0 [--dist D] [--version V] [--output F.svg] [--json FILE] [--columns osname|platform]\n";
if ($opt{help}) {
	print "usage: $0 [--dist D] [--version V] [--output F.svg] [--json FILE] [--columns osname|platform]\n";
	exit 0;
}
die "--columns must be osname or platform\n" unless $opt{columns} =~ /\A(?:osname|platform)\z/;
$opt{dist} //= dist_ini_name(File::Spec->catfile($FindBin::RealBin, 'dist.ini'));

my $reports = load_reports(\%opt);
die "no reports for $opt{dist}\n" unless @$reports;
my $want = $opt{version} // newest_version($reports);
my @these = grep { $_->{version} eq $want } @$reports;
die "no reports for $opt{dist} $want\n" unless @these;
my $output = $opt{output} // "cpantesters.$opt{dist}-$want.svg";
die "--output must end in .svg: only SVG keeps the links\n" unless $output =~ /\.svg\z/i;

# One cell per perl and column, holding that cell's reports.
my %cell;	# {perl label}{column label} => [reports]
my %row_key;	# perl label => sort key, so the newest perl is the top row
for my $r (@these) {
	my $perl = label($r->{perl});
	my $col  = label($r->{ $opt{columns} });
	push @{ $cell{$perl}{$col} }, $r;
	$row_key{$perl} //= perl_sort_key($perl);
}

# colored_table draws its rows in the sorted order of the keys of "data" and
# takes its row text from "row_labels", so the data is keyed by the sort key
# and labelled with the perl.
my (%data, %links, %grades);
for my $perl (keys %cell) {
	for my $col (keys %{ $cell{$perl} }) {
		my @r = @{ $cell{$perl}{$col} };
		my %n = (pass => 0, fail => 0, na => 0, unknown => 0);
		$n{ $_->{grade} // 'unknown' }++ for @r;
		$grades{$_} += $n{$_} for keys %n;
		my $graded = $n{pass} + $n{fail};
		# undef leaves a cell with only NA or UNKNOWN reports out of the color
		# scale, where it is drawn in undef_color like an empty cell; its text
		# and link still say what is there.
		$data{ $row_key{$perl} }{$col} = $graded ? $n{pass} / $graded : undef;
		my $text = "$n{pass}/$graded";
		my $other = $n{na} + $n{unknown};
		$text .= " +${other}?" if $other;	# "?" = NA or UNKNOWN reports, outside the pass/graded count
		$links{$perl}{$col} = [ report_url(pick_report(@r)), $text ];
	}
}
my @rows = sort { $row_key{$a} cmp $row_key{$b} } keys %row_key;
my %columns = map { %$_ } values %cell;
my $ncols = keys %columns;

my $fh = File::Temp->new(TEMPLATE => 'cpantesters-matrix-XXXXXX', SUFFIX => '.py', TMPDIR => 1, UNLINK => 0);
binmode $fh, ':encoding(UTF-8)';
print {$fh} python_prelude(\%links);
my $total = @these;
colored_table(
	output_file => $output,
	fh          => $fh,
	data        => \%data,
	row_labels  => \@rows,
	cmap        => 'RdYlGn',
	cb_min      => 0,
	cb_max      => 1,
	cblabel     => 'share of PASS among PASS + FAIL',
	undef_color => 'white',
	set_title   => sprintf('%s %s: %d reports, %d PASS, %d FAIL, %d NA, %d UNKNOWN',
		$opt{dist}, $want, $total, @grades{qw(pass fail na unknown)}),
	# Sizes in inches, chosen by eye on the 0.3418 matrix (29 perls x 7 OS):
	# about 0.3 in. a row and 1.1 in. a column keep "12/12 +1?" on one line.
	set_figheight => 1.5 + 0.3 * @rows,
	set_figwidth  => 2.5 + 1.1 * $ncols,
);
say "wrote $output: $total reports, " . @rows . " perls x $ncols $opt{columns} values";

# The top-level "name = ..." of a dist.ini: the key before its first [section].
sub dist_ini_name ($ini) {
	open my $in, '<', $ini or die "no --dist given, and cannot read \"$ini\": $!\n";
	while (my $line = <$in>) {
		last if $line =~ /^\s*\[/;
		if ($line =~ /^\s*name\s*=\s*(\S+)\s*$/) {
			close $in;
			return $1;
		}
	}
	close $in;
	die "no --dist given, and \"$ini\" has no top-level name\n";
}

# The summary array, from --json or from the API.
sub load_reports ($opt) {
	my $body;
	if (defined $opt->{json}) {
		open my $in, '<:raw', $opt->{json} or die "cannot read \"$opt->{json}\": $!\n";
		$body = do { local $/; <$in> };
		close $in;
	} else {
		# The whole history when no version is named, so the newest can be found.
		my $url = "https://api.cpantesters.org/v3/summary/$opt->{dist}"
			. (defined $opt->{version} ? "/$opt->{version}" : '');
		my $res = HTTP::Tiny->new(agent => 'cpantesters.matrix.pl/1', timeout => 120)->get($url);
		die "GET $url: $res->{status} $res->{reason}\n" unless $res->{success};
		$body = $res->{content};
	}
	my $all = JSON::PP->new->decode($body);
	die "the summary is not a JSON array\n" unless ref $all eq 'ARRAY';
	return [ grep { $_->{dist} eq $opt->{dist} } @$all ];
}

sub newest_version ($reports) {
	my %seen = map { $_->{version} => 1 } @$reports;
	my ($newest) = sort { version->parse($b) <=> version->parse($a) } keys %seen;
	return $newest;
}

# The report a cell links to: its newest FAIL, else its newest NA or UNKNOWN,
# else its newest PASS -- the one worth reading first.
sub pick_report (@r) {
	my %rank = (fail => 0, unknown => 1, na => 1, pass => 2);
	my @by = sort {
		($rank{ $a->{grade} } // 1) <=> ($rank{ $b->{grade} } // 1)
			or $b->{date} cmp $a->{date}
	} @r;
	return $by[0];
}

sub report_url ($r) {
	die "a report has no guid\n" unless ($r->{guid} // '') =~ /\A[0-9a-f-]+\z/;
	return "https://www.cpantesters.org/cpan/report/$r->{guid}";
}

# The API's strings go into Python as labels, so anything outside the
# characters a perl version, OS name or archname uses is replaced.
sub label ($s) {
	$s //= 'unknown';
	$s =~ s/[^\w.+-]+/ /g;
	return $s;
}

# "5.42.3" => "0994.0957.0996": descending by revision, version, subversion, so
# that a plain string sort puts the newest perl first. A perl string that is
# not three numbers sorts last, by its text.
sub perl_sort_key ($perl) {
	my @v = $perl =~ /\A(\d+)\.(\d+)\.(\d+)/ or return "zzzz.$perl";
	return join '.', (map { sprintf '%04d', 9999 - $_ } @v), $perl;
}

# Python run before Matplotlib::Simple's own script: wraps Axes.table() so that
# each data cell, found by its row and column label, gets its URL and text.
# The cell map travels as base64 JSON, so no label is ever Python source.
sub python_prelude ($links) {
	my $b64 = encode_base64(JSON::PP->new->canonical->utf8->encode($links), '');
	return <<"PY";
import base64, json
import matplotlib.axes
_links = json.loads(base64.b64decode('$b64').decode('utf-8'))
_axes_table = matplotlib.axes.Axes.table
def _linked_table(self, *args, **kwargs):
	table = _axes_table(self, *args, **kwargs)
	# colored_table's hidden imshow makes the axes square-celled, and the
	# table's own font sizing shrinks the text to fit; the figure size set
	# from the grid's shape is what should decide both
	self.set_aspect('auto')
	table.auto_set_font_size(False)
	table.set_fontsize(9)
	cells = table.get_celld()
	# (0, c) is the column-label row and (r, -1) the row-label column
	rows = {r: cell.get_text().get_text() for (r, c), cell in cells.items() if c == -1}
	cols = {c: cell.get_text().get_text() for (r, c), cell in cells.items() if r == 0}
	for (r, c), cell in cells.items():
		if r == 0 or c == -1:
			continue
		link = _links.get(rows.get(r), {}).get(cols.get(c))
		if link is None:
			cell.get_text().set_text('')
			continue
		cell.set_url(link[0])
		cell.get_text().set_url(link[0])
		cell.get_text().set_text(link[1])
	return table
matplotlib.axes.Axes.table = _linked_table
PY
}
