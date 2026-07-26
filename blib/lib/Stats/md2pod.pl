#!/usr/bin/env perl

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
#use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Markdown::To::POD 'markdown_to_pod';
use List::MoreUtils 'first_index';
use Test::More;
use Test::Pod;
use Test::CPAN::Changes;

sub file2string {
	my $file = shift;
	open my $fh, '<', $file;
	return do { local $/; <$fh> };
}

# Split one GFM table row into its cells.
#
# GFM lets a cell contain a literal pipe as "\\|", so splitting on every pipe
# tears such a cell in two and shifts every later cell one column left. Split
# on unescaped pipes only, then unescape. The -1 limit keeps a trailing empty
# cell, which plain split() would discard (leaving the row a column short).
sub split_row {
	my ($row) = @_;
	my @cells = split /(?<!\\)\|/, $row, -1;
	# a row wrapped in pipes yields an empty field at each end
	shift @cells if @cells && $cells[0]  =~ /^\s*$/ && $row =~ /^\s*\|/;
	pop   @cells if @cells && $cells[-1] =~ /^\s*$/ && $row =~ /\|\s*$/;
	for my $c (@cells) {
		$c =~ s/\\\|/|/g;
		$c =~ s/^\s+|\s+$//g;
	}
	return @cells;
}

# Make text safe to drop inside an HTML element. Without this a cell such as
# "`>` `<` `>=`" (filter's operator table) turns into "<code><</code>", which
# every HTML parser reads as an unclosed tag and swallows -- two operators went
# missing from the rendered documentation that way.
sub html_escape {
	my ($text) = @_;
	$text =~ s/&/&amp;/g;
	$text =~ s/</&lt;/g;
	$text =~ s/>/&gt;/g;
	return $text;
}

# Markdown inline formatting -> HTML, for the inside of a table cell. Escaping
# has to come first, so that the tags this adds are the only markup present.
sub md_inline_to_html {
	my ($cell) = @_;
	$cell = html_escape($cell);
	$cell =~ s/`([^`]+)`/<code>$1<\/code>/g;
	$cell =~ s/\*\*([^\*]+)\*\*/<b>$1<\/b>/g;
	$cell =~ s/\*([^\*]+)\*/<i>$1<\/i>/g;
	return $cell;
}

# Helper to build an HTML table from extracted Markdown rows
sub table_to_html {
	my ($header, $sep, $body_ref) = @_;
	my $html = "<table>\n";

	# Process header
	$html .= "<thead>\n<tr>\n";
	my @headers = split_row($header);
	$html .= '  <th>' . md_inline_to_html($_) . "</th>\n" for @headers;
	$html .= "</tr>\n</thead>\n<tbody>\n";

	# Process body
	for my $row (@$body_ref) {
		my @cells = split_row($row);

		# Pad a row that was short of cells, and fold a row that has too many
		# back into its last column rather than widening the whole table (which
		# would leave every later cell under the wrong heading).
		push @cells, '' while @cells < @headers;
		if (@headers && @cells > @headers) {
			warn "md2pod: table row has " . scalar(@cells) . ' cells for '
			   . scalar(@headers) . " headings, folding the surplus into the "
			   . "last column:\n  $row\n";
			my @tail = splice @cells, $#headers;
			$cells[$#headers] = join ' ', grep { length } @tail;
		}

		$html .= "<tr>\n";
		$html .= '  <td>' . md_inline_to_html($_) . "</td>\n" for @cells;
		$html .= "</tr>\n";
	}
	$html .= "</tbody>\n</table>\n";

	# Ensure blank lines around the =begin and =end directives for valid POD
	return "\n\n=begin html\n\n$html\n=end html\n\n";
}

# Pre-processor to extract GFM tables and replace them with alphanumeric placeholders
sub extract_and_convert_tables {
	my ($text) = @_;
	my @lines = split /\n/, $text;
	my @out;
	my @saved_tables;
	my $i = 0;

	while ($i < @lines) {
	  # Look for a table header followed by a standard GFM separator row
	  if ($i + 1 < @lines &&
		   $lines[$i] =~ /\|/ &&
		   $lines[$i+1] =~ /^[ \t]*\|?[ \t]*:?-+[-: \t]*\|/) {

		   my $header = $lines[$i];
		   my $sep = $lines[$i+1];
		   my @body;
		   $i += 2;

		   # Consume consecutive data rows (must contain at least one pipe)
		   while ($i < @lines && $lines[$i] =~ /\|/) {
		       push @body, $lines[$i];
		       $i++;
		   }

		   my $html = table_to_html($header, $sep, \@body);
		   push @saved_tables, $html;
		   # Use an alphanumeric placeholder to prevent Markdown parser interference
		   push @out, "\n\nHTMLTABLEPLACEHOLDER" . ($#saved_tables) . "\n\n";
	  } else {
		   push @out, $lines[$i];
		   $i++;
	  }
	}
	return (join("\n", @out), \@saved_tables);
}

# Markdown::To::POD applies "_..._" emphasis inside a word, which GFM does not.
# A bare Perl identifier is one word, so "_rename_inplace" came out as
# "I<rename>inplace" -- the heading named after that helper lost its name.
# Backslash-escaping the underscores is honoured by the converter and produces
# the identifier verbatim, so do that to every identifier-shaped token that is
# not already inside an inline code span (whose contents are passed through as
# written, and where backslashes would show up literally).
sub protect_underscore_identifiers {
	my ($text) = @_;
	my @out;
	for my $chunk (split /(`[^`\n]*`)/, $text) {
		if ($chunk =~ /\A`/) { push @out, $chunk; next }
		$chunk =~ s{(?<![\\\w])(_*[A-Za-z0-9]+(?:_[A-Za-z0-9]+)+_*)(?![\w])}
		           {my $t = $1; $t =~ s/_/\\_/g; $t}gex;
		push @out, $chunk;
	}
	return join '', @out;
}

# The GitHub anchor for a heading: lowercased, formatting and punctuation
# dropped, spaces turned into hyphens. Used to resolve "[text](#anchor)" back
# to the heading it points at.
sub gh_anchor {
	my ($text) = @_;
	$text = lc $text;
	$text =~ s/`//g;
	$text =~ s/\*//g;
	$text =~ s/[^\w\- ]//g;
	$text =~ s/ /-/g;
	return $text;
}

# Markdown inline formatting -> POD formatting codes, for link text.
sub md_inline_to_pod {
	my ($text) = @_;
	$text =~ s/`([^`]+)`/C<$1>/g;
	$text =~ s/\*\*([^\*]+)\*\*/B<$1>/g;
	$text =~ s/\*([^\*]+)\*/I<$1>/g;
	return $text;
}

# Markdown::To::POD mangles links: "[text](#anchor)" becomes "L<#anchor>",
# which is not POD link syntax, and the link text is thrown away entirely --
# "[`read_table`](#)" came out as "L<#>", with the words gone. An external
# "[CPAN](https://...)" fares worse and comes out as "LL<https://...>".
#
# Convert them here instead. An in-document link is resolved against the
# document's own headings and becomes a real POD section link; one that resolves
# to nothing keeps its text and loses only the link. Placeholders keep the
# result away from the converter, the same way the tables and code blocks do.
sub extract_links {
	my ($text, $headings) = @_;
	my @saved;
	my $unresolved = 0;

	$text =~ s{(?<!\\)\[([^\]\[]*)\]\(([^()\s]*)\)}{
		my ($label, $target) = ($1, $2);
		my $pod = md_inline_to_pod($label);
		if ($target =~ /\A#(.*)\z/s) {
			my $section = $headings->{$1};
			if (!defined $section) {
				$unresolved++;
				warn "md2pod: link [$label](#$1) matches no heading; "
				   . "keeping the text, dropping the link\n" if length $1;
			}
			elsif ($section eq $pod) { $pod = qq{L</"$section">} }
			else                     { $pod = qq{L<$pod|/"$section">} }
		}
		elsif (length $target) {
			$pod = qq{L<$pod|$target>};
		}
		push @saved, $pod;
		'PODLINKPLACEHOLDER' . $#saved . 'Z';
	}gex;

	return ($text, \@saved);
}

# Every ATX heading, keyed by its GitHub anchor, with the text as it will read
# in the POD. Fenced and indented code is skipped so a "# comment" line is
# never mistaken for a heading.
sub collect_headings {
	my ($text) = @_;
	my %by_anchor;
	my $in_fence = 0;
	for my $ln (split /\n/, $text, -1) {
		$in_fence = !$in_fence if $ln =~ /^[ \t]*(?:```|~~~)/;
		next if $in_fence;
		next if $ln =~ /^(?:\t| {4,})/;
		next unless $ln =~ /^\#{1,6}[ \t]+(\S.*?)[ \t]*\#*[ \t]*$/;
		my $heading = $1;
		my $anchor  = gh_anchor($heading);
		next unless length $anchor;
		# GitHub disambiguates repeats with -1, -2, ...; first one wins here
		$by_anchor{$anchor} = md_inline_to_pod($heading)
			unless exists $by_anchor{$anchor};
	}
	return \%by_anchor;
}

# Pull indented (4-space) code blocks out before conversion, and put them back
# as POD verbatim paragraphs afterwards.
#
# Handing them to the Markdown converter is not safe. It reads them as prose
# whenever the block follows a list, and prose is then reformatted: the
# indentation is flattened, a leading "#" on a comment line becomes a heading
# (README.md's dropna example produced "=head1 { A => [1, 2], ... }" that way,
# which truncates the section), "_name_" turns into italics, ">" turns into
# "E<gt>", and backticks turn into C<>. None of that belongs in a code sample.
#
# A block starts at a line indented by >= 4 columns that follows a blank line
# and is not a list marker (a nested "  - item" continuation can also be
# indented, and must be left to the converter). Interior blank lines are kept
# as long as indented content resumes after them.
sub extract_code_blocks {
	my ($text) = @_;
	my @lines = split /\n/, $text, -1;
	my @out;
	my @saved;
	my $blank = 1;                       # start of file counts as a blank line
	my $i = 0;

	while ($i < @lines) {
		my $is_indented = $lines[$i] =~ /^(?:\t| {4,})\S/
		               && $lines[$i] !~ /^\s+(?:[-*+]|[0-9]+[.)])\s/;
		if (!($blank && $is_indented)) {
			$blank = $lines[$i] =~ /^\s*$/ ? 1 : 0;
			push @out, $lines[$i++];
			next;
		}

		my @body;
		while ($i < @lines) {
			if ($lines[$i] =~ /^\s*$/) {           # keep an interior blank line
				my $j = $i;
				$j++ while $j < @lines && $lines[$j] =~ /^\s*$/;
				last if $j >= @lines || $lines[$j] !~ /^(?:\t| {4,})\S/;
				push @body, '' for $i .. $j - 1;
				$i = $j;
				next;
			}
			last unless $lines[$i] =~ /^(?:\t| {4,})/;
			push @body, $lines[$i++];
		}

		push @saved, \@body;
		push @out, '', 'PODVERBATIMPLACEHOLDER' . $#saved . 'Z', '';
		$blank = 1;
	}
	return (join("\n", @out), \@saved);
}

# One saved block as a POD verbatim paragraph: strip the indentation the whole
# block shares, then re-indent by a single space, so relative indentation
# inside the sample survives exactly.
sub code_block_to_pod {
	my ($body) = @_;
	my $common;
	for my $line (@$body) {
		next unless $line =~ /\S/;
		my ($lead) = $line =~ /^([ \t]*)/;
		$lead =~ s/\t/    /g;
		my $n = length $lead;
		$common = $n if !defined($common) || $n < $common;
	}
	$common = 0 unless defined $common;

	my @out;
	for my $line (@$body) {
		if ($line !~ /\S/) { push @out, ''; next }
		$line =~ s/\t/    /g;
		push @out, ' ' . substr($line, $common);
	}
	return join "\n", @out;
}

# Markdown::To::POD's list-detection regex has no blank line requirement before
# a following heading, so a heading glued straight onto a list (no blank line
# between them) gets swallowed into the final list item: the "=headN" is then
# emitted *inside* the "=over", and the list's "=back" lands after it
# ("You forgot a '=back' before '=headN'", "=back without =over"). This is a
# parse-time failure, so it must be repaired in the Markdown before conversion.
#
# Guarantee a blank line before every heading (ATX "# X" and Setext "X" over a
# row of "=" or "-"). The "\S" mirrors the converter's own header regex (header
# text required, so a bare "###" is left alone). Fenced code regions are
# skipped so "# comment" lines inside them are untouched. GFM table separators
# never match the Setext underline test because they contain pipes.
sub ensure_blank_before_headings {
	my ($text) = @_;
	my @lines = split /\n/, $text, -1;
	my @out;
	my $in_fence = 0;
	for my $j (0 .. $#lines) {
		my $ln = $lines[$j];
		$in_fence = !$in_fence if $ln =~ /^[ \t]*(?:```|~~~)/;
		my $is_atx = !$in_fence && $ln =~ /^\#{1,6}[ \t]*\S/;
		my $is_setext = !$in_fence && $ln =~ /\S/
			&& $j < $#lines && $lines[$j+1] =~ /^[ \t]*(?:=+|-+)[ \t]*$/;
		push @out, ''
			if ($is_atx || $is_setext) && @out && $out[-1] ne '';
		push @out, $ln;
	}
	return join "\n", @out;
}

# Markdown::To::POD emits a nested list's "=over" directly after the parent
# "=item" line with no blank line between them. POD requires a blank line
# before every command paragraph, so without it the "=over" is absorbed into
# the item's text rather than opening a list; the matching inner "=back" then
# closes the *outer* list, orphaning later "=item"/"=back" directives
# ("'=item' outside of any '=over'", "=back without =over").
#
# Repair by guaranteeing a blank line before every POD command paragraph.
# "=begin X" / "=end X" data blocks (the HTML tables) are copied verbatim so
# their raw contents are never rewritten.
sub fix_pod_command_spacing {
	my ($pod) = @_;
	my @in = split /\n/, $pod, -1;
	my @out;
	my $in_data = 0;
	for my $cmd_line (@in) {
		if ($in_data) {
			push @out, $cmd_line;
			$in_data = 0 if $cmd_line =~ /^=end\b/;
			next;
		}
		if ($cmd_line =~ /^=\w+/) {
			# a command paragraph must be preceded by a blank line
			push @out, '' if @out && $out[-1] ne '';
			push @out, $cmd_line;
			$in_data = 1 if $cmd_line =~ /^=begin\b/;
		}
		else {
			push @out, $cmd_line;
		}
	}
	return join "\n", @out;
}

# Guarantee balanced =over/=back. Even with the blank-line repairs above, the
# converter can emit a heading while a list is still open (e.g. a Setext-
# underlined heading, or any heading the pre-processor's normalization missed),
# leaving the matching =back stranded after the heading. Close any list still
# open when a heading / =cut / end-of-file is reached, and drop any =back that
# has no open =over. "=begin X" / "=end X" data blocks are passed through
# verbatim so their contents are never miscounted.
sub balance_pod_over_back {
	my ($pod) = @_;
	my @in = split /\n/, $pod, -1;
	my @out;
	my $depth = 0;
	my $in_data = 0;
	for my $bal_line (@in) {
		if ($in_data) {
			push @out, $bal_line;
			$in_data = 0 if $bal_line =~ /^=end\b/;
			next;
		}
		if ($bal_line =~ /^=begin\b/) {
			push @out, $bal_line;
			$in_data = 1;
			next;
		}
		if ($bal_line =~ /^=over\b/) {
			$depth++;
			push @out, $bal_line;
			next;
		}
		if ($bal_line =~ /^=back\b/) {
			# drop a =back that has no matching open =over
			if ($depth > 0) {
				$depth--;
				push @out, $bal_line;
			}
			next;
		}
		if ($bal_line =~ /^=(?:head\d+|cut|pod|encoding)\b/) {
			while ($depth > 0) {
				push @out, '', '=back';
				$depth--;
			}
			push @out, '' if @out && $out[-1] ne '';
			push @out, $bal_line;
			next;
		}
		push @out, $bal_line;
	}
	while ($depth > 0) {
		push @out, '', '=back';
		$depth--;
	}
	return join "\n", @out;
}

my $md = file2string('README.md');

# 0. Ensure headings are separated from preceding blocks so the converter's
#    list detection terminates correctly before them
$md = ensure_blank_before_headings($md);
my $md_later = $md;

# 1. Pre-process the Markdown to convert GFM tables into POD HTML blocks
my ($md_processed, $tables_ref) = extract_and_convert_tables($md);

# 1b. Set indented code blocks aside so the converter cannot reformat them
my ($code_ref);
($md_processed, $code_ref) = extract_code_blocks($md_processed);

# 1c. Resolve links against the document's headings, and keep identifiers with
#     a leading underscore from being read as emphasis. Both run after the code
#     blocks are out of the way, so code samples are never rewritten.
my $links_ref;
($md_processed, $links_ref) = extract_links($md_processed, collect_headings($md));
$md_processed = protect_underscore_identifiers($md_processed);

# 2. Convert standard markdown to POD
my $pod = markdown_to_pod($md_processed);

# 3. Restore the HTML tables back into the generated POD
for my $idx (0 .. $#$tables_ref) {
	my $table_html = $tables_ref->[$idx];
	# Anchor the end of the number with \b: without it the /g replace for a
	# short index (e.g. 1) also matches the prefix of longer placeholders
	# (HTMLTABLEPLACEHOLDER10, ...11), dropping the wrong table there and
	# leaving a stray leftover digit. \b stops after the last digit, so
	# HTMLTABLEPLACEHOLDER1 no longer matches inside HTMLTABLEPLACEHOLDER10.
	$pod =~ s/HTMLTABLEPLACEHOLDER${idx}\b/$table_html/g;
}

# 3b. Restore the code blocks as verbatim paragraphs. Anchored with \b like
#     the tables above, so PLACEHOLDER1 cannot match inside PLACEHOLDER10.
for my $idx (0 .. $#$code_ref) {
	my $verbatim = code_block_to_pod($code_ref->[$idx]);
	$pod =~ s/^[ \t]*PODVERBATIMPLACEHOLDER${idx}Z[ \t]*$/$verbatim/mg;
}
if ($pod =~ /(PODVERBATIMPLACEHOLDER\d+Z)/) {
	die "md2pod: code-block placeholder $1 survived conversion\n";
}

# 3c. Restore the links
for my $idx (0 .. $#$links_ref) {
	my $link = $links_ref->[$idx];
	$pod =~ s/PODLINKPLACEHOLDER${idx}Z/$link/g;
}
if ($pod =~ /(PODLINKPLACEHOLDER\d+Z)/) {
	die "md2pod: link placeholder $1 survived conversion\n";
}

# 4. Repair command-paragraph spacing so nested lists stay balanced POD
$pod = fix_pod_command_spacing($pod);

# 5. Close any list left open across a heading and drop stray =back directives
$pod = balance_pod_over_back($pod);

my @pod = split /\n/, $pod;
unshift @pod, "=encoding utf8\n";

say 'Writing read.me.pod from README.md, which must be copied into lib/Stats/LikeR.pm';
open my $fh, '>', 'read.me.pod';
say $fh join ("\n", @pod);
close $fh;

my $lib = file2string('lib/Stats/LikeR.pm');
my @lib = split /\n/, $lib;
my @marks = grep { $lib[$_] eq '1;' } 0 .. $#lib;
if (!@marks) {
	die "md2pod: no '1;' line in lib/Stats/LikeR.pm to append the POD after\n";
}
if (@marks > 1) {
	die 'md2pod: lib/Stats/LikeR.pm has ' . scalar(@marks) . " lines that are "
	  . "just '1;' (at " . join(', ', map { $_ + 1 } @marks) . '); cannot tell '
	  . "which one ends the code\n";
}

# Trim everything after `1;` to prep for new POD insertion
splice @lib, $marks[0] + 1;
push @lib, @pod; 

open my $out_fh, '>', 'lib/Stats/LikeR.pm';
say $out_fh join ("\n", @lib);
close $out_fh;

pod_file_ok( 'lib/Stats/LikeR.pm' );

my $outfile = 'Changes';
my $dist    = 'Stats-LikeR'; # Inferred from your documentation

open my $out, '>', $outfile or die "Cannot write $outfile: $!\n";

# Write the mandatory CPAN::Changes::Spec header
say $out "Revision history for $dist\n";

my ($needs_bullet, $in_code_block) = (0, 0);
my @md_later = split /\n/, $md_later;
my $fi = first_index {$_ eq '# Changes'} @md_later;
if ($fi == -1) {
	# first_index returns -1, and "-1 .. $#md_later" would quietly walk the
	# whole array from its last element, emitting a Changes file made of one
	# stray line. Say so instead.
	die "md2pod: README.md has no '# Changes' section\n";
}
# start past the heading itself, which is not a release
foreach my $i ($fi + 1 .. $#md_later) {
	my $line = $md_later[$i];
	# Stop at the copyright footer
	last if $line =~ /^#\s+COPYRIGHT AND LICENSE\s*$/i;
	# Toggle markdown code blocks (```)
	if ($line =~ /^```/) {
	  $in_code_block = !$in_code_block;
	  next;
	}
	# Handle Versions ("## 0.21 2026-07-07 CDT", or "## 0.27" with no date yet).
	# The date part has to be optional: without that, an undated version heading
	# fell through to the plain-text branch and was written out as " - ## 0.27",
	# so the release had no version line at all.
	if ($line =~ /^##\h+([\d._]+)(?:\h+(\S.*?))?\h*$/) {
	  my ($version, $date) = ($1, $2);
	  # CPAN::Changes wants a date on every release, and changes_file_ok() at
	  # the end of this script enforces it. Say which heading to fix rather
	  # than leaving the reader to work it out from the line number.
	  warn "md2pod: README.md's '## $version' heading has no release date, so "
	     . "Changes will not satisfy the CPAN::Changes spec. Write it as "
	     . "'## $version YYYY-MM-DD TZ', the way the older entries are.\n"
		unless defined $date;
	  say $out defined($date) ? "$version $date" : $version;
	  $needs_bullet = 1;
	} elsif ($line =~ /^###\s+(.+)/) {# Handle Groups (e.g., "### read_table")
	  print $out " [$1]\n";
	  $needs_bullet = 1;
	} elsif ($line =~ /^####\s+(.+)/) {
	# Handle Sub-Groups (e.g., "#### Bug fixes")
	  # CPAN Spec doesn't formally have sub-groups, so we format it as a distinct bulleted header
	  say $out " - $1:";
	  $needs_bullet = 1;
	} elsif ($line =~ /^\s*[-*]\s+(.+)/) {	# Handle explicit Markdown bullets
	  say $out " - $1";
	  $needs_bullet = 0;
	} elsif ($line =~ /^\s*$/) {# Handle empty lines
	  say $out '';
	  $needs_bullet = 1; # Reset so the next text block gets a bullet
	} else {# Handle normal text or indented code
		# 4-space indented code from Markdown, or a fenced block: keep it
		# indented for CPAN. Only the four columns Markdown uses to mark the
		# block are removed, so indentation *inside* the sample survives --
		# stripping all leading whitespace flattened nested code to one column.
		# ($1 was also read here when the fence branch matched, which left it
		# holding a capture from some earlier line entirely.)
		if ($in_code_block || $line =~ /^\h{4,}\S/) {
			my $code_line = $line;
			$code_line =~ s/\A\h{1,4}// unless $in_code_block;
			$code_line =~ s/\h+\z//;
			say $out "     $code_line";
		} else {
			# Strip leading/trailing formatting like **bold** just in case it breaks flow, 
			# though CPAN::Changes technically allows it as raw text.
			$line =~ s/\*\*(.+?)\*\*/$1/g; 
			if ($needs_bullet) {
				 print $out " - $line\n";
				 $needs_bullet = 0;
			} else {
				 # Continuation of the previous bullet
				 print $out "   $line\n";
			}
		}
	}
}

close $out;

say "Successfully generated '$outfile' from 'README.md'";
changes_file_ok('Changes');
done_testing();

