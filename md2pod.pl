#!/usr/bin/env perl

require 5.010;
use feature 'say';
use warnings FATAL => 'all';
use autodie ':default';
#use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Markdown::To::POD 'markdown_to_pod';
use List::MoreUtils 'first_index';
use Test::More;
use Test::Pod;

sub file2string {
	my $file = shift;
	open my $fh, '<', $file;
	return do { local $/; <$fh> };
}

# Helper to build an HTML table from extracted Markdown rows
sub table_to_html {
	my ($header, $sep, $body_ref) = @_;
	my $html = "<table>\n";

	# Process header
	$html .= "<thead>\n<tr>\n";
	my @headers = split /\|/, $header;
	# Clean empty elements if the row was wrapped in leading/trailing pipes
	shift @headers if @headers && $headers[0] =~ /^\s*$/ && $header =~ /^\s*\|/;
	pop @headers   if @headers && $headers[-1] =~ /^\s*$/ && $header =~ /\|\s*$/;
	for my $h (@headers) {
	  $h =~ s/^\s+|\s+$//g;
	  $html .= "  <th>$h</th>\n";
	}
	$html .= "</tr>\n</thead>\n<tbody>\n";
	# Process body
	for my $row (@$body_ref) {
		$html .= "<tr>\n";
		my @cells = split /\|/, $row;
		shift @cells if @cells && $cells[0] =~ /^\s*$/ && $row =~ /^\s*\|/;
		pop @cells   if @cells && $cells[-1] =~ /^\s*$/ && $row =~ /\|\s*$/;

		for my $c (@cells) {
			$c =~ s/^\s+|\s+$//g;
			# Convert Markdown inline formatting so it renders correctly inside the HTML block
			$c =~ s/`([^`]+)`/<code>$1<\/code>/g;
			$c =~ s/\*\*([^\*]+)\*\*/<b>$1<\/b>/g;
			$c =~ s/\*([^\*]+)\*/<i>$1<\/i>/g;
			$html .= "  <td>$c</td>\n";
		}

		# Pad with empty cells if a row was missing trailing pipes
		while (@cells < @headers) {
			$html .= "  <td></td>\n";
			push @cells, "";
		}
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

# 1. Pre-process the Markdown to convert GFM tables into POD HTML blocks
my ($md_processed, $tables_ref) = extract_and_convert_tables($md);

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
my $line = first_index {$_ eq '1;'} @lib;
if ($line == -1) {
	die 'Could not find correct line index';
}

# Trim everything after `1;` to prep for new POD insertion
splice @lib, 1-(scalar @lib - $line);
push @lib, @pod; 

open my $out_fh, '>', 'lib/Stats/LikeR.pm';
say $out_fh join ("\n", @lib);
close $out_fh;

pod_file_ok( 'lib/Stats/LikeR.pm' );
done_testing();
