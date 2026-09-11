#!/usr/bin/env perl
require 5.010;
use strict;
use warnings;
use Test::More;
use Test::Exception;
use File::Temp 'tempdir';
use Stats::LikeR 'read_table';

# The .xlsx worksheet tokenizer (xlsx_ws_scan in LikeR.xs, reached through
# Stats::LikeR::_parse_xlsx_sheet).  t/read_table.xlsx.t covers a workbook that
# looks the way Excel writes one; this file covers the shapes the tokenizer has
# to survive and the ones read_table's own fast path turns on, which the
# well-formed workbook never exercises:
#
#   * a cell with no r= at all, and an r= that is not the first attribute
#   * self-closing <c/> and <row/>, blank rows, and gaps mid-row
#   * every entity form _xml_unescape used to handle, decimal and hex numeric
#     character references included, and the ones that must NOT be decoded
#   * a shared-string index that is out of range or not a number
#   * t="str" / t="b" / t="e", which take the raw <v> like a number
#   * a column reference too long to be one
#   * the fast path (aoh/hoa, assembled in XS) and the callback path (a filter,
#     or hoh) agreeing cell for cell -- the same rows reach both
#
# Fixtures are built here with core IO::Compress::Zip, as t/read_table.xlsx.t
# builds its own, so the test needs no binary fixture and no CPAN reader.

no warnings 'once';	# $IO::Compress::Zip::ZipError is read, never assigned
my $have_zip = eval { require IO::Compress::Zip; 1 };
plan skip_all => 'IO::Compress::Zip (core) not available' unless $have_zip;

my $HAVE_LEAKTRACE;
BEGIN {
	$HAVE_LEAKTRACE = eval {
		require Test::LeakTrace;
		Test::LeakTrace->import('no_leaks_ok');
		1;
	} ? 1 : 0;
}

my $dir = tempdir(CLEANUP => 1);
my $seq = 0;

# Wrap a <sheetData> body (and an optional list of <si> elements) into a
# minimal one-worksheet workbook, and return its path.
sub mk {
	my ($sheetdata, $shared) = @_;
	my $path = "$dir/x" . $seq++ . '.xlsx';
	my $ns   = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main';
	my $rns  = 'http://schemas.openxmlformats.org/package/2006/relationships';
	my $ons  = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';
	my %part = (
		'[Content_Types].xml' => '<?xml version="1.0"?>'
			. '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
			. '<Default Extension="xml" ContentType="application/xml"/>'
			. '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
			. '</Types>',
		'_rels/.rels' => qq{<?xml version="1.0"?><Relationships xmlns="$rns">}
			. qq{<Relationship Id="rId1" Type="$ons/officeDocument" Target="xl/workbook.xml"/>}
			. '</Relationships>',
		'xl/workbook.xml' => qq{<?xml version="1.0"?><workbook xmlns="$ns" xmlns:r="$ons">}
			. '<sheets><sheet name="Data" sheetId="1" r:id="rId1"/></sheets></workbook>',
		'xl/_rels/workbook.xml.rels' => qq{<?xml version="1.0"?><Relationships xmlns="$rns">}
			. qq{<Relationship Id="rId1" Type="$ons/worksheet" Target="worksheets/sheet1.xml"/>}
			. '</Relationships>',
		'xl/worksheets/sheet1.xml' => qq{<?xml version="1.0"?><worksheet xmlns="$ns">}
			. "<sheetData>$sheetdata</sheetData></worksheet>",
	);
	my @order = ('[Content_Types].xml', '_rels/.rels', 'xl/workbook.xml',
	             'xl/_rels/workbook.xml.rels', 'xl/worksheets/sheet1.xml');
	if (defined $shared) {
		$part{'xl/sharedStrings.xml'} =
			qq{<?xml version="1.0"?><sst xmlns="$ns">$shared</sst>};
		push @order, 'xl/sharedStrings.xml';
	}
	my $first = 1;
	for my $name (@order) {
		IO::Compress::Zip::zip(\$part{$name}, $path,
			Name   => $name,
			Method => IO::Compress::Zip::ZIP_CM_DEFLATE(),
			Append => !$first)
			or die "Cannot add $name to $path: $IO::Compress::Zip::ZipError";
		$first = 0;
	}
	return $path;
}

# cells with no r= at all: each one lands in the column after the last
{
	my $f = mk('<row><c t="s"><v>0</v></c><c t="s"><v>1</v></c></row>'
	         . '<row><c t="s"><v>2</v></c><c><v>7</v></c></row>',
	           '<si><t>a</t></si><si><t>b</t></si><si><t>x</t></si>');
	is_deeply( read_table($f), [ { a => 'x', b => '7' } ],
		'cells with no r= fall into sequential columns' );
}

# r= present but not the first attribute
{
	my $f = mk('<row r="1"><c s="0" r="A1" t="s"><v>0</v></c><c s="0" r="B1" t="s"><v>1</v></c></row>'
	         . '<row r="2"><c s="0" r="B2"><v>9</v></c><c s="0" r="A2" t="s"><v>2</v></c></row>',
	           '<si><t>a</t></si><si><t>b</t></si><si><t>x</t></si>');
	is_deeply( read_table($f), [ { a => 'x', b => '9' } ],
		'r= is honoured wherever it sits, and out-of-order cells still land right' );
}

# self-closing cells and rows, a mid-row gap, and a blank row
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c><c r="C1" t="s"><v>2</v></c></row>'
	         . '<row r="2"><c r="A2"><v>1</v></c><c r="C2"><v>3</v></c></row>'
	         . '<row r="3"></row>'
	         . '<row r="4"/>'
	         . '<row r="5"><c r="A5"/><c r="B5"><v>5</v></c><c r="C5"><v>6</v></c></row>',
	           '<si><t>a</t></si><si><t>b</t></si><si><t>c</t></si>');
	is_deeply( read_table($f),
		[ { a => '1',   b => undef, c => '3' },
		  { a => undef, b => '5',   c => '6' } ],
		'self-closing cell/row, blank rows and a mid-row gap' );
}

# a self-closing <row/> BETWEEN two data rows does not swallow the next one
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
	         . '<row r="2"/>'
	         . '<row r="3"><c r="A3"><v>1</v></c><c r="B3"><v>2</v></c></row>',
	           '<si><t>a</t></si><si><t>b</t></si>');
	is_deeply( read_table($f), [ { a => '1', b => '2' } ],
		'an empty <row/> between data rows is dropped, not merged' );
}

# entity decoding, in inline strings and in shared strings
{
	my $f = mk('<row r="1"><c r="A1" t="inlineStr"><is><t>h&amp;1</t></is></c>'
	         .            '<c r="B1" t="inlineStr"><is><t>h2</t></is></c></row>'
	         . '<row r="2"><c r="A2" t="inlineStr"><is><t xml:space="preserve"> a &lt;b&gt; </t></is></c>'
	         .            '<c r="B2" t="inlineStr"><is><t>&#233;&#x4E2D;&#38;</t></is></c></row>'
	         . '<row r="3"><c r="A3" t="inlineStr"><is><t>&amp;lt;</t></is></c>'
	         .            '<c r="B3" t="inlineStr"><is><t>q&quot;&apos;x</t></is></c></row>'
	         . '<row r="4"><c r="A4" t="inlineStr"><is><t>&nbsp;&bad &amp;</t></is></c>'
	         .            '<c r="B4" t="inlineStr"><is><t/></is></c></row>'
	         . '<row r="5"><c r="A5" t="inlineStr"><is><t>&#999999999999;</t></is></c>'
	         .            '<c r="B5" t="inlineStr"><is><t>&#x7F;</t></is></c></row>');
	my $rows = read_table($f);
	is( $rows->[0]{'h&1'}, " a <b> ", 'lt/gt decoded, xml:space="preserve" kept' );
	# UTF-8 bytes, because read_table hands the file back as bytes and does not
	# decode it: e-acute is C3 A9 and U+4E2D is E4 B8 AD.
	is( $rows->[0]{h2}, "\xc3\xa9\xe4\xb8\xad&",
		'decimal and hex numeric character references become UTF-8 bytes' );
	is( $rows->[1]{'h&1'}, '&lt;',
		'&amp;lt; decodes once, to &lt; and not to <' );
	is( $rows->[1]{h2}, 'q"\'x', 'quot and apos decoded' );
	is( $rows->[2]{'h&1'}, '&nbsp;&bad &',
		'an unknown entity and a bare & are copied through' );
	# A character reference past 0x7FFFFFFF is left alone. XML 1.0 forbids one
	# above #x10FFFF anyway, and the perl parser's unguarded chr() was not the
	# same answer on every build: "&#999999999999;" came back as thirteen bytes
	# of perl's extended UTF-8 on ivsize=8 and died on the 32-bit perl with
	# "Use of code point 0xFFFFFFFF is not allowed".
	is( $rows->[3]{'h&1'}, '&#999999999999;',
		'a character reference past 0x7FFFFFFF is left in the text' );
	is( $rows->[3]{h2}, "\x{7f}",
		'0x7F, the largest one-byte reference, still decodes' );
	is( $rows->[2]{h2}, undef, 'an empty <t/> is an empty cell, so undef' );
	is( scalar @$rows, 4, 'four data rows' );
}

# rich text: the runs of one <si> are concatenated
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c></row>'
	         . '<row r="2"><c r="A2" t="s"><v>1</v></c></row>',
	           '<si><t>h1</t></si><si><r><rPr/><t>A &amp; </t></r><r><t>B</t></r></si>');
	is_deeply( read_table($f), [ { h1 => 'A & B' } ],
		"a shared string's rich-text runs are concatenated" );
}

# a shared-string index that is out of range, or not a number, is an empty cell;
# t="str" / t="b" / t="e" take the raw <v>, as a number does
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c><c r="C1" t="s"><v>2</v></c></row>'
	         . '<row r="2"><c r="A2" t="s"><v>99</v></c><c r="B2" t="str"><f>X()</f><v>calc</v></c>'
	         .            '<c r="C2" t="b"><v>1</v></c></row>'
	         . '<row r="3"><c r="A3" t="s"><v>x</v></c><c r="B3" t="e"><v>#N/A</v></c>'
	         .            '<c r="C3"><v>0</v></c></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si><si><t>h3</t></si>');
	is_deeply( read_table($f),
		[ { h1 => undef, h2 => 'calc', h3 => '1' },
		  { h1 => undef, h2 => '#N/A', h3 => '0' } ],
		'bad shared-string index -> empty; str/b/e take the raw <v>' );
}

# a '>' inside an attribute value must not end the tag early
{
	my $f = mk('<row r="1" spans="1:2"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
	         . '<row r="2" customFormat="a>b"><c r="A2" t="inlineStr"><is><t>x</t></is></c>'
	         .            '<c r="B2"><v xml:space="preserve">4</v></c></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si>');
	is_deeply( read_table($f), [ { h1 => 'x', h2 => '4' } ],
		"a '>' inside an attribute value, and an attributed <v>" );
}

# A reference past XFD, the last column ECMA-376 allows, is not a reference: the
# cell goes in the next column instead. Every row is padded to the widest column
# the sheet mentions, so this is what bounds a row at the 16,384 cells the format
# allows rather than the twelve million "ZZZZZ" would ask for.
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
	         . '<row r="2"><c r="A2"><v>1</v></c><c r="ZZZZZ2"><v>2</v></c></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si>');
	is_deeply( read_table($f), [ { h1 => '1', h2 => '2' } ],
		'an over-long column reference falls back to the next column' );
	is( Stats::LikeR::_xlsx_col_idx('XFD'), 16383, 'XFD is the last column, 16383' );
	is( Stats::LikeR::_xlsx_col_idx('XFE'), -1,    'XFE is one past it, so -1' );
	is( Stats::LikeR::_xlsx_col_idx('ZZZ'), -1,    'ZZZ is past it too' );
	is( Stats::LikeR::_xlsx_col_idx('AAAA'), -1,   'and a fourth letter is -1' );
	# XFD itself still places a cell -- the cap is on a fourth LETTER, not on a
	# fourth character, so "XFD1" must not be mistaken for one -- and a row that
	# uses it is legitimately 16,384 wide. Asserted on the parser directly,
	# because read_table would fold the 16,382 unnamed columns between into one
	# key and hide the width.
	my @w;
	Stats::LikeR::_parse_xlsx_sheet_xs(
		'<sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="XFD1" t="s"><v>1</v></c></row>'
	  . '<row r="2"><c r="A2"><v>1</v></c><c r="XFD2"><v>2</v></c></row></sheetData>',
		['h1', 'h2'], sub { push @w, scalar @{ $_[0] } });
	is_deeply( \@w, [ 16384, 16384 ], 'a cell at XFD makes the row 16,384 wide' );
}

# Two cells claiming the same column in one row: the last one wins, and the
# first is released rather than leaked or freed twice. Only malformed input
# gets here, and it is what a mutation fuzz of the worksheet XML turned up.
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
	         . '<row r="2"><c r="A2"><v>1</v></c><c r="A2" t="inlineStr"><is><t>dup</t></is></c>'
	         .            '<c r="B2"><v>2</v></c></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si>');
	is_deeply( read_table($f), [ { h1 => 'dup', h2 => '2' } ],
		'a repeated column reference in one row keeps the last cell' );
	is_deeply( read_table($f, filter => { 0 => sub { 1 } }), [ { h1 => 'dup', h2 => '2' } ],
		'and the callback path agrees' );
}

# A cell whose </c> is missing. The parser reads the sheet twice -- once for the
# width, once for the rows -- and both passes have to see the same cells, so the
# skip to </c> has to happen on both even though the first wants nothing out of
# the body. When it did not, the two passes disagreed here, rows came out at
# different widths, and the alignment croak that followed freed the row buffer
# twice and segfaulted. Found by a mutation fuzz of the worksheet XML.
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v><9c><c r="B1" t="s"><v>1</v></c>'
	         .            '<c r="C1" t="s"><v>2</v></c></row>'
	         . '<row r="2"><c r="A2"><v>1</v></c><c r="B2"><v>2</v></c><c r="C2"><v>3</v></c></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si><si><t>h3</t></si>');
	my $aoh = eval { read_table($f) };
	is( $@, '', 'a cell with no </c> does not crash the parser' );
	# The unterminated cell's body runs on to the next cell's </c>, so B1 is
	# swallowed and the header's middle column comes out unnamed. That is the
	# right answer for a broken file; what matters is that both passes reach it.
	is_deeply( $aoh, [ { h1 => '1', '' => '2', h3 => '3' } ],
		'the swallowed cell leaves an unnamed column, at the sheet width' );
	is_deeply( read_table($f, filter => { 0 => sub { 1 } }), $aoh,
		'the two passes agree, so both output paths do too' );
}

# lower-case references
{
	my $f = mk('<row r="1"><c r="a1" t="s"><v>0</v></c><c r="b1" t="s"><v>1</v></c></row>'
	         . '<row r="2"><c r="a2"><v>1</v></c><c r="b2"><v>2</v></c></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si>');
	is_deeply( read_table($f), [ { h1 => '1', h2 => '2' } ],
		'lower-case column letters are read the same as upper-case' );
}

# an empty sheet, and a sheet with a header and no data rows
{
	is_deeply( read_table(mk('', '<si><t>h1</t></si>')), [],
		'a worksheet with no rows at all is an empty table' );
	my $hdr = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>',
	             '<si><t>h1</t></si><si><t>h2</t></si>');
	is_deeply( read_table($hdr), [], 'a header with no data rows is an empty table' );
	# A worksheet with no data rows comes back as an empty hash rather than one
	# empty array per column, which is what it did before the parser moved into
	# XS. (A header-only CSV does keep the columns; the two have never agreed,
	# and making them agree is a change to read_table, not to this parser.)
	is_deeply( read_table($hdr, 'output.type' => 'hoa'), {},
		'a header with no data rows is an empty hash as a hoa' );
}

# a data row wider than the header is an alignment error naming the row
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
	         . '<row r="2"><c r="A2"><v>1</v></c><c r="B2"><v>2</v></c></row>'
	         . '<row r="3"><c r="A3"><v>3</v></c><c r="B3"><v>4</v></c><c r="D3"><v>5</v></c></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si>');
	# Every row is padded to the widest row in the sheet, the header included,
	# so a wide row does not make the table ragged: it widens the header with
	# unnamed columns instead. Those collapse to one '' key, which is where the
	# duplicate-name warning comes from.
	my @warn;
	my $rows = do { local $SIG{__WARN__} = sub { push @warn, $_[0] }; read_table($f) };
	is_deeply( $rows,
		[ { h1 => '1', h2 => '2', '' => undef },
		  { h1 => '3', h2 => '4', '' => '5'   } ],
		'a row wider than the header widens every row, header included' );
	is( scalar( grep { /duplicate column name/ } @warn ), 1,
		'and the unnamed columns it adds warn once as duplicates' );
}

# The fast path and the callback path must see the same cells. aoh and hoa are
# assembled in XS from the plan; a filter and hoh are not, so the same file read
# four ways is the check that both paths agree.
{
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c><c r="C1" t="s"><v>2</v></c></row>'
	         . '<row r="2"><c r="A2" t="s"><v>3</v></c><c r="B2"><v>1</v></c></row>'
	         . '<row r="3"><c r="A3" t="s"><v>4</v></c><c r="B3"><v>2</v></c><c r="C3" t="s"><v>3</v></c></row>',
	           '<si><t>name</t></si><si><t>n</t></si><si><t>tag</t></si>'
	         . '<si><t>alpha</t></si><si><t>beta</t></si>');
	my $aoh = read_table($f);
	is_deeply( $aoh,
		[ { name => 'alpha', n => '1', tag => undef },
		  { name => 'beta',  n => '2', tag => 'alpha' } ],
		'aoh through the XS fast path' );
	is_deeply( read_table($f, 'output.type' => 'hoa'),
		{ name => ['alpha','beta'], n => ['1','2'], tag => [undef,'alpha'] },
		'hoa through the XS fast path' );
	is_deeply( read_table($f, filter => { 0 => sub { 1 } }), $aoh,
		'a filter keeps every row through the perl callback path' );
	is_deeply( read_table($f, 'output.type' => 'hoh', 'row.names' => 'name'),
		{ alpha => { n => '1', tag => undef },
		  beta  => { n => '2', tag => 'alpha' } },
		'hoh through the perl callback path' );
	is_deeply( read_table($f, 'na.strings' => 'alpha'),
		[ { name => undef, n => '1', tag => undef },
		  { name => 'beta', n => '2', tag => undef } ],
		'na.strings is applied by the fast path too' );
}

SKIP: {
	skip 'Test::LeakTrace not installed', 8 unless $HAVE_LEAKTRACE;
	skip 'running under Devel::Cover', 8 if $INC{'Devel/Cover.pm'};

	# The parser is measured with the worksheet part already decompressed, so
	# that what is counted is the code this file is about. Every path through
	# xlsx_ws_row_end is here: the row handed to a callback that drops it, one
	# that keeps it, one that dies part-way (the row buffer is the parser's
	# until the call, and the callback's after it, and exactly one of them has
	# to release it), and a blank row, a gap and an entity along the way.
	my $sst = ['h1', 'h2', 'kept'];
	my $xml = '<sheetData>'
	        . '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
	        . '<row r="2"><c r="A2" t="s"><v>2</v></c><c r="B2"><v>1</v></c></row>'
	        . '<row r="3"/>'
	        . '<row r="4"><c r="B4" t="inlineStr"><is><t>x&amp;y</t></is></c></row>'
	        . '</sheetData>';
	my $parse = \&Stats::LikeR::_parse_xlsx_sheet_xs;
	my @kept;
	no_leaks_ok { $parse->($xml, $sst, sub { }) }
		'no leaks: the parser, with the callback dropping each row';
	no_leaks_ok { @kept = (); $parse->($xml, $sst, sub { push @kept, $_[0] }); @kept = () }
		'no leaks: the parser, with the callback keeping each row';
	no_leaks_ok { eval { $parse->($xml, $sst, sub { die "stop\n" }) } }
		'no leaks: a callback that dies on the first row';
	my $n = 0;
	no_leaks_ok { $n = 0; eval { $parse->($xml, $sst, sub { die "stop\n" if ++$n == 2 }) } }
		'no leaks: a callback that dies part-way through';

	# read_table's own paths need a real file, and decompressing one is not
	# clean on every perl: 5.10.1's bundled Compress::Raw::Zlib leaks 18 SVs per
	# read in crc32(undef) (IO/Uncompress/Unzip.pm line 608), which a leak check
	# around read_table measures along with this module. Ask this perl whether
	# it does rather than naming versions -- and warm the module first, so that
	# loading it is not counted as the leak.
	my $f = mk('<row r="1"><c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c></row>'
	         . '<row r="2"><c r="A2" t="s"><v>2</v></c><c r="B2"><v>1</v></c></row>'
	         . '<row r="3"><c r="A3" t="inlineStr"><is><t>x&amp;y</t></is></c><c r="B3"/></row>',
	           '<si><t>h1</t></si><si><t>h2</t></si><si><t>v</t></si>');
	Stats::LikeR::_unzip_member($f, 'xl/worksheets/sheet1.xml');
	my $unzip_leaks = Test::LeakTrace::leaked_count(
		sub { Stats::LikeR::_unzip_member($f, 'xl/worksheets/sheet1.xml') });
	skip "IO::Uncompress::Unzip leaks $unzip_leaks SV(s) per member on this perl",
		4 if $unzip_leaks;

	my $dup = mk('<row r="1"><c r="A1" t="s"><v>0</v></c></row>'
	           . '<row r="2"><c r="A2"><v>1</v></c><c r="A2" t="inlineStr"><is><t>d</t></is></c></row>',
	             '<si><t>h1</t></si>');
	# The filter is a plain sub, and no qr// is evaluated inside these blocks:
	# 5.10.0's pp_qr() leaks an SV every time one is, and CPAN smokers run it.
	no_leaks_ok { read_table($f) }                               'no leaks: aoh fast path';
	no_leaks_ok { read_table($f, 'output.type' => 'hoa') }       'no leaks: hoa fast path';
	no_leaks_ok { read_table($f, filter => { 0 => sub { 1 } }) } 'no leaks: callback path';
	no_leaks_ok { read_table($dup) }                             'no leaks: repeated column reference';
}

done_testing;
