#!/usr/bin/env perl

require 5.010;
use warnings FATAL => 'all';
use File::Temp;
use Stats::LikeR;
use Test::Exception; # dies_ok, throws_ok
use Test::More;
use Test::LeakTrace 'no_leaks_ok';

my $dir = File::Temp->newdir;
my $seq = 0;
# Return a fresh (csv, tex) path pair under the temp dir.
sub paths {
	$seq++;
	return ("$dir/t$seq.csv", "$dir/t$seq.tex");
}
sub slurp {
	my $file = shift;
	open my $fh, '<', $file or die "cannot read $file: $!";
	local $/;
	return <$fh>;
}
# index()-based containment checks avoid backslash-escaping headaches.
sub has {
	my ($hay, $needle, $name) = @_;
	ok(index($hay, $needle) != -1, $name) or diag("    missing substring: $needle");
}
sub lacks {
	my ($hay, $needle, $name) = @_;
	ok(index($hay, $needle) == -1, $name) or diag("    unexpected substring: $needle");
}
# Expected column spec for n columns at the default 'c' alignment.
sub spec {
	my $n = shift;
	return '\begin{tabular}{|' . ('c|' x $n) . '} \hline';
}

#--------
# byte-exact output (locks the whole format) + delimited file still written
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1], ['y', 2]], $csv,
		'tex.tab.file' => $tex, 'row.names' => 0);
	my $expected = <<'END_TEX';
% written by Stats::LikeR write_table
\begin{tabular}{|c|c|} \hline
\textbf{k} & \textbf{v} \\ \hline
\textbf{x} & 1\\
\textbf{y} & 2\\
\hline \end{tabular}
END_TEX
	ok(-e $tex, 'tex.tab.file: LaTeX file is created');
	is(slurp($tex), $expected, 'tex.tab.file: AoA output is byte-exact');
	is(slurp($csv), "k,v\nx,1\ny,2\n",
		'tex.tab.file: the delimited file is still written alongside');
}

#--------
# provenance comment + structural scaffolding
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv, 'tex.tab.file' => $tex, 'row.names' => 0);
	my $t = slurp($tex);
	has($t, '% written by Stats::LikeR write_table', 'provenance comment is present');
	has($t, spec(2), 'default alignment is c, one cell per column');
	has($t, '\hline \end{tabular}', 'table is closed with \hline \end{tabular}');
	has($t, '\textbf{k}', 'header cells are bold');
}

#--------
# tex.col.align
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv,
		'tex.tab.file' => $tex, 'tex.col.align' => 'l', 'row.names' => 0);
	has(slurp($tex), '\begin{tabular}{|l|l|} \hline', "tex.col.align => 'l'");
}
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv,
		'tex.tab.file' => $tex, 'tex.col.align' => 'r', 'row.names' => 0);
	has(slurp($tex), '\begin{tabular}{|r|r|} \hline', "tex.col.align => 'r'");
}

#--------
# tex.bold.1st.col
#--------
{
	my ($csv, $tex) = paths();
	# default (on): first data cell is wrapped in \textbf{}
	write_table([[qw(k v)], ['xx', 1]], $csv, 'tex.tab.file' => $tex, 'row.names' => 0);
	has(slurp($tex), '\textbf{xx}', 'tex.bold.1st.col defaults on: first column bold');
}
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['xx', 1]], $csv,
		'tex.tab.file' => $tex, 'tex.bold.1st.col' => 0, 'row.names' => 0);
	my $t = slurp($tex);
	lacks($t, '\textbf{xx}', 'tex.bold.1st.col => 0: first data cell not bold');
	has($t, '\textbf{k}', 'tex.bold.1st.col => 0: header still bold');
}

#--------
# tex.format (%.4g on numeric cells only)
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(num txt)], [10.12345, 'a_b']], $csv,
		'tex.tab.file' => $tex, 'tex.format' => 1, 'row.names' => 0);
	my $t = slurp($tex);
	has($t, '10.12', 'tex.format => 1: numeric cell rendered with %.4g');
	lacks($t, '10.12345', 'tex.format => 1: full-precision value replaced');
	has($t, 'a\_b', 'tex.format => 1: non-numeric cell still escaped, not formatted');
}
{
	my ($csv, $tex) = paths();
	write_table([[qw(num)], [10.12345]], $csv, 'tex.tab.file' => $tex, 'row.names' => 0);
	has(slurp($tex), '10.12345', 'tex.format off (default): numeric cell left as-is');
}

#--------
# tex.size (directive emitted after \begin{tabular})
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv,
		'tex.tab.file' => $tex, 'tex.size' => '\small', 'row.names' => 0);
	my $t = slurp($tex);
	has($t, '\small', 'tex.size: directive present');
	ok(index($t, '\small') > index($t, '\begin{tabular}'),
		'tex.size: directive follows \begin{tabular}');
}
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv, 'tex.tab.file' => $tex, 'row.names' => 0);
	lacks(slurp($tex), '\small', 'tex.size absent: no stray size directive');
}

#--------
# tex.comment (string and array ref), placed before \begin{tabular}
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv,
		'tex.tab.file' => $tex, 'tex.comment' => 'single note', 'row.names' => 0);
	my $t = slurp($tex);
	has($t, "% single note", 'tex.comment scalar: emitted as a % line');
	ok(index($t, '% single note') < index($t, '\begin{tabular}'),
		'tex.comment scalar: appears before the tabular');
}
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv, 'tex.tab.file' => $tex,
		'tex.comment' => ['line one', 'line two'], 'row.names' => 0);
	my $t = slurp($tex);
	has($t, "% line one", 'tex.comment arrayref: first line');
	has($t, "% line two", 'tex.comment arrayref: second line');
}
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], ['x', 1]], $csv, 'tex.tab.file' => $tex, 'row.names' => 0);
	# only the provenance comment line should be present
	my @comments = grep { /^%/ } split /\n/, slurp($tex);
	is(scalar(@comments), 1, 'tex.comment absent: only the provenance comment');
}

#--------
# cell escaping: # _ % &  and  > -> \textgreater{}
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(v)], ['a_b'], ['x>y'], ['50%'], ['p#q'], ['r&s']], $csv,
		'tex.tab.file' => $tex, 'row.names' => 0);
	my $t = slurp($tex);
	has($t, 'a\_b',              'escape: underscore');
	has($t, 'x\textgreater{}y',  'escape: > becomes \textgreater{}');
	has($t, '50\%',              'escape: percent');
	has($t, 'p\#q',              'escape: hash');
	has($t, 'r\&s',              'escape: ampersand');
}

#--------
# \includesvg{...svg} cells pass through unescaped
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(fig)], ['\includesvg{a_b.svg}']], $csv,
		'tex.tab.file' => $tex, 'row.names' => 0);
	my $t = slurp($tex);
	has($t, '\includesvg{a_b.svg}', 'includesvg: passed through verbatim');
	lacks($t, 'a\_b.svg', 'includesvg: underscore inside is NOT escaped');
}

#--------
# no tex.tab.file => no LaTeX file, delimited file unaffected
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(k v)], [1, 2]], $csv, 'row.names' => 0);
	ok(-e $csv,  'without tex.tab.file: delimited file still written');
	ok(!-e $tex, 'without tex.tab.file: no LaTeX file appears');
}

#--------
# LaTeX output works for every data-frame shape (column counts + files)
#--------
{
	my ($csv, $tex) = paths();
	write_table({ a => 1, b => 2, c => 3 }, $csv, 'tex.tab.file' => $tex);
	# flat hash + row.names(default on) => empty label col + a,b,c = 4 columns
	has(slurp($tex), spec(4), 'flat hash: 4 columns (label + a,b,c)');
	is(slurp($csv), ",a,b,c\n1,1,2,3\n", 'flat hash: delimited file correct');
}
{
	my ($csv, $tex) = paths();
	write_table({ p => [1, 2], q => [3, 4] }, $csv, 'tex.tab.file' => $tex, 'row.names' => 0);
	has(slurp($tex), spec(2), 'HoA: 2 columns');
	ok(-e $tex, 'HoA: LaTeX file created');
}
{
	my ($csv, $tex) = paths();
	write_table({ r1 => { x => 1, y => 2 }, r2 => { x => 3, y => 4 } }, $csv,
		'tex.tab.file' => $tex);
	# HoH + row.names(default on) => empty label col + x,y = 3 columns
	has(slurp($tex), spec(3), 'HoH: 3 columns (label + x,y)');
}
{
	my ($csv, $tex) = paths();
	write_table([{ name => 'A', age => 1 }, { name => 'B', age => 2 }], $csv,
		'tex.tab.file' => $tex, 'row.names' => 0);
	has(slurp($tex), spec(2), 'AoH: 2 columns');
}

#--------
# AoA input: first inner array is the header unless col.names is given
#--------
{
	my ($csv, $tex) = paths();
	write_table([[qw(gene score)], ['TP53', 0.9], ['BRCA1', 0.7]], $csv,
		'tex.tab.file' => $tex, 'row.names' => 0);
	is(slurp($csv), "gene,score\nTP53,0.9\nBRCA1,0.7\n",
		'AoA: first inner array used as header row');
	has(slurp($tex), '\textbf{gene} & \textbf{score}', 'AoA: header from first inner array');
}
{
	my ($csv, $tex) = paths();
	write_table([['TP53', 0.9], ['BRCA1', 0.7]], $csv,
		'col.names' => [qw(gene score)], 'tex.tab.file' => $tex);
	# col.names supplied => every inner array is data; row.names default => index col
	is(slurp($csv), ",gene,score\n1,TP53,0.9\n2,BRCA1,0.7\n",
		'AoA: col.names names columns; all inner arrays are data');
	has(slurp($tex), '\textbf{gene}', 'AoA: col.names header appears in LaTeX');
}

#--------
# error paths
#--------
{
	my ($csv, $tex) = paths();
	throws_ok {
		write_table([['h1', 'h2'], { a => 1 }], $csv,
			'tex.tab.file' => $tex, 'row.names' => 0)
	} qr/Array of Arrays/, 'AoA with a non-array element croaks';
}
{
	my ($csv, $tex) = paths();
	throws_ok {
		write_table([[qw(h)], [[1]]], $csv, 'tex.tab.file' => $tex, 'row.names' => 0)
	} qr/nested reference/, 'AoA with a nested-reference cell croaks';
}
{
	my ($csv, undef) = paths();
	throws_ok {
		write_table([[qw(k v)], [1, 2]], $csv,
			'tex.tab.file' => "$dir/no_such_subdir/out.tex", 'row.names' => 0)
	} qr/tex\.tab\.file/, 'unwritable tex.tab.file path croaks';
}

#--------
# leak safety: success and croak paths (mortal collector must be reclaimed)
#--------
no_leaks_ok {
	my ($csv, $tex) = paths();
	eval {
		write_table([[qw(gene score)], ['TP53', 0.912345], ['BRCA1', 0.7]], $csv,
			'tex.tab.file' => $tex, 'tex.format' => 1, 'tex.size' => '\small',
			'tex.comment' => ['run 3', 'q < 0.05'], 'row.names' => 0)
	}
} 'write_table(tex.tab.file): no memory leaks on success' unless $INC{'Devel/Cover.pm'};

no_leaks_ok {
	my ($csv, $tex) = paths();
	eval {
		write_table([{ x => 1, y => [1, 2] }], $csv, 'tex.tab.file' => $tex)
	}
} 'write_table(tex.tab.file): no leaks on AoH nested-ref croak' unless $INC{'Devel/Cover.pm'};

no_leaks_ok {
	my ($csv, $tex) = paths();
	eval {
		write_table([[qw(h)], [[1]]], $csv, 'tex.tab.file' => $tex, 'row.names' => 0)
	}
} 'write_table(tex.tab.file): no leaks on AoA nested-ref croak' unless $INC{'Devel/Cover.pm'};

no_leaks_ok {
	my ($csv, undef) = paths();
	eval {
		write_table([[qw(k v)], [1, 2]], $csv,
			'tex.tab.file' => "$dir/no_such_subdir/out.tex", 'row.names' => 0)
	}
} 'write_table(tex.tab.file): no leaks when the tex file cannot be opened' unless $INC{'Devel/Cover.pm'};

done_testing();
