#!/usr/bin/env perl
#
# Guard against building the module with the optimizer switched off.
#
# ExtUtils::MakeMaker's OPTIMIZE= *replaces* $Config{optimize} rather than
# adding to it.  compile.sh passed OPTIMIZE='-Wall', which is not a warning
# flag appended to perl's own -O2 -- it is the whole optimize line, so every
# build made by compile.sh, and every build `make install` then installed, was
# compiled at -O0.  Measured on perl-5.44.0 over 1e6 NVs: sum() 4.10 ms at -O0
# against 2.63 ms at -O2, cor() 25.5 ms against 12.9 ms.  Nothing failed; the
# module was simply half its speed, and had been for as long as compile.sh had
# said that.
#
# The same default lived in test.all.perls.pl, so the whole support matrix was
# built and timed at -O0 too.
#
# This is a property of the build scripts, not of the built module: there is
# no way to ask a loaded .so what -O it was compiled at.  So the check is on
# the scripts, and it can only run in the distribution root -- neither script
# ships (see MANIFEST.SKIP), so on a CPAN smoker there is nothing here to
# check and this file skips.  That is not a cross-validation test quietly
# skipping; it is a repo hygiene check that has no subject outside the repo.

require 5.010;
use strict;
use warnings FATAL => 'all';
use Test::More;
use File::Spec;

# The scripts sit next to t/, i.e. in the distribution root.
my ( $vol, $dir ) = File::Spec->splitpath( File::Spec->rel2abs(__FILE__) );
my $root = File::Spec->catpath( $vol, File::Spec->catdir( $dir, File::Spec->updir ), '' );

# What to look at in each script, and where its effective OPTIMIZE comes from.
#
# These are deliberately narrow patterns rather than one generic OPTIMIZE=
# sweep: test.all.perls.pl mentions the word twice more, once in its --help
# text and once interpolating "OPTIMIZE=$optimize" into the Makefile.PL command
# line, and neither of those carries the value.  The value is the default the
# variable is initialised to, so that is what gets read.
#
# test.all.perls.pl is not executed to ask it -- it says `use 5.044`, and this
# file is run under every perl in the matrix, 5.10.1 included.
my @SCRIPTS = (
	{
		file  => 'compile.sh',
		what  => "the OPTIMIZE= handed to Makefile.PL",
		regex => qr/(?<!\S)OPTIMIZE=(?|'([^']*)'|"([^"]*)"|(\S+))/,
	},
	{
		file  => 'test.all.perls.pl',
		what  => "the default \$optimize",
		regex => qr/\$optimize\s*=\s*(?|'([^']*)'|"([^"]*)")/,
	},
);

my @present = grep { -f File::Spec->catfile( $root, $_->{file} ) } @SCRIPTS;

plan skip_all => 'author build scripts are not part of the distribution'
	unless @present;

plan tests => scalar @present;

for my $spec (@present) {
	my $path = File::Spec->catfile( $root, $spec->{file} );
	open my $fh, '<', $path or die "cannot read $path: $!";
	my @lines = <$fh>;
	close $fh;

	my ( @found, @bad );
	for my $i ( 0 .. $#lines ) {
		my $line = $lines[$i];
		# Comment lines are skipped: compile.sh carries this file's own
		# rationale, which names the flag it warns against.
		next if $line =~ /^\s*#/;
		next unless $line =~ $spec->{regex};
		my $flags = $1;
		push @found, sprintf( 'line %d: %s', $i + 1, $flags );
		push @bad,   sprintf( 'line %d: %s', $i + 1, $flags )
			unless $flags =~ /(?:^|\s)-O/;
	}

	ok( @found && !@bad, "$spec->{file}: $spec->{what} enables the optimizer" )
		or diag(
		@bad
		? join( "\n",
			"OPTIMIZE= replaces perl's \$Config{optimize}; without an -O flag the",
			"module builds at -O0.  Found with no -O:", @bad )
		: "no match for $spec->{what} in $spec->{file} -- did it move?"
		);
}
