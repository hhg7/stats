#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use XS::Check;

my $check = XS::Check->new ();
$check->check_file ('LikeR.xs');

# Perl_isnan/Perl_isinf/Perl_isfinite must never be called from LikeR.xs: on
# every perl before 5.22 they can expand to perl.h's Perl_fp_class() block,
# which is written with an empty parameter list and compares against FP_CLASS_*
# names <ieeefp.h> does not define.  Configure only reaches that block where it
# fails to find isinf(), so the breakage is invisible here and fatal on
# illumos/Solaris -- it is what stopped 0.298 building there.  nv_isnan(),
# nv_isinf() and nv_isfinite() in LikeR.xs are the replacements.  Comments may
# still name the macros; only a call is an error.
{
	open my $xs, '<', 'LikeR.xs';
	my @bad;
	while (my $line = <$xs>) {
		push @bad, "LikeR.xs:$.: $line"
			if $line =~ m/\bPerl_is(?:nan|inf|finite)\s*\(/;
	}
	close $xs;
	if (@bad) {
		print "Perl_is*() must be nv_is*() -- see the nv_isnan comment in LikeR.xs:\n";
		print for @bad;
		exit 1;
	}
	say 'no Perl_isnan/Perl_isinf/Perl_isfinite calls in LikeR.xs';
}

