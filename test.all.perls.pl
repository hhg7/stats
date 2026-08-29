#!/usr/bin/env perl
# test_all_perls.pl - build & test Stats::LikeR against every perlbrew perl.
#
# Automates the manual loop of:
#     perlbrew use perl-5.10.1 && ./compile.sh
#     perlbrew use perl-5.12.5 && ./compile.sh
#     perlbrew use perl-5.44.0 && ./compile.sh
#
# Each version is built in its own environment (no `perlbrew use` needed: the
# perl binary is invoked directly and PATH is rewritten for the child), the
# tree is cleaned between versions, and a pass/fail summary is printed at the
# end.  Exit status is non-zero if any version failed.
#
# By default four perls are built at once (-P 4).  They cannot share the
# distribution root: there is one Makefile, one LikeR.c, one LikeR.o and one
# blib, and a LikeR.o built by another perl links without complaint and then
# dies at load.  So each parallel child builds in a private copy of the tree
# under the log directory and reports its result back to the parent, which
# leaves the root's Makefile and blib exactly as they were -- the same reason
# the NV-width compile checks avoid reconfiguring the tree.  Nothing forces
# Parallel::ForkManager here: the fork/throttle/reap loop is a few dozen lines
# and this script already forks in run_cmd, so the helper keeps needing nothing
# but core modules.
#
# -P 1 is the old behaviour: one perl at a time, in the distribution root,
# which is left built and installed against the perl that runs last -- the
# newest plain one, see the run order below.

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use Getopt::Long 'GetOptions';
# -P (how many perls at once) and -p (which perl) are different options, so the
# default case-folding of single-letter aliases has to go.
Getopt::Long::Configure('no_ignore_case');
use File::Spec;
use File::Temp ();
use File::Copy 'copy';
use File::Path 'remove_tree';
use IO::Handle;
use Cwd 'getcwd';
use Data::Dumper;
use POSIX 'strftime';
use Time::HiRes 'time';

my $PERLBREW_ROOT = $ENV{PERLBREW_ROOT} || File::Spec->catdir($ENV{HOME}, 'perl5', 'perlbrew');

# How many cores this machine has, or 0 when that cannot be determined.  Only
# an author-side helper runs this, but the probes stay portable so a run on a
# BSD or Solaris box degrades to a sensible default instead of dying.
sub cpu_count {
	if (open my $fh, '<', '/proc/cpuinfo') {           # linux
		my $n = grep { /^processor\s*:/ } <$fh>;
		close $fh;
		return $n if $n;
	}
	for my $probe (['getconf', '_NPROCESSORS_ONLN'],   # glibc, solaris, aix
		['sysctl', '-n', 'hw.ncpu']) {             # *bsd, darwin
		my $out = `@$probe 2>/dev/null`;
		return $1 if defined $out && $out =~ /(\d+)/ && $1 > 0;
	}
	return 0;
}

my ($help, $list, $install, $deps, $clean, $stop, $quiet, $jobs, $log_dir, $optimize, @only);
my ($par, $keep_work);
# On by default: this is the axis that let 0.31 reach CPAN broken, and an
# opt-in check is one nobody remembers to opt in to.  --no-x87 skips it.
my $x87 = 1;
my $ncpu  = cpu_count();
$install  = 1; # make install, as compile.sh does
$clean    = 1;
# OPTIMIZE= replaces perl's own $Config{optimize} rather than adding to it, so
# a bare '-Wall' here would build and time every perl in the matrix at -O0.
$optimize = '-O2 -Wall';
$log_dir  = File::Spec->catdir('.build', 'multiperl');
# Parallel by default, leaving 4 cores for the rest of the machine: a full
# matrix run is minutes of `make` and a test harness, not something that should
# need remembering to ask for.
my $budget = $ncpu > 5 ? $ncpu - 4 : 1;
# That budget is now split two ways: several perls at once, each with its own
# parallel make and harness.  Four perls is the whole local matrix bar one, and
# -P x -j multiply, so -j follows from -P rather than claiming the budget twice.
$par = 4;
my $par_default  = $par;
my $jobs_default = jobs_for($par);

sub jobs_for {
	my $p = shift;
	my $j = int($budget / ($p || 1));
	return $j > 1 ? $j : 1;
}

GetOptions(
	'perl|p=s@'     => \@only,
	'install!'      => \$install,
	'deps!'         => \$deps,
	'clean!'        => \$clean,
	'stop-on-fail!' => \$stop,
	'jobs|j=i'      => \$jobs,
	'parallel|P:i'  => \$par,
	'keep-work!'    => \$keep_work,
	'optimize=s'    => \$optimize,
	'x87!'          => \$x87,
	'log-dir=s'     => \$log_dir,
	'quiet|q'       => \$quiet,
	'list|l'        => \$list,
	'help|h'        => \$help,
) or usage(1);
usage(0) if $help;

sub usage {
	my $rc = shift;
	my $cpus = $ncpu ? "$ncpu CPUs - 4" : 'CPU count unknown';
	print STDERR <<"END";
usage: $0 [options]

Builds and tests the distribution in the current directory against each
perl installed under $PERLBREW_ROOT.

Run order: the quadmath build first, since its test suite is minutes where the
others are seconds and it must not sit in the queue behind them; then the rest
oldest first; then the newest plain perl (double NV, no ithreads) last, so a
serial run leaves the tree built and installed against the reference build
rather than an outlier.

One extra row, "<perl>+x87", rebuilds the newest plain perl with -mfpmath=387.
The perls above differ in how wide an NV is stored; they are all x86-64, so they
agree on doing the arithmetic in SSE at exactly that width.  32-bit x86 does not:
it computes in 80-bit x87 registers and rounds on a spill, which is how 0.31
reached CPAN with 73 subtests failing on an i686 smoker that nothing here could
reproduce.  The row is skipped where there is nothing to catch (a long-double or
quadmath perl, or a cc with no -mfpmath) and by --no-x87.

options:
  -p, --perl VERSION   only this perl (repeatable); accepts "5.10.1",
                       "perl-5.10.1" or an exact directory name such as
                       "5.44.0-quadmath".  default: every installed perl
  -l, --list           list the perls that would be tested, in run order, with
                       each one's NV width and threading, then exit
      --no-install     skip "make install" (build + test only)
      --no-clean       skip "make clean" before each version
      --deps           cpanm any missing PREREQ_PM for that perl first
      --stop-on-fail   abort at the first version that fails (with -P: launch
                       no more perls; the running ones finish)
  -P, --parallel [N]   build and test N perls at once, each in a private copy of
                       the tree under --log-dir/work, leaving the distribution
                       root untouched.  bare -P, or -P 0, uses one child per
                       perl, capped at the CPU count; -P 1 builds serially in
                       the distribution root and leaves it built against the
                       newest plain perl.  default: $par_default
      --keep-work      keep the private tree of a perl that passed (the tree of
                       one that failed is always kept, for inspection)
  -j, --jobs N         parallel make, and HARNESS_OPTIONS=j<N> for the tests;
                       0 builds and tests serially.  -P and -j multiply, so the
                       default is the CPU budget ($cpus) divided by -P
                       (default: $jobs_default at -P $par_default, $budget at -P 1)
      --optimize STR   OPTIMIZE= passed to Makefile.PL (default: $optimize)
      --no-x87         skip the extra x87 row (see below)
      --log-dir DIR    where per-version logs go (default: $log_dir)
  -q, --quiet          only write logs; do not echo build output.  implied by
                       -P > 1, where interleaved output would be unreadable
  -h, --help           this message

exit status: 0 if every perl built, tested and installed cleanly, else 1.
END
	exit $rc;
}

# ---------------------------------------------------------------- discovery --

my $perls_dir = File::Spec->catdir($PERLBREW_ROOT, 'perls');
die "$0: no perlbrew perls directory at $perls_dir\n" unless -d $perls_dir;

opendir my $dh, $perls_dir or die "$0: cannot read $perls_dir: $!\n";
my @installed = grep { -x File::Spec->catfile($perls_dir, $_, 'bin', 'perl') }
	grep { !/^\.\.?$/ } readdir $dh;
closedir $dh;

# numeric sort, oldest first, which is where the run order below starts from.
sub vkey {
	my $v = shift;
	my $variant = tgt_variant($v);
	($v = tgt_perl($v)) =~ s/^perl-//;
	my @p = ($v =~ /(\d+)/g);
	push @p, 0 while @p < 3;
	# the variant's key trails its base perl's, so it sorts next to the perl it
	# varies instead of landing between two versions
	return sprintf('%05d%05d%05d', @p[0 .. 2]) . ($variant ? "1$variant" : '0');
}
@installed = sort { vkey($a) cmp vkey($b) } @installed;

my @targets = @installed;
if (@only) {
	my %have = map { $_ => 1 } @installed;
	my (@want, @missing);
	for my $arg (map { split /,/ } @only) {
		# an exact directory name wins, so the quadmath build -- which is
		# installed as "5.44.0-quadmath", with no perl- prefix, and is one
		# of the three NV widths the matrix has to cover -- is selectable.
		my ($match) = grep { $have{$_} } $arg, "perl-$arg";
		if (defined $match) { push @want, $match }
		else                { push @missing, $arg }
	}
	die "$0: not installed under perlbrew: @missing\n(installed: @installed)\n" if @missing;
	my %seen;
	@targets = sort { vkey($a) cmp vkey($b) } grep { !$seen{$_}++ } @want;
}
die "$0: no perls found in $perls_dir\n" unless @targets;

# ------------------------------------------------------------------ run order --

# What each perl actually is.  The NV width and the threading model decide the
# order below, and the answer has to come from the interpreter: "5.44.0-quadmath"
# is a local naming habit rather than a promise, and the long-double build is
# only identifiable as "perl-5.12.5" plus an archname nobody parses.  A perl that
# will not answer is described as unknown, never as plain, so a broken
# interpreter cannot become the reference build the tree is left standing on.
my %facts;
sub facts {
	my $version = tgt_perl(shift);   # a variant is the same interpreter
	return $facts{$version} if $facts{$version};
	my $perl = File::Spec->catfile($perls_dir, $version, 'bin', 'perl');
	my %f = (nv => '?', threads => 0, known => 0);
	if (open my $fh, '-|', $perl, '-MConfig', '-e',
			'print "$Config{nvtype}\t", ($Config{useithreads} ? 1 : 0)') {
		my $line = <$fh>;
		close $fh;
		if (defined $line && $line =~ /^(\S[^\t]*)\t([01])/) {
			%f = (nv => $1, threads => $2, known => 1);
		}
	}
	$f{quadmath} = $f{nv} eq '__float128';
	# plain: the double-NV, unthreaded build that the long-double, quadmath and
	# threaded configurations are each a variation on.  Unknown does not count.
	$f{plain} = $f{known} && $f{nv} eq 'double' && !$f{threads};
	return $facts{$version} = \%f;
}

sub nv_label {
	my $f = shift;
	my %short = ('double' => 'double', 'long double' => 'long-double',
		'__float128' => 'quadmath');
	return ($short{ $f->{nv} } || $f->{nv}) . ($f->{threads} ? '-thr' : '');
}

# A target is a perl, optionally with a build variant after a '+'.  Everything
# that names a target -- its log, its private tree, its result file, its row in
# the summary -- uses the whole label, so a variant is a row of its own rather
# than a footnote on somebody else's; only the paths into perlbrew and the
# OPTIMIZE= string have to look past the suffix.
sub tgt_perl    { my $t = shift; $t =~ s/\+.*//; return $t }
sub tgt_variant { my $t = shift; return $t =~ /\+(.*)/ ? $1 : '' }
sub tgt_optimize {
	my $t = shift;
	# appended, not replacing: a --optimize the caller asked for still applies,
	# and -mfpmath=387 only decides which register file the arithmetic uses.
	return tgt_variant($t) eq 'x87' ? "$optimize -mfpmath=387" : $optimize;
}

# Run order, given @targets oldest first:
#
#   1. the quadmath builds, because their test suite is far and away the
#      slowest -- 279s against 8-22s for every other perl here, since a
#      __float128 libm in an inner loop costs what it costs -- and the matrix is
#      only as fast as the slowest perl that had to wait for a free slot;
#   2. everything else, still oldest first;
#   3. the newest plain perl last, because a serial run leaves the distribution
#      root built and installed against whichever perl went last, and that
#      should be the ordinary double-NV reference build rather than a
#      long-double, threaded or quadmath outlier.
sub order_targets {
	my @t = @_;
	my (@quad, @rest);
	push @{ facts($_)->{quadmath} ? \@quad : \@rest }, $_ for @t;
	# !tgt_variant: a serial run leaves the root built against whatever went
	# last, and that must be an ordinary build, not one forced onto the x87
	# stack by the row below.
	my ($plain) = grep { facts($_)->{plain} && !tgt_variant($_) } reverse @rest;
	@rest = grep { $_ ne $plain } @rest if defined $plain;
	return (@quad, @rest, defined $plain ? $plain : ());
}
# ----------------------------------------------------------- x87 variant --

# The matrix above varies the width an NV is stored in.  It does not vary the
# width the arithmetic is *done* in, and every perl in it is x86-64, where a
# double is computed in an SSE register and FLT_EVAL_METHOD is 0.  On 32-bit
# x86 it is 2: doubles are computed in the x87 register file at 80 bits and
# rounded only when they are spilled, so ceil(0.1 * 100) is 11 rather than 10
# and any product split into an integer part and a remainder can land in a
# different bucket.  That is what took Stats-LikeR 0.31 to 73 failing subtests
# on a CPAN smoker (i686-linux-64int, perl 5.32.1) across bedroc(), quantile(),
# var() and wilcox_test(), with nothing here able to see it.  So -mfpmath=387
# gets a row: it reproduced 63 of those 73, by test number, before the fix.
#
# One row is enough, and it goes on the newest plain perl.  The axis is the
# compiler's evaluation mode, not anything perl chooses, so a second x87 row on
# another double perl would hand the same source to the same cc with the same
# flags.  It is skipped on the long-double and quadmath perls because there is
# nothing there to catch -- an x87 register is exactly as wide as a long-double
# NV, and __float128 is not computed on the x87 stack at all -- and skipping
# quadmath in particular keeps a 290s no-op out of every run.
#
# What this row cannot see is the other half of 32-bit x86: i386 returns a
# double in st(0), so a callee hands its caller 80 bits, where x86-64 returns
# in xmm0 and narrows for free.  That is what moved wilcox_test()'s pseudomedian
# on the smoker, and only a real 32-bit perl will show it.

# Does this perl's cc take -mfpmath=387?  It is an x86-only gcc/clang option, so
# on any other machine or toolchain the row is simply absent -- there is no x87
# register file to force the arithmetic onto.  Both compilers exit non-zero on
# an -m switch they do not implement.
sub cc_takes_x87 {
	my $version = shift;
	my $perl = File::Spec->catfile($perls_dir, tgt_perl($version), 'bin', 'perl');
	my $cc;
	if (open my $fh, '-|', $perl, '-MConfig', '-e', 'print $Config{cc}') {
		$cc = <$fh>;
		close $fh;
	}
	return 0 unless defined $cc && length $cc;
	my $dir = File::Temp::tempdir(CLEANUP => 1);
	my $src = File::Spec->catfile($dir, 'x87probe.c');
	my $obj = File::Spec->catfile($dir, 'x87probe.o');
	open my $out, '>', $src or return 0;
	print {$out} "int main(void) { return 0; }\n";
	close $out;
	my $devnull = File::Spec->devnull;
	return system("$cc -mfpmath=387 -c $src -o $obj > $devnull 2> $devnull") == 0;
}

if ($x87) {
	my ($base) = grep { facts($_)->{plain} } reverse @targets;   # newest first
	if    (!defined $base)        { print "-- no plain (double NV) perl selected; no x87 row\n" }
	elsif (!cc_takes_x87($base))  { print "-- cc does not take -mfpmath=387 (not x86?); no x87 row\n" }
	else                          { push @targets, "$base+x87" }
}

@targets = order_targets(@targets);

if ($list) {
	printf "%-20s %-16s %s\n", $_,
		nv_label(facts($_)) . (tgt_variant($_) ? ' +' . tgt_variant($_) : ''),
		File::Spec->catfile($perls_dir, tgt_perl($_), 'bin', 'perl') for @targets;
	exit 0;
}

die "$0: no Makefile.PL in " . getcwd() . " - run this from the distribution root\n"
	unless -f 'Makefile.PL';

# ------------------------------------------------------------------ logging --

mkdirp($log_dir);
my $stamp = strftime '%Y%m%d-%H%M%S', localtime;

# ------------------------------------------------------------- parallelism --

# bare -P (or -P 0) means "all of them", within reason: more compilers than
# cores only makes every perl slower.
$par = 0 if $par < 0;
$par ||= do {
	my $cap = $ncpu || 4;   # an unknown CPU count stays conservative rather
	                        # than fanning out over however many perls exist
	@targets < $cap ? scalar @targets : $cap;
};
$par = @targets if $par > @targets;
$jobs = jobs_for($par) unless defined $jobs;   # -j follows -P unless asked for

my $work_root   = File::Spec->catdir($log_dir, 'work');
my $in_parallel = 0;   # set by the dispatch below, once it knows there is
                       # more than one perl to run
if ($par > 1) {
	# children chdir into their own tree, so log paths must not be relative
	$log_dir   = File::Spec->rel2abs($log_dir);
	$work_root = File::Spec->rel2abs($work_root);
}

sub mkdirp {
	my @parts = File::Spec->splitdir(shift);
	my $path;   # undef, not '': splitdir gives an absolute path a leading '',
	            # and catdir('', 'home') is '/home' where 'home' would be a
	            # directory of that name in the cwd.
	for my $p (@parts) {
		$path = defined $path ? File::Spec->catdir($path, $p) : $p;
		next if !length($path) || -d $path;
		mkdir $path or die "$0: mkdir $path: $!\n";
	}
}

# ------------------------------------------------------------- child runner --

# Run @cmd with STDERR folded into STDOUT, echoing to the terminal and to
# $logfh.  Returns ($exit_code, \@lines).
sub run_cmd {
	my ($cmd, $logfh) = @_;
	print $logfh "\n\$ @$cmd\n";
	print "\$ @$cmd\n" unless $quiet;

	my $pid = open my $fh, '-|';
	die "$0: fork failed: $!\n" unless defined $pid;
	if (!$pid) {                              # child
		open STDERR, '>&', \*STDOUT or die "$0: dup STDERR: $!\n";
		$| = 1;
		{ exec { $cmd->[0] } @$cmd; }
		print "exec @$cmd failed: $!\n";
		exit 127;
	}

	my @lines;
	while (defined(my $line = <$fh>)) {
		push @lines, $line;
		print $logfh $line;
		print $line unless $quiet;
	}
	close $fh;
	my $status = $?;
	my $code = $status == -1  ? -1
		: ($status & 127) ? 128 + ($status & 127)
		: ($status >> 8);
	return ($code, \@lines);
}

# PREREQ_PM names, scraped from Makefile.PL so this stays in sync with dist.ini.
sub prereqs {
	open my $fh, '<', 'Makefile.PL' or return ();
	local $/;
	my $src = <$fh>;
	close $fh;
	return () unless $src =~ /"PREREQ_PM"\s*=>\s*\{(.*?)\}/s;
	my $block = $1;
	my %seen;
	return grep { !$seen{$_}++ } ($block =~ /"([\w:]+)"\s*=>/g);
}
my @prereqs = prereqs();

# ------------------------------------------------------------ build one perl --

# Build, test and install the distribution in the current directory with one
# perl.  Returns the result hashref, and prints its own progress unless
# $silent: a parallel child is silent and the parent reports for it.
sub build_one {
	my ($version, $silent) = @_;
	my $root = File::Spec->catdir($perls_dir, tgt_perl($version));
	my $bin  = File::Spec->catdir($root, 'bin');
	my $perl = File::Spec->catfile($bin, 'perl');

	my $log = File::Spec->catfile($log_dir, "$version.$stamp.log");
	open my $logfh, '>', $log or die "$0: cannot write $log: $!\n";
	$logfh->autoflush(1);
	STDOUT->autoflush(1);

	unless ($silent) {
		print "\n", '=' x 72, "\n";
		printf "== %s   (log: %s)\n", $version, $log;
		print '=' x 72, "\n";
	}
	print $logfh "== $version at " . strftime('%F %T', localtime)
		. ' in ' . getcwd() . "\n";

	# Emulate `perlbrew use $version` for the children: this perl's bin first,
	# every other perlbrew perl stripped out, and local::lib / PERL5LIB
	# leftovers from the calling shell removed so nothing bleeds across
	# versions.
	local %ENV = %ENV;
	my @path = grep { index($_, File::Spec->catdir($perls_dir, '')) != 0 }
		split /:/, ($ENV{PATH} || '/usr/bin:/bin');
	$ENV{PATH}          = join ':', $bin, @path;
	$ENV{PERLBREW_ROOT} = $PERLBREW_ROOT;
	$ENV{PERLBREW_PERL} = tgt_perl($version);
	$ENV{PERLBREW_PATH} = $bin;
	delete @ENV{qw(PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
		PERLBREW_LIB PERL_MM_USE_DEFAULT)};
	$ENV{HARNESS_OPTIONS} = "j$jobs" if $jobs;

	my $t0 = time;
	my %r = (version => $version, log => $log, steps => [], warnings => 0);

	# sanity: the interpreter really is the version we think it is, and say so
	# when it is a threaded build -- an XS bug can be invisible on an
	# unthreaded perl and a hard compile error under MULTIPLICITY, so the
	# summary has to make that coverage visible rather than implied by a name.
	my ($vc, $vout) = run_cmd([$perl, '-MConfig', '-e',
		'printf "%vd%s\n", $^V, $Config{useithreads} ? "-thr" : ""'], $logfh);
	$r{reported} = $vc == 0 && @$vout ? do { my $s = $vout->[0]; chomp $s; $s } : '?';

	my @steps;
	if ($deps && @prereqs) {
		my @missing;
		for my $mod (@prereqs) {
			my ($c) = run_cmd([$perl, "-M$mod", '-e', '1'], $logfh);
			push @missing, $mod if $c != 0;
		}
		if (@missing) {
			my $cpanm = -x File::Spec->catfile($bin, 'cpanm')
				? File::Spec->catfile($bin, 'cpanm')
				: File::Spec->catfile($PERLBREW_ROOT, 'bin', 'cpanm');
			push @steps, ['deps', [$perl, $cpanm, '--notest', @missing]];
		}
	}

	if ($clean && -f 'Makefile') {
		push @steps, ['clean', ['make', 'clean'], 1];   # 1 = failure tolerated
	}

	push @steps, ['Makefile.PL',
		[$perl, 'Makefile.PL', 'OPTIMIZE=' . tgt_optimize($version)]];
	push @steps, ['make',        ['make', $jobs ? ("-j$jobs") : ()]];
	push @steps, ['make test',   ['make', 'test']];
	push @steps, ['make install',['make', 'install']] if $install;

	my $failed;
	for my $step (@steps) {
		my ($label, $cmd, $soft) = @$step;
		my $t_step = time;
		my ($code, $lines) = run_cmd($cmd, $logfh);
		my $secs = time - $t_step;
		printf $logfh "-- step '%s' exited %d after %.1fs\n", $label, $code, $secs;

		if ($label eq 'make test') {
			for my $l (@$lines) {
				$r{files}  = $1 if $l =~ /^(Files=\d+.*)/;
				$r{result} = $1 if $l =~ /^Result:\s*(\S+)/;
				$r{passed} = 1  if $l =~ /^All tests successful/;
			}
		}
		# residual build warnings are worth surfacing even on success
		$r{warnings} += grep { /: warning:/ } @$lines if $label eq 'make';

		push @{ $r{steps} }, { label => $label, code => $code, seconds => $secs };
		next if $code == 0 || $soft;
		$failed = $label;
		last;
	}

	if ($clean && !$failed && !$in_parallel) {
		# leave a clean tree behind before the next perl takes over.  in
		# parallel mode the whole private tree goes instead, so this would
		# only be a `make clean` nobody reads.
		run_cmd(['make', 'clean'], $logfh) unless $version eq $targets[-1];
	}

	$r{seconds} = time - $t0;
	$r{failed}  = $failed;
	close $logfh;

	report_one(\%r) unless $silent;
	return \%r;
}

# the two progress lines a finished perl prints, from the parent in either mode
sub report_one {
	my $r = shift;
	printf "-- %s: %s in %.1fs%s\n", $r->{version},
		($r->{failed} ? "FAILED at '$r->{failed}'" : 'ok'),
		$r->{seconds},
		(defined $r->{files} ? " ($r->{files})" : '');
	print '-- ', join('  ', map { sprintf '%s %.1fs', $_->{label}, $_->{seconds} }
		@{ $r->{steps} }), "\n";
}

# ------------------------------------------------------------ private trees --

# What must not be copied into a private build tree: the build products of
# whichever perl last used the source tree (copying them would recreate exactly
# the stale-LikeR.o hazard the private trees exist to avoid), the repository,
# the coverage database, a dzil build directory, and the log/work tree itself.
# The file patterns are MANIFEST.SKIP's, plus the gcov and dzil leftovers that
# a release archive never sees.
my %SKIP_DIR  = map { $_ => 1 } qw(.git .build blib cover_db);
my @SKIP_FILE = (qr/\.(?:o|a|so|bs|dylib|dll)$/, qr/^LikeR\.(?:c|def)$/,
	qr/^Makefile(?:\.old)?$/, qr/^pm_to_blib$/, qr/^blibdirs$/,
	qr/^MYMETA\./, qr/\.gc(?:da|no|ov)$/, qr/^Stats-LikeR-.*\.tar\.gz$/);

sub copy_tree {
	my ($src, $dst) = @_;
	mkdirp($dst);
	opendir my $dh, $src or die "$0: cannot read $src: $!\n";
	my @entries = grep { !/^\.\.?$/ } readdir $dh;
	closedir $dh;
	for my $e (@entries) {
		my $s = File::Spec->catfile($src, $e);
		my $d = File::Spec->catfile($dst, $e);
		if (-d $s) {
			next if $SKIP_DIR{$e} || $e =~ /^Stats-LikeR-[\d.]+$/;
			# a --log-dir inside the tree would copy itself forever
			my $abs = File::Spec->rel2abs($s);
			next if $abs eq $work_root || $abs eq $log_dir;
			copy_tree($s, $d);
			next;
		}
		next if grep { $e =~ $_ } @SKIP_FILE;
		copy($s, $d) or die "$0: copy $s -> $d: $!\n";
		chmod((stat $s)[2] & 07777, $d);
	}
}

# The child's result has to cross a fork, and %r is plain data, so Dumper out /
# eval in beats any IPC here: it also leaves the numbers next to the log when
# something needs explaining afterwards.
sub write_result {
	my ($file, $r) = @_;
	open my $fh, '>', $file or die "$0: cannot write $file: $!\n";
	local $Data::Dumper::Indent   = 0;
	local $Data::Dumper::Sortkeys = 1;
	print $fh Data::Dumper->Dump([$r], ['R']);
	close $fh or die "$0: close $file: $!\n";
}

sub read_result {
	my $file = shift;
	open my $fh, '<', $file or return undef;
	local $/;
	my $src = <$fh>;
	close $fh;
	my $R;
	eval $src;                      # our own Dumper output, nobody else's
	return ref $R eq 'HASH' ? $R : undef;
}

# ------------------------------------------------------------------ drivers --

sub run_serial {
	my @queue = @_;
	my @out;
	for my $version (@queue) {
		my $r = build_one($version);
		push @out, $r;
		next unless $r->{failed} && $stop;
		my %done = map { $_->{version} => 1 } @out;
		my @rest = grep { !$done{$_} } @queue;
		print "-- --stop-on-fail: skipping @rest\n" if @rest;
		last;
	}
	return @out;
}

# Fork up to $par children, each in its own copy of the tree, reaping them as
# they finish and starting the next perl in the freed slot.
sub run_parallel {
	my @order = @_;
	my @queue = @order;
	my $src   = getcwd();
	my %kid;                        # pid => { version, work, result, log }
	my (@out, $halt);

	my $reaper = sub {
		my ($pid, $status) = @_;
		my $kid = delete $kid{$pid} or return;
		my $r   = read_result($kid->{result});
		if (!$r) {
			# the child died without reporting: exec failure, signal, OOM
			$r = { version => $kid->{version}, log => $kid->{log},
				reported => '?', steps => [], warnings => 0, seconds => 0,
				failed => sprintf('child exited %d%s', $status >> 8,
					($status & 127) ? ' on signal ' . ($status & 127) : '') };
		}
		push @out, $r;
		report_one($r);
		$halt = 1 if $r->{failed} && $stop;
		if ($r->{failed} || $keep_work) {
			print "-- $kid->{version}: build tree kept at $kid->{work}\n";
		}
		else {
			remove_tree($kid->{work});
		}
	};

	local $SIG{INT} = local $SIG{TERM} = sub {
		my $sig = shift;
		print "\n-- $sig: stopping " . keys(%kid) . " running build(s)\n";
		# each child leads its own process group (see the fork below), so one
		# signal per group takes its make, its compilers and its harness with
		# it; signalling the child perl alone would orphan those.
		kill 'TERM', map { -$_ } keys %kid;
		sleep 1;
		kill 'KILL', map { -$_ } keys %kid;
		exit 130;
	};

	while (@queue || %kid) {
		while (@queue && keys(%kid) < $par && !$halt) {
			my $version = shift @queue;
			my $work    = File::Spec->catdir($work_root, "$version.$stamp");
			my $result  = File::Spec->catfile($log_dir, "$version.$stamp.result");
			my $log     = File::Spec->catfile($log_dir, "$version.$stamp.log");

			remove_tree($work) if -d $work;
			printf "-- %-16s starting (tree: %s)\n", $version, $work;
			copy_tree($src, $work);

			my $pid = fork;
			die "$0: fork failed: $!\n" unless defined $pid;
			if (!$pid) {                              # child
				$SIG{$_} = 'DEFAULT' for qw(INT TERM);
				setpgrp 0, 0;   # so ^C reaches this whole build, once, via
						# the parent's handler and not the terminal
				chdir $work or die "$0: chdir $work: $!\n";
				my $r = build_one($version, 1);
				write_result($result, $r);
				exit 0;
			}
			$kid{$pid} = { version => $version, work => $work,
				result => $result, log => $log };
		}
		last unless %kid;
		my $pid = waitpid -1, 0;
		last if $pid <= 0;
		$reaper->($pid, $?);
	}

	print "-- --stop-on-fail: skipping @queue\n" if $halt && @queue;
	# report in version order, not in the order they happened to finish
	my %by = map { $_->{version} => $_ } @out;
	return map { $by{$_} } grep { $by{$_} } @order;
}

# --------------------------------------------------------------- main loop --

my @results;
my $t_all = time;

if ($par > 1 && @targets > 1) {
	$in_parallel = 1;
	$quiet       = 1;   # N interleaved build logs on one terminal is noise
	printf "-- %d perl(s), %d at a time%s; per-version output goes to the logs\n",
		scalar @targets, $par, ($jobs ? " with make -j$jobs each" : '');
	print "-- note: --deps has the children sharing one ~/.cpanm; if a "
		. "prerequisite install misbehaves, run once with -P 1 --deps first\n"
		if $deps;
	@results = run_parallel(@targets);
}
else {
	@results = run_serial(@targets);
}

# ----------------------------------------------------------------- summary --

my $bad = grep { $_->{failed} } @results;
print "\n", '=' x 72, "\n";
printf "%-16s %-11s %-8s %-7s %-6s %s\n",
	qw(PERL REPORTED STATUS TIME WARN TESTS);
print '-' x 72, "\n";
for my $r (@results) {
	printf "%-16s %-11s %-8s %6.1fs %-6s %s\n",
		$r->{version},
		$r->{reported},
		($r->{failed} ? 'FAIL' : 'PASS'),
		$r->{seconds},
		($r->{warnings} || 0),
		($r->{failed} ? "failed at '$r->{failed}' - see $r->{log}"
			: ($r->{result} ? "Result: $r->{result}" : 'no test summary')),
		;
}
print '-' x 72, "\n";
my $skipped = @targets - @results;
printf "%d/%d perl(s) passed%s in %.1fs.  Logs in %s\n",
	scalar(@results) - $bad, scalar(@targets),
	($skipped ? " ($skipped not run)" : ''),
	time - $t_all, $log_dir;
# every perl installed into its own site_perl, but nothing was built here: say
# so, because a serial run leaves the root built against the newest plain perl
# and the
# NV-width checks assume they know what state the tree is in.
print "-- built in private trees; this directory is untouched (-P 1 to build here)\n"
	if $in_parallel;

exit($bad || $skipped ? 1 : 0);
