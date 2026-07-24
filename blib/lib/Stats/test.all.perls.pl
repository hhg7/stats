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

use 5.044;
no source::encoding;
use warnings FATAL => 'all';
use Getopt::Long 'GetOptions';
use File::Spec;
use IO::Handle;
use Cwd 'getcwd';
use POSIX 'strftime';
use Time::HiRes 'time';

my $PERLBREW_ROOT = $ENV{PERLBREW_ROOT} || File::Spec->catdir($ENV{HOME}, 'perl5', 'perlbrew');

my ($help, $list, $install, $deps, $clean, $stop, $quiet, $jobs, $log_dir, $optimize, @only);
$install  = 1; # make install, as compile.sh does
$clean    = 1;
$optimize = '-Wall';
$log_dir  = File::Spec->catdir('.build', 'multiperl');

GetOptions(
    'perl|p=s@'     => \@only,
    'install!'      => \$install,
    'deps!'         => \$deps,
    'clean!'        => \$clean,
    'stop-on-fail!' => \$stop,
    'jobs|j=i'      => \$jobs,
    'optimize=s'    => \$optimize,
    'log-dir=s'     => \$log_dir,
    'quiet|q'       => \$quiet,
    'list|l'        => \$list,
    'help|h'        => \$help,
) or usage(1);
usage(0) if $help;

sub usage {
    my $rc = shift;
    print STDERR <<"END";
usage: $0 [options]

Builds and tests the distribution in the current directory against each
perl installed under $PERLBREW_ROOT (oldest first).

options:
  -p, --perl VERSION   only this perl (repeatable); accepts "5.10.1" or
                       "perl-5.10.1".  default: every installed perl
  -l, --list           list the perls that would be tested, then exit
      --no-install     skip "make install" (build + test only)
      --no-clean       skip "make clean" before each version
      --deps           cpanm any missing PREREQ_PM for that perl first
      --stop-on-fail   abort at the first version that fails
  -j, --jobs N         parallel make, and HARNESS_OPTIONS=j<N> for the tests
      --optimize STR   OPTIMIZE= passed to Makefile.PL (default: $optimize)
      --log-dir DIR    where per-version logs go (default: $log_dir)
  -q, --quiet          only write logs; do not echo build output
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

# numeric sort so the oldest perl is exercised first and the tree is left
# built against the newest one.
sub vkey {
    my $v = shift;
    $v =~ s/^perl-//;
    my @p = ($v =~ /(\d+)/g);
    push @p, 0 while @p < 3;
    return sprintf '%05d%05d%05d', @p[0 .. 2];
}
@installed = sort { vkey($a) cmp vkey($b) } @installed;

my @targets = @installed;
if (@only) {
    my @want = map { my $x = $_; $x =~ s/^(?!perl-)/perl-/; $x } map { split /,/ } @only;
    my %have = map { $_ => 1 } @installed;
    my @missing = grep { !$have{$_} } @want;
    die "$0: not installed under perlbrew: @missing\n(installed: @installed)\n" if @missing;
    @targets = sort { vkey($a) cmp vkey($b) } @want;
}
die "$0: no perls found in $perls_dir\n" unless @targets;

if ($list) {
    printf "%-14s %s\n", $_, File::Spec->catfile($perls_dir, $_, 'bin', 'perl') for @targets;
    exit 0;
}

die "$0: no Makefile.PL in " . getcwd() . " - run this from the distribution root\n"
    unless -f 'Makefile.PL';

# ------------------------------------------------------------------ logging --

mkdirp($log_dir);
my $stamp = strftime '%Y%m%d-%H%M%S', localtime;

sub mkdirp {
    my @parts = File::Spec->splitdir(shift);
    my $path = '';
    for my $p (@parts) {
        $path = length($path) ? File::Spec->catdir($path, $p) : $p;
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
	my $code = $status == -1              ? -1
		    : ($status & 127)            ? 128 + ($status & 127)
		    :                              ($status >> 8);
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

# --------------------------------------------------------------- main loop --

my @results;
my $t_all = time;

for my $version (@targets) {
    my $root = File::Spec->catdir($perls_dir, $version);
    my $bin  = File::Spec->catdir($root, 'bin');
    my $perl = File::Spec->catfile($bin, 'perl');

    my $log = File::Spec->catfile($log_dir, "$version.$stamp.log");
    open my $logfh, '>', $log or die "$0: cannot write $log: $!\n";
    $logfh->autoflush(1);
    STDOUT->autoflush(1);

    print "\n", '=' x 72, "\n";
    printf "== %s   (log: %s)\n", $version, $log;
    print '=' x 72, "\n";
    print $logfh "== $version at " . strftime('%F %T', localtime) . "\n";

    # Emulate `perlbrew use $version` for the children: this perl's bin first,
    # every other perlbrew perl stripped out, and local::lib / PERL5LIB
    # leftovers from the calling shell removed so nothing bleeds across
    # versions.
    local %ENV = %ENV;
    my @path = grep { index($_, File::Spec->catdir($perls_dir, '')) != 0 }
               split /:/, ($ENV{PATH} || '/usr/bin:/bin');
    $ENV{PATH}          = join ':', $bin, @path;
    $ENV{PERLBREW_ROOT} = $PERLBREW_ROOT;
    $ENV{PERLBREW_PERL} = $version;
    $ENV{PERLBREW_PATH} = $bin;
    delete @ENV{qw(PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
                   PERLBREW_LIB PERL_MM_USE_DEFAULT)};
    $ENV{HARNESS_OPTIONS} = "j$jobs" if $jobs;

    my $t0 = time;
    my %r = (version => $version, log => $log, steps => [], warnings => 0);

    # sanity: the interpreter really is the version we think it is
    my ($vc, $vout) = run_cmd([$perl, '-e', 'printf "%vd\n", $^V'], $logfh);
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

    push @steps, ['Makefile.PL', [$perl, 'Makefile.PL', "OPTIMIZE=$optimize"]];
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

    if ($clean && !$failed) {
        # leave a clean tree behind before the next perl takes over
        run_cmd(['make', 'clean'], $logfh) unless $version eq $targets[-1];
    }

    $r{seconds} = time - $t0;
    $r{failed}  = $failed;
    close $logfh;

    printf "-- %s: %s in %.1fs%s\n", $version,
        ($failed ? "FAILED at '$failed'" : 'ok'),
        $r{seconds},
        (defined $r{files} ? " ($r{files})" : '');
    print '-- ', join('  ', map { sprintf '%s %.1fs', $_->{label}, $_->{seconds} }
                            @{ $r{steps} }), "\n";
    push @results, \%r;

    if ($failed && $stop) {
        my @rest = @targets[ scalar(@results) .. $#targets ];
        print "-- --stop-on-fail: skipping @rest\n" if @rest;
        last;
    }
}

# ----------------------------------------------------------------- summary --

my $bad = grep { $_->{failed} } @results;
print "\n", '=' x 72, "\n";
printf "%-14s %-9s %-8s %-7s %-6s %s\n",
    qw(PERL REPORTED STATUS TIME WARN TESTS);
print '-' x 72, "\n";
for my $r (@results) {
    printf "%-14s %-9s %-8s %6.1fs %-6s %s\n",
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

exit($bad || $skipped ? 1 : 0);
