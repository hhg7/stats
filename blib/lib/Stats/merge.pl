#!/usr/bin/env perl
# Exercises Stats::LikeR::merge across join types, key specs, and shapes,
# comparing the XS result against a pure-Perl reference implementation.
use 5.010;
use warnings FATAL => 'all';
use autodie ':default';
use Devel::Confess 'color';
use DDP {output => 'STDOUT', array_max => 20};
use Stats::LikeR;
use List::Util qw(any);

# ---- pure-Perl reference merge (AoH in, AoH out) ------------------------
sub ref_merge {
	my ($L, $R, %o) = @_;
	my $how = $o{how} // 'inner';
	my (@lk, @rk);
	if ($o{'left.on'} || $o{'right.on'}) {
		@lk = ref $o{'left.on'}  ? @{$o{'left.on'}}  : ($o{'left.on'});
		@rk = ref $o{'right.on'} ? @{$o{'right.on'}} : ($o{'right.on'});
	} elsif (defined $o{on} || defined $o{by}) {
		my $on = $o{on} // $o{by};
		@lk = ref $on ? @$on : ($on);
		@rk = @lk;
	} elsif ($how ne 'cross') {
		my %lseen; $lseen{$_}++ for map { keys %$_ } @$L;
		my %rseen; $rseen{$_}++ for map { keys %$_ } @$R;
		@lk = sort grep { $rseen{$_} } keys %lseen;
		@rk = @lk;
		die "ref_merge: no common columns\n" unless @lk;
	}
	my $nk = @lk;
	$nk = 0 if $how eq 'cross';

	my %lkset = map { $_ => 1 } @lk[0 .. $nk-1] if $nk;
	my %rkset = map { $_ => 1 } @rk[0 .. $nk-1] if $nk;
	my (@lc, %lc_seen); for my $row (@$L) { for my $c (sort keys %$row) { next if $lkset{$c}; push @lc, $c unless $lc_seen{$c}++; } }
	my (@rc, %rc_seen); for my $row (@$R) { for my $c (sort keys %$row) { next if $rkset{$c}; push @rc, $c unless $rc_seen{$c}++; } }
	my %lcset = map { $_ => 1 } @lc;
	my %rcset = map { $_ => 1 } @rc;
	my ($s0, $s1) = $o{suffixes} ? @{$o{suffixes}} : ('.x', '.y');
	my %lout = map { $_ => ($rcset{$_} ? "$_$s0" : $_) } @lc;
	my %rout = map { $_ => ($lcset{$_} ? "$_$s1" : $_) } @rc;

	my $keyf = sub {
		my $row = shift;
		my @p;
		for my $k (@_) { return undef unless defined $row->{$k}; push @p, length($row->{$k}) . "\x1e" . $row->{$k}; }
		join("\x1e", @p);
	};
	my $emit = sub {
		my ($li, $ri) = @_;
		my %row;
		for my $j (0 .. $nk-1) {
			my $v = (defined $li && defined $li->{$lk[$j]}) ? $li->{$lk[$j]} : (defined $ri ? $ri->{$rk[$j]} : undef);
			$row{$lk[$j]} = $v;
		}
		$row{$lout{$_}} = defined $li ? $li->{$_} : undef for @lc;
		$row{$rout{$_}} = defined $ri ? $ri->{$_} : undef for @rc;
		return \%row;
	};

	my @out;
	if ($how eq 'cross') {
		for my $li (@$L) { for my $ri (@$R) { push @out, $emit->($li, $ri); } }
		return \@out;
	}
	my %ridx;
	for my $j (0 .. $#$R) {
		my $k = $keyf->($R->[$j], @rk);
		next unless defined $k;
		push @{ $ridx{$k} }, $j;
	}
	my @matched;
	for my $li (@$L) {
		my $k = $keyf->($li, @lk);
		if (defined $k && $ridx{$k}) {
			for my $j (@{ $ridx{$k} }) { push @out, $emit->($li, $R->[$j]); $matched[$j] = 1; }
		} elsif ($how eq 'left' || $how eq 'outer' || $how eq 'full') {
			push @out, $emit->($li, undef);
		}
	}
	if ($how eq 'right' || $how eq 'outer' || $how eq 'full') {
		for my $j (0 .. $#$R) { push @out, $emit->(undef, $R->[$j]) unless $matched[$j]; }
	}
	return \@out;
}

# canonical signature of a result: multiset of sorted key=val strings per row.
# Accepts AoH or HoA (the two shapes merge can return).
sub sig {
	my $aoh = shift;
	if (ref $aoh eq 'HASH') {           # HoA -> AoH
		my @cols = keys %$aoh;
		my $n = @cols ? scalar @{ $aoh->{$cols[0]} } : 0;
		my @rows;
		for my $i (0 .. $n-1) { push @rows, { map { $_ => $aoh->{$_}[$i] } @cols }; }
		$aoh = \@rows;
	}
	my @sigs;
	for my $r (@$aoh) {
		push @sigs, join('|', map { "$_=" . (defined $r->{$_} ? $r->{$_} : 'UNDEF') } sort keys %$r);
	}
	return join("\n", sort @sigs);
}

my $pass = 0;
my $fail = 0;
sub check {
	my ($name, $got, $want) = @_;
	if (sig($got) eq sig($want)) { $pass++; say "ok   - $name"; }
	else {
		$fail++;
		say "FAIL - $name";
		say "  got : " . sig($got);
		say "  want: " . sig($want);
	}
}

# ---- test data ----------------------------------------------------------
my $emp = [
	{ id => 1, name => 'Alice',   dept => 10 },
	{ id => 2, name => 'Bob',     dept => 20 },
	{ id => 3, name => 'Carol',   dept => 30 },
	{ id => 4, name => 'Dave',    dept => 10 },
	{ id => 5, name => 'Eve',     dept => undef },
];
my $dept = [
	{ dept => 10, dname => 'Sales',       name => 'HQ' },
	{ dept => 20, dname => 'Engineering', name => 'Lab' },
	{ dept => 40, dname => 'Legal',       name => 'Annex' },
];

for my $how (qw(inner left right outer cross)) {
	my %o = $how eq 'cross' ? (how => 'cross') : (how => $how, on => 'dept');
	check("how=$how on=dept",
		Stats::LikeR::merge($emp, $dept, %o),
		ref_merge($emp, $dept, %o));
}

# natural join (no on) -> common cols {dept, name}
check("natural join",
	Stats::LikeR::merge($emp, $dept, how => 'inner'),
	ref_merge($emp, $dept, how => 'inner'));

# multi-key join
my $a = [ {k1=>1,k2=>'x',v=>'a'}, {k1=>1,k2=>'y',v=>'b'}, {k1=>2,k2=>'x',v=>'c'} ];
my $b = [ {k1=>1,k2=>'x',w=>'A'}, {k1=>2,k2=>'x',w=>'C'}, {k1=>2,k2=>'z',w=>'Z'} ];
for my $how (qw(inner left right outer)) {
	my %o = (how => $how, on => ['k1','k2']);
	check("multikey how=$how", Stats::LikeR::merge($a,$b,%o), ref_merge($a,$b,%o));
}

# left.on / right.on with different names
my $orders = [ {oid=>1, cust=>'c1'}, {oid=>2, cust=>'c2'}, {oid=>3, cust=>'c9'} ];
my $cust   = [ {cid=>'c1', city=>'NYC'}, {cid=>'c2', city=>'LA'} ];
for my $how (qw(inner left right outer)) {
	my %o = (how=>$how, 'left.on'=>'cust', 'right.on'=>'cid');
	check("left.on/right.on how=$how", Stats::LikeR::merge($orders,$cust,%o), ref_merge($orders,$cust,%o));
}

# suffixes
check("custom suffixes",
	Stats::LikeR::merge($emp,$dept, on=>'dept', suffixes=>['_emp','_dept']),
	ref_merge($emp,$dept, on=>'dept', suffixes=>['_emp','_dept']));

# many-to-many
my $m1 = [ {k=>1,l=>'a'}, {k=>1,l=>'b'} ];
my $m2 = [ {k=>1,r=>'X'}, {k=>1,r=>'Y'} ];
check("many-to-many inner", Stats::LikeR::merge($m1,$m2,on=>'k'), ref_merge($m1,$m2,on=>'k'));

# HoA input
my $hoaL = { id=>[1,2,3], grp=>['a','b','a'] };
my $hoaR = { grp=>['a','b'], score=>[100,200] };
check("HoA inputs -> HoA out",
	Stats::LikeR::merge($hoaL,$hoaR, on=>'grp', how=>'left'),
	ref_merge([map { {id=>$hoaL->{id}[$_], grp=>$hoaL->{grp}[$_]} } 0..2],
	          [map { {grp=>$hoaR->{grp}[$_], score=>$hoaR->{score}[$_]} } 0..1], on=>'grp', how=>'left'));

# output.type override
check("output.type hoa",
	Stats::LikeR::merge($emp,$dept, on=>'dept', how=>'inner', 'output.type'=>'hoa'),
	ref_merge($emp,$dept, on=>'dept', how=>'inner'));

say "";
say "pass=$pass fail=$fail";
exit($fail ? 1 : 0);
