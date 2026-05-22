#!/usr/bin/perl
use strict;
use warnings;
# use MyStats;   # your XS module

# ══════════════════════════════════════════════════════════════════════════
# MODE 1 – hash of groups (original behaviour)
#
# Keys   = group labels
# Values = array refs of numeric observations
# ══════════════════════════════════════════════════════════════════════════

my %groups = (
    yield => [5.5, 5.4, 5.8],
    ctrl  => [4.5, 4.8, 4.2],
);

my $r1 = oneway_test(\%groups);                    # Welch (default)
my $r2 = oneway_test(\%groups, var_equal => 1);    # classic ANOVA

print "── Mode 1: hash of groups ──────────────────────────────\n";
printf "Welch  F=%.4f  df=(%.4f, %.4f)  p=%.6f\n",
    $r1->{statistic}, $r1->{num_df}, $r1->{denom_df}, $r1->{p_value};
printf "Classic F=%.4f  df=(%.0f, %.0f)  p=%.6f\n",
    $r2->{statistic}, $r2->{num_df}, $r2->{denom_df}, $r2->{p_value};

# ══════════════════════════════════════════════════════════════════════════
# MODE 2 – formula "response ~ factor"
#
# Mirrors R's:
#   my_data <- stack(list(yield = yield, ctrl = ctrl))
#   oneway.test(Value ~ Group, data = my_data)
#
# The hash holds two parallel arrays of the same length.
# The factor array may contain numbers or strings — they are compared
# as strings so  0 / "0" / 0.0  are all the same group.
# ══════════════════════════════════════════════════════════════════════════

my %data = (
    yield => [5.5, 5.4, 5.8, 4.5, 4.8, 4.2],
    ctrl  => [1,   1,   1,   0,   0,   0  ],
);

my $r3 = oneway_test(\%data, formula => "yield ~ ctrl");
my $r4 = oneway_test(\%data, formula => "yield ~ ctrl", var_equal => 1);

print "\n── Mode 2: formula ─────────────────────────────────────\n";
printf "Formula : %s\n",  $r3->{formula};
printf "Method  : %s\n",  $r3->{method};
printf "Welch  F=%.4f  df=(%.4f, %.4f)  p=%.6f\n",
    $r3->{statistic}, $r3->{num_df}, $r3->{denom_df}, $r3->{p_value};
printf "Classic F=%.4f  df=(%.0f, %.0f)  p=%.6f\n",
    $r4->{statistic}, $r4->{num_df}, $r4->{denom_df}, $r4->{p_value};

# R reference output for the Welch case:
#   F = 7.5938, num df = 1, denom df = 2.6667, p-value = 0.08428

# ══════════════════════════════════════════════════════════════════════════
# Three groups with string labels — Mode 2
# ══════════════════════════════════════════════════════════════════════════

my %survey = (
    score => [72, 68, 75, 80, 85, 79, 91, 88, 94, 70],
    dept  => [qw(HR HR HR  IT IT  IT  RD RD  RD HR)],
);

my $r5 = oneway_test(\%survey, formula => "score ~ dept");
print "\n── Mode 2: string factor labels ───────────────────────\n";
printf "k=%d groups, n=%d obs\n", $r5->{k}, $r5->{n};
printf "Welch  F=%.4f  df=(%.4f, %.4f)  p=%.6f\n",
    $r5->{statistic}, $r5->{num_df}, $r5->{denom_df}, $r5->{p_value};

# ══════════════════════════════════════════════════════════════════════════
# Error cases
# ══════════════════════════════════════════════════════════════════════════

# Mismatched array lengths
eval { oneway_test({ y => [1,2,3], g => [1,2] }, formula => "y ~ g") };
print "\nCaught (length mismatch): $@" if $@;

# Malformed formula
eval { oneway_test(\%data, formula => "yield") };
print "Caught (bad formula)    : $@" if $@;

# Unknown variable name
eval { oneway_test(\%data, formula => "weight ~ ctrl") };
print "Caught (missing key)    : $@" if $@;
