#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

# -------------------------------------------------------------------
# Command-line argument handling
# -------------------------------------------------------------------
my $target = shift;
if (!$target) {
    die "Usage: perl xs.sxn.pl <function_name> [source_file.xs]\n";
}

my $filename = shift;
if (!$filename) {
    my @xs_files = glob("*.xs");
    if (@xs_files) {
        $filename = $xs_files[0];
    } else {
        die "Error: No .xs file found in the current directory. Please specify one manually.\n";
    }
}

if (!-e $filename) {
    die "Error: File '$filename' does not exist.\n";
}

open my $fh, '<', $filename or die "Can't open $filename: $!\n";
my @lines = <$fh>;
close $fh;

# -------------------------------------------------------------------
# Step 1: Chunk the XS file into contextual blocks
# -------------------------------------------------------------------
my @blocks;
my $id_counter = 0;

my $current_block = { lines => [], name => 'HEADER', is_function => 0, id => $id_counter++ };
push @blocks, $current_block;

for (my $i = 0; $i < @lines; $i++) {
    my $line = $lines[$i];
    chomp $line;

    if ($line =~ /^(?:MODULE|PACKAGE|BOOT|INCLUDE|TYPEMAP)\b/) {
        $current_block = { lines => [], name => 'BOILERPLATE', is_function => 0, id => $id_counter++ };
        push @blocks, $current_block;
    }
    elsif ($line =~ /^[A-Za-z_]/ && $line !~ /^#/) {
        my $next_line = $lines[$i + 1] // '';
        my $combined  = $line . " " . $next_line;
        
        if ($combined =~ /\b([A-Za-z_]\w*)\s*\(/) {
            my $possible_name = $1;
            if ($possible_name !~ /^(?:if|while|for|switch|return)$/) {
                $current_block = { lines => [], name => $possible_name, is_function => 1, id => $id_counter++ };
                push @blocks, $current_block;
            }
        }
    }

    push @{$current_block->{lines}}, $line;
}

# -------------------------------------------------------------------
# Step 2: Build the dependency tree mapping (Stricter matching)
# -------------------------------------------------------------------
my %function_map = map { $_->{name} => $_ } grep { $_->{is_function} } @blocks;

if (!exists $function_map{$target}) {
    die "Error: Target function '$target' was not found in '$filename'.\n";
}

for my $name (keys %function_map) {
    my $block = $function_map{$name};
    my $body_text = join("\n", @{$block->{lines}});
    
    my %deps;
    # FIX: Require an opening parenthesis down the line to confirm it's a call
    while ($body_text =~ /\b([A-Za-z_]\w*)\s*\(/g) {
        my $word = $1;
        if (exists $function_map{$word} && $word ne $name) {
            $deps{$word} = 1;
        }
    }
    $block->{deps} = [ keys %deps ];
}

# -------------------------------------------------------------------
# Step 3: Trace recursive dependencies (BFS)
# -------------------------------------------------------------------
my %needed;
my @queue = ($target);

while (@queue) {
    my $curr = shift @queue;
    next if $needed{$curr};
    
    if (exists $function_map{$curr}) {
        $needed{$curr} = 1;
        push @queue, @{$function_map{$curr}->{deps}};
    }
}

# -------------------------------------------------------------------
# Step 4: Write output while maintaining top-to-bottom compilation order
# -------------------------------------------------------------------
my $output_file = "$target.section.xs";
open my $out_fh, '>', $output_file or die "Can't write to $output_file: $!\n";

for my $block (sort { $a->{id} <=> $b->{id} } @blocks) {
    if (!$block->{is_function} || $needed{$block->{name}}) {
        next if @{$block->{lines}} == 0;
        print $out_fh join("\n", @{$block->{lines}}), "\n";
    }
}

close $out_fh;
print "Extraction complete: '$output_file' created.\n";