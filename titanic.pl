#!/usr/bin/env perl

use 5.042.2;
no source::encoding;
use warnings FATAL => 'all';
use autodie ':default';
use DDP {output => 'STDOUT', array_max => 10, show_memsize => 1};
use Devel::Confess 'color';
use Stats::LikeR;

my $titanic = read_table('titanic.more.complete.csv');
$titanic = cfilter($titanic, remove => ['ticketno', 'embarked']);
$titanic = assign($titanic, survived => sub { $_->{survived} eq 'yes' ? 1 : 0});
view($titanic);
my $glm = glm(formula => 'survived ~ age + class + gender + fare', data => $titanic, family => 'binomial');
p $glm;
my @predict;
for (my $age = 10; $age <= 90; $age += 10) {
	push @predict, {
		age => $age, class => '3rd', gender => 'male', fare => 7.11
	};
}
#my $p = predict($glm, \@predict );
#p $p;
#$titanic = aoh2hoa( $titanic );
my $vc = value_counts($titanic, 'country');
p $vc;
