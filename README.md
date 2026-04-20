# Synopsis

Get basic R statistical functions working in Perl as if they were part of List::Util, like `min`, `max`, `sum`, etc.
I've used Artificial Intelligence tools such as Claude, Gemini, and Grok to write this.
There are other similar tools on CPAN, but I want speed and a form like List::Util, which I've gotten here with the help of AI, which often required many attempts to do correctly.
This is meant to call subroutines directly through eXternal Subroutines (XS) for speed.

There **are** other modules on CPAN that can do **PARTS** of this, but this works the way that I **want** it to.

# Functions/Subroutines

## aov

    aov(
    {
        yield => [5.5, 5.4, 5.8, 4.5, 4.8, 4.2],
        ctrl  => [1,     1,   1,   0,   0,   0]
    },
    'yield ~ ctrl');

which returns

    {
        ctrl        {
            Df          1,
            "F value"   25.6000000000001,
            "Mean Sq"   1.70666666666667,
            Pr(>F)      0.00718232855871859,
            "Sum Sq"    1.70666666666667
        },
        Residuals   {
            Df          4,
            "Mean Sq"   0.0666666666666665,
            "Sum Sq"    0.266666666666666
       }
    }

#chisq_test

    my @test_data = ([762, 327, 468], [484, 239, 477]);
    my $test_data = chisq_test(\@test_data);

which outputs:

    {
    data.name   "Perl ArrayRef",
    expected    [
        [0] [
                [0] 703.671381936888,
                [1] 319.645266594124,
                [2] 533.683351468988
            ],
        [1] [
                [0] 542.328618063112,
                [1] 246.354733405876,
                [2] 411.316648531012
            ]
    ],
    method      "Pearson's Chi-squared test",
    observed    [
        [0] [
                [0] 762,
                [1] 327,
                [2] 468
            ],
        [1] [
                [0] 484,
                [1] 239,
                [2] 477
            ]
    ],
    p.value     2.95358918321176e-07,
    parameter   {
        df   2
    },
    statistic   {
        X-squared   30.0701490957547
    }
    }

## cor

    cor($array1, $array2, $method = 'pearson'),

that is, `pearson` is the default and will be used if `$method` is not specified.

Just like R, `pearson`, `spearman`, and `kendall` are available

## cor_test

    my $result = cor_test(
    		'x'         => $x,
    		'y'         => $y,
    		alternative => 'two.sided'
    		method      => 'pearson',
    		continuity  => 1
    	);

## cov

    cov($array1, $array2, 'pearson')

or

    cov($array1, $array2, 'spearman')

or

    cov($array1, $array2, 'kendall')

## fisher_test

    my $array_data = [
    	[10, 2],
    	[3, 15]
    ];
    my $res1 = fisher_test($array_data);

I have the p-value calculated very precisely, but there are some inexactness (approximately 1% for the confidence intervals) which I couldn't rectify.  The answers are very close to R besides the p-value, where they are identical.

## glm

takes a hash of an array as input

    my %tooth_growth = (
    	dose => [qw(0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
    0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
    2.0 2.0 2.0)],
    	len  => [qw(4.2 11.5  7.3  5.8  6.4 10.0 11.2 11.2  5.2  7.0 16.5 16.5 15.2 17.3 22.5
    17.3 13.6 14.5 18.8 15.5 23.6 18.5 33.9 25.5 26.4 32.5 26.7 21.5 23.3 29.5
    15.2 21.5 17.6  9.7 14.5 10.0  8.2  9.4 16.5  9.7 19.7 23.3 23.6 26.4 20.0
    25.2 25.8 21.2 14.5 27.3 25.5 26.4 22.4 24.5 24.8 30.9 26.4 27.3 29.4 23.0)],
    	supp => [qw(VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC VC
    VC VC VC VC VC OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ
    OJ OJ OJ OJ OJ OJ OJ OJ OJ OJ)]
    );

    my $glm_teeth = glm(
    	data    => \%tooth_growth,
    	formula => 'len ~ dose + supp',
    	family  => 'gaussian'
    );

I'm not completely confident that this is working perfectly, I've only gotten this subroutine to work for simple cases

## lm

This is the linear models function.

    $lm = lm(formula =>  'mpg ~ wt + hp', data => $mtcars);

where `$mtcars` is a hash of hashes

## matrix

    my $mat1 = matrix(
    	data => [1..6],
    	nrow => 2
    );

## mean

    mean(1,2,3);
    
or

    my @arr = 1..8;
    mean(@arr, 4, 5)

or

    mean([1,1], [2,2]) # 1.5

## median

works like mean, taking array references and arrays:

    median( $test_data[$i][0] )

## p_adjust

Returns array of false-discovery-rate-corrected p-values, where methods available are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"

    my @q = p_adjust(\@pvalues, $method);

## rbinom

Create a binomial distribution of numbers

    my $binom = rbinom( n => $n, prob => 0.5, size => 9);

## read_table

I've tried to make this as simple as possible, trying to follow from R:

    my $test_data = read_table('t/HepatitisCdata.csv');

output types can be AOH, HOA, HOH

    read_table($filename, 'output.type' => 'aoh');

    read_table($filename, 'output.type' => 'hoa');

## rnorm

Make a normal distribution of numbers, with pre-set mean `mean`, standard deviation `sd`, and number `n`.

    my ($rmean, $sd, $n) = (10, 2, 9999);
    my $normals = rnorm( n => $n, mean => $rmean, sd => $sd);

## runif

Make a distribution of approximately uniform distribution

    my $unif = runif( n => $n, min => 0, max => 1);

where `n` is the number of items, the values are between `min` and `max`

## scale

    my @scaled_results = scale(1..5);

## sd

    my $stdev = sd(2,4,4,4,5,5,7,9);

Correct answer is 2.1380899352994;

## seq

Works as closely as I can to R's seq, which is very similar to Perl's `for` loops.  Returns an array, not an array reference.

### Example 1: Standard integer sequence

    say 'seq(1, 5):';
    my @seq = seq(1, 5);
    say join(', ', @seq), "\n";

    say 'seq(1, 2, 0.25):';
    @seq = seq(1, 2, 0.25);

### Example 2: Fractional steps

    say 'seq(1, 2, 0.25):';
    @seq = seq(1, 2, 0.25);
    say join(", ", @seq), "\n";
    for (my $idx = 2; $idx >= 1; $idx -= 0.25) { # count down to pop
    	is_approx(pop @seq, $idx, "seq item $idx with fractional step");
    }

### Example 3: Negative steps

    say 'seq(10, 5, -1):';
    @seq = seq(10, 5, -1);
    say join(", ", @seq), "\n";
    for (my $idx = 5; $idx <= 10; $idx++) { # count down to pop
        is_approx(pop @seq, $idx, "seq item $idx with negative step");
    }
}

## shapiro_test

tests to see if an array reference is normally distributed, returns a p-value and a statistic

    my $shapiro = shapiro_test(
    	[1..5]
    );

and returns the hash reference:
    {
    p.value     0.589650577093106,
    p_value     0.589650577093106,
    statistic   0.960870680168535,
    W           0.960870680168535
    }

## t_test

There are 1-sample and 2-sample t-tests:

    my $t_test = t_test( 'x' => $test_data[$i][$j], mu => mean( $test_data[$i][$j] ));

or 2-sample:

    $t_test = t_test(
    	'x'    => $test_data[3][0],
    	'y'    => $test_data[3][1],
	    paired => true
    );


returns a hash reference, which looks like:

    conf_int     => [
        -0.06672889, 0.25672889
    ],
    df        => 5,
    estimate  => 0.095,
    p_value   => 0.19143688433660,
    statistic => 1.50996688705414

## var

as simple as possible:

    var(2, 4, 5, 8, 9)

## write_table

mimics R's "write.table", with data as first argument to subroutine, and output file as second

    write_table(\@data_aoh, $tmp_file, sep => "\t", 'row.names' => true);


