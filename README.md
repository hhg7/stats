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

## fisher_test

    my $array_data = [
    	[10, 2],
    	[3, 15]
    ];
    my $res1 = fisher_test($array_data);

I have the p-value calculated very precisely, but there are some inexactness (approximately 1% for the confidence intervals) which I couldn't rectify.  The answers are very close to R besides the p-value, where they are identical.

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

Returns array

    my @q = p_adjust(\@pvalues, $method);

## rbinom

    my $binom = rbinom( n => $n, prob => 0.5, size => 9);

## rnorm

    my ($rmean, $sd, $n) = (10, 2, 9999);
    my $normals = rnorm( n => $n, mean => $rmean, sd => $sd);

## runif

    my $unif = runif( n => $n, min => 0, max => 1);

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

## Shapiro Test

tests to see if an array reference is normally distributed, returns a p-value and a statistic

    my $shapiro = shapiro_test(
    	[1..5]
    );

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
