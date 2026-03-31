# Synopsis

Get basic R statistical functions working in Perl as if they were part of List::Util, like `min`, `max`, `sum`, etc.
I've used Artificial Intelligence tools such as Claude, Gemini, and Grok to write this.
There are other similar tools on CPAN, but I want speed and a form like List::Util, which I've gotten here with the help of AI, which often required many attempts to do correctly.
This is meant to call subroutines directly through eXternal Subroutines (XS) for speed.

# Functions/Subroutines

## cor

    cor($array1, $array2, $method = 'pearson'),

that is, `pearson` is the default and will be used if `$method` is not specified.

Just like R, `pearson`, `spearman`, and `kendall` are available

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

## median

## p_adjust

## pearson_r

## rnorm

## scale

## sd

## t_test

## var
