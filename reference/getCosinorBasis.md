# Basis matrix for cosinor

Generate basis matrix for cosinor regression.

## Usage

``` r
getCosinorBasis(x, period, intercept)
```

## Arguments

- x:

  Values of the predictor variable.

- period:

  Period for the predictor variable.

- intercept:

  If `TRUE`, a column of ones will be included in the basis.

## Value

A matrix with a row for each value of `x` and a column for each
component of the decomposition.

## Examples

``` r
b = getCosinorBasis(seq(0, 20, 4), period = 24, intercept = FALSE)
```
