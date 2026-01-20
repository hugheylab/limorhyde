# Basis matrix for periodic splines

Generate basis matrix for a periodic B-spline using
[`pbs::pbs()`](https://rdrr.io/pkg/pbs/man/pbs.html).

## Usage

``` r
getSplineBasis(x, period, nKnots, intercept)
```

## Arguments

- x:

  Values of the predictor variable.

- period:

  Period for the predictor variable.

- nKnots:

  Number of internal knots.

- intercept:

  If `TRUE`, a column of ones will be included in the basis.

## Value

A matrix with a row for each value of `x` and a column for each
component of the decomposition.

## Examples

``` r
b = getSplineBasis(seq(0, 20, 4), period = 24, nKnots = 3, intercept = FALSE)
```
