# Turn matrix 90 degrees counter-clockwise.

Help function to turn matrix 90 degrees counter-clockwise. `turnmat` is
used as an internal function in some of the plotting functions.
Depending on how the data is stored, it can be necessary to turn the
matrices 90 degrees counter-clockwise for being able to plot them
correctly.

## Usage

``` r
turnmat(x)
```

## Arguments

- x:

  Matrix to be turned.

## Value

Matrix `x`, turned by 90 degrees counter-clockwise.

## Examples

``` r
set.seed(987)
sampleMat <- matrix(stats::rnorm(100), nrow = 10)
sampleMatTurn <- turnmat(x = sampleMat)
```
