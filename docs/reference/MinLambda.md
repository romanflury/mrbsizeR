# Numerical optimization for finding appropriate smoothing levels.

Numerical optimization of an objective function \\G\\ is carried out to
find appropriate signal-dependent smoothing levels (\\\lambda\\'s). This
is easier than visual inspection via the signal-dependent tapering
function in
[`TaperingPlot`](https://romanflury.github.io/mrbsizeR/reference/TaperingPlot.md).

## Usage

``` r
MinLambda(Xmu, mm, nn, nGrid, nLambda = 2, lambda, sphere = FALSE)
```

## Arguments

- Xmu:

  Posterior mean of the input object as a vector.

- mm:

  Number of rows of the original input object.

- nn:

  Number of columns of the original input object.

- nGrid:

  Size of grid where objective function is evaluated (nGrid-by-nGrid).
  This argument is ignorded if a sequence `lambda` is specified.

- nLambda:

  Number of lambdas to minimize over. Possible arguments: 2 (default) or
  3.

- lambda:

  \\\lambda\\-sequence which is used for optimization. If nothing is
  provided,  
  `lambda <- 10^seq(-3, 10, len = nGrid)` is used for data on a grid
  and  
  `lambda <- 10^seq(-6, 1, len = nGrid)` is used for spherical data.

- sphere:

  `TRUE` or `FALSE`: Is the input object defined on a sphere?

## Value

A list with 3 objects:

`G` Value of objective function \\G\\.

`lambda` Evaluated smoothing parameters \\\lambda\\.

`minind` Index of minimal \\\lambda\\'s. `lambda`\[`minind`\] gives the
minimal values.

## Details

As signal-dependent tapering functions are quiet irregular, it is hard
to find appropriate smoothing values only by visual inspection of the
tapering function plot. A more formal approach is the numerical
optimization of an objective function.

Optimization can be carried out with 2 or 3 smoothing parameters. As the
smoothing parameters 0 and \\\infty\\ are always added, this results in
a mrbsizeR analysis with 4 or 5 smoothing parameters.

Sometimes, not all features of the input object can be extracted using
the smoothing levels proposed by `MinLambda`. It might then be necessary
to include additional smoothing levels.

[`plot.minLambda`](https://romanflury.github.io/mrbsizeR/reference/plot.minLambda.md)
creates a plot of the objective function \\G\\ on a grid. The minimum is
indicated with a white point. The minimum values of the \\\lambda\\'s
can be extracted from the output of `MinLambda`, see examples.

## Examples

``` r
# Artificial sample data
set.seed(987)
sampleData <- matrix(stats::rnorm(100), nrow = 10)
sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5

# Minimization of two lambdas on a 20-by-20-grid
minlamOut <- MinLambda(Xmu = c(sampleData), mm = 10, nn = 10, 
                       nGrid = 20, nLambda = 2)

# Minimal lambda values
minlamOut$lambda[minlamOut$minind]
#> [1]   0.1128838 297.6351442
```
