# Computation of simultaneous credible intervals.

Simultaneous credible intervals for all differences of smooths at
neighboring scales \\z\_{i}\\ are computed.

## Usage

``` r
CImap(smoothVec, mm, nn, prob = 0.95)
```

## Arguments

- smoothVec:

  Differences of smooths at neighboring scales.

- mm:

  Number of rows of the original input object.

- nn:

  Number of columns of the original input object.

- prob:

  Credibility level for the posterior credibility analysis. By default
  `prob = 0.95`.

## Value

An array with simultaneous credible intervals `VmapCI` and the
dimensions of the original input object, `mm` and `nn`.

## Details

`CImap` is an internal function of
[`mrbsizeRgrid`](https://romanflury.github.io/mrbsizeR/reference/mrbsizeRgrid.md)
and is usually not used independently. The output can be analyzed with
the plotting function
[`plot.CImapGrid`](https://romanflury.github.io/mrbsizeR/reference/plot.CImapGrid.md).

## Examples

``` r
# Artificial sample data: 10 observations (5-by-2 object), 10 samples
set.seed(987)
sampleData  <- matrix(stats::rnorm(100), nrow = 10)
sampleData [4:6, ] <- sampleData [4:6, ] + 5

# Calculation of the simultaneous credible intervals
CImap(smoothVec = sampleData , mm = 5, nn = 2, prob = 0.95)
#>      [,1] [,2]
#> [1,]  0.5  1.0
#> [2,]  0.5  0.5
#> [3,]  0.5  0.5
#> [4,]  1.0  0.5
#> [5,]  1.0  0.5
```
