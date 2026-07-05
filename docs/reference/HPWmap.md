# Computation of pointwise and highest pointwise probabilities.

Pointwise (PW) probabilities and highest pointwise (HPW) probabilities
of all differences of smooths at neighboring scales are computed.

## Usage

``` r
HPWmap(smoothVec, mm, nn, prob = 0.95)
```

## Arguments

- smoothVec:

  Differences of smooths at neighboring scales.

- mm:

  Number of rows of the original input image.

- nn:

  Number of columns of the original input image.

- prob:

  Credibility level for the posterior credibility analysis

## Value

List with two arrays:

- `pw`: Pointwise probabilities (`VmapPW`) including the dimensions of
  the original input image, `mm` and `nn`.

- `hpw`: Highest pointwise probabilities (`VmapHPW`) including the
  dimensions of the original input image, `mm` and `nn`.

## Details

`HPWmap` is an internal function of
[`mrbsizeRgrid`](https://romanflury.github.io/mrbsizeR/reference/mrbsizeRgrid.md)
and is usually not used independently. The output can be analyzed with
the plotting function
[`plot.HPWmapGrid`](https://romanflury.github.io/mrbsizeR/reference/plot.HPWmapGrid.md).

## Examples

``` r
# Artificial sample data: 10 observations (5-by-2 object), 10 samples
set.seed(987)
sampleData <- matrix(stats::rnorm(100), nrow = 10)
sampleData[4:6, ] <- sampleData[4:6, ] + 5

# Calculation of the simultaneous credible intervals
HPWmap(smoothVec = sampleData, mm = 5, nn = 2, prob = 0.95)
#> $pw
#>      [,1] [,2]
#> [1,]  0.5  1.0
#> [2,]  0.5  0.5
#> [3,]  0.5  0.5
#> [4,]  1.0  0.5
#> [5,]  1.0  0.5
#> 
#> $hpw
#>      [,1] [,2]
#> [1,]  0.5  1.0
#> [2,]  0.5  0.5
#> [3,]  0.5  0.5
#> [4,]  1.0  0.5
#> [5,]  1.0  0.5
#> 
```
