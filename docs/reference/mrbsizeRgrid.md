# Multiresolution analysis of random signals.

`mrbsizeRgrid` is the interface of the scale space multiresolution
method for data on a regular grid. Here, the differences of smooths as
well as the posterior credibility analysis are computed. The output can
be analyzed with the plotting functions
[`plot.smMeanGrid`](https://romanflury.github.io/mrbsizeR/reference/plot.smMeanGrid.md),
[`plot.CImapGrid`](https://romanflury.github.io/mrbsizeR/reference/plot.CImapGrid.md)
and
[`plot.HPWmapGrid`](https://romanflury.github.io/mrbsizeR/reference/plot.HPWmapGrid.md).

## Usage

``` r
mrbsizeRgrid(
  posteriorFile,
  mm,
  nn,
  lambdaSmoother,
  prob = 0.95,
  smoothOut = FALSE
)
```

## Arguments

- posteriorFile:

  Matrix with posterior samples as column vectors.

- mm:

  Number of rows of the original object.

- nn:

  Number of columns of the original object.

- lambdaSmoother:

  Vector consisting of the smoothing levels to be used.

- prob:

  Credibility level for the posterior credibility analysis.

- smoothOut:

  Should the differences of smooths at neighboring scales be returned as
  output (FALSE by default)?

## Value

A list containing the following sublists:

`smMean` Posterior mean of all differences of smooths created.

`hpout` Pointwise (PW) and highest pointwise (HPW) probabilities of all
differences of smooths created.

`ciout` Simultaneous credible intervals (CI) of all differences of
smooths created.

`smoothSamples` Samples of differences of smooths at neighboring scales,
as column vectors.

## Details

`mrbsizeRgrid` conducts two steps of the scale space multiresolution
analysis:

1.  Extraction of scale-dependent features from the reconstructed
    signal. This is done by smoothing at different smoothing levels and
    taking the difference of smooths at neighboring scales.

2.  Posterior credibility analysis of the differences of smooths
    created. Three different methods are applied: Pointwise
    probabilities (see
    [`HPWmap`](https://romanflury.github.io/mrbsizeR/reference/HPWmap.md)),
    highest pointwise probabilities (see
    [`HPWmap`](https://romanflury.github.io/mrbsizeR/reference/HPWmap.md))
    and simultaneous credible intervals (see
    [`CImap`](https://romanflury.github.io/mrbsizeR/reference/CImap.md)).

The signal can be reconstructed using the build-in multivariate
t-distribution sampling
[`rmvtDCT`](https://romanflury.github.io/mrbsizeR/reference/rmvtDCT.md).
It is also possible to provide samples generated with other methods, see
the parameter `posteriorFile` and the examples.

For further information and examples, see the vignette.

## Examples

``` r
# Artificial sample data
set.seed(987)
sampleData <- matrix(stats::rnorm(100), nrow = 10)
sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5

# Generate samples from multivariate t-distribution
tSamp <- rmvtDCT(object = sampleData, lambda = 0.2, sigma = 6, nu0 = 15,
                  ns = 1000)  
 
# mrbsizeRgrid analysis
mrbOut <- mrbsizeRgrid(posteriorFile = tSamp$sample, mm = 10, nn = 10, 
                       lambdaSmoother = c(1, 1000), prob = 0.95)
```
