# Multiresolution analysis of random signals for spherical data.

`mrbsizeRSphere` is the interface of the scale space multiresolution
method for spherical data. Here, the differences of smooths as well as
the posterior credibility analysis are computed. The output can be
analyzed with the plotting functions
[`plot.smMeanSphere`](https://romanflury.github.io/mrbsizeR/reference/plot.smMeanSphere.md),
[`plot.CImapSphere`](https://romanflury.github.io/mrbsizeR/reference/plot.CImapSphere.md)
and
[`plot.HPWmapSphere`](https://romanflury.github.io/mrbsizeR/reference/plot.HPWmapSphere.md).

## Usage

``` r
mrbsizeRsphere(
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

In contrast to `mrbsizeRgrid`, `mrbsizeRsphere` does not conduct
Bayesian signal reconstruction via sampling from a posterior
distribution. Samples of the posterior distribution have to be provided
instead.

For further information and examples, see
[`mrbsizeRgrid`](https://romanflury.github.io/mrbsizeR/reference/mrbsizeRgrid.md)
and the vignette.

## Examples

``` r
# Artificial spherical sample data
set.seed(987)
sampleData <- matrix(stats::rnorm(2000), nrow = 200)
sampleData[50:65, ] <- sampleData[50:65, ] + 5

# mrbsizeRsphere analysis
mrbOut <- mrbsizeRsphere(posteriorFile = sampleData, mm = 10, nn = 20,  
                          lambdaSmoother = c(1, 1000), prob = 0.95)
```
