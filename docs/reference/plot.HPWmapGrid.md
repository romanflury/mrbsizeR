# Plotting of pointwise and highest pointwise probabilities.

Maps with pointwise (PW) probabilities and/or highest pointwise (HPW)
probabilities of all differences of smooths at neighboring scales are
plotted.

## Usage

``` r
# S3 method for class 'HPWmapGrid'
plot(
  x,
  plotWhich = "Both",
  color = c("firebrick1", "gainsboro", "dodgerblue3"),
  turnOut = TRUE,
  title,
  aspRatio = 1,
  ...
)
```

## Arguments

- x:

  List containing the pointwise (PW) and highest pointwise (HPW)
  probabilities of all differences of smooths.

- plotWhich:

  Which probabilities shall be plotted? `HPW`, `PW` or `Both`?

- color:

  Vector of length 3 containing the colors to be used in the credibility
  maps. The first color represents the credibly negative pixels, the
  second color the pixels that are not credibly different from zero and
  the third color the credibly positive pixels.

- turnOut:

  Logical. Should the output images be turned 90 degrees
  counter-clockwise?

- title:

  Vector containing one string per plot. The required number of titles
  is equal to `length(mrbOut$hpout)`. If no `title` is passed, defaults
  are used.

- aspRatio:

  Adjust the aspect ratio of the plots. The default `aspRatio = 1`
  produces square plots.

- ...:

  Further graphical parameters can be passed.

## Value

Plots of pointwise and/or highest pointwise probabilities for all
differences of smooths are created.

## Details

The default colors of the maps have the following meaning:

- **Blue**: Credibly positive pixels.

- **Red**: Credibly negative pixels.

- **Grey**: Pixels that are not credibly different from zero.

`x` corresponds to the `hpout`-part of the output of
[`mrbsizeRgrid`](https://romanflury.github.io/mrbsizeR/reference/mrbsizeRgrid.md).

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

# Posterior mean of the differences of smooths
plot(x = mrbOut$smMean, turnOut = TRUE)


# Credibility analysis using pointwise (PW) maps
plot(x = mrbOut$hpout, plotWhich = "PW", turnOut = TRUE)


# Credibility analysis using highest pointwise probability (HPW) maps
plot(x = mrbOut$hpout, plotWhich = "HPW", turnOut = TRUE)

```
