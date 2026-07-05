# Sampling from marginal posterior multivariate t-distribution.

Samples from a marginal posterior multivariate t-distribution with
normal-inverse-chi-squared-prior are generated.

## Usage

``` r
rmvtDCT(object, lambda, sigma, nu0, ns)
```

## Arguments

- object:

  Observed object, as `matrix`.

- lambda:

  Scaling parameter (\\\lambda\\) of the
  normal-inverse-chi-squared-prior.

- sigma:

  Square root of the \\\sigma\_{0}^{2}\\ parameter of the
  normal-inverse-chi-squared-prior.

- nu0:

  Degrees of freedom (\\\nu\_{0}\\) of the
  normal-inverse-chi-square-prior.

- ns:

  Number of samples that should be generated.

## Value

A list containing the following elements:

`sample` Samples of the marginal posterior of the input.

`mu` Mean of the marginal posterior of the input.

## Details

An eigenvalue decomposition is used for sampling. To speed up
computations, a 2D discrete cosine transform (DCT) has been implemented,
see
[`dctMatrix`](https://romanflury.github.io/mrbsizeR/reference/dctMatrix.md).
The output is a list containing

1.  Samples of the marginal posterior of the input as column vectors.

2.  The mean of the marginal posterior of the input as a vector.

## Examples

``` r
# Artificial sample data
set.seed(987)
sampleData <- matrix(stats::rnorm(100), nrow = 10)
sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5

# Sampling from a multivariate t-distribution
t_dist_samp <- rmvtDCT(object = sampleData, lambda = 1, sigma = 10,
                       nu0 = 50, ns = 1000)
```
