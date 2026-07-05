# Generate eigenvalues of precision matrix Q on the surface of a sphere.

The eigenvalues of the precision matrix Q with dimension (`mm`, `nn`)
and polar angle limits `phimin`, `phimax` are calculated.

## Usage

``` r
eigenQsphere(phimin, phimax, mm, nn)
```

## Arguments

- phimin:

  Polar angle minimum.

- phimax:

  Polar angle maximum.

- mm:

  Number of rows of precision matrix Q.

- nn:

  Number of columns of precision matrix Q.

## Value

A list containing 2 elements:

- eigval Row vector containing the eigenvalues of Q.

- eigvec Matrix containing the eigenvectors of Q as columns.

## Details

The corresponding function for data on a grid is
` `[`eigenLaplace`](https://romanflury.github.io/mrbsizeR/reference/eigenLaplace.md).

## Examples

``` r
eig_out <- eigenQsphere(phimin = 180/10, phimax = 180 - 180/10, mm = 10, nn = 20)
```
