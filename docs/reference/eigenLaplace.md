# Generate eigenvalues of discrete Laplace matrix.

The eigenvalues of a discrete Laplace matrix with dimension (`mm`, `nn`)
are calculated.

## Usage

``` r
eigenLaplace(mm, nn)
```

## Arguments

- mm:

  Number of rows of the discrete Laplace matrix.

- nn:

  Number of columns of the discrete Laplace matrix.

## Value

A row vector containing the eigenvalues of the discrete laplace matrix.

## Examples

``` r
eigval <- eigenLaplace(5, 5)
```
