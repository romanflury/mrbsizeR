# Create a n-by-n discrete cosine transform matrix.

The discrete cosine transform (DCT) matrix for a given dimension n is
calculated.

## Usage

``` r
dctMatrix(n)
```

## Arguments

- n:

  Dimension for the DCT matrix.

## Value

The n-by-n DCT matrix.

## Details

The function can be used for 1D- or 2D-DCT transforms of data.

- **1D:** Let `Q` be a m-by-n matrix with some data. `D` is a m-by-m DCT
  matrix created by `dctMatrix(m)`. Then `D %*% Q` returns the discrete
  cosine transform of the columns of Q. `t(D) %*% Q` returns the inverse
  DCT of the columns of Q. As D is orthogonal, `solve(D) = t(D)`.

- **2D:** Let `Q` be a m-by-n matrix with some data. `D_m` is a m-by-m
  DCT matrix created by `dctMatrix(m)`, `D_n` a n-by-n DCT matrix
  created by `dctMatrix(n)`. `D_m %*% Q %*% t(D_n)` computes the 2D-DCT
  of Q. The inverse 2D-DCT of Q can be computed via  
  `t(D_mm) %*% DCT_Q %*% D_n`. D_m transforms along columns, D_n along
  rows. Since D is orthogonal, `solve(D) = t(D)`.

It can be faster to use `dctMatrix` than using a direct transformation,
especially when calculating several DCT's.

## Examples

``` r
D <- dctMatrix(5)
```
