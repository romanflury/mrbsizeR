# Inverse FFT shift of a 2d matrix.

`ifftshift` is an R equivalent to the Matlab function `ifftshift`
applied on matrices. For more information about `ifftshift` see the
Matlab documentation.

## Usage

``` r
ifftshift(inputMatrix, dimension = -1)
```

## Arguments

- inputMatrix:

  Matrix to be swapped.

- dimension:

  Which swap should be performed?

  - `1`: swap halves along the rows.

  - `2`: swap halves along the columns.

  - `-1`: swap first quadrant with third and second quadrant with
    fourth.

## Value

Swapped matrix.

## Details

`ifftshift` is the inverse function to
[`fftshift`](https://romanflury.github.io/mrbsizeR/reference/fftshift.md).
For more information see the details of
[`fftshift`](https://romanflury.github.io/mrbsizeR/reference/fftshift.md)

## Examples

``` r
set.seed(987)
sampleMat <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5)

# Swap halves along the rows:
ifftshift(sampleMat, dimension = 1)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    7    3    2    3    3
#> [2,]    2    2    5   10   10
#> [3,]    8    8   10    8    8
#> [4,]    9    1    4    3    9
#> [5,]    9    1    2    9    9

# Swap halves along the columns:
ifftshift(sampleMat, dimension = 2)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    4    3    9    9    1
#> [2,]    2    9    9    9    1
#> [3,]    2    3    3    7    3
#> [4,]    5   10   10    2    2
#> [5,]   10    8    8    8    8

# Swap first quadrant with third and second quadrant with fourth:
ifftshift(sampleMat, dimension = -1)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    2    3    3    7    3
#> [2,]    5   10   10    2    2
#> [3,]   10    8    8    8    8
#> [4,]    4    3    9    9    1
#> [5,]    2    9    9    9    1
```
