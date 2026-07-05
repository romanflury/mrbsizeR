# Swap the quadrants or halves of a 2d matrix.

`fftshift` is an R equivalent to the Matlab function `fftshift` applied
on matrices. For more information about `fftshift` see the Matlab
documentation.

## Usage

``` r
fftshift(inputMatrix, dimension = -1)
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

It is possible to swap the halves or the quadrants of the input matrix.
Halves can be swapped along the rows (`dimension = 1`) or along the
columns (`dimension = 2`). When swapping the quadrants, `fftshift` swaps
the first quadrant with the third and the second quadrant with the
fourth (`dimension = -1`).

## Examples

``` r
set.seed(987) 
sampleMat <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5)

# Swap halves along the rows:
fftshift(sampleMat, dimension = 1)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    2    2    5   10   10
#> [2,]    8    8   10    8    8
#> [3,]    9    1    4    3    9
#> [4,]    9    1    2    9    9
#> [5,]    7    3    2    3    3

# Swap halves along the columns:
fftshift(sampleMat, dimension = 2)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    3    9    9    1    4
#> [2,]    9    9    9    1    2
#> [3,]    3    3    7    3    2
#> [4,]   10   10    2    2    5
#> [5,]    8    8    8    8   10

# Swap first quadrant with third and second quadrant with fourth:
fftshift(sampleMat, dimension = -1)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   10   10    2    2    5
#> [2,]    8    8    8    8   10
#> [3,]    3    9    9    1    4
#> [4,]    9    9    9    1    2
#> [5,]    3    3    7    3    2
```
