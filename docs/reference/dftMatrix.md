# Create a n-by-n discrete Fourier transform matrix.

The discrete Fourier transform (DFT) matrix for a given dimension n is
calculated.

## Usage

``` r
dftMatrix(n)
```

## Arguments

- n:

  Dimension for the DFT matrix.

## Value

The n-by-n DFT matrix.

## Details

The DFT matrix can be used for computing the discrete Fourier transform
of a matrix or vector. `dftMatrix(n) %*% testMatrix` is the same as
`apply(testMatrix, MARGIN = 2, FUN = fft)`.

## Examples

``` r
set.seed(987)
testMatrix <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5)
D <- dftMatrix(5)

# Discrete Fourier transform with matrix multiplication:
D %*% testMatrix 
#>                     [,1]                [,2]                [,3]
#> [1,] 35.000000+0.000000i 15.000000+0.000000i 23.000000+0.000000i
#> [2,]  6.972136-3.889983i -0.263932+6.069610i  2.045085+9.371808i
#> [3,] -1.972136+4.167497i -4.736068+5.065553i -3.545085+1.849112i
#> [4,] -1.972136-4.167497i -4.736068-5.065553i -3.545085-1.849112i
#> [5,]  6.972136+3.889983i -0.263932-6.069610i  2.045085-9.371808i
#>                     [,4]                [,5]
#> [1,] 33.000000+0.000000i 39.000000+0.000000i
#> [2,] -2.263932+3.163440i  3.736068+3.163440i
#> [3,] -6.736068-7.245181i -0.736068-7.245181i
#> [4,] -6.736068+7.245181i -0.736068+7.245181i
#> [5,] -2.263932-3.163440i  3.736068-3.163440i

# Discrete Fourier transform with function fft: 
apply(testMatrix, MARGIN = 2, FUN = fft)
#>                     [,1]                [,2]                [,3]
#> [1,] 35.000000+0.000000i 15.000000+0.000000i 23.000000+0.000000i
#> [2,]  6.972136-3.889983i -0.263932+6.069610i  2.045085+9.371808i
#> [3,] -1.972136+4.167497i -4.736068+5.065553i -3.545085+1.849112i
#> [4,] -1.972136-4.167497i -4.736068-5.065553i -3.545085-1.849112i
#> [5,]  6.972136+3.889983i -0.263932-6.069610i  2.045085-9.371808i
#>                     [,4]                [,5]
#> [1,] 33.000000+0.000000i 39.000000+0.000000i
#> [2,] -2.263932+3.163440i  3.736068+3.163440i
#> [3,] -6.736068-7.245181i -0.736068-7.245181i
#> [4,] -6.736068+7.245181i -0.736068+7.245181i
#> [5,] -2.263932-3.163440i  3.736068-3.163440i
```
