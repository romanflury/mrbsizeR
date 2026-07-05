# Generate a tridiagonal matrix.

Generate a tridiagonal matrix with `upperDiag` as superdiagonal,
`lowerDiag` as subdiagonal and `mainDiag` as diagonal.

## Usage

``` r
tridiag(mainDiag, upperDiag, lowerDiag)
```

## Arguments

- mainDiag:

  Diagonal of tridiagonal matrix.

- upperDiag:

  Superdiagonal of tridiagonal matrix. Must have length
  `length(mainDiag) - 1`.

- lowerDiag:

  Subdiagonal of tridiagonal matrix. Must have length
  `length(mainDiag) - 1`.

## Value

Tridiagonal matrix.

## Examples

``` r
set.seed(987)
mainDiag <- sample(100:110, size = 6, replace = TRUE)
upperDiag <- sample(10:20, size = 5, replace = TRUE)
lowerDiag <- sample(1:10, size = 5, replace = TRUE)

tridiag(mainDiag, upperDiag, lowerDiag)  
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  108   10    0    0    0    0
#> [2,]    2  108   12    0    0    0
#> [3,]    0    2  106   11    0    0
#> [4,]    0    0    5  101   17    0
#> [5,]    0    0    0   10  107   13
#> [6,]    0    0    0    0    3  100
 
```
