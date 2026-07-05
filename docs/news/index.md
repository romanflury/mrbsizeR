# Changelog

## mrbsizeR 1.3.1

### Bug fixes

- new maintainer email.

## mrbsizeR 1.3

CRAN release: 2024-02-13

### Bug fixes

- fixed LazyData issue.
- removed old urls from DESCRIPTION and README.
- addressed NOTE from editable Rd contents.

### Bug fixes

- fixed type in rcpp code of
  [`MinLambda()`](https://romanflury.github.io/mrbsizeR/reference/MinLambda.md)
  (clang-UBSAN issue).
- fixed deprecated recycling of arrays in
  [`eigenQsphere()`](https://romanflury.github.io/mrbsizeR/reference/eigenQsphere.md).

## mrbsizeR 1.2

CRAN release: 2019-12-06

### Significant user-visible changes

- function
  [`TaperingPlot()`](https://romanflury.github.io/mrbsizeR/reference/TaperingPlot.md)
  with new argument returnseq(=FALSE), to return tapering sequences.

### Internal

- [`MinLambda()`](https://romanflury.github.io/mrbsizeR/reference/MinLambda.md):
  loops in rcpp.
- Modifications for pkgdown website.

### Minor bug fixes and improvements

- adjusting mrbsizeR vignette to adress pandoc2.8 warning.
- [`MinLambda()`](https://romanflury.github.io/mrbsizeR/reference/MinLambda.md):
  adjusted grid size, for a user defined sequence of lambdas.

## mrbsizeR 1.1.1

CRAN release: 2018-05-02

- Adressing Solaris compiler error in rcpp_dctmatrix.cpp.

## mrbsizeR 1.1.0

CRAN release: 2018-04-30

- New maintainer Roman Flury.
- [`dctMatrix()`](https://romanflury.github.io/mrbsizeR/reference/dctMatrix.md)
  and
  [`eigenLaplace()`](https://romanflury.github.io/mrbsizeR/reference/eigenLaplace.md)
  via rcpp.
