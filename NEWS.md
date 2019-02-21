# mrbsizeR 1.1.2.9000

## Significant user-visible changes

* function `TaperingPlot()` with new argument returnseq(=FALSE), to return tapering sequences.

## Internal

* `MinLambda()`: loops in rcpp.
* Modifications for pkgdown website.

## Minor bug fixes and improvements

* `MinLambda()`: adjusted grid size, for a user defined sequence of lambdas.

# mrbsizeR 1.1.1

* Adressing Solaris compiler error in rcpp_dctmatrix.cpp.

# mrbsizeR 1.1.0

* New maintainer Roman Flury.
* `dctMatrix()` and `eigenLaplace()` via rcpp.
