# mrbsizeR 1.3.1

## Bug fixes

* new maintainer email.


# mrbsizeR 1.3

## Bug fixes

* fixed LazyData issue.
* removed old urls from DESCRIPTION and README.
* addressed NOTE from editable Rd contents.

## Bug fixes

* fixed type in rcpp code of `MinLambda()` (clang-UBSAN issue).
* fixed deprecated recycling of arrays in `eigenQsphere()`.


# mrbsizeR 1.2

## Significant user-visible changes

* function `TaperingPlot()` with new argument returnseq(=FALSE), to return tapering sequences.

## Internal

* `MinLambda()`: loops in rcpp.
* Modifications for pkgdown website.

## Minor bug fixes and improvements

* adjusting mrbsizeR vignette to adress pandoc2.8 warning.
* `MinLambda()`: adjusted grid size, for a user defined sequence of lambdas.


# mrbsizeR 1.1.1

* Adressing Solaris compiler error in rcpp_dctmatrix.cpp.


# mrbsizeR 1.1.0

* New maintainer Roman Flury.
* `dctMatrix()` and `eigenLaplace()` via rcpp.
