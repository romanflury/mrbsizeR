#' Create a n-by-n discrete Fourier transform matrix.
#' 
#' The discrete Fourier transform (DFT) matrix for a given dimension n is
#' calculated.
#' 
#' The DFT matrix can be used for computing the discrete Fourier transform of 
#' a matrix or vector. \code{dftMatrix(n) \%*\% testMatrix} is the same as
#' \code{apply(testMatrix, MARGIN = 2, FUN = fft)}.
#' 
#' @param n Dimension for the DFT matrix.
#' @return The n-by-n DFT matrix.
#' @examples
#' set.seed(987)
#' testMatrix <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5)
#' D <- dftMatrix(5)
#' 
#' # Discrete Fourier transform with matrix multiplication:
#' D %*% testMatrix 
#' 
#' # Discrete Fourier transform with function fft: 
#' apply(testMatrix, MARGIN = 2, FUN = fft)
#' 
dftMatrix <- function(n) {
  omega <- exp(-2 * pi * (0 + 1i) / n)
  sq <- rep(0:(n - 1), n) * rep(0:(n - 1), each = n)
  AA <- sapply(sq, function(nn, omega) omega^nn, omega)
  dft.mat <- matrix(AA, n, n)
  return(dft.mat)
}

