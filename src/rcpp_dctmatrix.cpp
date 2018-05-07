#include <Rcpp.h>
using namespace Rcpp;

//' Create a n-by-n discrete cosine transform matrix.
//'
//' The discrete cosine transform (DCT) matrix for a given dimension n is
//' calculated.
//'
//' The function can be used for 1D- or 2D-DCT transforms of data.
//' \itemize{
//' \item \strong{1D:} Let \code{Q} be a m-by-n matrix with some data. \code{D} is a
//' m-by-m DCT matrix created by \code{dctMatrix(m)}. Then \code{D \%*\% Q} returns the
//' discrete cosine transform of the columns of Q. \code{t(D) \%*\% Q} returns the
//' inverse DCT of the columns of Q. As D is orthogonal, \code{solve(D) = t(D)}.
//' \item \strong{2D:} Let \code{Q} be a m-by-n matrix with some data. \code{D_m} is a
//' m-by-m DCT matrix created by \code{dctMatrix(m)}, \code{D_n} a n-by-n DCT matrix
//' created by \code{dctMatrix(n)}. \code{D_m \%*\% Q \%*\% t(D_n)} computes the 2D-DCT
//' of Q. The inverse 2D-DCT of Q can be computed via \cr \code{t(D_mm) \%*\% DCT_Q \%*\% D_n}.
//' D_m transforms along columns, D_n along rows. Since D is orthogonal, \code{solve(D) = t(D)}.
//' }
//' It can be faster to use \code{dctMatrix} than using a direct transformation,
//' especially when calculating several DCT's.
//'
//' @param n Dimension for the DCT matrix.
//' @return The n-by-n DCT matrix.
//' @examples
//' D <- dctMatrix(5)
// [[Rcpp::export]]
NumericMatrix dctMatrix(int n) {
  NumericMatrix dctMatrix(n, n);

  for(int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == 0) {
        dctMatrix(i, j) = 1/ std::sqrt(static_cast<double>(n));
      }
      if (i > 0) {
        dctMatrix(i, j) = std::sqrt(2/static_cast<double>(n)) * cos((2*j + 1) * i * M_PI / (2 * n));
      }
    }
  }

  return(dctMatrix);
}
