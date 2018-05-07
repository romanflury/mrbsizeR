#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector for_eigenLaplace(NumericVector mu, NumericVector lambda, int mm, int nn) {
  NumericVector Mu(nn * mm);

  for (int k = 0; k < nn; k++) {
    for (int s = 0; s < mm; s++) {
      Mu(k * mm + s) = lambda(s) + mu(k);
    }
  }

  return(Mu);
}

