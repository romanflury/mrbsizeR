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

// [[Rcpp::export]]
NumericMatrix initLambdaMat(int nGrid, int N, NumericVector lambda, NumericVector fac, NumericVector D) {
  NumericMatrix LambdaMat(N, nGrid);
  int l1 = 0;
  NumericVector lambda2(N);

  for (int i = 0; i < nGrid; i++) {
    l1 = lambda[i];
    lambda2 = fac / (1.0 + l1 * D);
    lambda2[0] = 0;
    LambdaMat(_, i) = lambda2;
  }
  return(LambdaMat);
}

// [[Rcpp::export]]
List min2Lambda(int nGrid, NumericMatrix LambdaMat, NumericVector lambda1, double minimum) {
  List ret;

  NumericVector lambda2, lambda3, diff12, diff23, diff34;
  int mini, minj;
  double val;
  NumericMatrix G(nGrid, nGrid);

  for (int i = 0; i < (nGrid-1); i++) {
    lambda2 = LambdaMat(_, i);

    for (int j = (i+1); j < nGrid; j++) {
      lambda3 = LambdaMat( _, j);
      diff12 = lambda1 - lambda2;
      diff12 = diff12 / std::sqrt(static_cast<double>(sum(diff12 * diff12)));

      diff23 = lambda2 - lambda3;
      diff23 = diff23 / std::sqrt(static_cast<double>(sum(diff23 * diff23)));

      diff34 = lambda3;
      diff34 = diff34 / std::sqrt(static_cast<double>(sum(diff34 * diff34)));

      val = std::abs(sum(diff12 * diff23)) + std::abs(sum(diff23 * diff34)) + std::abs(sum(diff12 * diff34));

      G(i, j) = val;

      if (val < minimum) {
        minimum = val;
        mini = i;
        minj = j;
      }
    }
  }

  ret["G"] = G;
  ret["mini"] = mini + 1;
  ret["minj"] = minj + 1;
  return ret;
}

// [[Rcpp::export]]
List min3Lambda(int i, int nGrid, NumericMatrix LambdaMat, NumericVector lambda1, NumericVector lambda2, double minimum) {
  List ret;
  double val;
  NumericVector lambda3, lambda4, diff12, diff23, diff34, diff45;

  int mini, minj, mink;
  NumericMatrix G(nGrid, nGrid);

  for (int j = (i+1); j < (nGrid - 1); j++) {
    lambda3 = LambdaMat(_,j);
    for (int k = (j + 1); k < nGrid; k++) {
      lambda4 = LambdaMat(_,k);

      diff12 = lambda1 - lambda2;
      diff12 = diff12 / std::sqrt(static_cast<double>(sum(diff12 * diff12)));

      diff23 = lambda2 - lambda3;
      diff23 = diff23 / std::sqrt(static_cast<double>(sum(diff23 * diff23)));

      diff34 = lambda3 - lambda4;
      diff34 = diff34 / std::sqrt(static_cast<double>(sum(diff34 * diff34)));

      diff45 = lambda4;
      diff45 = diff45 / std::sqrt(static_cast<double>(sum(diff45 * diff45)));

      val = std::abs(sum(diff12 * diff23)) + std::abs(sum(diff12 * diff34)) +
        std::abs(sum(diff12 * diff45)) + std::abs(sum(diff23 * diff34)) +
        std::abs(sum(diff23 * diff45)) + std::abs(sum(diff34 * diff45));

      G(j, k) = val;

      if (val < minimum) {
        minimum = val;
        mini = i;
        minj = j;
        mink = k;
      }
    }
  }

  ret["G"] = G;
  ret["mini"] = mini + 1;
  ret["minj"] = minj + 1;
  ret["mink"] = mink + 1;
  ret["minimum"] = minimum;
  return ret;
}







































