#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix costMatrix(NumericMatrix muhat, NumericMatrix nuhat, float q, float p) {
  int muhatSize = muhat.ncol();
  int nuhatSize = nuhat.ncol();
  NumericMatrix out(muhatSize, nuhatSize);
  for (int i = 0; i < muhatSize; i++) {
    for (int j = 0; j < nuhatSize; j++) {
      out(i, j) = pow(pow(mean(abs(pow(muhat(_, i) - nuhat(_, j), q))), 1/q), p);
    }
  }
  return out;
}

