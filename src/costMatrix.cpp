#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix costMatrix(NumericMatrix muhat, NumericMatrix nuhat, int q, int p) {
  int muhatSize = muhat.ncol();
  int nuhatSize = nuhat.ncol();
  NumericMatrix out(muhatSize, nuhatSize);
  for (int i = 0; i < muhatSize; i++) {
    for (int j = 0; j < nuhatSize; j++) {
      double temp = sum(pow(abs(muhat(_, i) - nuhat(_, j)), q));
      out(i, j) = pow(pow(temp, 1.0/q), p);
    }
  }
  return out;
}
