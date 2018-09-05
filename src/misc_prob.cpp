#include <Rcpp.h>
#include <cmath>
using namespace std;

// [[Rcpp::export]]
double meanC(Rcpp::NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

// // [[Rcpp::export]]
// //------------------------------------------------
// // draw from negative binomial distribution with mean lambda and variance gamma*lambda (gamma must be >1)
// int rnbinom1(double lambda, double gamma) {
//   double ret = R::rnbinom(lambda/(gamma-1), 1/gamma);
//   return(ret);
// }
