#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double f(double x) {
  return exp(-abs(x));
}

//[[Rcpp::export]]
NumericVector rwC_Metropolis (double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0; 
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= (f(y[0]) / f(x[i-1]))){
      x[i] = y[0];
    }
    else { 
      x[i] = x[i-1]; 
    }
  }
  return(x);
} 