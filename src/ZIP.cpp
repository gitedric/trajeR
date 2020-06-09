#include "CommonFunction.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::interfaces(r, cpp)]]

// ----------------------------------------------------------------------------
// gk ZIP
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double gkZIP_cpp(List beta,
                 List nu,
                   int i,
                   int k,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   double ymin,
                   double ymax, 
                   Nullable<NumericMatrix> TCOV,
                   Nullable<List> delta,
                   int nw){
  int period = A.ncol();
  k = k-1;
  i = i-1;
  NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
  double res = 1;
  for (int ind = 0; ind < period; ++ind){
    if (Y(i, ind) <= ymin){
      res = res*(R::pnorm((Y(i, ind)-muikt[ind])/sigma[k], 0.0, 1.0, true, false));
    }else if (Y(i, ind) >= ymax){
      res = res*(R::pnorm(-(Y(i, ind)-muikt[ind])/sigma[k], 0.0, 1.0, true, false));
    }else{
      res = res*(R::dnorm((Y(i, ind)-muikt[ind])/sigma[k], 0.0, 1.0, false)/sigma[k]);
    }
  }
  return(res);
}