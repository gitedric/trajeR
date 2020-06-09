#ifndef __LOGITUtils__
#define __LOGITUtils__
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double gkLOGIT_cpp(List beta,
                   int i,
                   int k,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   Nullable<NumericMatrix> TCOV,
                   Nullable<List> delta,
                   int nw);

NumericMatrix ftauxLOGIT_cpp(NumericVector pi,
                             NumericVector beta,
                             int ng,
                             IntegerVector nbeta,
                             int n,
                             NumericMatrix A,
                             NumericMatrix Y,
                             Nullable<NumericMatrix> TCOV,
                             Nullable<NumericVector> delta,
                             int nw,
                             int nx,
                             NumericMatrix X);

#endif //  __LOGITUtils__