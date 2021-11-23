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

double likelihoodLOGIT_cpp(NumericVector param,
                           int ng, 
                           int nx,
                           int n,
                           IntegerVector nbeta,
                           NumericMatrix A,
                           NumericMatrix Y,
                           NumericMatrix X,
                           Nullable<NumericMatrix> TCOV,
                           int nw);

NumericVector difLLOGIT_cpp(NumericVector param,
                            int ng, 
                            int nx,
                            int n,
                            IntegerVector nbeta,
                            NumericMatrix A,
                            NumericMatrix Y,
                            NumericMatrix X,
                            Nullable<NumericMatrix> TCOV,
                            int nw);

#endif //  __LOGITUtils__