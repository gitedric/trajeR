#ifndef __CensoredNormalUtils__
#define __CensoredNormalUtils__
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double gkCNORM_cpp(List beta,
                   NumericVector sigma,
                   int i,
                   int k,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   double ymin,
                   double ymax, 
                   Nullable<NumericMatrix> TCOV,
                   Nullable<List> delta,
                   int nw);

double gkalpha_cpp(List beta,
                   NumericVector alpha,
                   int i,
                   int k,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   double ymin,
                   double ymax, 
                   Nullable<NumericMatrix> TCOV,
                   Nullable<List> delta,
                   int nw);

NumericMatrix ftauxCNORM_cpp(NumericVector pi,
                             NumericVector beta,
                             NumericVector sigma,
                             int ng,
                             IntegerVector nbeta,
                             int n,
                             NumericMatrix A,
                             NumericMatrix Y,
                             double ymin,
                             double ymax, 
                             Nullable<NumericMatrix> TCOV,
                             Nullable<NumericVector> delta,
                             int nw,
                             int nx,
                             NumericMatrix X);
  
#endif //  __CensoredNormalUtils__