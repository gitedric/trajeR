#ifndef __ZIPUtils__
#define __ZIPUtils__
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double fSikt_cpp(NumericVector pi,
                 NumericVector beta,
                 NumericVector nu,
                 int k,
                 int i, 
                 int t,
                 IntegerVector nbeta,
                 IntegerVector nnu,
                 int n,
                 NumericMatrix A,
                 NumericMatrix Y,
                 Nullable<NumericMatrix> TCOV,
                 Nullable<NumericVector> delta,
                 int nw,
                 Nullable<IntegerVector> ndeltacum,
                 int period,
                 IntegerVector nbetacum,
                 IntegerVector nnucum);

NumericMatrix ftauxZIP_cpp(NumericVector pi,
                           NumericVector beta,
                           NumericVector nu,
                           int ng,
                           IntegerVector nbeta,
                           IntegerVector nnu,
                           int n,
                           NumericMatrix A,
                           NumericMatrix Y,
                           Nullable<NumericMatrix> TCOV,
                           Nullable<NumericVector> delta,
                           int nw,
                           int nx,
                           NumericMatrix X);


#endif //  __ZIPUtils__