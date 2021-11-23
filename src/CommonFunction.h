#ifndef __CommonFunction__
#define __CommonFunction__
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double facto(double nb);

NumericVector muikt_cpp(NumericVector beta,
                        int nbeta,
                        int i,
                        int period,
                        NumericMatrix A, 
                        Nullable<NumericMatrix> TCOV,
                        Nullable<List> delta,
                        int nw,
                        int k);

double piikIntern_cpp(NumericVector theta,
                      int i,
                      int k,
                      int ng,
                      NumericMatrix X);

double prodvect(NumericVector vec);

NumericMatrix submat_cpp(NumericMatrix X, LogicalVector condition);

NumericVector findtheta_cpp(NumericVector theta, 
                            NumericMatrix taux, 
                            NumericMatrix X, 
                            int n, 
                            int ng, 
                            int nx, 
                            int period, 
                            bool EMIRLS, 
                            int refgr);

double Wit_cpp(Nullable<NumericMatrix> TCOV,
               int period,
               Nullable<List> delta, 
               int nw,
               int i,
               int t,
               int k);

double WitEM_cpp(Nullable<NumericMatrix> TCOV,
               int period,
               Nullable<NumericVector> delta, 
               int nw,
               int i,
               int t,
               int k);
  
#endif //  __CommonFunction__
