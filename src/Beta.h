#ifndef __BETAUtils__
#define __BETAUtils__
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double gkBETA_cpp(List beta,
                  List phi,
                  int i,
                  int k,
                  IntegerVector nbeta,
                  IntegerVector nphi,
                  NumericMatrix A,
                  NumericMatrix Y,
                  Nullable<NumericMatrix> TCOV,
                  Nullable<List> delta,
                  int nw);

#endif //  __BETAUtils__
