#include "CommonFunction.h"
#include "Logit.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// matrix Pi
mat mPiLOGIT_cpp(int n,
                 int ng,
                 int nx,
                 NumericVector pi,                
                 NumericMatrix taux,
                 NumericMatrix X,
                 int refgr){
  // matrix  differential of Pi**2 or theta**2
  mat mPi;
  if (nx==1){
    mat mPitmp(ng-1, ng-1);
    for (int k = 0; k < ng-1; ++k){
      for (int l = 0; l < ng-1; ++l){
        if (k == l){
          mPitmp(k, l) = -sum(taux(_, k)/pow(pi[k], 2)+taux(_, ng-1)/pow(pi[ng-1], 2));
        }else{
          mPitmp(k, l) = -sum(taux(_, ng-1)/pow(pi[ng-1], 2));
        }
      }
    }
    mPi = mPitmp;
  }else{
    mat mPitmp((ng-1)*nx, (ng-1)*nx);
    for (int k = 1; k < ng; ++k){
      for (int kp = 1; kp < ng; ++kp){
        for (int l = 0; l < nx; ++l){
          for (int lp = 0; lp < nx; ++lp){
            double vtmp = 0;
            if (k == kp){
              for (int i = 0; i < n; ++i){
                double tmpPiik = piikIntern_cpp(pi, i, k, ng, X); 
                vtmp -= tmpPiik*(1-tmpPiik)*X(i, l)*X(i, lp); 
              }
            }else{
              for (int i = 0; i < n; ++i){
                vtmp += piikIntern_cpp(pi, i, k, ng, X)*piikIntern_cpp(pi, i, kp, ng, X)*X(i, l)*X(i, lp);
              }
            }
            mPitmp((k - 1)*nx+l, (kp - 1)*nx+lp) = vtmp;
          }
        }
      }
    }
  mPi= mPitmp;
  }
  return(mPi);
}

// -------------------------------------------------------------------
// Function exp
// -------------------------------------------------------------------
double fexp_cpp(int k,
                int i,
                int t,
                IntegerVector nbeta, 
                IntegerVector nbetacum, 
                NumericMatrix A, 
                NumericVector beta, 
                Nullable<NumericMatrix> TCOV,
                int period, 
                Nullable<NumericVector> deltainit, 
                Nullable<IntegerVector> ndeltacuminit, 
                int nw){
  
  double tmp;
  NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
  NumericVector deltak;
  if (nw!=0){
    NumericVector delta(deltainit.get());
    IntegerVector ndeltacum(ndeltacuminit.get());
    deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
  }
  
  //tmp = exp(sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,t]**(0:(nbeta[k]-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum))
  //double muikt = 0;
  NumericVector vtmp2;
  for (int po = 0; po < nbeta[k]; ++po){
    vtmp2.push_back(pow(A(i,t), po));
  }
  tmp = exp(sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, t, k));
  
  return(tmp/pow(1+tmp, 2));
}
// matrix beta
mat mbetaLOGIT_cpp(int i,
                   int t,
                   int ng,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericVector beta,
                   NumericMatrix taux,
                   IntegerVector nbetacum,
                   Nullable<NumericMatrix> TCOV,
                   int period,
                   Nullable<NumericVector> delta,
                   Nullable<IntegerVector> ndeltacum,
                   int nw){
  NumericMatrix res(sum(nbeta), sum(nbeta));
  for (int k = 0; k < ng; ++k){
    for (int l =  nbetacum[k]; l <  nbetacum[k+1]; ++l){
      for (int lp =  nbetacum[k]; lp <  nbetacum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*pow(A(i, t), l- nbetacum[k])*pow(A(i, t), lp-nbetacum[k])*fexp_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw);
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix delta
mat mdeltaLOGIT_cpp(int i,
                    int t,
                    int ng,
                    IntegerVector nbeta,
                    NumericMatrix A,
                    NumericVector beta,
                    NumericMatrix taux,
                    IntegerVector nbetacum,
                    NumericMatrix TCOV,
                    int period,
                    NumericVector delta,
                    IntegerVector ndeltacum,
                    int nw){
  NumericMatrix res(nw*ng, nw*ng);
  for (int k = 0; k < ng; ++k){
    for (int l =  ndeltacum[k]; l <  ndeltacum[k+1]; ++l){
      for (int lp =  ndeltacum[k]; lp <  ndeltacum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*TCOV(i, t + (l - ndeltacum[k])*period)*TCOV(i, t + (lp - ndeltacum[k])*period)*fexp_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw);
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix beta delta
mat mbetadeltaLOGIT_cpp(int i,
                        int t,
                        int ng,
                        IntegerVector nbeta,
                        NumericMatrix A,
                        NumericVector beta,
                        NumericMatrix taux,
                        IntegerVector nbetacum,
                        NumericMatrix TCOV,
                        int period,
                        NumericVector delta,
                        IntegerVector ndeltacum,
                        int nw){
  NumericMatrix res(sum(nbeta), ng*nw);
  for (int k = 0; k < ng; ++k){
    for (int l =  nbetacum[k]; l <  nbetacum[k+1]; ++l){
      for (int lp =  ndeltacum[k]; lp <  ndeltacum[k+1]; ++lp){
        res(l , k) = -taux(i, k)*TCOV(i, t + (lp - ndeltacum[k])*period)*pow(A(i, t), l-nbetacum[k])*fexp_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw);
      }
    }
  }
  return(as<arma::mat>(res));
}

// function Bikl
double BLiklLOGIT_cpp(int i,
                     int k,
                     int l,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     int period,
                     NumericVector beta,
                     NumericMatrix taux,
                     IntegerVector nbetacum,
                     Nullable<NumericMatrix> TCOVinit,
                     Nullable<NumericVector> deltainit,
                     Nullable<IntegerVector> ndeltacuminit,
                     int nw){
  double res = 0;
  NumericMatrix TCOV;
  IntegerVector ndeltacum;
  NumericVector delta;
  if (TCOVinit.isNotNull()){
    NumericMatrix tmp1(TCOVinit.get());
    IntegerVector tmp2(ndeltacuminit.get());
    NumericVector tmp3(deltainit.get());
    TCOV = tmp1;
    ndeltacum = tmp2;
    delta = tmp3;
  }
  NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
  NumericVector deltak(nw);
  if (TCOVinit.isNotNull()){
    deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
  }
  for (int t = 0; t < period; ++t){
    double tmp = 0;
    for (int po = 0; po < nbeta[k]; ++po){
      tmp += pow(A(i,t), po)*betak[po];
    }
    double tmp2 = exp(tmp + WitEM_cpp(TCOV, period, deltak, nw, i, t, k));
    res += pow(A(i, t), l)*(Y(i,t)-tmp2/(1+tmp2));
  }
  return(res);
}

// function Dikl
double DLiklLOGIT_cpp(int i,
                     int k,
                     int l,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     int period,
                     NumericVector beta,
                     NumericMatrix taux,
                     IntegerVector nbetacum,
                     Nullable<NumericMatrix> TCOVinit,
                     Nullable<NumericVector> deltainit,
                     Nullable<IntegerVector> ndeltacuminit,
                     int nw){
  double res = 0;
  NumericMatrix TCOV;
  IntegerVector ndeltacum;
  NumericVector delta;
  if (TCOVinit.isNotNull()){
    NumericMatrix tmp1(TCOVinit.get());
    IntegerVector tmp2(ndeltacuminit.get());
    NumericVector tmp3(deltainit.get());
    TCOV = tmp1;
    ndeltacum = tmp2;
    delta = tmp3;
  }
  NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
  NumericVector deltak(nw);
  if (TCOVinit.isNotNull()){
    deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
  }
  for (int t = 0; t < period; ++t){
    double muikt = 0;
    for (int po = 0; po < nbeta[k]; ++po){
      muikt += pow(A(i,t), po)*betak[po];
    }
    double tmp2 = exp(muikt + WitEM_cpp(TCOV, period, deltak, nw, i, t, k));
    res += TCOV(i, t + l*period)*(Y(i,t)-tmp2/(1+tmp2));
  }
  return(res);
}

// matrix cov Pi Betak
mat covPiBetakLOGIT_cpp(int k,
                        int ng,
                        int n,
                        IntegerVector nbeta,
                        NumericMatrix A,
                        NumericMatrix Y,
                        int period,
                        NumericVector beta,
                        NumericMatrix taux,
                        IntegerVector nbetacum,
                        Nullable<NumericMatrix> TCOVinit,
                        Nullable<NumericVector> deltainit,
                        Nullable<IntegerVector> ndeltacuminit,
                        int nw,
                        NumericVector pi){
  NumericMatrix res(ng-1, nbeta[k]);
  for (int kp = 0; kp < (ng-1); ++kp){
    for (int l = 0; l < nbeta[k]; ++l){
      double tmp = 0;
      if (kp == k){
        for (int i = 0; i < n; ++i){
          tmp += BLiklLOGIT_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,kp)*((1-taux(i,kp))/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp += BLiklLOGIT_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,k)*(-taux(i,kp)/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}

// matrix cov Pi Deltak
mat covPiDeltakLOGIT_cpp(int k,
                         int ng,
                         int n,
                         IntegerVector nbeta,
                         NumericMatrix A,
                         NumericMatrix Y,
                         int period,
                         NumericVector beta,
                         NumericMatrix taux,
                         IntegerVector nbetacum,
                         Nullable<NumericMatrix> TCOVinit,
                         Nullable<NumericVector> deltainit,
                         Nullable<IntegerVector> ndeltacuminit,
                         int nw,
                         NumericVector pi){
  NumericMatrix res(ng-1, nw);
  for (int kp = 0; kp < (ng-1); ++kp){
    for (int l = 0; l < nw; ++l){
      double tmp = 0;
      if (kp == k){
        for (int i = 0; i < n; ++i){
          tmp += DLiklLOGIT_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,kp)*((1-taux(i,kp))/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp += DLiklLOGIT_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,k)*(-taux(i,kp)/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}

// matrix cov Betak Betal
mat covBetakBetalLOGIT_cpp(int k,
                           int l,
                           int n,
                           IntegerVector nbeta,
                           NumericMatrix A,
                           NumericMatrix Y,
                           int period,
                           NumericVector beta,
                           NumericMatrix taux,
                           IntegerVector nbetacum,
                           Nullable<NumericMatrix> TCOVinit,
                           Nullable<NumericVector> deltainit,
                           Nullable<IntegerVector> ndeltacuminit,
                           int nw){
  NumericMatrix res(nbeta[k], nbeta[l]);
  if (k ==l){
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nbeta[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += BLiklLOGIT_cpp(i, k, p, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*BLiklLOGIT_cpp(i, k, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i, k));
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nbeta[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp -= BLiklLOGIT_cpp(i, k, p, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*BLiklLOGIT_cpp(i, l, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*taux(i, l);
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Betak Deltal
mat covBetakDeltalLOGIT_cpp(int k,
                            int l,
                            int n,
                            IntegerVector nbeta,
                            NumericMatrix A,
                            NumericMatrix Y,
                            int period,
                            NumericVector beta,
                            NumericMatrix taux,
                            IntegerVector nbetacum,
                            Nullable<NumericMatrix> TCOVinit,
                            Nullable<NumericVector> deltainit,
                            Nullable<IntegerVector> ndeltacuminit,
                            int nw){
  NumericMatrix res(nbeta[k], nw);
  if (k == l){
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += BLiklLOGIT_cpp(i, k, p, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DLiklLOGIT_cpp(i, k, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i, k));
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp -= BLiklLOGIT_cpp(i, k, p, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DLiklLOGIT_cpp(i, l, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*taux(i, l);
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}


// matrix cov Deltak Deltal
mat covDeltakDeltalLOGIT_cpp(int k,
                             int l,
                             int n,
                             IntegerVector nbeta,
                             NumericMatrix A,
                             NumericMatrix Y,
                             int period,
                             NumericVector beta,
                             NumericMatrix taux,
                             IntegerVector nbetacum,
                             Nullable<NumericMatrix> TCOVinit,
                             Nullable<NumericVector> deltainit,
                             Nullable<IntegerVector> ndeltacuminit,
                             int nw){
  NumericMatrix res(nw, nw);
  if (k ==l){
    for (int p = 0; p < nw; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += DLiklLOGIT_cpp(i, k, p, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DLiklLOGIT_cpp(i, k, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i, k));
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nw; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp -= DLiklLOGIT_cpp(i, k, p, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DLiklLOGIT_cpp(i, l, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*taux(i, l);
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}

// matrix cov Pi
mat covPiLOGIT_cpp(int n,
                   int ng,
                   int nx,
                   NumericVector pi,
                   NumericMatrix X,
                   NumericMatrix taux){
  mat res;
  if (nx == 1){
    mat covPi(ng-1, ng-1);
    for (int k = 0; k < ng-1; ++k){
      for (int l = 0; l < ng-1; ++l){
      double tmp =0;
        if (k == l){
          for (int i = 0; i < n; ++i){
            tmp += taux(i, k)*(1-taux(i, k))/pow(pi[k], 2)+taux(i, ng-1)*(1-taux(i, ng-1))/pow(pi[ng-1], 2)-2*taux(i, k)*taux(i, ng-1)/(pi[k]*pi[ng-1]);
          }
        }else{
          for (int i = 0; i < n; ++i){
            tmp += -taux(i, k)*taux(i, l)/(pi[k]*pi[l])+taux(i, k)*taux(i, ng-1)/(pi[ng-1]*pi[k])+taux(i, l)*taux(i, ng-1)/(pi[l]*pi[ng-1])+taux(i, ng-1)*(1-taux(i, ng-1))/pow(pi[ng-1],2);
          }
        }
        covPi(k, l) = tmp;
      }
    }
    res= covPi;
  }else{
    mat covPi((ng-1)*nx, (ng-1)*nx);
    for (int k = 1; k < ng; ++k){
      for (int l = 1; l < ng; ++l){
        for(int p = 0; p < nx; ++p){
          for (int q = 0; q < nx; ++q){
            double tmp =0;
            if (k == l){
              for (int i = 0; i < n; ++i){
                tmp += taux(i, k)*(1-taux(i, k))*X(i, p)*X(i, q);
              }
            }else{
              for (int i = 0; i < n; ++i){
                tmp -= taux(i, k)*taux(i, l)*X(i, p)*X(i, q);
              }
            }
            covPi((k - 1)*nx+p, (l - 1)*nx+q) = tmp;
          }
        }
      }
    }
    res= covPi;
  }

  return(res);
}
// matrix cov Pi Beta
mat covPiBetaLOGIT_cpp(int n,
                       int ng,
                       int nx,
                       NumericVector pi,
                       NumericVector beta,
                       Nullable<NumericVector> delta,
                       NumericMatrix X,
                       NumericMatrix taux,
                       IntegerVector nbeta,
                       NumericMatrix A,
                       NumericMatrix Y,
                       int period,
                       IntegerVector nbetacum,
                       Nullable<NumericMatrix> TCOV,
                       Nullable<IntegerVector> ndeltacum,
                       int nw){
  mat res;
  if (nx == 1){
    mat covPiBeta;
    for (int k = 0; k < ng-1; ++k){
      covPiBeta = join_rows(covPiBeta, covPiBetakLOGIT_cpp(k, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi));
    }
    mat mtmp(ng-1, nbeta[ng-1]);
    for (int kp = 0; kp < ng-1; ++kp){
      for (int l = 0; l < nbeta[ng-1]; ++l){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += BLiklLOGIT_cpp(i, ng-1, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux(i, ng-1)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, kp)/pi[kp]);
        }
        mtmp(kp, l) = tmp;
      }
    }
    res = join_rows(covPiBeta, mtmp);
  }else{
    mat covPiBeta((ng-1)*nx, sum(nbeta));
    for (int k = 0; k < ng-1; ++k){
      for (int l = 0; l < ng; ++l){
        for (int p = 0; p < nx; ++p){
          for (int q = 0; q < nbeta[l]; ++q){
            double tmp = 0;
            if (k == l){
              for (int i = 0; i < n; ++i){
                tmp += taux(i, k)*(1-taux(i, k))*X(i, p)*BLiklLOGIT_cpp(i, k, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw);
              }
            }else{
              for (int i = 0; i < n; ++i){
                tmp -= taux(i, k)*taux(i, l)*X(i, p)*BLiklLOGIT_cpp(i, l, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw);
              }
            }
            covPiBeta(k*nx + p, nbetacum[l] + q) = tmp;
          }
        }
      }
    }
    res = covPiBeta;
  }

  return(res);
}

// matrix cov Pi Delta
mat covPiDeltaLOGIT_cpp(int n,
                       int ng,
                       int nx,
                       NumericVector pi,
                       NumericVector beta,
                       Nullable<NumericVector> delta,
                       NumericMatrix X,
                       NumericMatrix taux,
                       IntegerVector nbeta,
                       NumericMatrix A,
                       NumericMatrix Y,
                       int period,
                       IntegerVector nbetacum,
                       Nullable<NumericMatrix> TCOV,
                       Nullable<IntegerVector> ndeltacuminit,
                       int nw){
  mat res;
  IntegerVector ndeltacum(ndeltacuminit);
  if (nx == 1){
    mat covPiDelta;
    for (int k = 0; k < ng-1; ++k){
      covPiDelta = join_rows(covPiDelta, covPiDeltakLOGIT_cpp(k, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi));
    }
    mat mtmp(ng-1, nw);
    for (int kp = 0; kp < ng-1; ++kp){
      for (int l = 0; l < nw; ++l){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += DLiklLOGIT_cpp(i, ng-1, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux(i, ng-1)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, kp)/pi[kp]);
        }
        mtmp(kp, l) = tmp;
      }
    }
    res = join_rows(covPiDelta, mtmp);
  }else{
    mat covPiDelta((ng-1)*nx*nw, nw*ng);
    for (int k = 0; k < ng-1; ++k){
      for (int l = 0; l < ng; ++l){
        for (int p = 0; p < nx; ++p){
          for (int q = 0; q < nw; ++q){
            double tmp = 0;
            if (k == l){
              for (int i = 0; i < n; ++i){
                tmp += taux(i, k)*(1-taux(i, k))*X(i, p)*DLiklLOGIT_cpp(i, k, q, nbeta, A, Y, period, beta,  taux, nbetacum, TCOV, delta, ndeltacum, nw);
              }
            }else{
              for (int i = 0; i < n; ++i){
                tmp -= taux(i, k)*taux(i, l)*X(i, p)*DLiklLOGIT_cpp(i, l, q, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw);
              }
            }
            covPiDelta(k*nx + p, ndeltacum[l] + q) = tmp;
          }
        }
      }
    }
    res = covPiDelta;
  }

  return(res);
}

// matrix cov Beta Delta
mat covBetaDeltaLOGIT_cpp(int n,
                        int ng,
                        int nx,
                        NumericVector pi,
                        NumericVector beta,
                        Nullable<NumericVector> delta,
                        NumericMatrix X,
                        NumericMatrix taux,
                        IntegerVector nbeta,
                        NumericMatrix A,
                        NumericMatrix Y,
                        int period,
                        IntegerVector nbetacum,
                        Nullable<NumericMatrix> TCOV,
                        Nullable<IntegerVector> ndeltacum,
                        int nw){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covBetakDeltalLOGIT_cpp(k, l, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}

// matrix cov Delta
mat covDeltaLOGIT_cpp(int n,
                          int ng,
                          int nx,
                          NumericVector pi,
                          NumericVector beta,
                          Nullable<NumericVector> delta,
                          NumericMatrix X,
                          NumericMatrix taux,
                          IntegerVector nbeta,
                          NumericMatrix A,
                          NumericMatrix Y,
                          int period,
                          IntegerVector nbetacum,
                          Nullable<NumericMatrix> TCOV,
                          Nullable<IntegerVector> ndeltacum,
                          int nw){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covDeltakDeltalLOGIT_cpp(k, l, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}

// matrix cov Beta
mat covBetaLOGIT_cpp(int n,
                      int ng,
                      int nx,
                      NumericVector pi,
                      NumericVector beta,
                      Nullable<NumericVector> delta,
                      NumericMatrix X,
                      NumericMatrix taux,
                      IntegerVector nbeta,
                      NumericMatrix A,
                      NumericMatrix Y,
                      int period,
                      IntegerVector nbetacum,
                      Nullable<NumericMatrix> TCOV,
                      Nullable<IntegerVector> ndeltacum,
                      int nw){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covBetakBetalLOGIT_cpp(k, l, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}

// Main function
// [[Rcpp::export]]
arma::vec IEMLOGIT_cpp(NumericVector param,
                  int ng,
                  int nx,
                  IntegerVector nbeta,
                  int n,
                  NumericMatrix A,
                  NumericMatrix Y,
                  NumericMatrix X,
                  Nullable<NumericMatrix> TCOV,
                  int nw,
                  int refgr){
  int period = A.ncol();

  // Compute Matirx information
  NumericVector pi(ng);
  NumericVector beta;
  NumericVector delta;
  IntegerVector ndeltacum(ng);
  IntegerVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-1)];
    beta = param[Range(ng,ng+sum(nbeta)-1)];
    if (param.length() > ng*nx+sum(nbeta)){
      delta = param[Range(ng+sum(nbeta), param.length() - 1)];
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
    }
  }else{
    pi = param[Range(0,ng*nx-1)];
    beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
    if (param.length() > ng*nx+sum(nbeta)){
      delta = param[Range(ng*nx+sum(nbeta), param.length() - 1)];
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
    }
  }
  NumericMatrix  taux = ftauxLOGIT_cpp(pi, beta, ng, nbeta, n, A, Y, TCOV, delta, nw, nx, X);

  // matrix -B
  mat B(sum(nbeta)+ng*nw, sum(nbeta)+ng*nw);
  B.fill(0);
  // we take of the dimension of mPi
  if (nw != 0){
    NumericMatrix TCOVv(TCOV.get());
    for (int i = 0; i < n; ++i){
      for (int t = 0; t < period; ++t){
        mat tmp2 = mbetadeltaLOGIT_cpp(i, t, ng, nbeta, A, beta, taux, nbetacum, TCOVv, period, delta, ndeltacum, nw);
        B += join_cols(join_rows(mbetaLOGIT_cpp(i, t, ng, nbeta, A, beta, taux, nbetacum, TCOVv, period, delta, ndeltacum, nw), tmp2),
                       join_rows(trans(tmp2), mdeltaLOGIT_cpp(i, t, ng, nbeta, A, beta, taux, nbetacum, TCOVv, period, delta, ndeltacum, nw))
        );
      }
    }
  }else{
    for (int i = 0; i < n; ++i){
      for (int t = 0; t < period; ++t){
        B += mbetaLOGIT_cpp(i, t, ng, nbeta, A, beta, taux, nbetacum, TCOV, period, delta, ndeltacum, nw);
      }
    }
  }
  B = join_rows(zeros(sum(nbeta)+nw*ng, (ng-1)*nx), B);
  B = join_cols(join_rows(mPiLOGIT_cpp(n, ng, nx, pi, taux , X, refgr), zeros((ng-1)*nx, sum(nbeta)+nw*ng)), B);

  // Matrix of covariance of score function
  mat cov ;
  if (nw != 0){
    mat tmp1 = covPiBetaLOGIT_cpp(n, ng, nx, pi, beta, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp2 = covPiDeltaLOGIT_cpp(n, ng, nx, pi, beta, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp4 = covBetaDeltaLOGIT_cpp(n, ng, nx, pi, beta, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);

    cov = join_rows(covPiLOGIT_cpp(n, ng, nx, pi, X, taux), tmp1, tmp2);
    cov = join_cols(cov, join_rows(trans(tmp1), covBetaLOGIT_cpp(n, ng, nx, pi, beta, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw), tmp4));
    cov = join_cols(cov, join_rows(trans(tmp2), trans(tmp4), covDeltaLOGIT_cpp(n, ng, nx, pi, beta, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw)));
   }else{
    mat tmp1 = covPiBetaLOGIT_cpp(n, ng, nx, pi, beta, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
 
    cov = join_rows(covPiLOGIT_cpp(n, ng, nx, pi, X, taux), tmp1);
    cov = join_cols(cov, join_rows(trans(tmp1), covBetaLOGIT_cpp(n, ng, nx, pi, beta, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw)));
 
  }
  // Information matrix of Fisher
  mat IEM = -B - cov;

  return(sqrt(diagvec(inv(IEM))));
}
