#include "CommonFunction.h"
#include "CensoredNormal.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// matrix Pi
mat mPiCNORM_cpp(int n,
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
    mat mPitmp(ng*nx, ng*nx);
    for (int k = 0; k < ng; ++k){
      for (int kp = 0; kp < ng; ++kp){
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
            mPitmp(k*nx+l, kp*nx+lp) = vtmp;
          }
        }
      }
    }
    mPitmp.shed_cols((refgr-1)*nx, (refgr-1)*nx+nx-1);
    mPitmp.shed_rows((refgr-1)*nx, (refgr-1)*nx+nx-1);
    mPi= mPitmp;
  }
  return(mPi);
}
// matrix beta
mat mbetaCNORM_cpp(int i,
                   int t,
                   int ng,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericVector sigma,
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
        res(l , lp) = -taux(i, k)*pow(A(i, t), l- nbetacum[k])*pow(A(i, t), lp-nbetacum[k])/pow(sigma[k], 2);
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix delta
mat mdeltaCNORM_cpp(int i,
                    int t,
                    int ng,
                    IntegerVector nbeta,
                    NumericMatrix A,
                    NumericVector sigma,
                    NumericMatrix taux,
                    IntegerVector nbetacum,
                    Nullable<NumericMatrix> TCOVinit,
                    int period,
                    Nullable<NumericVector> delta,
                    Nullable<IntegerVector> ndeltacuminit,
                    int nw){
  NumericMatrix TCOV;  
  IntegerVector ndeltacum;
  if (TCOVinit.isNotNull()){
    NumericMatrix tmp1(TCOVinit.get());  
    IntegerVector tmp2(ndeltacuminit.get());
    TCOV = tmp1;
    ndeltacum = tmp2;
  }
  NumericMatrix res(nw*ng, nw*ng);
  for (int k = 0; k < ng; ++k){
    for (int l =  ndeltacum[k]; l <  ndeltacum[k+1]; ++l){
      for (int lp =  ndeltacum[k]; lp <  ndeltacum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*TCOV(i, t + (l - ndeltacum[k])*period)*TCOV(i, t + (lp - ndeltacum[k])*period)/pow(sigma[k], 2);
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix beta delta
mat mbetadeltaCNORM_cpp(int i,
                        int t,
                        int ng,
                        IntegerVector nbeta,
                        NumericMatrix A,
                        NumericMatrix Y,
                        NumericVector beta,
                        NumericVector sigma,
                        NumericMatrix taux,
                        IntegerVector nbetacum,
                        Nullable<NumericMatrix> TCOVinit,
                        int period,
                        Nullable<NumericVector> deltainit,
                        Nullable<IntegerVector> ndeltacuminit,
                        int nw){
  NumericMatrix TCOV;  
  IntegerVector ndeltacum;
  if (TCOVinit.isNotNull()){
    NumericMatrix tmp1(TCOVinit.get());  
    IntegerVector tmp2(ndeltacuminit.get());
    TCOV = tmp1;
    ndeltacum = tmp2;
  }
  NumericMatrix res(sum(nbeta), ng*nw);
  for (int k = 0; k < ng; ++k){
    for (int l =  nbetacum[k]; l <  nbetacum[k+1]; ++l){
      for (int lp =  ndeltacum[k]; lp <  ndeltacum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*TCOV(i, t + (lp - ndeltacum[k])*period)*pow(A(i, t), l-nbetacum[k])/pow(sigma[k], 2);
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix beta sigma
mat mbetasigmaCNORM_cpp(int i,
                        int t,
                        int ng,
                        IntegerVector nbeta,
                        NumericMatrix A,
                        NumericMatrix Y,
                        NumericVector beta,
                        NumericVector sigma,
                        NumericMatrix taux,
                        IntegerVector nbetacum,
                        Nullable<NumericMatrix> TCOVinit,
                        int period,
                        Nullable<NumericVector> deltainit,
                        Nullable<IntegerVector> ndeltacuminit,
                        int nw){
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
  NumericMatrix res(sum(nbeta), ng);
  for (int k = 0; k < ng; ++k){
    for (int l =  nbetacum[k]; l <  nbetacum[k+1]; ++l){
      NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
      NumericVector deltak;
      if (TCOVinit.isNotNull()){
         deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
      }
      double muikt = 0;
      for (int po = 0; po < nbeta[k]; ++po){
        muikt += pow(A(i,t), po)*betak[po];
      }
      muikt += WitEM_cpp(TCOV, period, deltak, nw, i, t, k);
      res(l , k) = -2*taux(i, k)*pow(A(i, t), l-nbetacum[k])*(Y(i, t)-muikt)/pow(sigma[k], 3);
    }
  }
  return(as<arma::mat>(res));
}
// matrix delta sigma
mat mdeltasigmaCNORM_cpp(int i,
                         int t,
                         int ng,
                         IntegerVector nbeta,
                         NumericMatrix A,
                         NumericMatrix Y,
                         NumericVector beta,
                         NumericVector sigma,
                         NumericMatrix taux,
                         IntegerVector nbetacum,
                         Nullable<NumericMatrix> TCOVinit,
                         int period,
                         Nullable<NumericVector> deltainit,
                         Nullable<IntegerVector> ndeltacuminit,
                         int nw){
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
  NumericMatrix res(nw*ng, ng);
  for (int k = 0; k < ng; ++k){
    for (int l =  ndeltacum[k]; l <  ndeltacum[k+1]; ++l){
      NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
      NumericVector deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
      double muikt = 0;
      for (int po = 0; po < nbeta[k]; ++po){
        muikt += pow(A(i,t), po)*betak[po];
      }
      muikt += WitEM_cpp(TCOV, period, deltak, nw, i, t, k);
      res(l , k) = -2*taux(i, k)*TCOV(i, t + (l - ndeltacum[k])*period)*(Y(i, t)-muikt)/pow(sigma[k], 3);
    }
  }
  return(as<arma::mat>(res));
}
// matrix sigma
mat msigmaCNORM_cpp(int i,
                    int t,
                    int ng,
                    IntegerVector nbeta,
                    NumericMatrix A,
                    NumericMatrix Y,
                    NumericVector beta,
                    NumericVector sigma,
                    NumericMatrix taux,
                    IntegerVector nbetacum,
                    Nullable<NumericMatrix> TCOVinit,
                    int period,
                    Nullable<NumericVector> deltainit,
                    Nullable<IntegerVector> ndeltacuminit,
                    int nw){
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
  NumericMatrix res(ng, ng);
  for (int k = 0; k < ng; ++k){
    NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
    NumericVector deltak(nw);
    if (TCOVinit.isNotNull()){
      deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
    }
    double muikt = 0;
    for (int po = 0; po < nbeta[k]; ++po){
      muikt += pow(A(i,t), po)*betak[po];
    }
    muikt += WitEM_cpp(TCOV, period, deltak, nw, i, t, k);
    res(k , k) = -taux(i, k)*(-pow(sigma[k], 2)+3*pow(Y(i, t)-muikt, 2))/pow(sigma[k], 4);
  }
  return(as<arma::mat>(res));
}
// function Bikl
double BiklCNORM_cpp(int i,
                     int k,
                     int l,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     int period,
                     NumericVector beta,
                     NumericVector sigma,
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
    muikt += WitEM_cpp(TCOV, period, deltak, nw, i, t, k);    
    res += pow(A(i, t), l)*(Y(i,t)-muikt)/pow(sigma[k], 2);
  }
  return(res);
}
// function Dikl
double DiklCNORM_cpp(int i,
                     int k,
                     int l,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     int period,
                     NumericVector beta,
                     NumericVector sigma,
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
    muikt += WitEM_cpp(TCOV, period, deltak, nw, i, t, k);    
    res += TCOV(i, t + l*period)*(Y(i,t)-muikt)/pow(sigma[k], 2);
  }
  return(res);
}
// function Sik
double SikCNORM_cpp(int i,
                     int k,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     int period,
                     NumericVector beta,
                     NumericVector sigma,
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
    muikt += WitEM_cpp(TCOV, period, deltak, nw, i, t, k);    
    res -= (pow(sigma[k], 2)-pow(Y(i,t)-muikt, 2))/pow(sigma[k], 3);
  }
  return(res);
}
// matrix cov Pi Betak
mat covPiBetakCNORM_cpp(int k,
                        int ng,
                        int n,
                        IntegerVector nbeta,
                        NumericMatrix A,
                        NumericMatrix Y,
                        int period,
                        NumericVector beta,
                        NumericVector sigma,
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
          tmp += BiklCNORM_cpp(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,kp)*((1-taux(i,kp))/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp += BiklCNORM_cpp(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,k)*(-taux(i,kp)/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Pi Deltak
mat covPiDeltakCNORM_cpp(int k,
                         int ng,
                         int n,
                         IntegerVector nbeta,
                         NumericMatrix A,
                         NumericMatrix Y,
                         int period,
                         NumericVector beta,
                         NumericVector sigma,
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
          tmp += DiklCNORM_cpp(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,kp)*((1-taux(i,kp))/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp += DiklCNORM_cpp(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,k)*(-taux(i,kp)/pi[kp]+taux(i,ng-1)/pi(ng-1));
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Beta Sigmak
mat covBetaSigmakCNORM_cpp(int k,
                           IntegerVector nbeta,
                           int n,
                           int ng, 
                           NumericMatrix A,
                           NumericMatrix Y,
                           int period,
                           NumericVector beta,
                           NumericVector sigma,
                           NumericMatrix taux,
                           IntegerVector nbetacum,
                           Nullable<NumericMatrix> TCOVinit,
                           Nullable<NumericVector> deltainit,
                           Nullable<IntegerVector> ndeltacuminit,
                           int nw){
  NumericMatrix res(nbeta[k], ng);
  for (int kp = 0; kp < nbeta[k]; ++kp){
    for (int l = 0; l < ng; ++l){
      double tmp = 0;
      for (int i = 0; i < n; ++i){
        if (k == l){
          tmp += BiklCNORM_cpp(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*SikCNORM_cpp(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i,k));
        }else{
          tmp -= BiklCNORM_cpp(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*SikCNORM_cpp(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,k)*taux(i, l);
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Betak Betal
mat covBetakBetalCNORM_cpp(int k,
                           int l,
                           int n,
                           IntegerVector nbeta,
                           NumericMatrix A,
                           NumericMatrix Y,
                           int period,
                           NumericVector beta,
                           NumericVector sigma,
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
          tmp += BiklCNORM_cpp(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*BiklCNORM_cpp(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i, k));  
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nbeta[l]; ++q){
        double tmp = 0;  
        for (int i = 0; i < n; ++i){
          tmp -= BiklCNORM_cpp(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*BiklCNORM_cpp(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*taux(i, l);  
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Betak Deltal
mat covBetakDeltalCNORM_cpp(int k,
                            int l,
                            int n,
                            IntegerVector nbeta,
                            NumericMatrix A,
                            NumericMatrix Y,
                            int period,
                            NumericVector beta,
                            NumericVector sigma,
                            NumericMatrix taux,
                            IntegerVector nbetacum,
                            Nullable<NumericMatrix> TCOVinit,
                            Nullable<NumericVector> deltainit,
                            Nullable<IntegerVector> ndeltacuminit,
                            int nw){
  NumericMatrix res(nbeta[k], nw);
  if (k ==l){
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;  
        for (int i = 0; i < n; ++i){
          tmp += BiklCNORM_cpp(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DiklCNORM_cpp(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i, k));  
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;  
        for (int i = 0; i < n; ++i){
          tmp -= BiklCNORM_cpp(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DiklCNORM_cpp(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*taux(i, l);  
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Deltak Deltal
mat covDeltakDeltalCNORM_cpp(int k,
                             int l,
                             int n,
                             IntegerVector nbeta,
                             NumericMatrix A,
                             NumericMatrix Y,
                             int period,
                             NumericVector beta,
                             NumericVector sigma,
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
          tmp += DiklCNORM_cpp(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DiklCNORM_cpp(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i, k));  
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nw; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;  
        for (int i = 0; i < n; ++i){
          tmp -= DiklCNORM_cpp(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*DiklCNORM_cpp(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*taux(i, l);  
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Delta Sigmak
mat covDeltaSigmakCNORM_cpp(int k,
                            IntegerVector nbeta,
                            int n,
                            int ng, 
                            NumericMatrix A,
                            NumericMatrix Y,
                            int period,
                            NumericVector beta,
                            NumericVector sigma,
                            NumericMatrix taux,
                            IntegerVector nbetacum,
                            Nullable<NumericMatrix> TCOVinit,
                            Nullable<NumericVector> deltainit,
                            Nullable<IntegerVector> ndeltacuminit,
                            int nw){
  NumericMatrix res(nw, ng);
  for (int kp = 0; kp < nw; ++kp){
    for (int l = 0; l < ng; ++l){
      double tmp = 0;
      if (k == l){
        for (int i = 0; i < n; ++i){
          tmp += DiklCNORM_cpp(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*SikCNORM_cpp(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i, k)*(1-taux(i,k));
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp -= DiklCNORM_cpp(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*SikCNORM_cpp(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw)*taux(i,k)*taux(i, l);
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Pi
mat covPiCNORM_cpp(int n,
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
    for (int k = 0; k < ng-1; ++k){
      for (int l = 0; l < ng-1; ++l){
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
            covPi(k*nx+p, l*nx+q) = tmp;
          }
        }
      }
    }
    res= covPi;
  }
  
  return(res);
}
// matrix cov Pi Beta
mat covPiBetaCNORM_cpp(int n,
                       int ng, 
                       int nx,
                       NumericVector pi,
                       NumericVector beta,
                       NumericVector sigma,
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
      covPiBeta = join_rows(covPiBeta, covPiBetakCNORM_cpp(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi));
    }
    mat mtmp(ng-1, nbeta[ng-1]);
    for (int kp = 0; kp < ng-1; ++kp){
      for (int l = 0; l < nbeta[ng-1]; ++l){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += BiklCNORM_cpp(i, ng-1, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux(i, ng-1)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, kp)/pi[kp]);
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
                tmp += taux(i, k)*(1-taux(i, k))*X(i, p)*BiklCNORM_cpp(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw);
              } 
            }else{
              for (int i = 0; i < n; ++i){
                tmp -= taux(i, k)*taux(i, l)*X(i, p)*BiklCNORM_cpp(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw);
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
mat covPiDeltaCNORM_cpp(int n,
                       int ng, 
                       int nx,
                       NumericVector pi,
                       NumericVector beta,
                       NumericVector sigma,
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
      covPiDelta = join_rows(covPiDelta, covPiDeltakCNORM_cpp(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi));
    }
    mat mtmp(ng-1, nw);
    for (int kp = 0; kp < ng-1; ++kp){
      for (int l = 0; l < nw; ++l){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += DiklCNORM_cpp(i, ng-1, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux(i, ng-1)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, kp)/pi[kp]);
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
                tmp += taux(i, k)*(1-taux(i, k))*X(i, p)*DiklCNORM_cpp(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw);
              } 
            }else{
              for (int i = 0; i < n; ++i){
                tmp -= taux(i, k)*taux(i, l)*X(i, p)*DiklCNORM_cpp(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw);
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
mat covBetaDeltaCNORM_cpp(int n,
                        int ng, 
                        int nx,
                        NumericVector pi,
                        NumericVector beta,
                        NumericVector sigma,
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
      mtmp = join_rows(mtmp, covBetakDeltalCNORM_cpp(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov Delta
mat covDeltaCNORM_cpp(int n,
                          int ng, 
                          int nx,
                          NumericVector pi,
                          NumericVector beta,
                          NumericVector sigma,
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
      mtmp = join_rows(mtmp, covDeltakDeltalCNORM_cpp(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov Delta Sigma
mat covDeltaSigmaCNORM_cpp(int n,
                      int ng, 
                      int nx,
                      NumericVector pi,
                      NumericVector beta,
                      NumericVector sigma,
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
      res = join_cols(res, covDeltaSigmakCNORM_cpp(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw));
  }
  return(res);
}
// matrix cov PI Sigma
mat covPiSigmaCNORM_cpp(int n,
                        int ng, 
                        int nx,
                        NumericVector pi,
                        NumericVector beta,
                        NumericVector sigma,
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
    mat covPiSigma(ng-1, ng);
    for (int k = 0; k < (ng-1); ++k){
      for (int l = 0; l < ng; ++l){
        double tmp = 0;
        if (l != (ng-1)){
          if (k == l){
            for (int i = 0; i < n; ++i){
              tmp += taux(i, k)*SikCNORM_cpp(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*((1-taux(i, k))/pi[k]+taux(i, ng-1)/pi[ng-1]);
            }
          }else{
            for (int i = 0; i < n; ++i){
              tmp += taux(i, l)*SikCNORM_cpp(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)/pi[l]*(-taux(i, k)/pi[k]+taux(i, ng-1)/pi[ng-1]);
            }
          }
        }else{
          for (int i = 0; i < n; ++i){
            tmp += taux(i, ng-1)*SikCNORM_cpp(i, ng-1, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, k)/pi[k]);
          }
        }
        covPiSigma(k, l) = tmp;  
      }
    }
    res= covPiSigma;
  }else{
    mat covPiSigma((ng-1)*nx, ng);
    for (int k = 0; k < (ng-1); ++k){
      for (int l = 0; l < ng; ++l){
        for (int p = 0; p < nx; ++p){
          double tmp = 0;
          if (k == l){
            for (int i = 0; i < n; ++i){
              tmp += taux(i, k)*(1-taux(i, k))*X(i, p)*SikCNORM_cpp(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw);
            }
          }else{
            for (int i = 0; i < n; ++i){
              tmp -= taux(i, k)*taux(i, l)*X(i, p)*SikCNORM_cpp(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw);
            }
          }
          covPiSigma(k*nx + p, l) = tmp;            
        } 
      }
    }
    res= covPiSigma;
  }
  return(res);
}
// matrix cov Beta
mat covBetaCNORM_cpp(int n,
                      int ng, 
                      int nx,
                      NumericVector pi,
                      NumericVector beta,
                      NumericVector sigma,
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
      mtmp = join_rows(mtmp, covBetakBetalCNORM_cpp(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov Beta Sigma
mat covBetaSigmaCNORM_cpp(int n,
                     int ng, 
                     int nx,
                     NumericVector pi,
                     NumericVector beta,
                     NumericVector sigma,
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
    res = join_cols(res, covBetaSigmakCNORM_cpp(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw));
  }
  return(res);
}
// matrix cov Sigma
mat covSigmaCNORM_cpp(int n,
                      int ng, 
                      int nx,
                      NumericVector pi,
                      NumericVector beta,
                      NumericVector sigma,
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
  mat res(ng, ng);
  for (int k = 0; k < ng; ++k){
    for (int l = 0; l < ng; ++l){
      double tmp = 0;
      if (k == l){
        for (int i = 0; i < n; ++i){
          tmp += pow(SikCNORM_cpp(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw), 2)*taux(i, k)*(1-taux(i, k));
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp -= SikCNORM_cpp(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*SikCNORM_cpp(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux(i, k)*taux(i, l);
        }
      }
      res(k, l) = tmp;
    }
  }
  return(res);
}
// Main function
// [[Rcpp::export]]
arma::vec IEM_cpp(NumericVector param,
                  int ng,
                  int nx,
                  IntegerVector nbeta,
                  int n,
                  NumericMatrix A,
                  NumericMatrix Y,
                  NumericMatrix X,
                  double ymin,
                  double ymax,
                  Nullable<NumericMatrix> TCOV,
                  int nw,
                  int refgr){
  int period = A.ncol();
  
  // Compute Matirx information
  NumericVector pi(ng);
  NumericVector beta;
  NumericVector sigma;
  NumericVector delta;
  IntegerVector ndeltacum(ng);
  IntegerVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-1)];
    beta = param[Range(ng,ng+sum(nbeta)-1)];
    sigma = param[Range(ng+sum(nbeta), ng+sum(nbeta)+ng-1)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng+sum(nbeta)+ng, param.length() - 1)];
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
    }
  }else{
    pi = param[Range(0,ng*nx-1)];
    beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
    sigma = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+ng-1)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng*nx+sum(nbeta)+ng, param.length() - 1)];
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
    }
  }
  NumericMatrix  taux = ftauxCNORM_cpp(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X);
  
  
  
  // matrix -B
  mat B(sum(nbeta)+ng*(nw+1), sum(nbeta)+ng*(nw+1)); 
  B.fill(0);
  // we take fof the dimension of mPi
  if (nw !=0){
    for (int i = 0; i < n; ++i){
      for (int t = 0; t < period; ++t){
        mat tmp1 = mbetasigmaCNORM_cpp(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw);
        mat tmp2 = mbetadeltaCNORM_cpp(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw);
        mat tmp3 = mdeltasigmaCNORM_cpp(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw);
        B += join_cols(join_rows(mbetaCNORM_cpp(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw), tmp2, tmp1),
                       join_rows(trans(tmp2), mdeltaCNORM_cpp(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw), tmp3),
                       join_rows(trans(tmp1), trans(tmp3), msigmaCNORM_cpp(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw))
        );
      }
    }
  }else{
    for (int i = 0; i < n; ++i){
      for (int t = 0; t < period; ++t){
        mat tmp1 = mbetasigmaCNORM_cpp(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw);
        B += join_cols(join_rows(mbetaCNORM_cpp(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw), tmp1),
                       join_rows(trans(tmp1), msigmaCNORM_cpp(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)));
      }
    }
  }
  B = join_rows(zeros(sum(nbeta)+ng+nw*ng, (ng-1)*nx), B);
  B = join_cols(join_rows(mPiCNORM_cpp(n, ng, nx, pi, taux , X, refgr), zeros((ng-1)*nx, sum(nbeta)+ng+nw*ng)), B);
  
  // Matrix of covariance of score function
  mat cov ;
  if (nw != 0){
    mat tmp1 = covPiBetaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp2 = covPiDeltaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp3 = covPiSigmaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp4 = covBetaDeltaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp5 = covBetaSigmaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp6 = covDeltaSigmaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    
    cov = join_rows(covPiCNORM_cpp(n, ng, nx, pi, X, taux), tmp1, tmp2, tmp3);
    cov = join_cols(cov, join_rows(trans(tmp1), covBetaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw), tmp4, tmp5));
    cov = join_cols(cov, join_rows(trans(tmp2), trans(tmp4), covDeltaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw), tmp6));
    cov = join_cols(cov, join_rows(trans(tmp3), trans(tmp5), trans(tmp6), covSigmaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw)));
  }else{
    mat tmp1 = covPiBetaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp3 = covPiSigmaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    mat tmp5 = covBetaSigmaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw);
    
    cov = join_rows(covPiCNORM_cpp(n, ng, nx, pi, X, taux), tmp1, tmp3);
    cov = join_cols(cov, join_rows(trans(tmp1), covBetaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw), tmp5));
    cov = join_cols(cov, join_rows(trans(tmp3), trans(tmp5), covSigmaCNORM_cpp(n, ng, nx, pi, beta, sigma, delta, X, taux, nbeta, A, Y, period, nbetacum, TCOV, ndeltacum, nw)));
    
  }
  // Information matrix of Fisher  
  mat IEM = -B - cov;
  
  return(sqrt(diagvec(inv(IEM))));
}
