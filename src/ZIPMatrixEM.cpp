#include "CommonFunction.h"
#include "ZIP.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// matrix Pi
mat mPiZIP_cpp(int n,
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
// -------------------------------------------------------------------
// Functions  lamda and nu
// -------------------------------------------------------------------
double lambdaikt_cpp(int k, 
                     int i, 
                     int t, 
                     IntegerVector nbeta, 
                     IntegerVector nbetacum, 
                     NumericMatrix A, 
                     NumericVector beta, 
                     Nullable<NumericMatrix> TCOV, 
                     int period, 
                     Nullable<NumericVector> delta, 
                     Nullable<IntegerVector> ndeltacum, 
                     int nw){
  NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
  NumericVector deltak;
  if (nw != 0){
    NumericVector vtmp(delta.get());
    NumericVector ntmp(ndeltacum.get());
    deltak = vtmp[Range(ntmp[k], ntmp[k+1]-1)]; 
  }  NumericVector vtmp2;
  double betaikt =0;
  for (int po = 0; po < nbeta[k]; ++po){
    betaikt += betak[po]*pow(A(i,t), po);
  }
  return(exp(betaikt + WitEM_cpp(TCOV, period, deltak, nw, i, t, k)));
}

double rhoikt_cpp(int k, 
                  int i,  
                  int t, 
                  IntegerVector nnu, 
                  IntegerVector nnucum,
                  NumericMatrix A, 
                  NumericVector nu){
  NumericVector nuk = nu[Range(nnucum[k], nnucum[k+1]-1)];
  double nuikt = 0;
  for (int po = 0; po < nnu[k]; ++po){
    nuikt += nuk[po]*pow(A(i,t), po);
  }
  return(exp(nuikt)/(1+exp(nuikt)));
}

// matrix beta

mat mbetaZIP_cpp(int i,
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
                 int nw,
                 IntegerVector nnucum,
                 IntegerVector nnu,
                 NumericVector nu, 
                 NumericVector pi,
                 int n,
                 NumericMatrix Y){
  NumericMatrix res(sum(nbeta), sum(nbeta));
  for (int k = 0; k < ng; ++k){
    double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
    for (int l =  nbetacum[k]; l <  nbetacum[k+1]; ++l){
      for (int lp =  nbetacum[k]; lp <  nbetacum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*(1-sikt)*pow(A(i, t), l- nbetacum[k])*pow(A(i, t), lp-nbetacum[k])*lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw);
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix delta

mat mdeltaZIP_cpp(int i,
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
                 int nw,
                 IntegerVector nnucum,
                 IntegerVector nnu,
                 NumericVector nu, 
                 NumericVector pi,
                 int n,
                 NumericMatrix Y){
  NumericMatrix res(nw*ng, nw*ng);
  for (int k = 0; k < ng; ++k){
    double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
    for (int l =  ndeltacum[k]; l <  ndeltacum[k+1]; ++l){
      for (int lp =  ndeltacum[k]; lp <  ndeltacum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*(1-sikt)*TCOV(i, t + (l - ndeltacum[k])*period)*TCOV(i, t + (lp - ndeltacum[k])*period)*lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw);
      }
    }
  }
  return(as<arma::mat>(res));
}

// matrix beta delta

mat mbetadeltaZIP_cpp(int i,
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
                      int nw,
                      IntegerVector nnucum,
                      IntegerVector nnu,
                      NumericVector nu, 
                      NumericVector pi,
                      int n,
                      NumericMatrix Y){
  NumericMatrix res(sum(nbeta), ng*nw);
  for (int k = 0; k < ng; ++k){
    double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
    for (int l =  nbetacum[k]; l <  nbetacum[k+1]; ++l){
      for (int lp =  ndeltacum[k]; lp <  ndeltacum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*(1-sikt)*TCOV(i, t + (lp - ndeltacum[k])*period)*pow(A(i, t), l-nbetacum[k])*lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw);
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix nu

mat mnuZIP_cpp(int i,
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
               int nw,
               IntegerVector nnucum,
               IntegerVector nnu,
               NumericVector nu, 
               NumericVector pi,
               int n,
               NumericMatrix Y){
  NumericMatrix res(sum(nnu), sum(nnu));
  for (int k = 0; k < ng; ++k){
    for (int l =  nnucum[k]; l <  nnucum[k+1]; ++l){
      for (int lp =  nnucum[k]; lp <  nnucum[k+1]; ++lp){
        res(l , lp) = -taux(i, k)*pow(A(i, t), l- nnucum[k])*pow(A(i, t), lp-nnucum[k])*rhoikt_cpp(k, i, t, nnu, nnucum, A, nu)*(1-rhoikt_cpp(k, i, t, nnu, nnucum, A, nu));
      }
    }
  }
  return(as<arma::mat>(res));
}

// function BPiikl

double BPiiklZIP_cpp(int i,
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
                     int nw,
                    IntegerVector nnucum,
                    IntegerVector nnu,
                    NumericVector nu, 
                    NumericVector pi,
                    int n){
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
    res += pow(A(i, t), l)*(Y(i,t)-lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))*(fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)-1);
  }
  return(res);
}
// function DPiiklZIP

double DPiiklZIP_cpp(int i,
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
                     int nw,
                     IntegerVector nnucum,
                     IntegerVector nnu,
                     NumericVector nu, 
                     NumericVector pi,
                     int n){
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
    res +=  TCOV(i, t + l*period)*(Y(i,t)-lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))*(fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)-1);
  }
  return(res);
}
// function NiklZIP 
 
double NiklZIP_cpp(int i,
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
                   int nw,
                   IntegerVector nnucum,
                   IntegerVector nnu,
                   NumericVector nu, 
                   NumericVector pi,
                   int n){
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
    res += pow(A(i, t), l)*(fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)-rhoikt_cpp(k, i, t, nnu, nnucum, A, nu));
  }
  return(res);
}



// matrix cov Pi Betak
 
mat covPiBetakZIP_cpp(int k,
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
                      NumericVector pi,
                      IntegerVector nnucum,
                      NumericVector nu,
                      IntegerVector nnu){
  NumericMatrix res(ng-1, nbeta[k]);
  for (int kp = 0; kp < (ng-1); ++kp){
    for (int l = 0; l < nbeta[k]; ++l){
      double tmp = 0;
      if (kp == k){
        for (int i = 0; i < n; ++i){
          tmp += BPiiklZIP_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw, nnucum, nnu, nu, pi, n)*taux(i,kp)*((1-taux(i,kp))/pi[kp]+taux(i,ng-1)/pi[ng-1]);
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp += BPiiklZIP_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw, nnucum, nnu, nu, pi, n)*taux(i,k)*(-taux(i,kp)/pi[kp]+taux(i,ng-1)/pi[ng-1]);
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Pi Betak
 
mat covPiDeltakZIP_cpp(int k,
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
                      NumericVector pi,
                      IntegerVector nnucum,
                      NumericVector nu,
                      IntegerVector nnu){
  NumericMatrix res(ng-1, nw);
  for (int kp = 0; kp < (ng-1); ++kp){
    for (int l = 0; l < nw; ++l){
      double tmp = 0;
      if (kp == k){
        for (int i = 0; i < n; ++i){
          tmp += DPiiklZIP_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw, nnucum, nnu, nu, pi, n)*taux(i,kp)*((1-taux(i,kp))/pi[kp]+taux(i,ng-1)/pi[ng-1]);
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp += DPiiklZIP_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw, nnucum, nnu, nu, pi, n)*taux(i,k)*(-taux(i,kp)/pi[kp]+taux(i,ng-1)/pi[ng-1]);
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Pi Betak
 
mat covPiNukZIP_cpp(int k,
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
                    NumericVector pi,
                    IntegerVector nnucum,
                    NumericVector nu,
                    IntegerVector nnu){
  NumericMatrix res(ng-1, nnu[k]);
  for (int kp = 0; kp < (ng-1); ++kp){
    for (int l = 0; l < nnu[k]; ++l){
      double tmp = 0;
      if (kp == k){
        for (int i = 0; i < n; ++i){
          tmp += NiklZIP_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw, nnucum, nnu, nu, pi, n)*taux(i,kp)*((1-taux(i,kp))/pi[kp]+taux(i,ng-1)/pi[ng-1]);
        }
      }else{
        for (int i = 0; i < n; ++i){
          tmp += NiklZIP_cpp(i, k, l, nbeta, A, Y, period, beta, taux, nbetacum, TCOVinit, deltainit, ndeltacuminit, nw, nnucum, nnu, nu, pi, n)*taux(i,k)*(-taux(i,kp)/pi[kp]+taux(i,ng-1)/pi[ng-1]);
        }
      }
      res(kp, l) = tmp;
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Betak Betal
 
mat covBetakBetalZIP_cpp(int k,
                         int l,
                         int ng,
                         int n,
                         IntegerVector nbeta,
                         NumericMatrix A,
                         NumericMatrix Y,
                         int period,
                         NumericVector beta,
                         NumericMatrix taux,
                         IntegerVector nbetacum,
                         Nullable<NumericMatrix> TCOV,
                         Nullable<NumericVector> delta,
                         Nullable<IntegerVector> ndeltacum,
                         int nw,
                         NumericVector pi,
                         IntegerVector nnucum,
                         NumericVector nu,
                         IntegerVector nnu){
  NumericMatrix res(nbeta[k], nbeta[l]);
  if (k == l){
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nbeta[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pikpt = pow(A(i,t), p)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            double pikqt = pow(A(i,t), q)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pikpt*pikqt*taux(i,k)*((1-taux(i,k))*(1-2*sikt)-sikt*(1-taux(i, k)*sikt));
          }
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nbeta[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pikpt = pow(A(i,t), p)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            double pilqt = pow(A(i,t), q)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pikpt*pilqt*taux(i,k)*taux(i,l)*(sikt-1)*(1-silt);
          }
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Betak Betal
 
mat covBetakNulZIP_cpp(int k,
                       int l,
                       int ng,
                       int n,
                       IntegerVector nbeta,
                       NumericMatrix A,
                       NumericMatrix Y,
                       int period,
                       NumericVector beta,
                       NumericMatrix taux,
                       IntegerVector nbetacum,
                       Nullable<NumericMatrix> TCOV,
                       Nullable<NumericVector> delta,
                       Nullable<IntegerVector> ndeltacum,
                       int nw,
                       NumericVector pi,
                       IntegerVector nnucum,
                       NumericVector nu,
                       IntegerVector nnu){
  NumericMatrix res(nbeta[k], nnu[l]);
  if (k == l){
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nnu[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pilpt = pow(A(i,t), p)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pow(A(i,t),q)*pilpt*taux(i,k)*(sikt-1)*(taux(i,k)*sikt+rhoikt_cpp(k, i, t, nnu, nnucum, A, nu)*(1-taux(i, k)));
          }
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nnu[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pilpt = pow(A(i,t), p)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pow(A(i,t), q)*pilpt*taux(i,k)*taux(i,l)*(silt-1)*(sikt-rhoikt_cpp(k, i, t, nnu, nnucum, A, nu));
          }
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Betak Deltal
 
mat covBetakDeltalZIP_cpp(int k,
                         int l,
                         int ng,
                         int n,
                         IntegerVector nbeta,
                         NumericMatrix A,
                         NumericMatrix Y,
                         int period,
                         NumericVector beta,
                         NumericMatrix taux,
                         IntegerVector nbetacum,
                         NumericMatrix TCOV,
                         NumericVector delta,
                         IntegerVector ndeltacum,
                         int nw,
                         NumericVector pi,
                         IntegerVector nnucum,
                         NumericVector nu,
                         IntegerVector nnu){
  NumericMatrix res(nbeta[k], nw);
  if (k == l){
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pikpt = pow(A(i,t), p)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            double pikqt =  TCOV(i, t + q*period)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pikpt*pikqt*taux(i,k)*((1-taux(i,k))*(1-2*sikt)-sikt*(1-taux(i, k)*sikt));
          }
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nbeta[k]; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pikpt = pow(A(i,t), p)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            double pilqt =  TCOV(i, t + q*period)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pikpt*pilqt*taux(i,k)*taux(i,l)*(sikt-1)*(1-silt);
          }
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Nuk Nul
 
mat covNukNulZIP_cpp(int k,
                       int l,
                       int ng,
                       int n,
                       IntegerVector nbeta,
                       NumericMatrix A,
                       NumericMatrix Y,
                       int period,
                       NumericVector beta,
                       NumericMatrix taux,
                       IntegerVector nbetacum,
                       Nullable<NumericMatrix> TCOV,
                       Nullable<NumericVector> delta,
                       Nullable<IntegerVector> ndeltacum,
                       int nw,
                       NumericVector pi,
                       IntegerVector nnucum,
                       NumericVector nu,
                       IntegerVector nnu){
  NumericMatrix res(nnu[k], nnu[l]);
  if (k == l){
    for (int p = 0; p < nnu[k]; ++p){
      for (int q = 0; q < nnu[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double roikt = rhoikt_cpp(k, i, t, nnu, nnucum, A, nu);
            tmp += pow(A(i,t), p)*pow(A(i,t), q)*taux(i,k)*(sikt*(1-taux(i,k)*sikt)+roikt*roikt*(1-taux(i,k))+2*roikt*sikt*(1-taux(i,k)));
          }
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nnu[k]; ++p){
      for (int q = 0; q < nnu[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            tmp += pow(A(i,t), p)*pow(A(i,t), q)*taux(i,k)*taux(i,l)*(-sikt*silt+rhoikt_cpp(l, i, t, nnu, nnucum, A, nu)*sikt+rhoikt_cpp(k, i, t, nnu, nnucum, A, nu)*silt-rhoikt_cpp(k, i, t, nnu, nnucum, A, nu)*rhoikt_cpp(l, i, t, nnu, nnucum, A, nu));
          }
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Deltak Nul
 
mat covDeltakNulZIP_cpp(int k,
                       int l,
                       int ng,
                       int n,
                       IntegerVector nbeta,
                       NumericMatrix A,
                       NumericMatrix Y,
                       int period,
                       NumericVector beta,
                       NumericMatrix taux,
                       IntegerVector nbetacum,
                       NumericMatrix TCOV,
                       NumericVector delta,
                       IntegerVector ndeltacum,
                       int nw,
                       NumericVector pi,
                       IntegerVector nnucum,
                       NumericVector nu,
                       IntegerVector nnu){
  NumericMatrix res(nw, nnu[l]);
  if (k == l){
    for (int p = 0; p < nw; ++p){
      for (int q = 0; q < nnu[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pilpt =  TCOV(i, t + p*period)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pow(A(i,t),q)*pilpt*taux(i,k)*(sikt-1)*(taux(i,k)*sikt+rhoikt_cpp(k, i, t, nnu, nnucum, A, nu)*(1-taux(i, k)));
          }
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nw; ++p){
      for (int q = 0; q < nnu[l]; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pilpt = TCOV(i, t + p*period)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pow(A(i,t), q)*pilpt*taux(i,k)*taux(i,l)*(silt-1)*(sikt-rhoikt_cpp(k, i, t, nnu, nnucum, A, nu));
          }
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Deltak Deltal
 
mat covDeltakDeltalZIP_cpp(int k,
                           int l,
                           int ng,
                           int n,
                           IntegerVector nbeta,
                           NumericMatrix A,
                           NumericMatrix Y,
                           int period,
                           NumericVector beta,
                           NumericMatrix taux,
                           IntegerVector nbetacum,
                           NumericMatrix TCOV,
                           NumericVector delta,
                           IntegerVector ndeltacum,
                           int nw,
                           NumericVector pi,
                           IntegerVector nnucum,
                           NumericVector nu,
                           IntegerVector nnu){
  NumericMatrix res(nw, nw);
  if (k == l){
    for (int p = 0; p < nw; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pikpt =TCOV(i, t + p*period)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            double pikqt = TCOV(i, t + q*period)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pikpt*pikqt*taux(i,k)*((1-taux(i,k))*(1-2*sikt)-sikt*(1-taux(i, k)*sikt));
          }
        }
        res(p, q) = tmp;
      }
    }
  }else{
    for (int p = 0; p < nw; ++p){
      for (int q = 0; q < nw; ++q){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
            double pikpt = TCOV(i, t + p*period)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            double pilqt = TCOV(i, t + q*period)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
            tmp += pikpt*pilqt*taux(i,k)*taux(i,l)*(sikt-1)*(1-silt);
          }
        }
        res(p, q) = tmp;
      }
    }
  }
  return(as<arma::mat>(res));
}
// matrix cov Pi
 
mat covPiZIP_cpp(int n,
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
 
mat covPiBetaZIP_cpp(int ng,
                     int n,
                     int nx,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     NumericMatrix X,
                     int period,
                     NumericVector beta,
                     NumericMatrix taux,
                     IntegerVector nbetacum,
                     Nullable<NumericMatrix> TCOV,
                     Nullable<NumericVector> delta,
                     Nullable<IntegerVector> ndeltacum,
                     int nw,
                     NumericVector pi,
                     IntegerVector nnucum,
                     NumericVector nu,
                     IntegerVector nnu){
  mat res;
  if (nx == 1){
    mat covPiBeta;
    for (int k = 0; k < ng-1; ++k){
      covPiBeta = join_rows(covPiBeta, covPiBetakZIP_cpp(k, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    mat mtmp(ng-1, nbeta[ng-1]);
    for (int kp = 0; kp < ng-1; ++kp){
      for (int l = 0; l < nbeta[ng-1]; ++l){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += BPiiklZIP_cpp(i, ng-1, l, nbeta, A, Y, period, beta,  taux, nbetacum, TCOV, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n)*taux(i, ng-1)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, kp)/pi[kp]);
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
                for (int t = 0 ; t < period; ++t){
                  double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
                  double pikqt = pow(A(i,t), q)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
                  tmp += X(i,p)*pikqt*taux(i,k)*taux(i,k)*(1-taux(i,k))*(1-sikt);
                }
              }
            }else{
              for (int i = 0; i < n; ++i){
                for (int t = 0 ; t < period; ++t){
                  double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
                  double pilqt = pow(A(i,t), q)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
                  tmp += X(i,p)*pilqt*taux(i,k)*taux(i,l)*(1+silt);
                }
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
// matrix cov Pi Beta
 
mat covPiNuZIP_cpp(int ng,
                     int n,
                     int nx,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     NumericMatrix X,
                     int period,
                     NumericVector beta,
                     NumericMatrix taux,
                     IntegerVector nbetacum,
                     Nullable<NumericMatrix> TCOV,
                     Nullable<NumericVector> delta,
                     Nullable<IntegerVector> ndeltacum,
                     int nw,
                     NumericVector pi,
                     IntegerVector nnucum,
                     NumericVector nu,
                     IntegerVector nnu){
  mat res;
  if (nx == 1){
    mat covPiNu;
    for (int k = 0; k < ng-1; ++k){
      covPiNu = join_rows(covPiNu, covPiNukZIP_cpp(k, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    mat mtmp(ng-1, nnu[ng-1]);
    for (int kp = 0; kp < ng-1; ++kp){
      for (int l = 0; l < nnu[ng-1]; ++l){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp -= NiklZIP_cpp(i, ng-1, l, nbeta, A, Y, period, beta,  taux, nbetacum, TCOV, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n)*taux(i, ng-1)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, kp)/pi[kp]);
        }
        mtmp(kp, l) = tmp;
      }
    }
    res = join_rows(covPiNu, mtmp);
  }else{
    mat covPiNu((ng-1)*nx, sum(nnu));
    for (int k = 0; k < ng-1; ++k){
      for (int l = 0; l < ng; ++l){
        for (int p = 0; p < nx; ++p){
          for (int q = 0; q < nnu[l]; ++q){
            double tmp = 0;
            if (k == l){
              for (int i = 0; i < n; ++i){
                for (int t = 0 ; t < period; ++t){
                  double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
                  tmp += X(i,p)*pow(A(i,t), q)*taux(i,k)*(1-taux(i,k))*(sikt-rhoikt_cpp(k, i, t, nnu, nnucum, A, nu));
                }
              }
            }else{
              for (int i = 0; i < n; ++i){
                for (int t = 0 ; t < period; ++t){
                  double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
                  tmp += X(i,p)*pow(A(i,t),q)*taux(i,k)*taux(i,l)*(rhoikt_cpp(l, i, t, nnu, nnucum, A, nu)-silt);
                }
              }
            }
            covPiNu(k*nx + p, nnucum[l] + q) = tmp;
          }
        }
      }
    }
    res = covPiNu;
  }
  return(res);
}
// matrix cov Pi Beta
 
mat covPiDeltaZIP_cpp(int ng,
                     int n,
                     int nx,
                     IntegerVector nbeta,
                     NumericMatrix A,
                     NumericMatrix Y,
                     NumericMatrix X,
                     int period,
                     NumericVector beta,
                     NumericMatrix taux,
                     IntegerVector nbetacum,
                     NumericMatrix TCOV,
                     NumericVector delta,
                     IntegerVector ndeltacum,
                     int nw,
                     NumericVector pi,
                     IntegerVector nnucum,
                     NumericVector nu,
                     IntegerVector nnu){
  mat res;
  if (nx == 1){
    mat covPiDelta;
    for (int k = 0; k < ng-1; ++k){
      covPiDelta = join_rows(covPiDelta, covPiDeltakZIP_cpp(k, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    mat mtmp(ng-1, nw);
    for (int kp = 0; kp < ng-1; ++kp){
      for (int l = 0; l < nw; ++l){
        double tmp = 0;
        for (int i = 0; i < n; ++i){
          tmp += DPiiklZIP_cpp(i, ng-1, l, nbeta, A, Y, period, beta,  taux, nbetacum, TCOV, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n)*taux(i, ng-1)*((1-taux(i, ng-1))/pi[ng-1]+taux(i, kp)/pi[kp]);
        }
        mtmp(kp, l) = tmp;
      }
    }
    res = join_rows(covPiDelta, mtmp);
  }else{
    mat covPiDelta((ng-1)*nx, nw*ng);
    for (int k = 0; k < ng-1; ++k){
      for (int l = 0; l < ng; ++l){
        for (int p = 0; p < nx; ++p){
          for (int q = 0; q < nw; ++q){
            double tmp = 0;
            if (k == l){
              for (int i = 0; i < n; ++i){
                for (int t = 0 ; t < period; ++t){
                  double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
                  double pikqt =  TCOV(i, t + q*period)*(Y(i, t) - lambdaikt_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
                  tmp += X(i,p)*pikqt*taux(i,k)*taux(i,k)*(1-taux(i,k))*(1-sikt);
                }
              }
            }else{
              for (int i = 0; i < n; ++i){
                for (int t = 0 ; t < period; ++t){
                  double silt = fSikt_cpp(pi, beta, nu, l, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
                  double pilqt = TCOV(i, t + q*period)*(Y(i, t) - lambdaikt_cpp(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw));
                  tmp += X(i,p)*pilqt*taux(i,k)*taux(i,l)*(1+silt);
                }
              }
            }
            covPiDelta(k*nx + p, nw + q) = tmp;
          }
        }
      }
    }
    res = covPiDelta;
  }
  return(res);
}
// matrix cov Beta
 
mat covBetaZIP_cpp(int ng,
                   int n,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   int period,
                   NumericVector beta,
                   NumericMatrix taux,
                   IntegerVector nbetacum,
                   Nullable<NumericMatrix> TCOV,
                   Nullable<NumericVector> delta,
                   Nullable<IntegerVector> ndeltacum,
                   int nw,
                   NumericVector pi,
                   IntegerVector nnucum,
                   NumericVector nu,
                   IntegerVector nnu){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covBetakBetalZIP_cpp(k, l, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov beta Nu
 
mat covBetaNuZIP_cpp(int ng,
                   int n,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   int period,
                   NumericVector beta,
                   NumericMatrix taux,
                   IntegerVector nbetacum,
                   Nullable<NumericMatrix> TCOV,
                   Nullable<NumericVector> delta,
                   Nullable<IntegerVector> ndeltacum,
                   int nw,
                   NumericVector pi,
                   IntegerVector nnucum,
                   NumericVector nu,
                   IntegerVector nnu){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covBetakNulZIP_cpp(k, l, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov Beta Delta
 
mat covBetaDeltaZIP_cpp(int ng,
                   int n,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   int period,
                   NumericVector beta,
                   NumericMatrix taux,
                   IntegerVector nbetacum,
                   NumericMatrix TCOV,
                   NumericVector delta,
                   IntegerVector ndeltacum,
                   int nw,
                   NumericVector pi,
                   IntegerVector nnucum,
                   NumericVector nu,
                   IntegerVector nnu){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covBetakDeltalZIP_cpp(k, l, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov Nu
mat covNuZIP_cpp(int ng,
                   int n,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   int period,
                   NumericVector beta,
                   NumericMatrix taux,
                   IntegerVector nbetacum,
                   Nullable<NumericMatrix> TCOV,
                   Nullable<NumericVector> delta,
                   Nullable<IntegerVector> ndeltacum,
                   int nw,
                   NumericVector pi,
                   IntegerVector nnucum,
                   NumericVector nu,
                   IntegerVector nnu){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covNukNulZIP_cpp(k, l, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov Delta
 
mat covDeltaZIP_cpp(int ng,
                   int n,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   int period,
                   NumericVector beta,
                   NumericMatrix taux,
                   IntegerVector nbetacum,
                   NumericMatrix TCOV,
                   NumericVector delta,
                   IntegerVector ndeltacum,
                   int nw,
                   NumericVector pi,
                   IntegerVector nnucum,
                   NumericVector nu,
                   IntegerVector nnu){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covDeltakDeltalZIP_cpp(k, l, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}
// matrix cov Delta Nu
 
mat covDeltaNuZIP_cpp(int ng,
                    int n,
                    IntegerVector nbeta,
                    NumericMatrix A,
                    NumericMatrix Y,
                    int period,
                    NumericVector beta,
                    NumericMatrix taux,
                    IntegerVector nbetacum,
                    NumericMatrix TCOV,
                    NumericVector delta,
                    IntegerVector ndeltacum,
                    int nw,
                    NumericVector pi,
                    IntegerVector nnucum,
                    NumericVector nu,
                    IntegerVector nnu){
  mat res;
  for (int k = 0; k < ng; ++k){
    mat mtmp;
    for (int l = 0; l < ng; ++l){
      mtmp = join_rows(mtmp, covDeltakNulZIP_cpp(k, l, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu));
    }
    res= join_cols(res, mtmp);
  }
  return(res);
}

// Main function
// [[Rcpp::export]]
arma::vec IEMZIP_cpp(NumericVector param,
                  int ng,
                  int nx,
                  IntegerVector nbeta,
                  IntegerVector nnu,
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
  NumericVector nu;
  NumericVector delta;
  
  IntegerVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  IntegerVector nnucum(nnu.size());
  std::partial_sum(nnu.begin(), nnu.end(), nnucum.begin());
  nnucum.push_front(0);
  IntegerVector ndeltacum;
  
  if (nx == 1){
    pi = param[Range(0,ng-1)];
    beta = param[Range(ng,ng+sum(nbeta)-1)];
    nu = param[Range(ng+sum(nbeta), ng+sum(nbeta)+sum(nnu)-1)];
    if (param.length() > ng*nx+sum(nbeta)+sum(nnu)){
      delta = param[Range(ng+sum(nbeta)+sum(nnu), param.length()-1)];
      NumericVector deltatmp(ng);
      IntegerVector ndeltacumtmp(nw*ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacumtmp.begin());
      ndeltacumtmp.push_front(0);
      ndeltacum = ndeltacumtmp;
    }
  }else{
    pi = param[Range(0,ng*nx-1)];
    beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
    nu = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+sum(nnu)-1)];
    if (param.length() > ng*nx+sum(nbeta)+sum(nnu)){
      delta = param[Range(ng*nx+sum(nbeta)+sum(nnu), param.length()-1)];
      NumericVector deltatmp(ng);
      IntegerVector ndeltacumtmp(nw*ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacumtmp.begin());
      ndeltacumtmp.push_front(0);
      ndeltacum = ndeltacumtmp;
    }
  }
  NumericMatrix  taux = ftauxZIP_cpp(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X);
  
  // matrix -B
  mat B(sum(nbeta)+sum(nnu)+ng*nw, sum(nbeta)+sum(nnu)+ng*nw);
  B.fill(0);
  // we take fof the dimension of mPi
  if (nw != 0){
    NumericMatrix TCOVv(TCOV.get());
    mat tmp0 = zeros(sum(nbeta), sum(nnu));
    mat tmp2 = zeros(sum(nw*ng), sum(nnu)); 
    for (int i = 0; i < n; ++i){
      for (int t = 0; t < period; ++t){
        mat tmp1 = mbetadeltaZIP_cpp(i, t, ng, nbeta, A, beta,  taux, nbetacum, TCOVv, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y);
        B += join_cols(join_rows(mbetaZIP_cpp(i, t, ng, nbeta, A, beta,  taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y), tmp0, tmp1),
                       join_rows(trans(tmp0), mnuZIP_cpp(i, t, ng, nbeta, A, beta,  taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y), trans(tmp2)),
                       join_rows(trans(tmp1), tmp2, mdeltaZIP_cpp(i, t, ng, nbeta, A, beta,  taux, nbetacum, TCOVv, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y))
        );
      }
    }
  }else{
    mat tmp0 = zeros(sum(nbeta), sum(nnu));
    for (int i = 0; i < n; ++i){
      for (int t = 0; t < period; ++t){
         B += join_cols(join_rows(mbetaZIP_cpp(i, t, ng, nbeta, A, beta,  taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y), tmp0),
                       join_rows(trans(tmp0), mnuZIP_cpp(i, t, ng, nbeta, A, beta,  taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y))
        );
      }
    }
  }
  B = join_rows(zeros(sum(nbeta)+sum(nnu)+nw*ng, (ng-1)*nx), B);
  B = join_cols(join_rows(mPiZIP_cpp(n, ng, nx, pi, taux , X, refgr), zeros((ng-1)*nx, sum(nbeta)+sum(nnu)+nw*ng)), B);

  // Matrix of covariance of score function
  mat cov ;
  if (nw != 0){
    NumericMatrix TCOVv(TCOV.get());
    mat tmp1 = covPiBetaZIP_cpp(ng, n, nx, nbeta, A, Y, X, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    mat tmp2 = covPiNuZIP_cpp(ng, n, nx, nbeta, A, Y, X, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    mat tmp3 = covPiDeltaZIP_cpp(ng, n, nx, nbeta, A, Y, X, period, beta, taux, nbetacum, TCOVv, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    mat tmp4 = covBetaNuZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    mat tmp5 = covBetaDeltaZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOVv, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    mat tmp6 = covDeltaNuZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOVv, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    
    cov = join_rows(covPiZIP_cpp(n, ng, nx, pi, X, taux), tmp1, tmp2, tmp3);
    cov = join_cols(cov, join_rows(trans(tmp1), covBetaZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu), tmp4, tmp5));
    cov = join_cols(cov, join_rows(trans(tmp2), trans(tmp4), covNuZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu), trans(tmp6)));
    cov = join_cols(cov, join_rows(trans(tmp3), trans(tmp5), tmp6, covDeltaZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOVv, delta, ndeltacum, nw, pi, nnucum, nu, nnu)));
   }else{
    mat tmp1 = covPiBetaZIP_cpp(ng, n, nx, nbeta, A, Y, X, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    mat tmp2 = covPiNuZIP_cpp(ng, n, nx, nbeta, A, Y, X, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    mat tmp3 = covBetaNuZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu);
    
    cov = join_rows(covPiZIP_cpp(n, ng, nx, pi, X, taux), tmp1, tmp2);
    cov = join_cols(cov, join_rows(trans(tmp1), covBetaZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu), tmp3));
    cov = join_cols(cov, join_rows(trans(tmp2), trans(tmp3), covNuZIP_cpp(ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)));
  }
  // Information matrix of Fisher
  mat IEM = -B - cov;

  return(sqrt(diagvec(inv(IEM))));
}

