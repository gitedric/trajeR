#include "CommonFunction.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::interfaces(r, cpp)]]

// ----------------------------------------------------------------------------
// gk Logit
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double gkLOGIT_cpp(List beta,
                   int i,
                   int k,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   Nullable<NumericMatrix> TCOV,
                   Nullable<List> delta,
                   int nw){
  int period = A.ncol();
  NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
  double res = 1;
  for (int t = 0; t < period; ++t){
    if (R_IsNA(Y(i, t)) == FALSE){
      res *= pow(1-1/(1+exp(muikt[t])), Y(i,t))*pow(1/(1+exp(muikt[t])), 1-Y(i,t));
    }
  }
  return(res);
}
// ----------------------------------------------------------------------------
// dif likelihood theta
// ----------------------------------------------------------------------------
NumericVector difLthetakLOGIT_cpp(NumericVector theta,
                                  List beta,
                                  Nullable<List> delta,
                                  int k,
                                  int ng,
                                  int nx,
                                  int n,
                                  IntegerVector nbeta,
                                  NumericMatrix A,
                                  NumericMatrix Y,
                                  NumericMatrix X,
                                  Nullable<NumericMatrix> TCOV,
                                  int nw){
  NumericVector thetas;
  for (int l = 0; l < nx; ++l){
    double a = 0;
    for (int i = 0 ; i < n; ++i){
      NumericVector vtmp;
      for (int s = 0; s < ng; ++s){
        double tmp = 0;
        for (int ind = 0; ind < nx; ++ind){
          tmp += theta[s*nx+ind]*X(i,ind); 
        }
        vtmp.push_back(exp(tmp));
      }
      double tmp1 = 0;
      for (int s = 0; s < ng; ++s){
        tmp1 += vtmp[s]*(gkLOGIT_cpp(beta, i, k, nbeta, A, Y, TCOV, delta, nw)-gkLOGIT_cpp(beta, i, s, nbeta, A, Y, TCOV, delta, nw));
      }
      tmp1 = tmp1*vtmp[k];
      double tmp2 = 0;
      for (int s = 0; s < ng; ++s){
        tmp2 += vtmp[s]*gkLOGIT_cpp(beta, i, s, nbeta, A, Y, TCOV, delta, nw);
      }
      a += X(i, l)*tmp1/(sum(vtmp)*tmp2);
    }
    thetas.push_back(a);
  }
  
  return(thetas);
}
// ----------------------------------------------------------------------------
// dif likelihood betak
// ----------------------------------------------------------------------------
NumericVector difLbetakLOGIT_cpp(NumericVector theta,
                                 List beta,
                                 Nullable<List> delta,
                                 int k,
                                 int ng,
                                 int nx,
                                 int n,
                                 IntegerVector nbeta,
                                 NumericMatrix A,
                                 NumericMatrix Y,
                                 NumericMatrix X,
                                 Nullable<NumericMatrix> TCOV,
                                 int nw){
  NumericVector betas;
  int period = A.ncol();
  
  for (int l = 0; l < nbeta[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      // initialize mean
      NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      NumericVector tmp1;
      for (int t = 0; t < period; ++t){
        tmp1.push_back(pow(1-1/(1+exp(muikt[t])), Y(i,t))*pow(1/(1+exp(muikt[t])), 1-Y(i,t)));
      }
      double tmp2 = 0;
      for (int t = 0; t < period; ++t){
        if (R_IsNA(Y(i, t)) == FALSE){
          NumericVector tmp1minus = tmp1;
          tmp1minus.erase(tmp1minus.begin() + t);
          double ytmp = -1;
          if (Y(i, t) == 1){ ytmp = 1;}
          tmp2 += ytmp*pow(A(i, t), l)/(1+exp(muikt[t]))*(1-1/(1+exp(muikt[t])))*prodvect(tmp1minus);
        }
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkLOGIT_cpp(beta, i, s, nbeta, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*tmp2;
    }
    betas.push_back(a);
  }
  return(betas); 
}
// ----------------------------------------------------------------------------
// dif likelihood deltak
// ----------------------------------------------------------------------------
NumericVector difLdeltakLOGIT_cpp(NumericVector theta,
                                  List beta,
                                  Nullable<List> delta,
                                  int k,
                                  int ng,
                                  int nx,
                                  int n,
                                  IntegerVector nbeta,
                                  NumericMatrix A,
                                  NumericMatrix Y,
                                  NumericMatrix X,
                                  Nullable<NumericMatrix> initTCOV,
                                  int nw){
  NumericVector deltas;
  int period = A.ncol();
  NumericMatrix TCOV(initTCOV);
  
  for (int l = 0; l < nw; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      // initialize mean
      NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      NumericVector tmp1;
      for (int t = 0; t < period; ++t){
        tmp1.push_back(pow(1-1/(1+exp(muikt[t])), Y(i,t))*pow(1/(1+exp(muikt[t])), 1-Y(i,t)));
      }
      double tmp2 = 0;
      for (int t = 0; t < period; ++t){
        if (R_IsNA(Y(i, t)) == FALSE){
          NumericVector tmp1minus = tmp1;
          tmp1minus.erase(tmp1minus.begin() + t);
          tmp2 += TCOV(i, t + l*period)*exp(Y(i, t)*muikt[t])/pow((1+exp(muikt[t])), 2)*(Y(i, t) - (1-Y(i, t))*exp(muikt[t]))*prodvect(tmp1minus);
        }
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkLOGIT_cpp(beta, i, s, nbeta, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*tmp2;
    }
    deltas.push_back(a);
  }
  return(deltas); 
}
// ----------------------------------------------------------------------------
// dif likelihood 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector difLLOGIT_cpp(NumericVector param,
                                  int ng, 
                                  int nx,
                                  int n,
                                  IntegerVector nbeta,
                                  NumericMatrix A,
                                  NumericMatrix Y,
                                  NumericMatrix X,
                                  Nullable<NumericMatrix> TCOV,
                                  int nw){
  NumericVector out;
  //NumericVector theta = param[Range(0,ng*nx-1)];
  NumericVector theta = param[Range(0,(ng-1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
    }
  //NumericVector beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
  NumericVector beta = param[Range((ng-1)*nx,(ng-1)*nx+sum(nbeta)-1)];
  // create a list for beta
  List betaL(ng);
  int ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nbeta[i]; j++){
      tmp.push_back(beta[ind + j]);
    }
    ind += nbeta[i];
    betaL[i] = tmp;
  }
  NumericVector delta;
  List deltaL(ng);
  if (param.length() > (ng-1)*nx+sum(nbeta)){
    delta = param[Range((ng-1)*nx+sum(nbeta), param.length())];
    if (nw != 0){
      int ind = 0;
      for (int i = 0; i < ng; i++){
        NumericVector tmp;
        for (int j = 0; j < nw; j++){
          tmp.push_back(delta[ind + j]);
        }
        ind += nw;
        deltaL[i] = tmp;
      }  
    } 
  }
  
  //for (int k = 0; k < ng; ++k){
  for (int k = 1; k < ng; ++k){
    NumericVector tmp = difLthetakLOGIT_cpp(theta, betaL, deltaL, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLbetakLOGIT_cpp(theta, betaL, deltaL, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  if (nw != 0){
    for (int k = 0; k < ng; ++k){
      NumericVector tmp =  difLdeltakLOGIT_cpp(theta, betaL, deltaL, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw);
      for (int i = 0; i < tmp.length(); ++i){
        out.push_back(tmp[i]);
      }
    }
  }
  return(out);
}
// ----------------------------------------------------------------------------
// Likelihood 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double likelihoodLOGIT_cpp(NumericVector param,
                           int ng, 
                           int nx,
                           int n,
                           IntegerVector nbeta,
                           NumericMatrix A,
                           NumericMatrix Y,
                           NumericMatrix X,
                           Nullable<NumericMatrix> TCOV,
                           int nw){
  //NumericVector theta = param[Range(0,ng*nx-1)];
  NumericVector theta = param[Range(0,(ng-1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  //NumericVector beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
  NumericVector beta = param[Range((ng-1)*nx,(ng-1)*nx+sum(nbeta)-1)];
  // create a list for beta
  List betaL(ng);
  int ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nbeta[i]; j++){
      tmp.push_back(beta[ind + j]);
    }
    ind += nbeta[i];
    betaL[i] = tmp;
  }
  NumericVector delta;
  List deltaL(ng);
  if (param.length() > (ng-1)*nx+sum(nbeta)){
    delta = param[Range((ng-1)*nx+sum(nbeta), param.length())];
    if (nw != 0){
      int ind = 0;
      for (int i = 0; i < ng; i++){
        NumericVector tmp;
        for (int j = 0; j < nw; j++){
          tmp.push_back(delta[ind + j]);
        }
        ind += nw;
        deltaL[i] = tmp;
      }  
    } 
  }
  double out = 0;
  for (int i =0; i < n; ++i){
    double a = 0;
    for (int s = 0; s < ng; ++s){
      a +=  piikIntern_cpp(theta, i, s, ng, X)*gkLOGIT_cpp(betaL, i, s, nbeta, A, Y, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}
// ----------------------------------------------------------------------------
// Function rate
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
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
                             NumericMatrix X){
  // create a list for beta
  List betaL(ng);
  int ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nbeta[i]; j++){
      tmp.push_back(beta[ind + j]);
    }
    ind += nbeta[i];
    betaL[i] = tmp;
  }
  List deltaL(ng);
  if (delta.isNotNull()){
    // create a list for delta
    NumericVector deltatmp(delta.get());
    if (nw != 0){
      int ind = 0;
      for (int i = 0; i < ng; i++){
        NumericVector tmp1;
        for (int j = 0; j < nw; j++){
          tmp1.push_back(deltatmp[ind + j]);
        }
        ind += nw;
        deltaL[i] = tmp1;
      }
    }  
  }
  NumericMatrix mtmp(n,ng);
  if (nx == 1){
    for (int i = 0; i < n; ++i){
      double s = 0;
      for (int k = 0; k < ng; ++k){
        mtmp(i,k) = pi[k]*gkLOGIT_cpp(betaL, i, k, nbeta, A, Y, TCOV, deltaL, nw);
        s += mtmp(i,k);
      } 
      mtmp(i, _) = 1/(1+ (s - mtmp(i, _))/mtmp(i, _));
    }
  }else{
    for (int i = 0; i < n; ++i){
      double s = 0;
      for (int k = 0; k < ng; ++k){
        mtmp(i,k) = piikIntern_cpp(pi, i, k, ng, X)*gkLOGIT_cpp(betaL, i, k, nbeta, A, Y, TCOV, deltaL, nw);
        s += mtmp(i,k);
      } 
      mtmp(i, _) = 1/(1+ (s - mtmp(i, _))/mtmp(i, _));
    }
  }
  return(mtmp);
}
// ----------------------------------------------------------------------------
// Q function for betak
// ----------------------------------------------------------------------------
double QbetakLOGIT_cpp(NumericVector betak,
                       NumericMatrix taux,
                       int k,
                       int n,
                       int ng,
                       IntegerVector nbeta,
                       NumericMatrix A,
                       NumericMatrix Y,
                       Nullable<NumericMatrix> TCOV,
                       Nullable<NumericVector> deltainit,
                       int nw){
  int period = A.ncol();
  NumericVector delta;
  
  double a = 0;
  NumericVector deltak;
  if (nw!=0){
    NumericVector tmp3(deltainit.get());
    delta = tmp3;
    NumericVector ndeltacum(ng);
    NumericVector deltatmp(ng);
    deltatmp.fill(nw);
    std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
    ndeltacum.push_front(0);
    deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
  }
 
  for (int i = 0; i < n; ++i){
    double b = 0;
    double muikt = 0;
    for (int s = 0; s < period; ++s){
      NumericVector vtmp2;
      for (int po = 0; po < nbeta[k]; ++po){
        vtmp2.push_back(pow(A(i,s), po));
      }
      muikt = sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, s, k);
      b += Y(i, s)*muikt-log(1+exp(muikt));
    }
    a += taux(i, k)*b;
  }
  return(a);
}
// ----------------------------------------------------------------------------
// Q function for deltak
// ----------------------------------------------------------------------------
double QdeltakLOGIT_cpp(NumericVector deltak,
                       NumericMatrix taux,
                       int k,
                       int n,
                       int ng,
                       IntegerVector nbeta,
                       NumericMatrix A,
                       NumericMatrix Y,
                       Nullable<NumericMatrix> TCOV,
                       NumericVector beta,
                       int nw){
  int period = A.ncol();
  NumericVector delta;
  
  double a = 0;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
  
  
  for (int i = 0; i < n; ++i){
    double b = 0;
    double muikt = 0;
    for (int s = 0; s < period; ++s){
      NumericVector vtmp2;
      for (int po = 0; po < nbeta[k]; ++po){
        vtmp2.push_back(pow(A(i,s), po));
      }
      muikt = sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, s, k);
      b += Y(i, s)*muikt-log(1+exp(muikt));
    }
    a += taux(i, k)*b;
  }
  return(a);
}
// ----------------------------------------------------------------------------
// Differential of betak Q function
// ----------------------------------------------------------------------------
NumericVector difQbetakLOGIT_cpp(NumericVector betak,
                                 NumericMatrix taux,
                                 int k,
                                 int n,
                                 int ng,
                                 IntegerVector nbeta,
                                 NumericMatrix A,
                                 NumericMatrix Y,
                                 Nullable<NumericMatrix> TCOV,
                                 Nullable<NumericVector> deltainit,
                                 int nw){
  int period = A.ncol();
  NumericVector delta;
  NumericVector deltak;
  if (nw!=0){
    NumericVector tmp3(deltainit.get());
    delta = tmp3;
    NumericVector ndeltacum(ng);
    NumericVector deltatmp(ng);
    deltatmp.fill(nw);
    std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
    ndeltacum.push_front(0);
    deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
  }
  NumericVector betas;
  for (int l = 0; l < nbeta[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; ++i){
      double muikt = 0;
      for (int s = 0; s < period; ++s){
        NumericVector vtmp2;
        for (int po = 0; po < nbeta[k]; ++po){
          vtmp2.push_back(pow(A(i,s), po));
        }
        muikt = exp(sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, s, k));
        a += taux(i, k)*pow(A(i, s), l)*(Y(i,s)-muikt/(1+muikt));
      }
    } 
    betas.push_back(a);
  }
  return(betas);
}
// ----------------------------------------------------------------------------
// Differential of deltak Q function
// ----------------------------------------------------------------------------
NumericVector difQdeltakLOGIT_cpp(NumericVector deltak,
                                 NumericMatrix taux,
                                 int k,
                                 int n,
                                 int ng,
                                 IntegerVector nbeta,
                                 NumericMatrix A,
                                 NumericMatrix Y,
                                 NumericMatrix TCOV,
                                 NumericVector beta,
                                 int nw){
  int period = A.ncol();
  NumericVector delta;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
  
  NumericVector deltas;
  for (int l = 0; l < nw; ++l){
    double a = 0;
    for (int i = 0; i < n; ++i){
      double muikt = 0;
      for (int s = 0; s < period; ++s){
        NumericVector vtmp2;
        for (int po = 0; po < nbeta[k]; ++po){
          vtmp2.push_back(pow(A(i,s), po));
        }
        muikt = exp(sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, s, k));
        a += taux(i, k)*TCOV(i, s + l*period)*(Y(i,s)-muikt/(1+muikt));
      }
    } 
    deltas.push_back(a);
  }
  return(deltas);
}
// ----------------------------------------------------------------------------
// Q function for beta
// ----------------------------------------------------------------------------
double QbetaLOGIT_cpp(NumericVector beta,
                      NumericMatrix taux,
                      int n,
                      int ng,
                      IntegerVector nbeta,
                      NumericMatrix A,
                      NumericMatrix Y,
                      Nullable<NumericMatrix> TCOV,
                      Nullable<NumericVector> delta,
                      int nw){
  double a = 0;
  for (int k = 0; k < ng; ++k){
    a += QbetakLOGIT_cpp(beta, taux, k, n, ng, nbeta, A, Y, TCOV, delta, nw);
  }
  return(a);
}
// ----------------------------------------------------------------------------
// differential of Q beta
// ----------------------------------------------------------------------------
NumericVector difQbetaLOGIT_cpp(NumericVector beta,
                      NumericMatrix taux,
                      int n,
                      int ng,
                      IntegerVector nbeta,
                      NumericMatrix A,
                      NumericMatrix Y,
                      Nullable<NumericMatrix> TCOV,
                      Nullable<NumericVector> delta,
                      int nw){
  vec res;
  for (int k = 0; k < ng; ++k){
    NumericVector tmp = difQbetakLOGIT_cpp(beta, taux, k, n, ng, nbeta, A, Y, TCOV, delta, nw);
    res = join_cols(res, as<arma::vec>(tmp));
  }
  return(Rcpp::NumericVector(res.begin(), res.end()));
}

// ----------------------------------------------------------------------------
// 
// Algorithm EM
//
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double likelihoodEMLOGIT_cpp(int n,
                             int ng,
                             IntegerVector nbeta,
                             NumericVector beta,
                             NumericVector pi,
                             NumericMatrix A,
                             NumericMatrix Y,
                             Nullable<NumericMatrix> TCOV,
                             Nullable<NumericVector> delta,
                             int nw){
  double out = 0;
  // create a list for beta
  List betaL(ng);
  int ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nbeta[i]; j++){
      tmp.push_back(beta[ind + j]);
    }
    ind += nbeta[i];
    betaL[i] = tmp;
  }
  // create a list for delta
  List deltaL(ng);
  NumericVector deltatmp(delta.get());
  if (nw != 0){
    int ind = 0;
    for (int i = 0; i < ng; i++){
      NumericVector tmp1;
      for (int j = 0; j < nw; j++){
        tmp1.push_back(deltatmp[ind + j]);
      }
      ind += nw;
      deltaL[i] = tmp1;
    }
  }
  for (int i = 0; i < n; ++i){
    double a = 0;
    for (int s = 0; s < ng; ++s){
      // add 1 top i and s because gkCNORM_cpp is an export function that may be used in R
      a +=  pi[s]*gkLOGIT_cpp(betaL, i, s, nbeta, A, Y, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}
// ----------------------------------------------------------------------------
// EM
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EMLOGIT_cpp(NumericVector param,
                          int ng,
                          int nx,
                          int n,
                          IntegerVector nbeta,
                          NumericMatrix A,
                          NumericMatrix Y,
                          NumericMatrix X,
                          Nullable<NumericMatrix> TCOV,
                          int nw,
                          int itermax,
                          bool EMIRLS,
                          int refgr){
  int period = A.ncol();
  double prec = 0.000001;
  NumericVector pi(ng*nx);
  NumericVector beta;
  NumericVector delta;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    if (param.length() > ng*nx+sum(nbeta)-1){
      delta = param[Range(ng+sum(nbeta)- 1, param.length() - 1)];
    }
    pi.push_back(1-sum(pi));
  }else{
    pi = param[Range(0,(ng-1)*nx-1)];
    for (int i = 0; i < nx; i++){
      pi.push_front(0);
    }
    beta = param[Range((ng-1)*nx,(ng-1)*nx+sum(nbeta)-1)];
    if (param.length() > (ng-1)*nx+sum(nbeta)){
      delta = param[Range((ng-1)*nx+sum(nbeta), param.length() - 1)];
    }
  }
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta), as<arma::rowvec>(delta));
  int tour = 1;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEMLOGIT_cpp(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodLOGIT_cpp(NumericVector(vparam.begin() + nx, vparam.end()), ng, nx, n, nbeta, A, Y, X, TCOV, nw));
    }
    // E-step
    NumericMatrix  taux = ftauxLOGIT_cpp(pi, beta, ng, nbeta, n, A, Y, TCOV, delta, nw, nx, X);
    NumericVector newbeta;
    NumericVector  newdelta;
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];
    
    vec vtmp;
    for (int k = 0; k < ng; ++k){
      NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
      List tmp = optim(Rcpp::_["par"] = betak,
                       Rcpp::_["fn"] = Rcpp::InternalFunction(&QbetakLOGIT_cpp),
                       Rcpp::_["gr"] = Rcpp::InternalFunction(&difQbetakLOGIT_cpp),
                       Rcpp::_["method"] = "BFGS",
                       Rcpp::_["taux"] = taux,
                       Rcpp::_["k"] = k,
                       Rcpp::_["n"] = n,
                       Rcpp::_["ng"] = ng,
                       Rcpp::_["nbeta"] = nbeta,
                       Rcpp::_["A"] = A,
                       Rcpp::_["Y"] = Y,
                       Rcpp::_["TCOV"] = TCOV,
                       Rcpp::_["delta"] = delta,
                       Rcpp::_["nw"] = nw,
                       Rcpp::_["hessian"] =  0,
                       Rcpp::_["control"] = List::create(Named("fnscale")=-1)
      );
      vtmp = join_cols(vtmp, as<arma::vec>(tmp[0]));
    }
    newbeta = NumericVector(vtmp.begin(), vtmp.end());
    
    if (nw!=0){
      vec vtmp;
      NumericVector ndeltacum(ng);
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
      for (int k = 0; k < ng; ++k){
        NumericVector deltak = beta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
        List tmp = optim(Rcpp::_["par"] = deltak,
                         Rcpp::_["fn"] = Rcpp::InternalFunction(&QdeltakLOGIT_cpp),
                         Rcpp::_["gr"] = Rcpp::InternalFunction(&difQdeltakLOGIT_cpp),
                         Rcpp::_["method"] = "BFGS",
                         Rcpp::_["taux"] = taux,
                         Rcpp::_["k"] = k,
                         Rcpp::_["n"] = n,
                         Rcpp::_["ng"] = ng,
                         Rcpp::_["nbeta"] = nbeta,
                         Rcpp::_["A"] = A,
                         Rcpp::_["Y"] = Y,
                         Rcpp::_["TCOV"] = TCOV,
                         Rcpp::_["beta"] = beta,
                         Rcpp::_["nw"] = nw,
                         Rcpp::_["hessian"] =  0,
                         Rcpp::_["control"] = List::create(Named("fnscale")=-1)
        );
        vtmp = join_cols(vtmp, as<arma::vec>(tmp[0]));
      }
      newdelta = NumericVector(vtmp.begin(), vtmp.end());
    }
    // calculus of pi
    NumericVector newpi;
    if (nx == 1){
      NumericVector tmp(ng);
      for (int i = 0; i < ng; ++i){
        tmp[i] = sum(taux(_, i));
      }
      newpi = tmp/n;
    }else{
      newpi = findtheta_cpp(pi, taux, X, n, ng, nx, period, EMIRLS, refgr);
    }
    // stop test
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), as<arma::rowvec>(newbeta), as<arma::rowvec>(newdelta));
    rowvec tmp1(newparam.size());
    tmp1.fill(prec);
    if (all(abs(newparam-vparam)<tmp1)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    delta = newdelta;
    pi = newpi;
  }

  return(NumericVector(vparam.begin(), vparam.end()));
}
// ----------------------------------------------------------------------------
// EM  IRLS
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EMLOGITIRLS_cpp(NumericVector param,
                              int ng, 
                              int nx,
                              int n,
                              IntegerVector nbeta,
                              NumericMatrix A,
                              NumericMatrix Y,
                              NumericMatrix X,
                              Nullable<NumericMatrix> TCOVinit,
                              int nw, 
                              int itermax, 
                              bool EMIRLS,
                              int refgr){
  int period = A.ncol();
  NumericMatrix TCOV;
  double prec = 0.000001;
  NumericVector pi;
  NumericVector beta;
  NumericVector delta;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    if (param.length() > ng*nx+sum(nbeta)-1){
      delta = param[Range(ng+sum(nbeta)- 1, param.length() - 1)];
      TCOV = TCOVinit.get();
    }
    pi.push_back(1-sum(pi));
  }else{
    pi = param[Range(0,(ng-1)*nx-1)];
    for (int i = 0; i < nx; i++){
      pi.push_front(0);
    }
    beta = param[Range((ng-1)*nx,(ng-1)*nx+sum(nbeta)-1)];
    if (param.length() > (ng-1)*nx+sum(nbeta)){
      delta = param[Range((ng-1)*nx+sum(nbeta), param.length() - 1)];
      TCOV = TCOVinit.get();
    }
  }
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta), as<arma::rowvec>(delta));
  int tour = 1;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEMLOGIT_cpp(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodLOGIT_cpp(NumericVector(vparam.begin() + nx, vparam.end()), ng, nx, n, nbeta, A, Y, X, TCOV, nw));
    }
    // E-step
    NumericMatrix  taux = ftauxLOGIT_cpp(pi, beta, ng, nbeta, n, A, Y, TCOV, delta, nw, nx, X);
    rowvec newbeta;
    rowvec newdelta;
    
    if (nw == 0){
      for (int k = 0; k < ng; ++k){
        vec newbetaIRLS;
        NumericVector betaIRLS = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
        arma::vec precIRLS(betaIRLS.length());
        precIRLS.fill(1);
        int stop = 0;
        
        while(all(abs(precIRLS) > 0.000001) &&  stop < 300){
          stop +=1;
          
          mat Aw(nbeta[k], n*period);
          NumericVector tmpr;
          vec S(n*period);
          vec tmp2;
          for (int i = 0; i < n; ++i){
            for (int t = 0; t < period; ++t){
              NumericVector vtmp2;
              for (int po = 0; po < nbeta[k]; ++po){
                vtmp2.push_back(pow(A(i, t), po));
              }
              double betaAit = sum(betaIRLS*vtmp2);
              double rhoikt = exp(betaAit)/(1+exp(betaAit));
              for (int kk = 0; kk < nbeta[k]; ++kk){
                Aw(kk, period*i + t) = pow(A(i, t), kk);
              }
              tmpr.push_back(rhoikt*(1-rhoikt));
              S[period*i + t] = betaAit + (Y(i, t) - rhoikt)/(rhoikt*(1-rhoikt));
            }
            NumericVector vtmp = rep(taux(i, k), period);
            tmp2 = join_cols(tmp2, as<arma::vec>(vtmp));
          }
          // mat W = diagmat(as<arma::vec>(tmpr));
          // mat Z = diagmat(tmp2);
          // newbetaIRLS = inv(Aw*Z*W*trans(Aw))*Aw*Z*W*S;
          
          mat ZW = diagmat(tmp2 % as<arma::vec>(tmpr));
          
          mat Q;
          mat R;
          qr(Q, R, trans(Aw*sqrt(ZW)));
          mat Rc = R.submat(0, 0, nbeta[k] - 1, nbeta[k] - 1);
          mat Qc = Q.submat(0, 0, Q.n_rows - 1, nbeta[k] - 1);
          newbetaIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*S);
          
          precIRLS = Rcpp::as<arma::vec>(betaIRLS)-newbetaIRLS;
          betaIRLS = Rcpp::NumericVector(newbetaIRLS.begin(), newbetaIRLS.end());
        }
        newbeta = join_rows(newbeta, as<arma::rowvec>(betaIRLS));
      }
    }else{
      IntegerVector ndeltacum(ng);
      IntegerVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
      for (int k = 0; k < ng; ++k){
        vec newbetaIRLS;
        vec newdeltaIRLS;
        NumericVector betaIRLS = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
        NumericVector deltaIRLS = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
        
        arma::vec precIRLS(betaIRLS.length());
        precIRLS.fill(1);
        int stop = 0;
        
        
        while(all(abs(precIRLS) > 0.000001) &&  stop < 300){
          stop +=1;
          mat Aw(nbeta[k], n*period);
          mat Ww(nw, n*period);
          NumericVector tmpr;
          vec S(n*period);
          vec Sw(n*period);
          vec tmp2;
          for (int i = 0; i < n; ++i){
            for (int t = 0; t < period; ++t){
              NumericVector vtmp2;
              for (int po = 0; po < nbeta[k]; ++po){
                vtmp2.push_back(pow(A(i, t), po));
              }
              double betaAit = sum(betaIRLS * vtmp2);
              NumericVector vtmp3;
              for (int po = 0; po < nbeta[k]; ++po){
                vtmp3.push_back(TCOV(i, t + po*period));
              }
              double deltaWit = WitEM_cpp(TCOV, period, deltaIRLS, nw, i, t, k);
              double rhoikt = exp(betaAit + deltaWit)/(1+exp(betaAit + deltaWit));
              for (int kk = 0; kk < nbeta[k]; ++kk){
                Aw(kk, period*i + t) = pow(A(i, t), kk);
              }
              for (int kk = 0; kk < nw; ++kk){
                Ww(kk, period*i + t) = TCOV(i, t + kk*period);
              }
              tmpr.push_back(rhoikt*(1-rhoikt));
              S[period*i + t] = betaAit + (Y(i, t) - rhoikt)/(rhoikt*(1-rhoikt));
              Sw[period*i + t] = deltaWit + (Y(i, t) - rhoikt)/(rhoikt*(1-rhoikt));
            }
            NumericVector vtmp = rep(taux(i, k), period);
            tmp2 = join_cols(tmp2, as<arma::vec>(vtmp));
          }
          mat ZW = diagmat(tmp2 % as<arma::vec>(tmpr));
          
          mat Q;
          mat R;
          
          qr(Q, R, trans(Aw*sqrt(ZW)));
          mat Rc = R.submat(0, 0, nbeta[k] - 1, nbeta[k] - 1);
          mat Qc = Q.submat(0, 0, Q.n_rows - 1, nbeta[k] - 1);
          newbetaIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*S);
          
          qr(Q, R, trans(Ww*sqrt(ZW)));
          Rc = R.submat(0, 0, nw - 1, nw - 1);
          Qc = Q.submat(0, 0, Q.n_rows - 1, nw - 1);
          newdeltaIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*Sw);

          precIRLS = join_cols(Rcpp::as<arma::vec>(betaIRLS),Rcpp::as<arma::vec>(deltaIRLS))-join_cols(newbetaIRLS, newdeltaIRLS);
          betaIRLS = Rcpp::NumericVector(newbetaIRLS.begin(), newbetaIRLS.end());
          deltaIRLS = Rcpp::NumericVector(newdeltaIRLS.begin(), newdeltaIRLS.end());
        }
        newbeta = join_rows(newbeta, as<arma::rowvec>(betaIRLS));
        newdelta = join_rows(newdelta, as<arma::rowvec>(deltaIRLS));
      }
    }
    // calculus of pi
    NumericVector newpi;
    if (nx == 1){
      NumericVector tmp(ng);
      for (int i = 0; i < ng; ++i){
        tmp[i] = sum(taux(_, i));
      }
      newpi = tmp/n;
    }else{
      newpi = findtheta_cpp(pi, taux, X, n, ng, nx, period, EMIRLS, refgr);
    }
    // stop test
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), newbeta, newdelta);
    rowvec tmp(newparam.size());
    tmp.fill(prec);
    if (all(abs(newparam-vparam)<tmp)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    delta = newdelta;
    pi = newpi;
  }
  
  return(NumericVector(vparam.begin(), vparam.end()));
}
