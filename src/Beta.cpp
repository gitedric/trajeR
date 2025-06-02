#include "CommonFunction.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::interfaces(r, cpp)]]

// ----------------------------------------------------------------------------
// gk Beta
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
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
                  int nw){
  List deltaNULL;
  double epsilon = std::numeric_limits<double>::epsilon();
  int period = A.ncol();

  NumericVector tmpbeta = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
  NumericVector tmpphi = muikt_cpp(phi[k], nphi[k], i, period, A, TCOV, deltaNULL, 0, k);
  NumericVector muikt = pmax(pmin(1/(1+exp(-tmpbeta)), 1 - epsilon), epsilon);
  NumericVector phikt = exp(tmpphi);
  
  double res = 1;
  for (int t = 0; t < period; ++t){
    if (R_IsNA(Y(i, t)) == FALSE){
     res *= R::dbeta(Y(i, t), muikt[t]*phikt[t], (1-muikt[t])*phikt[t], false);
    }
  }
  return(res);
}
// ----------------------------------------------------------------------------
// dif likelihood theta
// ----------------------------------------------------------------------------
NumericVector difLthetakBETA_cpp(NumericVector theta,
                                  List beta,
                                  List phi,
                                  Nullable<List> delta,
                                  int k,
                                  int ng,
                                  int nx,
                                  IntegerVector nbeta,
                                  IntegerVector nphi,
                                  int n,
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
        tmp1 += vtmp[s]*(gkBETA_cpp(beta, phi, i, k, nbeta, nphi, A, Y, TCOV, delta, nw)-gkBETA_cpp(beta, phi, i, s, nbeta, nphi, A, Y, TCOV, delta, nw));
      }
      tmp1 = tmp1*vtmp[k];
      double tmp2 = 0;
      for (int s = 0; s < ng; ++s){
        tmp2 += vtmp[s]*gkBETA_cpp(beta, phi, i, s, nbeta, nphi, A, Y, TCOV, delta, nw);
      }
      a += X(i, l)*tmp1/(sum(vtmp)*tmp2);
    }
    thetas.push_back(a);
  }
  
  return(thetas);
}
// ----------------------------------------------------------------------------
// dif likelihood beta
// ----------------------------------------------------------------------------
NumericVector difLbetakBETA_cpp(NumericVector theta,
                                 List beta,
                                 List phi,
                                 Nullable<List> delta,
                                 int k,
                                 int ng,
                                 int nx,
                                 IntegerVector nbeta,
                                 IntegerVector nphi,
                                 int n,
                                 NumericMatrix A,
                                 NumericMatrix Y,
                                 NumericMatrix X,
                                 Nullable<NumericMatrix> TCOV,
                                 int nw){
  NumericVector betas;
  List deltaNULL;
  double epsilon = std::numeric_limits<double>::epsilon();
  int period = A.ncol();
  
  for (int l = 0; l < nbeta[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      
      NumericVector tmpbeta = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      NumericVector tmpphi = muikt_cpp(phi[k], nphi[k], i, period, A, TCOV, deltaNULL, 0, k);
      NumericVector muikt = pmax(pmin(1/(1+exp(-tmpbeta)), 1 - epsilon), epsilon);
      NumericVector phikt = exp(tmpphi);
      double tmp2 = 0;
      NumericVector tmp1;
      
      for (int t = 0; t < period; ++t){
        if (R_IsNA(Y(i, t)) == FALSE){
          tmp1.push_back(R::dbeta(Y(i, t), muikt[t]*phikt[t], (1-muikt[t])*phikt[t], false));
      }
      }
      for (int t = 0; t < period; ++t){
        if (R_IsNA(Y(i, t)) == FALSE){
          tmp2 += exp(tmpbeta[t])/pow(1 + exp(tmpbeta[t]), 2)*pow(A(i,t),l)*phikt[t]*(log(Y(i,t)/(1-Y(i,t))) - R::digamma(muikt[t]*phikt[t]) + R::digamma((1-muikt[t])*phikt[t]))*prodvect(tmp1);
        }
      }
      
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkBETA_cpp(beta, phi, i, s, nbeta, nphi, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*tmp2;
    }
    betas.push_back(a);
  }
  return(betas); 
}
// ----------------------------------------------------------------------------
// dif likelihood delta
// ----------------------------------------------------------------------------
NumericVector difLdeltakBETA_cpp(NumericVector theta,
                                List beta,
                                List phi,
                                Nullable<List> delta,
                                int k,
                                int ng,
                                int nx,
                                IntegerVector nbeta,
                                IntegerVector nphi,
                                int n,
                                NumericMatrix A,
                                NumericMatrix Y,
                                NumericMatrix X,
                                Nullable<NumericMatrix> initTCOV,
                                int nw){
  NumericVector deltas;
  List deltaNULL;
  NumericMatrix TCOV(initTCOV);
  double epsilon = std::numeric_limits<double>::epsilon();
  int period = A.ncol();
  
  for (int l = 0; l < nw; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      
      NumericVector tmpbeta = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      NumericVector tmpphi = muikt_cpp(phi[k], nphi[k], i, period, A, TCOV, deltaNULL, 0, k);
      NumericVector muikt = pmax(pmin(1/(1+exp(-tmpbeta)), 1 - epsilon), epsilon);
      NumericVector phikt = exp(tmpphi);
      double tmp2 = 0;
      NumericVector tmp1;
      
      for (int t = 0; t < period; ++t){
        if (R_IsNA(Y(i, t)) == FALSE){
          tmp1.push_back(R::dbeta(Y(i, t), muikt[t]*phikt[t], (1-muikt[t])*phikt[t], false));
        }
      }
      for (int t = 0; t < period; ++t){
        if ((R_IsNA(Y(i, t)) == FALSE) && (R_IsNA(TCOV(i, t + l*period)) == FALSE)){
          tmp2 += exp(tmpbeta[t])/pow(1 + exp(tmpbeta[t]), 2)*TCOV(i, t + l*period)*phikt[t]*(log(Y(i,t)/(1-Y(i,t))) - R::digamma(muikt[t]*phikt[t]) + R::digamma((1-muikt[t])*phikt[t]))*prodvect(tmp1);
        }
      }
      
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkBETA_cpp(beta, phi, i, s, nbeta, nphi, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*tmp2;
    }
    deltas.push_back(a);
  }
  return(deltas); 
}
// ----------------------------------------------------------------------------
// dif likelihood phi
// ----------------------------------------------------------------------------
NumericVector difLphikBETA_cpp(NumericVector theta,
                                List beta,
                                List phi,
                                Nullable<List> delta,
                                int k,
                                int ng,
                                int nx,
                                IntegerVector nbeta,
                                IntegerVector nphi,
                                int n,
                                NumericMatrix A,
                                NumericMatrix Y,
                                NumericMatrix X,
                                Nullable<NumericMatrix> TCOV,
                                int nw){
  NumericVector phis;
  List deltaNULL;
  double epsilon = std::numeric_limits<double>::epsilon();
  int period = A.ncol();
  
  for (int l = 0; l < nphi[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      
      NumericVector tmpbeta = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      NumericVector tmpphi = muikt_cpp(phi[k], nphi[k], i, period, A, TCOV, deltaNULL, 0, k);
      NumericVector muikt = pmax(pmin(1/(1+exp(-tmpbeta)), 1 - epsilon), epsilon);
      NumericVector phikt = exp(tmpphi);
      double tmp2 = 0;
      NumericVector tmp1;
      
      for (int t = 0; t < period; ++t){
        if (R_IsNA(Y(i, t)) == FALSE){
          tmp1.push_back(R::dbeta(Y(i, t), muikt[t]*phikt[t], (1-muikt[t])*phikt[t], false));
        }
      }
      for (int t = 0; t < period; ++t){
        if (R_IsNA(Y(i, t)) == FALSE){
          tmp2 += phikt[t]*pow(A(i,t),l)*(muikt[t]*(log(Y(i,t)/(1-Y(i,t)))-R::digamma(muikt[t]*phikt[t]) + R::digamma((1-muikt[t])*phikt[t]))+log(1-Y(i,t))-R::digamma((1-muikt[t])*phikt[t])+R::digamma(phikt[t]))*prodvect(tmp1);
        }
      }
      
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkBETA_cpp(beta, phi, i, s, nbeta, nphi, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*tmp2;
    }
    phis.push_back(a);
  }
  return(phis); 
}
// ----------------------------------------------------------------------------
// dif likelihood BETA
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector difLBETA_cpp(NumericVector param,
                            int ng, 
                            int nx,
                            IntegerVector nbeta,
                            IntegerVector nphi,
                            int n,
                            NumericMatrix A,
                            NumericMatrix Y,
                            NumericMatrix X,
                            Nullable<NumericMatrix> TCOV,
                            int nw){
  NumericVector out;
  NumericVector theta = param[Range(0,(ng - 1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  NumericVector beta = param[Range((ng - 1)*nx, (ng - 1)*nx+sum(nbeta)-1)];
  NumericVector phi = param[Range((ng - 1)*nx+sum(nbeta), (ng - 1)*nx+sum(nbeta)+sum(nphi)-1)];
  // create a list 
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
  List phiL(ng);
  ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nphi[i]; j++){
      tmp.push_back(phi[ind + j]);
    }
    ind += nphi[i];
    phiL[i] = tmp;
  }
  NumericVector delta;
  List deltaL(ng);
  if (param.length() > (ng - 1)*nx+sum(nbeta)+sum(nphi)){
    delta = param[Range((ng - 1)*nx+sum(nbeta)+sum(nphi), param.length())];
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

  for (int k = 1; k < ng; ++k){
    NumericVector tmp =  difLthetakBETA_cpp(theta, betaL, phiL, deltaL, k, ng, nx, nbeta, nphi, n, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLbetakBETA_cpp(theta, betaL, phiL, deltaL, k, ng, nx, nbeta, nphi, n, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLphikBETA_cpp(theta, betaL, phiL, deltaL, k, ng, nx, nbeta, nphi, n, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  if (nw != 0){
    for (int k = 0; k < ng; ++k){
      NumericVector tmp =  difLdeltakBETA_cpp(theta, betaL, phiL, deltaL, k, ng, nx, nbeta, nphi, n, A, Y, X, TCOV, nw);
      for (int i = 0; i < tmp.length(); ++i){
        out.push_back(tmp[i]);
      }
    }
  }
  return(out);
}
// ----------------------------------------------------------------------------
// Likelihood sigmawith same sigma and exponential reparmetrization
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double LikelihoodBETA_cpp(NumericVector param,
                          int ng, 
                          int nx,
                          IntegerVector nbeta,
                          IntegerVector nphi,
                          int n,
                          NumericMatrix A,
                          NumericMatrix Y,
                          NumericMatrix X,
                          Nullable<NumericMatrix> TCOV,
                          int nw){
  NumericVector theta = param[Range(0,(ng - 1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  NumericVector beta = param[Range((ng - 1)*nx, (ng - 1)*nx+sum(nbeta)-1)];
  NumericVector phi = param[Range((ng - 1)*nx+sum(nbeta), (ng - 1)*nx+sum(nbeta)+sum(nphi)-1)];
  // create a list 
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
  List phiL(ng);
  ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nphi[i]; j++){
      tmp.push_back(phi[ind + j]);
    }
    ind += nphi[i];
    phiL[i] = tmp;
  }
  NumericVector delta;
  List deltaL(ng);
  if (param.length() > (ng - 1)*nx+sum(nbeta)+sum(nphi)){
    delta = param[Range((ng - 1)*nx+sum(nbeta)+sum(nphi), param.length() - 1)];
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

  double out = 0;
  for (int i =0; i < n; ++i){
    double a = 0;
    for (int s = 0; s < ng; ++s){
      a +=  piikIntern_cpp(theta, i, s, ng, X)*gkBETA_cpp(betaL, phiL, i, s, nbeta, nphi, A, Y, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}


