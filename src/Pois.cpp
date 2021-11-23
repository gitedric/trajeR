#include "CommonFunction.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include "ZIP.h"
using namespace Rcpp;
using namespace arma;
// [[Rcpp::interfaces(r, cpp)]]


// ----------------------------------------------------------------------------
// gk ZIP
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double gkPois_cpp(List beta,
                  int i,
                  int k,
                  IntegerVector nbeta,
                  NumericMatrix A,
                  NumericMatrix Y,
                  Nullable<NumericMatrix> TCOV,
                  Nullable<List> delta,
                  int nw){
  int period = A.ncol();
  NumericVector lambdaikt = exp(muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k));
  
  double res = 1;
  for (int t = 0; t < period; ++t){
    if (R_IsNA(Y(i, t)) == FALSE){
      if (lambdaikt[t]>20){
        res *= normpdf(Y(i, t), lambdaikt[t], sqrt(lambdaikt[t]));}
      else{
        res *= pow(lambdaikt[t], Y(i, t))*exp(-lambdaikt[t])/facto(Y(i, t));
      }
    }
  }
  return(res);
}

// ----------------------------------------------------------------------------
// dif likelihood theta
// ----------------------------------------------------------------------------
NumericVector difLthetakPois_cpp(NumericVector theta,
                                List beta,
                                Nullable<List> delta,
                                int k,
                                int ng,
                                int nx,
                                IntegerVector nbeta,
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
        tmp1 += vtmp[s]*(gkPois_cpp(beta, i, k, nbeta, A, Y, TCOV, delta, nw)-gkPois_cpp(beta, i, s, nbeta, A, Y, TCOV, delta, nw));
      }
      tmp1 = tmp1*vtmp[k];
      double tmp2 = 0;
      for (int s = 0; s < ng; ++s){
        tmp2 += vtmp[s]*gkPois_cpp(beta, i, s, nbeta, A, Y, TCOV, delta, nw);
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
// [[Rcpp::export]]
NumericVector difLbetakPois_cpp(NumericVector theta,
                                List beta,
                                Nullable<List> delta,
                                int k,
                                int ng,
                                int nx,
                                IntegerVector nbeta,
                                int n,
                                NumericMatrix A,
                                NumericMatrix Y,
                                NumericMatrix X,
                                Nullable<NumericMatrix> TCOV,
                                int nw){
  NumericVector betas;
  int period = A.ncol();
  NumericVector zero(period);
  NumericVector rhoikt(period);
  
  for (int l = 0; l < nbeta[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      double difbkl = 0;
      
      // initialize values
      NumericVector lambdaikt = exp(muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k));

      LogicalVector ind(period);
      for (int t = 0; t < period; t++){
        if (R_IsNA(Y(i, t))){
          ind[t] = FALSE;
        }else{
          ind[t] = TRUE;
        }
      }
      // create a matrix with Yit, muikt ans Ait
      NumericMatrix YimuiA(3, period);
      YimuiA(0, _) = Y(i, _);
      YimuiA(1, _) = lambdaikt;
      YimuiA(2, _) = A(i, _);
      // construct a matrix with element of Yimui for which ind is true
      NumericMatrix mtmp = submat_cpp(YimuiA, ind);
      NumericVector vytmp = mtmp(0, _);
      NumericVector vlambdaikt = mtmp(1, _);
      NumericVector vAtmp = mtmp(2, _);
      for (int t = 0; t < vytmp.size(); ++t){
        
          // vector with all elements except the element t
          NumericVector vlambdaiktminus = vlambdaikt;
          vlambdaiktminus.erase(vlambdaiktminus.begin() + t);
          
          NumericVector vyminus = vytmp;
          vyminus.erase(vyminus.begin() + t);
          NumericVector vlambdaiktminuspower;
          for (int t = 0; t < vlambdaiktminus.size(); ++t){
            vlambdaiktminuspower.push_back(pow(vlambdaiktminus[t], vyminus[t]));
          }
          difbkl += pow(vAtmp[t], l)*pow(vlambdaikt[t], vytmp[t])*exp(-vlambdaikt[t])*(vytmp[t]-vlambdaikt[t])/facto(vytmp[t])*prodvect(vlambdaiktminuspower*exp(-vlambdaiktminus)/factorial(vyminus));
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkPois_cpp(beta, i, s, nbeta, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*difbkl;
    }
    betas.push_back(a);
  }
  return(betas); 
}

// ----------------------------------------------------------------------------
// dif likelihood for POIS
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector difLPois_cpp(NumericVector param,
                          int ng, 
                          int nx,
                          IntegerVector nbeta,
                          int n,
                          NumericMatrix A,
                          NumericMatrix Y,
                          NumericMatrix X,
                          Nullable<NumericMatrix> TCOV,
                          int nw){
  NumericVector out;
  NumericVector theta = param[Range(0,(ng-1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
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
    NumericVector tmp =  difLthetakPois_cpp(theta, betaL, deltaL, k, ng, nx, nbeta, n, A, Y, X, TCOV, nw);
    for (int i =0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLbetakPois_cpp(theta, betaL, deltaL, k, ng, nx, nbeta, n, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  // if (nw != 0){
  //   for (int k = 0; k < ng; ++k){
  //     NumericVector tmp =  difLdeltakZIP_cpp(theta, betaL, nuL, deltaL, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw);
  //     for (int i = 0; i < tmp.length(); ++i){
  //       out.push_back(tmp[i]);
  //     }
  //   }
  // }
  return(out);
}

// ----------------------------------------------------------------------------
// Likelihood 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double likelihoodPois_cpp(NumericVector param,
                           int ng, 
                           int nx,
                           IntegerVector nbeta,
                           int n,
                           NumericMatrix A,
                           NumericMatrix Y,
                           NumericMatrix X,
                           Nullable<NumericMatrix> TCOV,
                           int nw){
  // create a list for beta
  NumericVector theta = param[Range(0,(ng-1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  NumericVector beta = param[Range((ng-1)*nx,(ng-1)*nx+sum(nbeta)-1)];
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
      a +=  piikIntern_cpp(theta, i, s, ng, X)*gkPois_cpp(betaL, i, s, nbeta, A, Y, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}
