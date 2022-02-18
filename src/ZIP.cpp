#include "CommonFunction.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::interfaces(r, cpp)]]

// compute sum of product nu by AI
NumericVector nuikt_cpp(NumericVector nu,
                        int nnu,
                        int i,
                        int period,
                        NumericMatrix A, 
                        int k){
  NumericVector nuikt;
  for (int s = 0; s < period; ++s){
    NumericVector vtmp2;
    for (int po = 0; po < nnu; ++po){
      vtmp2.push_back(pow(A(i,s), po));
    }
    nuikt.push_back(sum(nu*vtmp2));
  }
  return(nuikt);
}
// ----------------------------------------------------------------------------
// gk ZIP
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double gkZIP_cpp(List beta,
                 List nu,
                 int i,
                 int k,
                 IntegerVector nbeta,
                 IntegerVector nnu,
                 NumericMatrix A,
                 NumericMatrix Y,
                 Nullable<NumericMatrix> TCOV,
                 Nullable<List> delta,
                 int nw){
  int period = A.ncol();
  NumericVector rhoikt = 1/(1+exp(-nuikt_cpp(nu[k], nnu[k], i, period, A, k)));
  NumericVector lambdaikt = exp(muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k));
  
  double res = 1;
  for (int t = 0; t < period; ++t){
    if (R_IsNA(Y(i, t)) == FALSE){
      if (Y(i, t) == 0){
        res *= rhoikt[t]+(1-rhoikt[t])*exp(-lambdaikt[t]);
      }else{
        res *= (1-rhoikt[t])*pow(lambdaikt[t], Y(i, t))*exp(-lambdaikt[t])/facto(Y(i, t));
      }
    }
  }
  
  return(res);
}
// ----------------------------------------------------------------------------
// dif likelihood theta
// ----------------------------------------------------------------------------
NumericVector difLthetakZIP_cpp(NumericVector theta,
                                List beta,
                                List nu,
                                Nullable<List> delta,
                                int k,
                                int ng,
                                int nx,
                                IntegerVector nbeta,
                                IntegerVector nnu,
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
        tmp1 += vtmp[s]*(gkZIP_cpp(beta, nu, i, k, nbeta, nnu, A, Y, TCOV, delta, nw)-gkZIP_cpp(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw));
      }
      tmp1 = tmp1*vtmp[k];
      double tmp2 = 0;
      for (int s = 0; s < ng; ++s){
        tmp2 += vtmp[s]*gkZIP_cpp(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw);
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
NumericVector difLbetakZIP_cpp(NumericVector theta,
                               List beta,
                               List nu,
                               Nullable<List> delta,
                               int k,
                               int ng,
                               int nx,
                               IntegerVector nbeta,
                               IntegerVector nnu,
                               int n,
                               NumericMatrix A,
                               NumericMatrix Y,
                               NumericMatrix X,
                               Nullable<NumericMatrix> TCOV,
                               int nw){
  NumericVector betas;
  int period = A.ncol();
  NumericVector zero(period);
  
  for (int l = 0; l < nbeta[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      double difbkl = 0;
      double difbkl0 = 0;
      double py = 1;
      double py0 = 1;
      // initialize values
      NumericVector rhoikt = 1/(1+exp(-nuikt_cpp(nu[k], nnu[k], i, period, A, k)));
      NumericVector lambdaikt = exp(muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k));
      // create a matrix with Yit, muikt ans Ait
      NumericMatrix YimuiA(4, period);
      YimuiA(0, _) = Y(i, _);
      YimuiA(1, _) = lambdaikt;
      YimuiA(2, _) = A(i, _);
      YimuiA(3, _) = rhoikt;
      // check if there exist some censored values
      LogicalVector ind0 = (Y(i, _) == zero);
      LogicalVector ind = !(ind0);
      for (int t = 0; t < period; t++){
        if (R_IsNA(Y(i, t))){
          ind0[t] = FALSE;
          ind[t] = FALSE;
        }
      }
      if (is_true(any(ind0 == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind0);
        NumericVector vytmp = mtmp(0, _);
        NumericVector vlambdaikt = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vrhoikt = mtmp(3, _);
        
        if (sum(ind0) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vlambdaiktminus = vlambdaikt;
            vlambdaiktminus.erase(vlambdaiktminus.begin() + t);
            NumericVector vrhoiktminus = vrhoikt;
            vrhoiktminus.erase(vrhoiktminus.begin() + t);
            difbkl0 -= pow(vAtmp[t], l)*vlambdaikt[t]*(1-vrhoikt[t])*exp(-vlambdaikt[t])*prodvect(vrhoiktminus + (1-vrhoiktminus)*exp(-vlambdaiktminus));
          }
        }else{
          difbkl0 = -pow(vAtmp[0], l)*vlambdaikt[0]*(1-vrhoikt[0])*exp(-vlambdaikt[0]);
        }
        py0 = prodvect(vrhoikt +(1-vrhoikt)*exp(-vlambdaikt));
      }
      if (is_true(any(ind == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind);
        NumericVector vytmp = mtmp(0, _);
        NumericVector vlambdaikt = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vrhoikt = mtmp(3, _);
        if (sum(ind) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vlambdaiktminus = vlambdaikt;
            vlambdaiktminus.erase(vlambdaiktminus.begin() + t);
            NumericVector vrhoiktminus = vrhoikt;
            vrhoiktminus.erase(vrhoiktminus.begin() + t);
            NumericVector vyminus = vytmp;
            vyminus.erase(vyminus.begin() + t);
            NumericVector vlambdaiktminuspower;
            for (int t = 0; t < vlambdaiktminus.size(); ++t){
              vlambdaiktminuspower.push_back(pow(vlambdaiktminus[t], vyminus[t]));
            }
            difbkl += pow(vAtmp[t], l)*pow(vlambdaikt[t], vytmp[t])*(1-vrhoikt[t])*exp(-vlambdaikt[t])*(vytmp[t]-vlambdaikt[t])/facto(vytmp[t])*prodvect((1-vrhoiktminus)*vlambdaiktminuspower*exp(-vlambdaiktminus)/factorial(vyminus));
          }
        }else{
          difbkl =  pow(vAtmp[0], l)*pow(vlambdaikt[0], vytmp[0])*(1-vrhoikt[0])*exp(-vlambdaikt[0])*(vytmp[0]-vlambdaikt[0])/facto(vytmp[0]);
        }
        NumericVector vlambdaiktpower;
        for (int t = 0; t < vlambdaikt.size(); ++t){
          vlambdaiktpower.push_back(pow(vlambdaikt[t], vytmp[t]));
        }
        py = prodvect((1-vrhoikt)*vlambdaiktpower*exp(-vlambdaikt)/factorial(vytmp));
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkZIP_cpp(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*(difbkl0*py+difbkl*py0);
    }
    betas.push_back(a);
  }
  return(betas); 
}
// ----------------------------------------------------------------------------
// dif likelihood nuk
// ----------------------------------------------------------------------------
NumericVector difLnukZIP_cpp(NumericVector theta,
                             List beta,
                             List nu,
                             Nullable<List> delta,
                             int k,
                             int ng,
                             int nx,
                             IntegerVector nbeta,
                             IntegerVector nnu,
                             int n,
                             NumericMatrix A,
                             NumericMatrix Y,
                             NumericMatrix X,
                             Nullable<NumericMatrix> TCOV,
                             int nw){
  NumericVector nus;
  int period = A.ncol();
  NumericVector zero(period);
  
  for (int l = 0; l < nnu[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      double difbkl = 0;
      double difbkl0 = 0;
      double py = 1;
      double py0 = 1;
      // initialize values
      NumericVector enuikt = exp(nuikt_cpp(nu[k], nnu[k], i, period, A, k));
      NumericVector rhoikt = enuikt/(1+enuikt);
      NumericVector lambdaikt = exp(muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k));
      // create a matrix with Yit, muikt ans Ait
      NumericMatrix YimuiA(5, period);
      YimuiA(0, _) = Y(i, _);
      YimuiA(1, _) = lambdaikt;
      YimuiA(2, _) = A(i, _);
      YimuiA(3, _) = rhoikt;
      YimuiA(4, _) = enuikt;
      // check if there exist some censored values
      LogicalVector ind0 = (Y(i, _) == zero);
      LogicalVector ind = !(ind0);
      for (int t = 0; t < period; t++){
        if (R_IsNA(Y(i, t))){
          ind0[t] = FALSE;
          ind[t] = FALSE;
        }
      }
      if (is_true(any(ind0 == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind0);
        NumericVector vytmp = mtmp(0, _);
        NumericVector vlambdaikt = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vrhoikt = mtmp(3, _);
        NumericVector venuikt = mtmp(4, _);
        
        if (sum(ind0) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vlambdaiktminus = vlambdaikt;
            vlambdaiktminus.erase(vlambdaiktminus.begin() + t);
            NumericVector vrhoiktminus = vrhoikt;
            vrhoiktminus.erase(vrhoiktminus.begin() + t);
            NumericVector venuiktminus = venuikt;
            venuiktminus.erase(venuiktminus.begin() + t);
            difbkl0 += pow(vAtmp[t], l)*vrhoikt[t]*(1-exp(-vlambdaikt[t]))/(1+venuikt[t])*prodvect(vrhoiktminus + (1-vrhoiktminus)*exp(-vlambdaiktminus));
          }
        }else{
          difbkl0 = pow(vAtmp[0], l)*vrhoikt[0]*(1-exp(-vlambdaikt[0]))/(1+venuikt[0]);
        }
        py0 = prodvect(vrhoikt +(1-vrhoikt)*exp(-vlambdaikt));
      }
      if (is_true(any(ind == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind);
        NumericVector vytmp = mtmp(0, _);
        NumericVector vlambdaikt = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vrhoikt = mtmp(3, _);
        NumericVector venuikt = mtmp(4, _);
        if (sum(ind) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vlambdaiktminus = vlambdaikt;
            vlambdaiktminus.erase(vlambdaiktminus.begin() + t);
            NumericVector vrhoiktminus = vrhoikt;
            vrhoiktminus.erase(vrhoiktminus.begin() + t);
            NumericVector vyminus = vytmp;
            vyminus.erase(vyminus.begin() + t);
            NumericVector vlambdaiktminuspower;
            NumericVector venuiktminus = venuikt;
            venuiktminus.erase(venuiktminus.begin() + t);
            for (int t = 0; t < vlambdaiktminus.size(); ++t){
              vlambdaiktminuspower.push_back(pow(vlambdaiktminus[t], vyminus[t]));
            }
            difbkl -= pow(vAtmp[t], l)*vrhoikt[t]*pow(vlambdaikt[t], vytmp[t])*exp(-vlambdaikt[t])/(facto(vytmp[t])*(1+venuikt[t]))*prodvect((1-vrhoiktminus)*vlambdaiktminuspower*exp(-vlambdaiktminus)/factorial(vyminus));
          }
        }else{
          difbkl =  -pow(vAtmp[0], l)*vrhoikt[0]*pow(vlambdaikt[0], vytmp[0])*exp(-vlambdaikt[0])/(facto(vytmp[0])*(1+venuikt[0]));
        }
        NumericVector vlambdaiktpower;
        for (int t = 0; t < vlambdaikt.size(); ++t){
          vlambdaiktpower.push_back(pow(vlambdaikt[t], vytmp[t]));
        }
        py = prodvect((1-vrhoikt)*vlambdaiktpower*exp(-vlambdaikt)/factorial(vytmp));
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkZIP_cpp(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*(difbkl0*py+difbkl*py0);
    }
    nus.push_back(a);
  }
  return(nus); 
}
// ----------------------------------------------------------------------------
// dif likelihood deltak
// ----------------------------------------------------------------------------
NumericVector difLdeltakZIP_cpp(NumericVector theta,
                                List beta,
                                List nu,
                                Nullable<List> delta,
                                int k,
                                int ng,
                                int nx,
                                IntegerVector nbeta,
                                IntegerVector nnu,
                                int n,
                                NumericMatrix A,
                                NumericMatrix Y,
                                NumericMatrix X,
                                Nullable<NumericMatrix> TCOVinit,
                                int nw){
  NumericVector betas;
  int period = A.ncol();
  NumericVector zero(period);
  NumericMatrix TCOV(TCOVinit);
  
  for (int l = 0; l < nw; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      double difbkl = 0;
      double difbkl0 = 0;
      double py = 1;
      double py0 = 1;
      // initialize values
      NumericVector rhoikt = 1/(1+exp(-nuikt_cpp(nu[k], nnu[k], i, period, A, k)));
      NumericVector lambdaikt = exp(muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k));
      // create a matrix with Yit, muikt ans Ait
      NumericMatrix YimuiA(4, period);
      YimuiA(0, _) = Y(i, _);
      YimuiA(1, _) = lambdaikt;
      for (int t = 0; t < period; t++){
        YimuiA(2, t) = TCOV(i, t + l*period);   
      }
      YimuiA(3, _) = rhoikt;
      // check if there exist some censored values
      LogicalVector ind0 = (Y(i, _) == zero);
      LogicalVector ind = !(ind0);
      for (int t = 0; t < period; t++){
        if (R_IsNA(Y(i, t))){
          ind0[t] = FALSE;
          ind[t] = FALSE;
        }
      }
      if (is_true(any(ind0 == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind0);
        NumericVector vytmp = mtmp(0, _);
        NumericVector vlambdaikt = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vrhoikt = mtmp(3, _);
        
        if (sum(ind0) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vlambdaiktminus = vlambdaikt;
            vlambdaiktminus.erase(vlambdaiktminus.begin() + t);
            NumericVector vrhoiktminus = vrhoikt;
            vrhoiktminus.erase(vrhoiktminus.begin() + t);
            difbkl0 -= vAtmp[t]*vlambdaikt[t]*(1-vrhoikt[t])*exp(-vlambdaikt[t])*prodvect(vrhoiktminus + (1-vrhoiktminus)*exp(-vlambdaiktminus));
          }
        }else{
          difbkl0 = -vAtmp[0]*vlambdaikt[0]*(1-vrhoikt[0])*exp(-vlambdaikt[0]);
        }
        py0 = prodvect(vrhoikt +(1-vrhoikt)*exp(-vlambdaikt));
      }
      if (is_true(any(ind == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind);
        NumericVector vytmp = mtmp(0, _);
        NumericVector vlambdaikt = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vrhoikt = mtmp(3, _);
        if (sum(ind) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vlambdaiktminus = vlambdaikt;
            vlambdaiktminus.erase(vlambdaiktminus.begin() + t);
            NumericVector vrhoiktminus = vrhoikt;
            vrhoiktminus.erase(vrhoiktminus.begin() + t);
            NumericVector vyminus = vytmp;
            vyminus.erase(vyminus.begin() + t);
            NumericVector vlambdaiktminuspower;
            for (int t = 0; t < vlambdaiktminus.size(); ++t){
              vlambdaiktminuspower.push_back(pow(vlambdaiktminus[t], vyminus[t]));
            }
            difbkl += vAtmp[t]*pow(vlambdaikt[t], vytmp[t])*(1-vrhoikt[t])*exp(-vlambdaikt[t])*(vytmp[t]-vlambdaikt[t])/facto(vytmp[t])*prodvect((1-vrhoiktminus)*vlambdaiktminuspower*exp(-vlambdaiktminus)/factorial(vyminus));
          }
        }else{
          difbkl =  vAtmp[0]*pow(vlambdaikt[0], vytmp[0])*(1-vrhoikt[0])*exp(-vlambdaikt[0])*(vytmp[0]-vlambdaikt[0])/facto(vytmp[0]);
        }
        NumericVector vlambdaiktpower;
        for (int t = 0; t < vlambdaikt.size(); ++t){
          vlambdaiktpower.push_back(pow(vlambdaikt[t], vytmp[t]));
        }
        py = prodvect((1-vrhoikt)*vlambdaiktpower*exp(-vlambdaikt)/factorial(vytmp));
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkZIP_cpp(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*(difbkl0*py+difbkl*py0);
    }
    betas.push_back(a);
  }
  return(betas); 
}
// ----------------------------------------------------------------------------
// dif likelihood for ZIP
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector difLZIP_cpp(NumericVector param,
                          int ng, 
                          int nx,
                          IntegerVector nbeta,
                          IntegerVector nnu,
                          int n,
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
  NumericVector beta = param[Range((ng-1)*nx,(ng-1)*nx+sum(nbeta)-1)];
  //NumericVector beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
  NumericVector nu = param[Range((ng-1)*nx+sum(nbeta), (ng-1)*nx+sum(nbeta)+sum(nnu) - 1)];
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
  List nuL(ng);
  ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nnu[i]; j++){
      tmp.push_back(nu[ind + j]);
    }
    ind += nnu[i];
    nuL[i] = tmp;
  }
  NumericVector delta;
  List deltaL(ng);
  if (param.length() > (ng-1)*nx+sum(nbeta)+sum(nnu)){
    delta = param[Range((ng-1)*nx+sum(nbeta)+sum(nnu), param.length() -1)];
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
  //for (int k = 0; k < ng; ++k){
  for (int k = 1; k < ng; ++k){
    NumericVector tmp =  difLthetakZIP_cpp(theta, betaL, nuL, deltaL, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw);
    for (int i =0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLbetakZIP_cpp(theta, betaL, nuL, deltaL, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLnukZIP_cpp(theta, betaL, nuL, deltaL, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  if (nw != 0){
    for (int k = 0; k < ng; ++k){
      NumericVector tmp =  difLdeltakZIP_cpp(theta, betaL, nuL, deltaL, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw);
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
double likelihoodZIP_cpp(NumericVector param,
                         int ng, 
                         int nx,
                         IntegerVector nbeta,
                         IntegerVector nnu,
                         int n,
                         NumericMatrix A,
                         NumericMatrix Y,
                         NumericMatrix X,
                         Nullable<NumericMatrix> TCOV,
                         int nw){
  // NumericVector theta = param[Range(0,ng*nx-1)];
  // NumericVector beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
  // NumericVector nu = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+sum(nnu))];
  // create a list for beta
  NumericVector theta = param[Range(0,(ng-1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  NumericVector beta = param[Range((ng-1)*nx,(ng-1)*nx+sum(nbeta)-1)];
  NumericVector nu = param[Range((ng-1)*nx+sum(nbeta), (ng-1)*nx+sum(nbeta)+sum(nnu) - 1)];
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
  List nuL(ng);
  ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nnu[i]; j++){
      tmp.push_back(nu[ind + j]);
    }
    ind += nnu[i];
    nuL[i] = tmp;
  }
  NumericVector delta;
  List deltaL(ng);
  if (param.length() > (ng-1)*nx+sum(nbeta)+sum(nnu)){
    delta = param[Range((ng-1)*nx+sum(nbeta)+sum(nnu), param.length() - 1)];
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
      a +=  piikIntern_cpp(theta, i, s, ng, X)*gkZIP_cpp(betaL, nuL, i, s, nbeta, nnu, A, Y, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}

// ----------------------------------------------------------------------------
// 
// Algorithm EM
//
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double likelihoodEMZIP_cpp(int n,
                           int ng,
                           IntegerVector nbeta,
                           IntegerVector nnu,
                           NumericVector beta,
                           NumericVector nu,
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
  // create a list for nu
  List nuL(ng);
  ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nnu[i]; j++){
      tmp.push_back(nu[ind + j]);
    }
    ind += nnu[i];
    nuL[i] = tmp;
  }
  // create a list for delta
  List deltaL(ng);
  if (nw != 0 ) {
    NumericVector deltatmp(delta.get());
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
      a +=  pi[s]*gkZIP_cpp(betaL, nuL, i, s, nbeta, nnu, A, Y, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}

// ----------------------------------------------------------------------------
// Function rate
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
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
  // create a list for nu
  List nuL(ng);
  ind = 0;
  for (int i = 0; i < ng; i++){
    NumericVector tmp;
    for (int j = 0; j < nnu[i]; j++){
      tmp.push_back(nu[ind + j]);
    }
    ind += nnu[i];
    nuL[i] = tmp;
  }
  // create a list for delta
  List deltaL(ng);
  if (nw != 0 ) {
    NumericVector deltatmp(delta.get());
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
  NumericMatrix mtmp(n,ng);
  if (nx == 1){
    for (int i = 0; i < n; ++i){
      double s = 0;
      for (int k = 0; k < ng; ++k){
        mtmp(i,k) = pi[k]*gkZIP_cpp(betaL, nuL, i, k, nbeta, nnu, A, Y, TCOV, deltaL, nw);
        s += mtmp(i,k);
      } 
      mtmp(i, _) = 1/(1+ (s - mtmp(i, _))/mtmp(i, _));
    }
  }else{
    for (int i = 0; i < n; ++i){
      double s = 0;
      for (int k = 0; k < ng; ++k){
        mtmp(i,k) = piikIntern_cpp(pi, i, k, ng, X)*gkZIP_cpp(betaL, nuL, i, k, nbeta, nnu, A, Y, TCOV, deltaL, nw);
        s += mtmp(i,k);
      } 
      mtmp(i, _) = 1/(1+ (s - mtmp(i, _))/mtmp(i, _));
    }
  }
  return(mtmp);
}
// ----------------------------------------------------------------------------
// Function rate
// ----------------------------------------------------------------------------
double fzkSikt_cpp(NumericVector pi,
                   NumericVector beta,
                   NumericVector nu,
                   NumericMatrix zk,
                   int k,
                   int i, 
                   int t,
                   int ng,
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
                   IntegerVector nnucum){
  double prob;
  if (Y(i, t) > 0){
    prob = 0;
  }else{
    NumericVector deltak;
    if (nw != 0){
      NumericVector vtmp(delta.get());
      NumericVector ntmp(ndeltacum.get());
      deltak = vtmp[Range(ntmp[k], ntmp[k+1]-1)]; 
    }
    NumericVector nuk = nu[Range(nnucum[k], nnucum[k+1]-1)];
    NumericVector vtmp1;
    for (int po = 0; po < nnu[k]; ++po){
      vtmp1.push_back(pow(A(i,t), po));
    }
    double nuikt = sum(nuk*vtmp1);
    NumericVector vtmp2;
    for (int po = 0; po < nbeta[k]; ++po){
      vtmp2.push_back(pow(A(i,t), po));
    }
    NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
    double lambdaikt =  exp(sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, t, k));
    prob = zk(i, k)/(1+exp(-nuikt-lambdaikt));
  }
  return(prob);
}
// ----------------------------------------------------------------------------
// Function rate
// ----------------------------------------------------------------------------
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
                 IntegerVector nnucum){
  double prob;
  if (Y(i, t) > 0){
    prob = 0;
  }else{
    NumericVector deltak;
    if (nw != 0){
      NumericVector vtmp(delta.get());
      NumericVector ntmp(ndeltacum.get());
      deltak = vtmp[Range(ntmp[k], ntmp[k+1]-1)]; 
    }
    NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
    NumericVector nuk = nu[Range(nnucum[k], nnucum[k+1]-1)];
    NumericVector vtmp1;
    for (int po = 0; po < nnu[k]; ++po){
      vtmp1.push_back(pow(A(i,t), po));
    }
    double nuikt = sum(nuk*vtmp1);
    NumericVector vtmp2;
    for (int po = 0; po < nbeta[k]; ++po){
      vtmp2.push_back(pow(A(i,t), po));
    }
    double lambdaikt =  exp(sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, t, k));
    prob = 1/(1+exp(-nuikt-lambdaikt));
  }
  return(prob);
}

// ----------------------------------------------------------------------------
//  Q function for betak, nuk and deltak given a value k.  The goal is to sum all Qbetaknuk after
// ----------------------------------------------------------------------------
double QbetakZIP_cpp(NumericVector beta,
                     NumericMatrix zk, 
                     NumericMatrix Sikt,
                     int k, 
                     int nbeta,
                     int nnu,
                     int n,
                     NumericMatrix A,
                     NumericMatrix Y,
                     Nullable<NumericMatrix> TCOV,
                     Nullable<NumericVector> delta,
                     int nw,
                     Nullable<IntegerVector> ndeltacum){
  int period = A.ncol();
  double a = 0;
  for (int i = 0; i < n; ++i){
    double zik = zk(i, k);
    for (int t = 0; t < period; ++t){
      double sikt = Sikt(i, t);
      NumericVector vtmp2;
      for (int po = 0; po < nbeta; ++po){
        vtmp2.push_back(pow(A(i,t), po));
      }
      double betaAit = sum(beta*vtmp2) + WitEM_cpp(TCOV, period, delta, nw, i, t, k); 
      a += zik*(1 - sikt)*(Y(i,t)*betaAit - exp(betaAit));
    }
  }
  return(a);
}

double QnukZIP_cpp(NumericVector nu,
                   NumericMatrix zk, 
                   NumericMatrix Sikt,
                   int k, 
                   int nbeta,
                   int nnu,
                   int n,
                   NumericMatrix A,
                   NumericMatrix Y){
  int period = A.ncol();
  double a = 0;
  for (int i = 0; i < n; ++i){
    double zik = zk(i, k);
    for (int t = 0; t < period; ++t){
      double sikt = Sikt(i, t);
      NumericVector vtmp1;
      for (int po = 0; po < nnu; ++po){
        vtmp1.push_back(pow(A(i,t), po));
      }
      double nuAit = sum(nu*vtmp1);
      a += zik*(sikt*nuAit - log(1+exp(nuAit)));
    }
  }
  return(a);
}

double QdeltakZIP_cpp(NumericVector delta,
                      NumericMatrix zk, 
                      NumericMatrix Sikt,
                      int k, 
                      int nbeta,
                      int nnu,
                      int n,
                      NumericMatrix A,
                      NumericMatrix Y,
                      Nullable<NumericMatrix> TCOV,
                      NumericVector beta,
                      int nw){
  int period = A.ncol();
  double a = 0;
  for (int i = 0; i < n; ++i){
    double zik = zk(i, k);
    for (int t = 0; t < period; ++t){
      double sikt = Sikt(i, t);
      NumericVector vtmp2;
      for (int po = 0; po < nbeta; ++po){
        vtmp2.push_back(pow(A(i,t), po));
      }
      double betaAit = sum(beta*vtmp2) + WitEM_cpp(TCOV, period, delta, nw, i, t, k); 
      a += zik*(1 - sikt)*(Y(i,t)*betaAit - exp(betaAit));
    }
  }
  return(a);
}

double QbetadeltakZIP_cpp(NumericVector betadelta,
                     NumericMatrix zk, 
                     NumericMatrix Sikt,
                     int k, 
                     int nbeta,
                     int nnu,
                     int n,
                     NumericMatrix A,
                     NumericMatrix Y,
                     NumericMatrix TCOV,
                     int nw){
  int period = A.ncol();
  
  
  NumericVector betak = betadelta[Range(0, nbeta-1)];
  NumericVector deltak = betadelta[Range(nbeta, betadelta.length()-1)];
  double a = 0;
  for (int i = 0; i < n; ++i){
    double zik = zk(i, k);
    for (int t = 0; t < period; ++t){
      double sikt = Sikt(i, t);
      NumericVector vtmp2;
      for (int po = 0; po < nbeta; ++po){
        vtmp2.push_back(pow(A(i,t), po));
      }
      double betaAit = sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, t, k); 
      a += zik*(1- sikt)*(Y(i,t)*betaAit - exp(betaAit));
    }
  }
  return(a);
}
// ----------------------------------------------------------------------------
//  Differential of betak, nuk and deltak Q function
// ----------------------------------------------------------------------------

NumericVector difQbetakZIP_cpp(NumericVector beta,
                               NumericMatrix zk, 
                               NumericMatrix Sikt,
                               int k, 
                               int nbeta,
                               int nnu,
                               int n,
                               NumericMatrix A,
                               NumericMatrix Y,
                               Nullable<NumericMatrix> TCOV,
                               Nullable<NumericVector> delta,
                               int nw,
                               Nullable<IntegerVector> ndeltacum){
  int period = A.ncol();
  NumericVector betas;
  for (int l = 0; l < nbeta; ++l){
    double a = 0;
    for (int i = 0; i < n; ++i){
      double zik = zk(i, k);
      for (int t = 0; t < period; ++t){
        double sikt = Sikt(i, t);
        NumericVector vtmp2;
        for (int po = 0; po < nbeta; ++po){
          vtmp2.push_back(pow(A(i,t), po));
        }
        double betaAit = sum(beta*vtmp2) + WitEM_cpp(TCOV, period, delta, nw, i, t, k); 
        a += pow(A(i,t), l)*zik*(1-sikt)*(Y(i,t)-exp(betaAit));
      }
    }  
    betas.push_back(a);
  }
  return(betas);
}

NumericVector difQnukZIP_cpp(NumericVector nu,
                             NumericMatrix zk, 
                             NumericMatrix Sikt,
                             int k, 
                             int nbeta,
                             int nnu,
                             int n,
                             NumericMatrix A,
                             NumericMatrix Y){
  int period = A.ncol();
  NumericVector nus;
  for (int l = 0; l < nnu; ++l){
    double a = 0;
    for (int i = 0; i < n; ++i){
      double zik = zk(i, k);
      for (int t = 0; t < period; ++t){
        double sikt = Sikt(i, t);
        NumericVector vtmp1;
        for (int po = 0; po < nnu; ++po){
          vtmp1.push_back(pow(A(i,t), po));
        }
        double nuAit = sum(nu*vtmp1);
        a += pow(A(i, t), l)*zik*(sikt-exp(nuAit)/(1+exp(nuAit)));
      }
    }
    nus.push_back(a);
  }
  return(nus);
}

NumericVector difQdeltakZIP_cpp(NumericVector delta,
                                NumericMatrix zk, 
                                NumericMatrix Sikt,
                                int k, 
                                int nbeta,
                                int nnu,
                                int n,
                                NumericMatrix A,
                                NumericMatrix Y,
                                NumericMatrix TCOV,
                                NumericVector beta,
                                int nw){
  int period = A.ncol();
  NumericVector deltas;
  for (int l = 0; l < nw; ++l){
    double a = 0;
    for (int i = 0; i < n; ++i){
      double zik = zk(i, k);
      for (int t = 0; t < period; ++t){
        double sikt = Sikt(i, t);
        NumericVector vtmp2;
        for (int po = 0; po < nbeta; ++po){
          vtmp2.push_back(pow(A(i,t), po));
        }
        double betaAit = sum(beta*vtmp2) + WitEM_cpp(TCOV, period, delta, nw, i, t, k); 
        a +=  TCOV(i, t + l*period)*zik*(1-sikt)*(Y(i,t)-exp(betaAit));
      }
    }  
    deltas.push_back(a);
  }
  return(deltas);
}

NumericVector difQbetadeltakZIP_cpp(NumericVector betadelta,
                               NumericMatrix zk, 
                               NumericMatrix Sikt,
                               int k, 
                               int nbeta,
                               int nnu,
                               int n,
                               NumericMatrix A,
                               NumericMatrix Y,
                               NumericMatrix TCOV,
                               int nw){
  int period = A.ncol();
  NumericVector betak = betadelta[Range(0, nbeta-1)];
  NumericVector deltak = betadelta[Range(nbeta, betadelta.length()-1)];
  
  NumericVector betadeltas;
  for (int l = 0; l < nbeta; ++l){
    double a = 0;
    for (int i = 0; i < n; ++i){
      double zik = zk(i, k);
      for (int t = 0; t < period; ++t){
        double sikt = Sikt(i, t);
        NumericVector vtmp2;
        for (int po = 0; po < nbeta; ++po){
          vtmp2.push_back(pow(A(i,t), po));
        }
        double betaAit = sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, t, k); 
        a += pow(A(i,t), l)*zik*(1-sikt)*(Y(i,t)-exp(betaAit));
      }
    }  
    betadeltas.push_back(a);
  }
  for (int l = 0; l < nw; ++l){
    double a = 0;
    for (int i = 0; i < n; ++i){
      double zik = zk(i, k);
      for (int t = 0; t < period; ++t){
        double sikt = Sikt(i, t);
        NumericVector vtmp2;
        for (int po = 0; po < nbeta; ++po){
          vtmp2.push_back(pow(A(i,t), po));
        }
        double betaAit = sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, t, k); 
        a += TCOV(i, t + l*period)*zik*(1-sikt)*(Y(i,t)-exp(betaAit));
      }
    }  
    betadeltas.push_back(a);
  }
  return(betadeltas);
}
// ----------------------------------------------------------------------------
//  EM 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EMZIP_cpp(NumericVector param,
                        int ng,
                        int nx,
                        int n,
                        IntegerVector nbeta,
                        IntegerVector nnu,
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
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    nu = param[Range(ng+sum(nbeta)-1, ng+sum(nbeta)+sum(nnu)-2)];
    if (param.length() > ng*nx+sum(nbeta)+sum(nnu)-1){
      delta = param[Range(ng+sum(nbeta)+sum(nnu)-1, param.length()-1)];
      NumericVector deltatmp(ng);
      IntegerVector ndeltacumtmp(nw*ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacumtmp.begin());
      ndeltacumtmp.push_front(0);
      ndeltacum = ndeltacumtmp;
    }
    pi.push_back(1-sum(pi));
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
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta),  as<arma::rowvec>(nu), as<arma::rowvec>(delta));
  int tour = 1;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEMZIP_cpp(n, ng, nbeta, nnu, beta, nu, pi, A, Y, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodZIP_cpp(NumericVector(vparam.begin() + nx, vparam.end()), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw));
    }
    // E-step
    NumericMatrix zk = ftauxZIP_cpp(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X);
    
    
    NumericVector newbeta;
    NumericVector newnu;
    NumericVector newdelta;
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];
    
    vec vtmp;
    vec vtmpdel;
    vec vtmpnu; 
    if (nw == 0){
      for (int k = 0; k < ng; ++k){
        NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
        NumericVector nuk = nu[Range(nnucum[k], nnucum[k+1]-1)];
        
        NumericMatrix Sikt(n, period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            Sikt(i, t) = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
          }
        } 
        
        List tmp = optim(Rcpp::_["par"] = betak,
                         Rcpp::_["fn"] = Rcpp::InternalFunction(&QbetakZIP_cpp),
                         Rcpp::_["gr"] = Rcpp::InternalFunction(&difQbetakZIP_cpp),
                         Rcpp::_["method"] = "BFGS",
                         Rcpp::_["zk"] = zk,
                         Rcpp::_["Sikt"] = Sikt,
                         Rcpp::_["k"] = k,
                         Rcpp::_["nbeta"] = nbeta[k],
                         Rcpp::_["nnu"] = nnu[k],
                         Rcpp::_["n"] = n,
                         Rcpp::_["A"] =  A,
                         Rcpp::_["Y"] =  Y,
                         Rcpp::_["TCOV"] = TCOV,
                         Rcpp::_["delta"] = delta,
                         Rcpp::_["nw"] = nw,
                         Rcpp::_["ndeltacum"] = ndeltacum,
                         Rcpp::_["hessian"] =  0,
                         Rcpp::_["control"] = List::create(Named("fnscale")=-1)
        );
        List tmpnu = optim(Rcpp::_["par"] = nuk,
                           Rcpp::_["fn"] = Rcpp::InternalFunction(&QnukZIP_cpp),
                           Rcpp::_["gr"] = Rcpp::InternalFunction(&difQnukZIP_cpp),
                           Rcpp::_["method"] = "BFGS",
                           Rcpp::_["zk"] = zk,
                           Rcpp::_["Sikt"] = Sikt,
                           Rcpp::_["k"] = k,
                           Rcpp::_["nbeta"] = nbeta[k],
                           Rcpp::_["nnu"] = nnu[k],
                           Rcpp::_["n"] = n,
                           Rcpp::_["A"] = A,
                           Rcpp::_["Y"] = Y,
                           Rcpp::_["hessian"] =  0,
                           Rcpp::_["control"] = List::create(Named("fnscale")=-1)
        );
        vtmp = join_cols(vtmp, as<arma::vec>(tmp[0]));
        vtmpnu = join_cols(vtmpnu, as<arma::vec>(tmpnu[0]));
      }
      newbeta = NumericVector(vtmp.begin(), vtmp.end());
      newnu = NumericVector(vtmpnu.begin(), vtmpnu.end());
    }else{
      IntegerVector ndeltacum(ng);
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);

      for (int k = 0; k < ng; ++k){
        NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
        NumericVector nuk = nu[Range(nnucum[k], nnucum[k+1]-1)];
        NumericVector deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];

        NumericMatrix Sikt(n, period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            Sikt(i, t) = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
          }
        } 

        List tmpnu = optim(Rcpp::_["par"] = nuk,
                           Rcpp::_["fn"] = Rcpp::InternalFunction(&QnukZIP_cpp),
                           Rcpp::_["gr"] = Rcpp::InternalFunction(&difQnukZIP_cpp),
                           Rcpp::_["method"] = "BFGS",
                           Rcpp::_["zk"] = zk,
                           Rcpp::_["Sikt"] = Sikt,
                           Rcpp::_["k"] = k,
                           Rcpp::_["nbeta"] = nbeta[k],
                           Rcpp::_["nnu"] = nnu[k],                    
                           Rcpp::_["n"] = n,
                           Rcpp::_["A"] = A,
                           Rcpp::_["Y"] = Y,
                           Rcpp::_["hessian"] =  0,
                           Rcpp::_["control"] = List::create(Named("fnscale")=-1)
        );
        List tmp = optim(Rcpp::_["par"] = betak,
                         Rcpp::_["fn"] = Rcpp::InternalFunction(&QbetakZIP_cpp),
                         Rcpp::_["gr"] = Rcpp::InternalFunction(&difQbetakZIP_cpp),
                         Rcpp::_["method"] = "BFGS",
                         Rcpp::_["zk"] = zk,
                         Rcpp::_["Sikt"] = Sikt,
                         Rcpp::_["k"] = k,
                         Rcpp::_["nbeta"] = nbeta[k],
                         Rcpp::_["nnu"] = nnu[k],
                         Rcpp::_["n"] = n,
                         Rcpp::_["A"] =  A,
                         Rcpp::_["Y"] =  Y,
                         Rcpp::_["TCOV"] = TCOV,
                         Rcpp::_["delta"] = deltak,
                         Rcpp::_["nw"] = nw,
                         Rcpp::_["ndeltacum"] = ndeltacum,
                         Rcpp::_["hessian"] =  0,
                         Rcpp::_["control"] = List::create(Named("fnscale")=-1)
        );
        List tmpdel = optim(Rcpp::_["par"] = deltak,
                         Rcpp::_["fn"] = Rcpp::InternalFunction(&QdeltakZIP_cpp),
                         Rcpp::_["gr"] = Rcpp::InternalFunction(&difQdeltakZIP_cpp),
                         Rcpp::_["method"] = "BFGS",
                         Rcpp::_["zk"] = zk,
                         Rcpp::_["Sikt"] = Sikt,
                         Rcpp::_["k"] = k,
                         Rcpp::_["nbeta"] = nbeta[k],
                         Rcpp::_["nnu"] = nnu[k],
                         Rcpp::_["n"] = n,
                         Rcpp::_["A"] =  A,
                         Rcpp::_["Y"] =  Y,
                         Rcpp::_["TCOV"] = TCOV,
                         Rcpp::_["beta"] = betak,
                         Rcpp::_["nw"] = nw,
                         Rcpp::_["hessian"] =  0,
                         Rcpp::_["control"] = List::create(Named("fnscale")=-1)
        );
        vtmpnu = join_cols(vtmpnu, as<arma::vec>(tmpnu[0]));
        vtmp = join_cols(vtmp, as<arma::vec>(tmp[0]));
        vtmpdel = join_cols(vtmpdel, as<arma::vec>(tmpdel[0]));
      }
      
      newbeta = NumericVector(vtmp.begin(), vtmp.end());
      newnu = NumericVector(vtmpnu.begin(), vtmpnu.end());
      newdelta = NumericVector(vtmpdel.begin(), vtmpdel.end());
    }
    
    // calculus of pi
    NumericVector newpi;
    if (nx == 1){
      NumericVector tmp(ng);
      for (int i = 0; i < ng; ++i){
        tmp[i] = sum(zk(_, i));
      }
      newpi = tmp/n;
    }else{
      newpi = findtheta_cpp(pi, zk, X, n, ng, nx, period, EMIRLS, refgr);
    }
    // stop test
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), as<arma::rowvec>(newbeta), as<arma::rowvec>(newnu), as<arma::rowvec>(newdelta));
    rowvec tmp1(newparam.size());
    tmp1.fill(prec);
    if (all(abs(newparam-vparam)<tmp1)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    nu = newnu;
    delta = newdelta;
    pi = newpi;
    
  }
  return(NumericVector(vparam.begin(), vparam.end()));
}

// ----------------------------------------------------------------------------
//  EM IRLS
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EMZIPIRLS_cpp(NumericVector param,
                        int ng,
                        int nx,
                        int n,
                        IntegerVector nbeta,
                        IntegerVector nnu,
                        NumericMatrix A,
                        NumericMatrix Y,
                        NumericMatrix X,
                        Nullable<NumericMatrix> TCOVinit,
                        int nw,
                        int itermax,
                        bool EMIRLS,
                        int refgr){

  // int period = A.ncol();
   NumericMatrix TCOV;
  // double prec = 0.000001;
  // NumericVector pi(ng);
  // NumericVector beta;
  // NumericVector nu;
  // NumericVector delta;
  // 
  // IntegerVector nbetacum(nbeta.size());
  // std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  // nbetacum.push_front(0);
  // IntegerVector nnucum(nnu.size());
  // std::partial_sum(nnu.begin(), nnu.end(), nnucum.begin());
  // nnucum.push_front(0);
  // IntegerVector ndeltacum(nw*ng);
  // 
  // if (nx == 1){
  //   pi = param[Range(0,ng-2)];
  //   beta = param[Range(ng-1,ng+sum(nbeta)-2)];
  //   nu = param[Range(ng+sum(nbeta)-1, ng+sum(nbeta)+sum(nnu)-2)];
  //   if (param.length() > ng*nx+sum(nbeta)+sum(nnu)-1){
  //     delta = param[Range(ng+sum(nbeta)+sum(nnu)- 1, param.length()-1)];
  //   }
  //   pi.push_back(1-sum(pi));
  // }else{
  //   pi = param[Range(0,ng*nx-1)];
  //   beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
  //   nu = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+sum(nnu)-1)];
  //   if (param.length() > ng*nx+sum(nbeta)+sum(nnu)){
  //     delta = param[Range(ng*nx+sum(nbeta)+sum(nnu), param.length()-1)];
  //   }
  // }
  // 
  // 
  // rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta),  as<arma::rowvec>(nu), as<arma::rowvec>(delta));
  // int tour = 1;
  int period = A.ncol();
  double prec = 0.000001;
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
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    nu = param[Range(ng+sum(nbeta)-1, ng+sum(nbeta)+sum(nnu)-2)];
    if (param.length() > ng*nx+sum(nbeta)+sum(nnu)-1){
      delta = param[Range(ng+sum(nbeta)+sum(nnu)-1, param.length()-1)];
      NumericVector deltatmp(ng);
      IntegerVector ndeltacumtmp(nw*ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacumtmp.begin());
      ndeltacumtmp.push_front(0);
      ndeltacum = ndeltacumtmp;
    }
    pi.push_back(1-sum(pi));
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
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta),  as<arma::rowvec>(nu), as<arma::rowvec>(delta));
  int tour = 1;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEMZIP_cpp(n, ng, nbeta, nnu, beta, nu, pi, A, Y, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodZIP_cpp(NumericVector(vparam.begin() + nx, vparam.end()), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw));
    }
    // E-step
    NumericMatrix zk = ftauxZIP_cpp(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X);
    rowvec  newbeta;
    rowvec  newnu;
    rowvec  newdelta;

    if (nw == 0){
      for (int k = 0; k < ng; ++k){
        vec newbetaIRLS;
        NumericVector betaIRLS = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
        vec newnuIRLS;
        NumericVector nuIRLS = nu[Range(nnucum[k], nnucum[k+1]-1)];
        arma::vec precIRLS(betaIRLS.length());
        precIRLS.fill(1);
        int stop = 0;
        
        while(all(abs(precIRLS) > 0.000001) &&  stop < 300){
          stop +=1;
          mat Aw(nnu[k], n*period);
          mat Awp(nbeta[k], n*period);
          vec Sn(n*period);
          vec Sw(n*period);
          vec Sp(n*period);
          vec W(n*period);
          vec Wp(n*period);
          vec Z(n*period);
          
          for (int i = 0; i < n; ++i){
            for (int t = 0; t < period; ++t){
              
              double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
              double nuAit = 0;
              for (int po = 0; po < nnu[k]; ++po){
                nuAit += pow(A(i, t), po)*nuIRLS[po];
                Aw(po, period*i + t) = pow(A(i, t), po);
              }
              
              double rhoikt = exp(nuAit)/(1+exp(nuAit));
              W[period*i + t] = rhoikt*(1-rhoikt);
              Sn[period*i + t] = nuAit + (sikt - rhoikt)/(rhoikt*(1-rhoikt));
              
              double betaAit = 0;
              for (int po = 0; po < nbeta[k]; ++po){
                betaAit += pow(A(i, t), po)*betaIRLS[po];
                Awp(po, period*i + t) = pow(A(i, t), po);
              }
              
              double lambdaikt = exp(betaAit);
              Sp[period*i + t] = 1-sikt;
              Sw[period*i + t] = betaAit+Y(i,t)/lambdaikt-1;
              Wp[period*i + t] = lambdaikt;
              
              Z[period*i + t] = zk(i, k);
            }
          }
          
          mat Q;
          mat R;
          mat ZW = diagmat(Z % W);
          qr(Q, R, trans(Aw*sqrt(ZW)));
          mat Rc = R.submat(0, 0, nnu[k]-1, nnu[k]-1);
          mat Qc = Q.submat(0, 0, Q.n_rows-1, nnu[k]-1);
          newnuIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*Sn);
       
          mat SpZWp = diagmat(Sp % Z % Wp);
          qr(Q, R, trans(Awp*sqrt(SpZWp)));
          Rc = R.submat(0, 0, nbeta[k]-1, nbeta[k]-1);
          Qc = Q.submat(0, 0, Q.n_rows-1, nbeta[k]-1);
          newbetaIRLS = solve(Rc, trans(Qc)*sqrt(SpZWp)*Sw);
        
          precIRLS = join_cols(Rcpp::as<arma::vec>(betaIRLS),Rcpp::as<arma::vec>(nuIRLS))-join_cols(newbetaIRLS, newnuIRLS);
          betaIRLS = Rcpp::NumericVector(newbetaIRLS.begin(), newbetaIRLS.end());
          nuIRLS = Rcpp::NumericVector(newnuIRLS.begin(), newnuIRLS.end());
        }
        newbeta = join_rows(newbeta, as<arma::rowvec>(betaIRLS));
        newnu = join_rows(newnu, as<arma::rowvec>(nuIRLS));
      }
    }else{
      IntegerVector ndeltacum(ng);
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
      NumericMatrix tmp1(TCOVinit.get());  
      TCOV = tmp1;
      
      for (int k = 0; k < ng; ++k){
        
        vec newbetaIRLS;
        NumericVector betaIRLS = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
        vec newnuIRLS;
        NumericVector nuIRLS = nu[Range(nnucum[k], nnucum[k+1]-1)];
        vec newdeltaIRLS;
        NumericVector deltaIRLS = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
        NumericVector deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
        arma::vec precIRLS(betaIRLS.length());
        precIRLS.fill(1);
        int stop = 0;
        
        while(all(abs(precIRLS) > 0.000001) &&  stop < 300){
          stop +=1;
          mat Aw(nnu[k], n*period);
          mat Awp(nbeta[k], n*period);
          mat Ww(nw, n*period);
          vec Sn(n*period);
          vec Sb(n*period);
          vec Sp(n*period);
          vec Sw(n*period);
          vec W(n*period);
          vec Wp(n*period);
          vec Z(n*period);
          
          for (int i = 0; i < n; ++i){
            for (int t = 0; t < period; ++t){
              
              double sikt = fSikt_cpp(pi, beta, nu, k, i, t, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum);
              double nuAit = 0;
              for (int po = 0; po < nnu[k]; ++po){
                nuAit += pow(A(i, t), po)*nuIRLS[po];
                Aw(po, period*i + t) = pow(A(i, t), po);
              }
              for (int kk = 0; kk < nw; ++kk){
                Ww(kk, period*i + t) = TCOV(i, t + kk*period);
              }
              
              double rhoikt = exp(nuAit)/(1+exp(nuAit));
              W[period*i + t] = rhoikt*(1-rhoikt);
              Sn[period*i + t] = nuAit + (sikt - rhoikt)/(rhoikt*(1-rhoikt));
              
              double betaAit = 0;
              for (int po = 0; po < nbeta[k]; ++po){
                betaAit += pow(A(i, t), po)*betaIRLS[po];
                Awp(po, period*i + t) = pow(A(i, t), po);
              }
              
              double lambdaikt = exp(betaAit +  WitEM_cpp(TCOV, period, deltak, nw, i, t, k));
              Sp[period*i + t] = 1-sikt;
              Sb[period*i + t] = betaAit+Y(i,t)/lambdaikt-1;
              Wp[period*i + t] = lambdaikt;
              
              Sw[period*i + t] = WitEM_cpp(TCOV, period, deltak, nw, i, t, k)+Y(i,t)/lambdaikt-1;
              
              Z[period*i + t] = zk(i, k);
            }
          }
          
          mat Q;
          mat R;
          mat ZW = diagmat(Z % W);
          qr(Q, R, trans(Aw*sqrt(ZW)));
          mat Rc = R.submat(0, 0, nnu[k]-1, nnu[k]-1);
          mat Qc = Q.submat(0, 0, ZW.n_rows-1, nnu[k]-1);
          newnuIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*Sn);
          
          mat SpZWp = diagmat(Sp % Z % Wp);
          qr(Q, R, trans(Awp*sqrt(SpZWp)));
          Rc = R.submat(0, 0, nbeta[k]-1, nbeta[k]-1);
          Qc = Q.submat(0, 0, SpZWp.n_rows-1, nbeta[k]-1);
          newbetaIRLS = solve(Rc, trans(Qc)*sqrt(SpZWp)*Sb);
         
          mat ZsWp = diagmat(Z % Sp % Wp);
          qr(Q, R, trans(Ww*sqrt(ZsWp)));
          
          Rc = R.submat(0, 0, nw-1, nw-1);
          Qc = Q.submat(0, 0, ZsWp.n_rows-1, nw-1);
          
          newdeltaIRLS = solve(R, trans(Q)*(sqrt(ZsWp)*Sw));
          
          precIRLS = join_cols(Rcpp::as<arma::vec>(betaIRLS),Rcpp::as<arma::vec>(nuIRLS), Rcpp::as<arma::vec>(deltaIRLS))-join_cols(newbetaIRLS, newnuIRLS, newdeltaIRLS);
          
          betaIRLS = Rcpp::NumericVector(newbetaIRLS.begin(), newbetaIRLS.end());
          nuIRLS = Rcpp::NumericVector(newnuIRLS.begin(), newnuIRLS.end());
          deltaIRLS = Rcpp::NumericVector(newdeltaIRLS.begin(), newdeltaIRLS.end());
        }
        newbeta = join_rows(newbeta, as<arma::rowvec>(betaIRLS));
        newnu = join_rows(newnu, as<arma::rowvec>(nuIRLS));
        newdelta = join_rows(newdelta, as<arma::rowvec>(deltaIRLS));
      }
    }

    // calculus of pi
    NumericVector newpi;
    if (nx == 1){
      NumericVector tmp(ng);
      for (int i = 0; i < ng; ++i){
        tmp[i] = sum(zk(_, i));
      }
      newpi = tmp/n;
    }else{
      newpi = findtheta_cpp(pi, zk, X, n, ng, nx, period, EMIRLS, refgr);
    }
    // stop test
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), newbeta, newnu, newdelta);
    rowvec tmp1(newparam.size());
    tmp1.fill(prec);
    if (all(abs(newparam-vparam)<tmp1)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    nu = newnu;
    delta = newdelta;
    pi = newpi;
    
  }
 return(NumericVector(vparam.begin(), vparam.end()));
}
