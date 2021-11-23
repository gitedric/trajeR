#include "CommonFunction.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::interfaces(r, cpp)]]


// ----------------------------------------------------------------------------
// gk 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double gkCNORM_cpp(List beta,
                   NumericVector sigma,
                   int i,
                   int k,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   double ymin,
                   double ymax, 
                   Nullable<NumericMatrix> TCOV,
                   Nullable<List> delta,
                   int nw){
  int period = A.ncol();
  k = k-1;
  i = i-1;
  NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
  double res = 1;
  for (int ind = 0; ind < period; ++ind){
    if (R_IsNA(Y(i, ind)) == FALSE){
      if (Y(i, ind) <= ymin){
        res = res*(R::pnorm((Y(i, ind)-muikt[ind])/sigma[k], 0.0, 1.0, true, false));
      }else if (Y(i, ind) >= ymax){
        res = res*(R::pnorm(-(Y(i, ind)-muikt[ind])/sigma[k], 0.0, 1.0, true, false));
      }else{
        res = res*(R::dnorm((Y(i, ind)-muikt[ind])/sigma[k], 0.0, 1.0, false)/sigma[k]);
      }
    }
  }
  return(res);
}
// ----------------------------------------------------------------------------
// gk with exponential parametrization
// ----------------------------------------------------------------------------
double gkalpha_cpp(List beta,
                   NumericVector alpha,
                   int i,
                   int k,
                   IntegerVector nbeta,
                   NumericMatrix A,
                   NumericMatrix Y,
                   double ymin,
                   double ymax, 
                   Nullable<NumericMatrix> TCOV,
                   Nullable<List> delta,
                   int nw){
  int period = A.ncol();
  NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
  NumericVector sigma = exp(alpha);
  double res = 1;
  for (int t = 0; t < period; ++t){
    if (R_IsNA(Y(i, t)) == FALSE){
      if (Y(i, t) <= ymin){
        res *= R::pnorm((Y(i, t)-muikt[t])/sigma[k], 0.0, 1.0, true, false);
      }else if (Y(i, t) >= ymax){
        res *= R::pnorm(-(Y(i, t)-muikt[t])/sigma[k], 0.0, 1.0, true, false);
      }else{
        res *= R::dnorm((Y(i, t)-muikt[t])/sigma[k], 0.0, 1.0, false)/sigma[k];
      }
    }
  }
  return(res);
}
// ----------------------------------------------------------------------------
// dif likelihood theta
// ----------------------------------------------------------------------------
NumericVector difLthetakalpha_cpp(NumericVector theta,
                                  List beta,
                                  NumericVector alpha,
                                  Nullable<List> delta,
                                  int k,
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
        tmp1 += vtmp[s]*(gkalpha_cpp(beta, alpha, i, k, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)-gkalpha_cpp(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw));
      }
      tmp1 = tmp1*vtmp[k];
      double tmp2 = 0;
      for (int s = 0; s < ng; ++s){
        tmp2 += vtmp[s]*gkalpha_cpp(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw);
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
NumericVector difLbetakalpha_cpp(NumericVector theta,
                                 List beta,
                                 NumericVector alpha,
                                 Nullable<List> delta,
                                 int k,
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
                                 int nw){
  NumericVector betas;
  int period = A.ncol();
  
  for (int l = 0; l < nbeta[k]; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      double difbkl = 0;
      double difbklmin = 0;
      double difbklmax = 0;
      double py = 1;
      double pymax = 1;
      double pymin = 1;
      // initialize mean
      NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      // create a matrix with Yit, muikt ans Ait
      NumericMatrix YimuiA(3, period);
      YimuiA(0, _) = Y(i, _);
      YimuiA(1, _) = muikt;
      YimuiA(2, _) = A(i, _);
      // check if there exist some censored values
      LogicalVector indmin = (Y(i, _) <= rep(ymin, period));
      LogicalVector indmax = (Y(i, _) >= rep(ymax, period));
      LogicalVector ind = !(indmin | indmax);
      for (int t = 0; t < period; t++){
        if (R_IsNA(Y(i, t))){
          indmin[t] = FALSE;
          indmax[t] = FALSE;
          ind[t] = FALSE;
          }
        }
      if (is_true(any(indmin == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, indmin); 
        NumericVector vytmp = rep(ymin, sum(indmin));
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(indmin) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbklmin -= pow(vAtmp[t], l)*exp(-alpha[k])*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(vmeantmpminus, 0.0, 1.0, true, false));   
          } 
        }else{
          difbklmin = -pow(vAtmp[0], l)*exp(-alpha[k])*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        pymin = prodvect(Rcpp::pnorm(vmean, 0.0, 1.0, true, false));
      }
      if (is_true(any(indmax == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, indmax); 
        NumericVector vytmp = rep(ymax, sum(indmax));
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(indmax) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbklmax += pow(vAtmp[t], l)*exp(-alpha[k])*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(-vmeantmpminus, 0.0, 1.0, true, false));   
          } 
        }else{
          difbklmax = pow(vAtmp[0], l)*exp(-alpha[k])*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        pymax = prodvect(Rcpp::pnorm(-vmean, 0.0, 1.0, true, false));
      }
      if (is_true(any(ind == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind); 
        NumericVector vytmp = mtmp(0, _);
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(ind) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbkl += pow(vAtmp[t], l)*exp(-2*alpha[k])*vmean[t]*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmeantmpminus, 0.0, 1.0, false)));   
          } 
        }else{
          difbkl = pow(vAtmp[0], l)*exp(-2*alpha[k])*vmean[0]*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        py = prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmean, 0.0, 1.0, false)));
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkalpha_cpp(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax);
    }
    betas.push_back(a);
  }
  return(betas); 
}
// ----------------------------------------------------------------------------
// dif likelihood deltak
// ----------------------------------------------------------------------------
NumericVector difLdeltakalpha_cpp(NumericVector theta,
                                 List beta,
                                 NumericVector alpha,
                                 Nullable<List> delta,
                                 int k,
                                 int ng,
                                 int nx,
                                 IntegerVector nbeta,
                                 int n,
                                 NumericMatrix A,
                                 NumericMatrix Y,
                                 NumericMatrix X,
                                 double ymin,
                                 double ymax, 
                                 Nullable<NumericMatrix> initTCOV,
                                 int nw){
  NumericVector deltas;
  int period = A.ncol();
  NumericMatrix TCOV(initTCOV);
  
  for (int l = 0; l < nw; ++l){
    double a = 0;
    for (int i = 0; i < n; i++){
      double difbkl = 0;
      double difbklmin = 0;
      double difbklmax = 0;
      double py = 1;
      double pymax = 1;
      double pymin = 1;
      // initialize mean
      NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      // create a matrix with Yit, muikt ans Ait
      NumericMatrix YimuiA(3, period);
      YimuiA(0, _) = Y(i, _);
      YimuiA(1, _) = muikt;
      for (int t = 0; t < period; t++){
        YimuiA(2, t) = TCOV(i, t + l*period);   
      }
      // check if there exist some censored values
      LogicalVector indmin = (Y(i, _) <= rep(ymin, period));
      LogicalVector indmax = (Y(i, _) >= rep(ymax, period));
      LogicalVector ind = !(indmin | indmax);
      for (int t = 0; t < period; t++){
        if (R_IsNA(Y(i, t))){
          indmin[t] = FALSE;
          indmax[t] = FALSE;
          ind[t] = FALSE;
        }
      }
      if (is_true(any(indmin == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, indmin); 
        NumericVector vytmp = rep(ymin, sum(indmin));
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(indmin) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbklmin -= vAtmp[t]*exp(-alpha[k])*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(vmeantmpminus, 0.0, 1.0, true, false));   
          } 
        }else{
          difbklmin = -vAtmp[0]*exp(-alpha[k])*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        pymin = prodvect(Rcpp::pnorm(vmean, 0.0, 1.0, true, false));
      }
      if (is_true(any(indmax == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, indmax); 
        NumericVector vytmp = rep(ymax, sum(indmax));
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(indmax) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbklmax += vAtmp[t]*exp(-alpha[k])*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(-vmeantmpminus, 0.0, 1.0, true, false));   
          } 
        }else{
          difbklmax = vAtmp[0]*exp(-alpha[k])*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        pymax = prodvect(Rcpp::pnorm(-vmean, 0.0, 1.0, true, false));
      }
      if (is_true(any(ind == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind); 
        NumericVector vytmp = mtmp(0, _);
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(ind) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbkl += vAtmp[t]*exp(-2*alpha[k])*vmean[t]*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmeantmpminus, 0.0, 1.0, false)));   
          } 
        }else{
          difbkl = vAtmp[0]*exp(-2*alpha[k])*vmean[0]*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        py = prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmean, 0.0, 1.0, false)));
      }
      double so = 0;
      for (int s = 0; s < ng; ++s){
        so += piikIntern_cpp(theta, i, s, ng, X)*gkalpha_cpp(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw);
      }
      a += piikIntern_cpp(theta, i, k, ng, X)/so*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax);
    }
    deltas.push_back(a);
  }
  return(deltas); 
}
// ----------------------------------------------------------------------------
// dif likelihood sigma with reparametrization alpha
// ----------------------------------------------------------------------------
double difLsigmakalpha_cpp(NumericVector theta,
                                  List beta,
                                  NumericVector alpha,
                                  Nullable<List> delta,
                                  int k,
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
                                  int nw){
  int period = A.ncol();
  double alphas = 0;
  for (int i = 0; i < n; i++){
    double difbkl = 0;
    double difbklmin = 0;
    double difbklmax = 0;
    double py = 1;
    double pymax = 1;
    double pymin = 1;
    // initialize mean
    NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
    // create a matrix with Yit, muikt ans Ait
    NumericMatrix YimuiA(3, period);
    YimuiA(0, _) = Y(i, _);
    YimuiA(1, _) = muikt;
    YimuiA(2, _) = A(i, _);
    // check if there exist some censored values
    LogicalVector indmin = (Y(i, _) <= rep(ymin, period));
    LogicalVector indmax = (Y(i, _) >= rep(ymax, period));
    LogicalVector ind = !(indmin | indmax);
    for (int t = 0; t < period; t++){
      if (R_IsNA(Y(i, t))){
        indmin[t] = FALSE;
        indmax[t] = FALSE;
        ind[t] = FALSE;
      }
    }
    if (is_true(any(indmin == TRUE))){
      // construct a matrix with element of Yimui for which ind is true
      NumericMatrix mtmp = submat_cpp(YimuiA, indmin); 
      NumericVector vytmp = rep(ymin, sum(indmin));
      NumericVector vmutmp = mtmp(1, _);
      NumericVector vAtmp = mtmp(2, _);
      NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
      if (sum(indmin) > 1){
        for (int t = 0; t < vytmp.size(); ++t){
          // vector with all elements except the element t
          NumericVector vmeantmpminus = vmean;
          vmeantmpminus.erase(vmeantmpminus.begin() + t);
          difbklmin -= vmean[t]*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(vmeantmpminus, 0.0, 1.0, true, false));   
        } 
      }else{
        difbklmin = -vmean[0]*(R::dnorm(vmean[0], 0.0, 1.0, false));   
      }
      pymin = prodvect(Rcpp::pnorm(vmean, 0.0, 1.0, true, false));
    }
    if (is_true(any(indmax == TRUE))){
      // construct a matrix with element of Yimui for which ind is true
      NumericMatrix mtmp = submat_cpp(YimuiA, indmax); 
      NumericVector vytmp = rep(ymax, sum(indmax));
      NumericVector vmutmp = mtmp(1, _);
      NumericVector vAtmp = mtmp(2, _);
      NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
      if (sum(indmax) > 1){
        for (int t = 0; t < vytmp.size(); ++t){
          // vector with all elements except the element t
          NumericVector vmeantmpminus = vmean;
          vmeantmpminus.erase(vmeantmpminus.begin() + t);
          difbklmax += vmean[t]*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(-vmeantmpminus, 0.0, 1.0, true, false));   
        } 
      }else{
        difbklmax = vmean[0]*(R::dnorm(vmean[0], 0.0, 1.0, false));   
      }
      pymax = prodvect(Rcpp::pnorm(-vmean, 0.0, 1.0, true, false));
    }
    if (is_true(any(ind == TRUE))){
      // construct a matrix with element of Yimui for which ind is true
      NumericMatrix mtmp = submat_cpp(YimuiA, ind); 
      NumericVector vytmp = mtmp(0, _);
      NumericVector vmutmp = mtmp(1, _);
      NumericVector vAtmp = mtmp(2, _);
      NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
      if (sum(ind) > 1){
        for (int t = 0; t < vytmp.size(); ++t){
          // vector with all elements except the element t
          NumericVector vmeantmpminus = vmean;
          vmeantmpminus.erase(vmeantmpminus.begin() + t);
          difbkl += exp(-alpha[k])*(-1+pow(vmean[t], 2))*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmeantmpminus, 0.0, 1.0, false)));   
        } 
      }else{
        difbkl = exp(-alpha[k])*(-1+pow(vmean[0], 2))*(R::dnorm(vmean[0], 0.0, 1.0, false));   
      }
      py = prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmean, 0.0, 1.0, false)));
    }
    double so = 0;
    for (int s = 0; s < ng; ++s){
      so += piikIntern_cpp(theta, i, s, ng, X)*gkalpha_cpp(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw);
    }
    alphas += piikIntern_cpp(theta, i, k, ng, X)/so*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax);
  }
  return(alphas); 
}
// ----------------------------------------------------------------------------
// dif likelihood sigma with the same value for each element
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double difLsigmaalphaunique_cpp(NumericVector theta,
                                       List beta,
                                       NumericVector alpha,
                                       Nullable<List> delta,
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
                                       int nw){
  int period = A.ncol();
  double a = 0;
  for (int i = 0; i < n; i++){
    double s1 = 0;
    double s2 = 0;
    for (int k = 0; k < ng; ++k){
      double difbkl = 0;
      double difbklmin = 0;
      double difbklmax = 0;
      double py = 1;
      double pymax = 1;
      double pymin = 1;
      // initialize mean
      NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
      // create a matrix with Yit, muikt ans Ait
      NumericMatrix YimuiA(3, period);
      YimuiA(0, _) = Y(i, _);
      YimuiA(1, _) = muikt;
      YimuiA(2, _) = A(i, _);
      // check if there exist some censored values
      LogicalVector indmin = (Y(i, _) <= rep(ymin, period));
      LogicalVector indmax = (Y(i, _) >= rep(ymax, period));
      LogicalVector ind = !(indmin | indmax);
      for (int t = 0; t < period; t++){
        if (R_IsNA(Y(i, t))){
          indmin[t] = FALSE;
          indmax[t] = FALSE;
          ind[t] = FALSE;
        }
      }
      if (is_true(any(indmin == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, indmin); 
        NumericVector vytmp = rep(ymin, sum(indmin));
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(indmin) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbklmin -= vmean[t]*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(vmeantmpminus, 0.0, 1.0, true, false));   
          } 
        }else{
          difbklmin = -vmean[0]*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        pymin = prodvect(Rcpp::pnorm(vmean, 0.0, 1.0, true, false));
      }
      if (is_true(any(indmax == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, indmax); 
        NumericVector vytmp = rep(ymax, sum(indmax));
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(indmax) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbklmax += vmean[t]*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(Rcpp::pnorm(-vmeantmpminus, 0.0, 1.0, true, false));   
          } 
        }else{
          difbklmax = vmean[0]*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        pymax = prodvect(Rcpp::pnorm(-vmean, 0.0, 1.0, true, false));
      }
      if (is_true(any(ind == TRUE))){
        // construct a matrix with element of Yimui for which ind is true
        NumericMatrix mtmp = submat_cpp(YimuiA, ind); 
        NumericVector vytmp = mtmp(0, _);
        NumericVector vmutmp = mtmp(1, _);
        NumericVector vAtmp = mtmp(2, _);
        NumericVector vmean = (vytmp - vmutmp)*exp(-alpha[k]);
        if (sum(ind) > 1){
          for (int t = 0; t < vytmp.size(); ++t){
            // vector with all elements except the element t
            NumericVector vmeantmpminus = vmean;
            vmeantmpminus.erase(vmeantmpminus.begin() + t);
            difbkl += exp(-alpha[k])*(-1+pow(vmean[t], 2))*(R::dnorm(vmean[t], 0.0, 1.0, false))*prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmeantmpminus, 0.0, 1.0, false)));   
          } 
        }else{
          difbkl = exp(-alpha[k])*(-1+pow(vmean[0], 2))*(R::dnorm(vmean[0], 0.0, 1.0, false));   
        }
        py = prodvect(exp(-alpha[k])*(Rcpp::dnorm(vmean, 0.0, 1.0, false)));
      }
      for (int s = 0; s < ng; ++s){
        s1 += piikIntern_cpp(theta, i, s, ng, X)*gkalpha_cpp(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw);
      }
      s2 += piikIntern_cpp(theta, i, k, ng, X)*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax);
    }
    a += s2/s1;
  }
  
  return(a); 
}
// ----------------------------------------------------------------------------
// dif likelihood with exponential reparmetrization
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector difLalpha_cpp(NumericVector param,
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
                                  bool ssigma){
  NumericVector out;
  NumericVector theta = param[Range(0,(ng - 1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  NumericVector beta = param[Range((ng - 1)*nx, (ng - 1)*nx+sum(nbeta)-1)];
  NumericVector alpha = param[Range((ng - 1)*nx+sum(nbeta), (ng - 1)*nx+sum(nbeta)+ng-1)];
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
  if (param.length() > (ng - 1)*nx+sum(nbeta)+ng){
    delta = param[Range((ng - 1)*nx+sum(nbeta)+ng, param.length())];
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
    NumericVector tmp =  difLthetakalpha_cpp(theta, betaL, alpha, deltaL, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLbetakalpha_cpp(theta, betaL, alpha, deltaL, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    double tmp =  difLsigmakalpha_cpp(theta, betaL, alpha, deltaL, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
    out.push_back(tmp);
  }
  if (nw != 0){
    for (int k = 0; k < ng; ++k){
      NumericVector tmp =  difLdeltakalpha_cpp(theta, betaL, alpha, deltaL, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
      for (int i = 0; i < tmp.length(); ++i){
        out.push_back(tmp[i]);
      }
    }
  }
  return(out);
}
// ----------------------------------------------------------------------------
// dif likelihood with same sigma and exponential reparmetrization
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector difLalphaunique_cpp(NumericVector param,
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
                                 bool ssigma){
  NumericVector out;
  NumericVector theta = param[Range(0,(ng - 1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  NumericVector beta = param[Range((ng - 1)*nx, (ng - 1)*nx+sum(nbeta)-1)];
  int indtmp;
  indtmp = (ng - 1)*nx+sum(nbeta) + 1;
  NumericVector alpha = param[Range((ng - 1)*nx+sum(nbeta), indtmp - 1)];
  NumericVector alpharep = rep(alpha, ng);
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
  if (param.length() > indtmp){
    delta = param[Range(indtmp, param.length())];
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

  for (int k = 1; k < ng; ++k){
    NumericVector tmp =  difLthetakalpha_cpp(theta, betaL, alpharep, deltaL, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
    for (int i = 0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  for (int k = 0; k < ng; ++k){
    NumericVector tmp =  difLbetakalpha_cpp(theta, betaL, alpharep, deltaL, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
    for (int i =0; i < tmp.length(); ++i){
      out.push_back(tmp[i]);
    }
  }
  double tmp =  difLsigmaalphaunique_cpp(theta, betaL, alpharep, deltaL, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
  out.push_back(tmp);
  if (nw != 0){
    for (int k = 0; k < ng; ++k){
      NumericVector tmp =  difLdeltakalpha_cpp(theta, betaL, alpharep, deltaL, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw);
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
double Likelihoodalpha_cpp(NumericVector param,
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
                           bool ssigma){
  NumericVector theta = param[Range(0,(ng - 1)*nx-1)];
  for (int i = 0; i < nx; i++){
    theta.push_front(0);
  }
  NumericVector beta = param[Range((ng - 1)*nx, (ng - 1)*nx+sum(nbeta)-1)];
  int indtmp;
  if (ssigma){
    indtmp = (ng - 1)*nx+sum(nbeta) + 1;
  }else{
    indtmp = (ng - 1)*nx+sum(nbeta) + ng;
    }
  NumericVector alpha = param[Range((ng - 1)*nx+sum(nbeta), indtmp - 1)];
  NumericVector alpharep;
  if (ssigma){
    alpharep = rep(alpha, ng);
  }else{
    alpharep = alpha;
  }
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
  if (param.length() > indtmp){
    delta = param[Range(indtmp, param.length())];
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
      a +=  piikIntern_cpp(theta, i, s, ng, X)*gkalpha_cpp(betaL, alpharep, i, s, nbeta, A, Y, ymin, ymax, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}
// ----------------------------------------------------------------------------
// Likelihood 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double likelihoodCNORM_cpp(NumericVector param,
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
                           int nw){
  NumericVector theta = param[Range(0, ng*nx-1)];
  // for (int i = 0; i < nx; i++){
  //   theta.push_front(0);
  // }
  NumericVector beta = param[Range(ng*nx, ng*nx+sum(nbeta)-1)];
  NumericVector sigma = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+ng-1)];
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
  if (param.length() > ng*nx+sum(nbeta)+ng){
    delta = param[Range(ng*nx+sum(nbeta)+ng, param.length())];
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
      a +=  piikIntern_cpp(theta, i, s, ng, X)*gkCNORM_cpp(betaL, sigma, i+1, s+1, nbeta, A, Y, ymin, ymax, TCOV, deltaL, nw);
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
double likelihoodEM_cpp(int n,
                        int ng,
                        IntegerVector nbeta,
                        NumericVector beta,
                        NumericVector sigma,
                        NumericVector pi,
                        NumericMatrix A,
                        NumericMatrix Y,
                        double ymin,
                        double ymax, 
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
      a +=  pi[s]*gkCNORM_cpp(betaL, sigma, i+1, s+1, nbeta, A, Y, ymin, ymax, TCOV, deltaL, nw);
    }
    out += log(a);
  }
  return(out);
}
// ----------------------------------------------------------------------------
// Function rate
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix ftauxCNORM_cpp(NumericVector pi,
                        NumericVector beta,
                        NumericVector sigma,
                        int ng,
                        IntegerVector nbeta,
                        int n,
                        NumericMatrix A,
                        NumericMatrix Y,
                        double ymin,
                        double ymax, 
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
        mtmp(i,k) = pi[k]*gkCNORM_cpp(betaL, sigma, i+1, k+1, nbeta, A, Y, ymin, ymax, TCOV, deltaL, nw);
        s += mtmp(i,k);
      } 
      mtmp(i, _) = 1/(1+ (s - mtmp(i, _))/mtmp(i, _));
    }
  }else{
      for (int i = 0; i < n; ++i){
        double s = 0;
        for (int k = 0; k < ng; ++k){
          mtmp(i,k) = piikIntern_cpp(pi, i, k, ng, X)*gkCNORM_cpp(betaL, sigma, i+1, k+1, nbeta, A, Y, ymin, ymax, TCOV, deltaL, nw);
          s += mtmp(i,k);
        } 
        mtmp(i, _) = 1/(1+ (s - mtmp(i, _))/mtmp(i, _));
      }
  }
  return(mtmp);
}
// ----------------------------------------------------------------------------
// EM not censored
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EM_cpp(NumericVector param,
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
                                  int itermax, 
                                  bool EMIRLS,
                                  int refgr){
  int period = A.ncol();
  double prec = 0.000001;
  NumericVector pi(ng);
  NumericVector beta;
  NumericVector sigma;
  NumericVector delta;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    sigma = param[Range(ng+sum(nbeta)-1, ng+sum(nbeta)+ng-2)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng+sum(nbeta)+ng - 1, param.length() - 1)];
    }
    pi.push_back(1-sum(pi));  
  }else{
    pi = param[Range(0,ng*nx-1)];
    beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
    sigma = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+ng-1)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng*nx+sum(nbeta)+ng, param.length() - 1)];
    }
  }
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta), as<arma::rowvec>(sigma), as<arma::rowvec>(delta));
  int tour = 1;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEM_cpp(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodCNORM_cpp(NumericVector(vparam.begin(), vparam.end()), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw));
    }
    // E-step
    NumericMatrix  taux = ftauxCNORM_cpp(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X);
    rowvec newbeta;
    rowvec newdelta;
    rowvec newsigma(ng);
    if (nw == 0){
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        double b = 0;
        mat Ai(nbeta[k], period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
          }
          rowvec vtmp = Y(i, _);
          a += taux(i, k)*vtmp*trans(Ai);
          rowvec mtmp = vtmp - trans((Rcpp::as<arma::vec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1))*Ai;
          b += taux(i, k)*as_scalar(mtmp*trans(mtmp));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
        newsigma[k]= sqrt(b/(period*sum(taux(_, k))));
      }
    }else{
      // TCOV is not NULL
      NumericVector ndeltacum(ng);
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
      mat mTCOV = as<arma::mat>(TCOV);
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        rowvec c(nw);
        c.fill(0);
        mat Sw(nw, nw);
        Sw.fill(0);
        double b = 0;
        mat Ai(nbeta[k], period);
        mat Wi(nw, period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
            for (int kk = 0; kk < nw; ++kk){
              Wi(kk, t) = mTCOV(i, t + kk*period);
            }
          }
          rowvec vtmp = Y(i, _);
          a += taux(i, k)*(vtmp*trans(Ai) - (Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi*trans(Ai));
          c += taux(i, k)*(vtmp*trans(Wi) - (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai*trans(Wi));
          Sw += taux(i, k)*Wi*trans(Wi);
          rowvec mtmp = vtmp - (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai-(Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi;
          b += taux(i, k)*as_scalar(mtmp*trans(mtmp));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
        newdelta = join_rows(newdelta, c*inv(Sw)); 
        newsigma[k]= sqrt(b/(period*sum(taux(_, k))));
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
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), newbeta, newsigma, newdelta);
    rowvec tmp(newparam.size());
    tmp.fill(prec);
    if (all(abs(newparam-vparam)<tmp)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    sigma = newsigma;
    delta = newdelta;
    pi = newpi;
  }
  return(NumericVector(vparam.begin(), vparam.end()));
}
// ----------------------------------------------------------------------------
// EM not censored sigma unique
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EMSigmaunique_cpp(NumericVector param,
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
                     int itermax, 
                     bool EMIRLS,
                     int refgr){
  int period = A.ncol();
  double prec = 0.000001;
  NumericVector pi(ng);
  NumericVector beta;
  NumericVector sigma;
  NumericVector delta;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    sigma = param[Range(ng+sum(nbeta)-1, ng+sum(nbeta)+ng-2)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng+sum(nbeta)+ng - 1, param.length() - 1)];
    }
    pi.push_back(1-sum(pi));  
  }else{
    pi = param[Range(0,ng*nx-1)];
    beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
    sigma = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+ng-1)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng*nx+sum(nbeta)+ng, param.length() - 1)];
    }
  }
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta), as<arma::rowvec>(sigma), as<arma::rowvec>(delta));
  int tour = 1;
  double b = 0;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEM_cpp(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodCNORM_cpp(NumericVector(vparam.begin(), vparam.end()), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw));
    }
    // E-step
    NumericMatrix  taux = ftauxCNORM_cpp(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X);
    rowvec newbeta;
    rowvec newdelta;
    if (nw == 0){
      b = 0;
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        mat Ai(nbeta[k], period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
          }
          rowvec vtmp = Y(i, _);
          a += taux(i, k)*vtmp*trans(Ai);
          rowvec mtmp = vtmp - trans((Rcpp::as<arma::vec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1))*Ai;
          b += taux(i, k)*as_scalar(mtmp*trans(mtmp));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
      }
    }else{
      // TCOV is not NULL
      NumericVector ndeltacum(ng);
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
      mat mTCOV = as<arma::mat>(TCOV);
      b = 0;
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        rowvec c(nw);
        c.fill(0);
        mat Sw(nw, nw);
        Sw.fill(0);
        mat Ai(nbeta[k], period);
        mat Wi(nw, period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
            for (int kk = 0; kk < nw; ++kk){
              Wi(kk, t) = mTCOV(i, t + kk*period);
            }
          }
          rowvec vtmp = Y(i, _);
          a += taux(i, k)*(vtmp*trans(Ai) - (Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi*trans(Ai));
          c += taux(i, k)*(vtmp*trans(Wi) - (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai*trans(Wi));
          Sw += taux(i, k)*Wi*trans(Wi);
          rowvec mtmp = vtmp - (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai-(Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi;
          b += taux(i, k)*as_scalar(mtmp*trans(mtmp));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
        newdelta = join_rows(newdelta, c*inv(Sw)); 
      }
    }
    NumericVector vstmp = rep(sqrt(b/(period*n)), ng);
    rowvec newsigma(ng);
    newsigma = Rcpp::as<arma::rowvec>(vstmp);
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
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), newbeta, newsigma, newdelta);
    rowvec tmp(newparam.size());
    tmp.fill(prec);
    if (all(abs(newparam-vparam)<tmp)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    sigma = newsigma;
    delta = newdelta;
    pi = newpi;
  }
  return(NumericVector(vparam.begin(), vparam.end()));
}
// ----------------------------------------------------------------------------
// EM censored
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EMCensored_cpp(NumericVector param,
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
                     int itermax, 
                     bool EMIRLS,
                     int refgr){
  int period = A.ncol();
  double prec = 0.000001;
  NumericVector pi(ng);
  NumericVector beta;
  NumericVector sigma;
  NumericVector delta;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    sigma = param[Range(ng+sum(nbeta)-1, ng+sum(nbeta)+ng-2)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng+sum(nbeta)+ng - 1, param.length() - 1)];
    }
    pi.push_back(1-sum(pi));  
  }else{
    pi = param[Range(0,ng*nx-1)];
    beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
    sigma = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+ng-1)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng*nx+sum(nbeta)+ng, param.length() - 1)];
    }
  }
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta), as<arma::rowvec>(sigma), as<arma::rowvec>(delta));
  int tour = 1;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEM_cpp(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodCNORM_cpp(NumericVector(vparam.begin(), vparam.end()), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw));
    }
    // E-step
    NumericMatrix  taux = ftauxCNORM_cpp(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X);
    rowvec newbeta;
    rowvec newdelta;
    rowvec newsigma(ng);
    if (nw == 0){
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        double b = 0;
        mat Ai(nbeta[k], period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
          }
          NumericVector muikt;
          for (int s = 0; s < period; ++s){
            NumericVector vtmp2;
            for (int po = 0; po < nbeta[k]; ++po){
              vtmp2.push_back(pow(A(i,s), po));
            }
            NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
            muikt.push_back(sum(betak*vtmp2));
          }
          NumericVector vymax = rep(ymax, period) - muikt;
          rowvec alphamax = as<arma::rowvec>(vymax)/sigma[k];
          NumericVector vymin = rep(ymin, period) - muikt;
          rowvec alphamin = as<arma::rowvec>(vymin)/sigma[k]; 
          rowvec qmax = normpdf(alphamax)/normcdf(-alphamax);
          rowvec qmin = normpdf(alphamin)/normcdf(alphamin);
          rowvec Ytildei(period);
          rowvec Ytilde2i(period);
          for (int t = 0; t < period; ++t){
            if (Y(i, t) <= ymin){
              Ytildei[t] = muikt[t]-sigma[k]*qmin[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1-alphamin[t]*qmin[t])+pow(muikt[t], 2)-2*muikt[t]*sigma[k]*qmin[t];
            }else if (Y(i, t) >= ymax){
              Ytildei[t] = muikt[t]+sigma[k]*qmax[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1+alphamax[t]*qmax[t])+pow(muikt[t], 2)+2*muikt[t]*sigma[k]*qmax[t];
            }else{
              Ytildei[t] = Y(i, t);
              Ytilde2i[t] = Y(i,t)*Y(i,t);
            }
          }
          a += taux(i, k)*Ytildei*trans(Ai);
          rowvec vtmp = trans((Rcpp::as<arma::vec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1))*Ai;
          mat un(period, 1);
          un.ones();
          b += taux(i, k)*as_scalar(Ytilde2i*un-2*Ytildei*trans(vtmp)+vtmp*trans(vtmp));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
        newsigma[k]= sqrt(b/(period*sum(taux(_, k))));
      }
    }else{
      // TCOV is not NULL
      NumericVector ndeltacum(ng);
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
      mat mTCOV = as<arma::mat>(TCOV);
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        rowvec c(nw);
        c.fill(0);
        mat Sw(nw, nw);
        Sw.fill(0);
        double b = 0;
        mat Ai(nbeta[k], period);
        mat Wi(nw, period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
            for (int kk = 0; kk < nw; ++kk){
              Wi(kk, t) = mTCOV(i, t + kk*period);
            }
          }
          NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
          NumericVector deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
          NumericVector muikt;
          for (int s = 0; s < period; ++s){
            NumericVector vtmp2;
            for (int po = 0; po < nbeta[k]; ++po){
              vtmp2.push_back(pow(A(i,s), po));
            }
            muikt.push_back(sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, s, k));
          }
          NumericVector vymax = rep(ymax, period) - muikt;
          rowvec alphamax = as<arma::rowvec>(vymax)/sigma[k];
          NumericVector vymin = rep(ymin, period) - muikt;
          rowvec alphamin = as<arma::rowvec>(vymin)/sigma[k]; 
          rowvec qmax = normpdf(alphamax)/normcdf(-alphamax);
          rowvec qmin = normpdf(alphamin)/normcdf(alphamin);
          rowvec Ytildei(period);
          rowvec Ytilde2i(period);
          for (int t = 0; t < period; ++t){
            if (Y(i, t) <= ymin){
              Ytildei[t] = muikt[t]-sigma[k]*qmin[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1-alphamin[t]*qmin[t])+pow(muikt[t], 2)-2*muikt[t]*sigma[k]*qmin[t];
            }else if (Y(i, t) >= ymax){
              Ytildei[t] = muikt[t]+sigma[k]*qmax[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1+alphamax[t]*qmax[t])+pow(muikt[t], 2)+2*muikt[t]*sigma[k]*qmax[t];
            }else{
              Ytildei[t] = Y(i, t);
              Ytilde2i[t] = Y(i,t)*Y(i,t);
            }
          }
          a += taux(i, k)*(Ytildei*trans(Ai) - (Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi*trans(Ai));
          c += taux(i, k)*(Ytildei*trans(Wi) - (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai*trans(Wi));
          Sw += taux(i, k)*Wi*trans(Wi);
          rowvec mtmp = (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai+(Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi;
          mat un(period, 1);
          un.ones();
          b += taux(i, k)*(as_scalar(Ytilde2i*un)-2*as_scalar(Ytildei*trans(mtmp))+as_scalar(mtmp*trans(mtmp)));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
        newdelta = join_rows(newdelta, c*inv(Sw)); 
        newsigma[k]= sqrt(b/(period*sum(taux(_, k))));
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
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), newbeta, newsigma, newdelta);
    rowvec tmp(newparam.size());
    tmp.fill(prec);
    if (all(abs(newparam-vparam)<tmp)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    sigma = newsigma;
    delta = newdelta;
    pi = newpi;
  }
  return(NumericVector(vparam.begin(), vparam.end()));
}
// ----------------------------------------------------------------------------
// EM censored same sigma
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector EMCensoredSigmaunique_cpp(NumericVector param,
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
                             int itermax, 
                             bool EMIRLS,
                             int refgr){
  int period = A.ncol();
  double prec = 0.000001;
  NumericVector pi(ng);
  NumericVector beta;
  NumericVector sigma;
  NumericVector delta;
  NumericVector nbetacum(nbeta.size());
  std::partial_sum(nbeta.begin(), nbeta.end(), nbetacum.begin());
  nbetacum.push_front(0);
  if (nx == 1){
    pi = param[Range(0,ng-2)];
    beta = param[Range(ng-1,ng+sum(nbeta)-2)];
    sigma = param[Range(ng+sum(nbeta)-1, ng+sum(nbeta)+ng-2)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng+sum(nbeta)+ng - 1, param.length() - 1)];
    }
    pi.push_back(1-sum(pi));  
  }else{
    pi = param[Range(0,ng*nx-1)];
    beta = param[Range(ng*nx,ng*nx+sum(nbeta)-1)];
    sigma = param[Range(ng*nx+sum(nbeta), ng*nx+sum(nbeta)+ng-1)]; 
    if (param.length() > ng*nx+sum(nbeta)+ng){
      delta = param[Range(ng*nx+sum(nbeta)+ng, param.length() - 1)];
    }
  }
  rowvec vparam = join_rows(as<arma::rowvec>(pi), as<arma::rowvec>(beta), as<arma::rowvec>(sigma), as<arma::rowvec>(delta));
  int tour = 1;
  double b = 0;
  while (tour < itermax){
    if (nx == 1){
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodEM_cpp(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw));
    }else{
      // a modifier
      Rprintf("iter %3d value ", tour);
      Rprintf("%.6f\n", -likelihoodCNORM_cpp(NumericVector(vparam.begin(), vparam.end()), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw));
    }
    // E-step
    NumericMatrix  taux = ftauxCNORM_cpp(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X);
    rowvec newbeta;
    rowvec newdelta;
    if (nw == 0){
      b = 0;
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        mat Ai(nbeta[k], period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
          }
          NumericVector muikt;
          for (int s = 0; s < period; ++s){
            NumericVector vtmp2;
            for (int po = 0; po < nbeta[k]; ++po){
              vtmp2.push_back(pow(A(i,s), po));
            }
            NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
            muikt.push_back(sum(betak*vtmp2));
          }
          NumericVector vymax = rep(ymax, period) - muikt;
          rowvec alphamax = as<arma::rowvec>(vymax)/sigma[k];
          NumericVector vymin = rep(ymin, period) - muikt;
          rowvec alphamin = as<arma::rowvec>(vymin)/sigma[k]; 
          rowvec qmax = normpdf(alphamax)/normcdf(-alphamax);
          rowvec qmin = normpdf(alphamin)/normcdf(alphamin);
          rowvec Ytildei(period);
          rowvec Ytilde2i(period);
          for (int t = 0; t < period; ++t){
            if (Y(i, t) <= ymin){
              Ytildei[t] = muikt[t]-sigma[k]*qmin[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1-alphamin[t]*qmin[t])+pow(muikt[t], 2)-2*muikt[t]*sigma[k]*qmin[t];
            }else if (Y(i, t) >= ymax){
              Ytildei[t] = muikt[t]+sigma[k]*qmax[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1+alphamax[t]*qmax[t])+pow(muikt[t], 2)+2*muikt[t]*sigma[k]*qmax[t];
            }else{
              Ytildei[t] = Y(i, t);
              Ytilde2i[t] = Y(i,t)*Y(i,t);
            }
          }
          a += taux(i, k)*Ytildei*trans(Ai);
          rowvec vtmp = trans((Rcpp::as<arma::vec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1))*Ai;
          mat un(period, 1);
          un.ones();
          b += taux(i, k)*as_scalar(Ytilde2i*un-2*Ytildei*trans(vtmp)+vtmp*trans(vtmp));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
      }
    }else{
      // TCOV is not NULL
      NumericVector ndeltacum(ng);
      NumericVector deltatmp(ng);
      deltatmp.fill(nw);
      std::partial_sum(deltatmp.begin(), deltatmp.end(), ndeltacum.begin());
      ndeltacum.push_front(0);
      mat mTCOV = as<arma::mat>(TCOV);
      b = 0;
      for (int k = 0; k < ng; ++k){
        rowvec a(nbeta[k]);
        a.fill(0);
        rowvec c(nw);
        c.fill(0);
        mat Sw(nw, nw);
        Sw.fill(0);
        mat Ai(nbeta[k], period);
        mat Wi(nw, period);
        for (int i = 0; i < n; ++i){
          for (int t = 0; t < period; ++t){
            for (int kk = 0; kk < nbeta[k]; ++kk){
              Ai(kk, t) = pow(A(i, t), kk);
            }
            for (int kk = 0; kk < nw; ++kk){
              Wi(kk, t) = mTCOV(i, t + kk*period);
            }
          }
          NumericVector betak = beta[Range(nbetacum[k], nbetacum[k+1]-1)];
          NumericVector deltak = delta[Range(ndeltacum[k], ndeltacum[k+1]-1)];
          NumericVector muikt;
          for (int s = 0; s < period; ++s){
            NumericVector vtmp2;
            for (int po = 0; po < nbeta[k]; ++po){
              vtmp2.push_back(pow(A(i,s), po));
            }
            muikt.push_back(sum(betak*vtmp2) + WitEM_cpp(TCOV, period, deltak, nw, i, s, k));
          }
          NumericVector vymax = rep(ymax, period) - muikt;
          rowvec alphamax = as<arma::rowvec>(vymax)/sigma[k];
          NumericVector vymin = rep(ymin, period) - muikt;
          rowvec alphamin = as<arma::rowvec>(vymin)/sigma[k]; 
          rowvec qmax = normpdf(alphamax)/normcdf(-alphamax);
          rowvec qmin = normpdf(alphamin)/normcdf(alphamin);
          rowvec Ytildei(period);
          rowvec Ytilde2i(period);
          for (int t = 0; t < period; ++t){
            if (Y(i, t) <= ymin){
              Ytildei[t] = muikt[t]-sigma[k]*qmin[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1-alphamin[t]*qmin[t])+pow(muikt[t], 2)-2*muikt[t]*sigma[k]*qmin[t];
            }else if (Y(i, t) >= ymax){
              Ytildei[t] = muikt[t]+sigma[k]*qmax[t];
              Ytilde2i[t] = pow(sigma[k], 2)*(1+alphamax[t]*qmax[t])+pow(muikt[t], 2)+2*muikt[t]*sigma[k]*qmax[t];
            }else{
              Ytildei[t] = Y(i, t);
              Ytilde2i[t] = Y(i,t)*Y(i,t);
            }
          }
          a += taux(i, k)*(Ytildei*trans(Ai) - (Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi*trans(Ai));
          c += taux(i, k)*(Ytildei*trans(Wi) - (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai*trans(Wi));
          Sw += taux(i, k)*Wi*trans(Wi);
          rowvec mtmp = (Rcpp::as<arma::rowvec>(beta)).subvec(nbetacum[k], nbetacum[k+1]-1)*Ai+(Rcpp::as<arma::rowvec>(delta)).subvec(ndeltacum[k], ndeltacum[k+1]-1)*Wi;
          mat un(period, 1);
          un.ones();
          b += taux(i, k)*(as_scalar(Ytilde2i*un)-2*as_scalar(Ytildei*trans(mtmp))+as_scalar(mtmp*trans(mtmp)));
        }
        newbeta = join_rows(newbeta, a*inv(Ai*trans(Ai))/sum(taux(_, k)));
        newdelta = join_rows(newdelta, c*inv(Sw)); 
      }
    }
    NumericVector vstmp = rep(sqrt(b/(period*n)), ng);
    rowvec newsigma(ng);
    newsigma = Rcpp::as<arma::rowvec>(vstmp);
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
    rowvec newparam = join_rows(as<arma::rowvec>(newpi), newbeta, newsigma, newdelta);
    rowvec tmp(newparam.size());
    tmp.fill(prec);
    if (all(abs(newparam-vparam)<tmp)){
      tour = itermax + 2;
    }
    ++tour;
    vparam = newparam;
    beta = newbeta;
    sigma = newsigma;
    delta = newdelta;
    pi = newpi;
  }
  return(NumericVector(vparam.begin(), vparam.end()));
}

