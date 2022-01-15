#include "CommonFunction.h"
#include "CensoredNormal.h"
#include "Logit.h"
#include "ZIP.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::interfaces(r, cpp)]]


// ----------------------------------------------------------------------------
// convtolist
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List convtolist_cpp(arma::vec v,
                    arma::vec nelt){
  List ltmp;
  vec  z = zeros<vec>(1);
  nelt = join_cols(z, nelt);
  vec eletcum = cumsum(nelt);
  int imax = nelt.size();
  
  for (int nl = 0; nl < imax - 1; ++nl){
    ltmp.push_back(v.subvec(eletcum[nl], eletcum[nl+1]-1));
  }
  return(ltmp);
}

// ----------------------------------------------------------------------------
// contruct the matrix in the proof of the number of psy's parameters 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat mPsi_cpp(arma::vec psi,
                   List lng){
  
  vec psitmp = psi;
  mat mres;
  
  int cpt = 0;
  for (int i = 0; i < (lng.size()-1); ++i){
    mat mres2;
    int row = lng[i];
    for (int j = 0; j < i+1; ++j){
      int col = lng[j];
      mat mrzero(row, col, arma::fill::zeros);
      mres2 = join_rows(mres2, mrzero);
    }
    for (int j = i+1;  j < (lng.size()); ++j){
      int col = lng[j];
      mat mtmp(row-1, col-1);
      for (int jm = 0; jm < col-1; ++jm){
        for (int im = 0; im < row-1; ++im){
          mtmp(im, jm) = psitmp[cpt];
          cpt += 1;
        }
      }
      mat mrzero(row-1, 1, arma::fill::zeros);
      mat mczero(1, col, arma::fill::zeros);
      
      mres2 = join_rows(mres2, join_cols(mczero, join_rows(mrzero, mtmp)));
    }
    mres = join_cols(mres, mres2);
  }
  
  return(mres);
}
// ----------------------------------------------------------------------------
// extract a submatrix of spi
// ----------------------------------------------------------------------------
arma::mat extmat_cpp(arma::mat m,
                     int i,
                     int j,
                     List lng){
  int si = 0;
  int sj = 0;
  
  if (i > 0){
    for (int ind = 0; ind < i; ++ind){
      int tmp = lng[ind];
      si += tmp;
    }
  }
  if (j > 0){
    for (int ind = 0; ind < j; ++ind){
      int tmp = lng[ind];
      sj += tmp;
    }
  }
  int itmp = lng[i];
  int jtmp = lng[j];
  
  return(m.submat(si, sj, si+itmp-1, sj+jtmp-1));
}


// ----------------------------------------------------------------------------
//  compute the value muk for a given vector (k1, ..., kJ)
// ----------------------------------------------------------------------------
double mukMult_cpp(List ltheta,
                   arma::mat mPsi,
                   int i,
                   arma::vec vk,
                   List lng,
                   List lX){
  double muk = 0;
  
  for( int j = 0; j < ltheta.size(); ++j){
    mat mtmp = lX[j];
    int ntheta = mtmp.n_cols;  
    double spsi = 0;
    if (j != 0){
      for (int h = 0; h < j; ++h){
        mat mext = extmat_cpp(mPsi, h, j, lng);
        spsi += mext(vk[h], vk[j]);
      } 
    }
    vec vtheta = ltheta[j];
    colvec vtmp = mtmp.row(i);
    
    muk += as_scalar(vtheta.subvec(vk[j]*ntheta, (vk[j]+1)*ntheta - 1)*vtmp) + spsi;
  }
  return(muk);
}

// ----------------------------------------------------------------------------
//  compute the joint probability 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double piikMult_cpp(List ltheta,
                    arma::mat mPsi,
                    int i,
                    arma::vec vk,
                    List lng,
                    List lX,
                    arma::mat mk){
  double s = 0;
  int itmp = mk.n_rows;
  for (int l = 0; l < itmp; ++l){
    s += exp(mukMult_cpp(ltheta, mPsi, i, conv_to<vec>::from(mk.row(l)), lng, lX));
  }
  return(exp(mukMult_cpp(ltheta, mPsi, i, vk, lng, lX))/s);
}

// ----------------------------------------------------------------------------
//  choice of the density in the calculus of the probability
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double densityChoice_cpp(List beta,
                         Nullable<NumericVector> alphainit,
                         int i,
                         int k,
                         IntegerVector nbeta,
                         NumericMatrix A,
                         NumericMatrix Y,
                         double ymin,
                         double ymax, 
                         Nullable<NumericMatrix> TCOV,
                         Nullable<List> delta,
                         int nw,
                         Nullable<List> nuinit,
                         Nullable<IntegerVector> nnuinit,
                         std::string model){
  double s;
  List nu(nuinit);
  IntegerVector nnu(nnuinit);
  NumericVector alpha(alphainit);

  if (model == "LOGIT"){
    s = gkLOGIT_cpp(beta, i, k, nbeta, A, Y, TCOV, delta, nw);
  }else if (model == "CNORM"){
    alpha = alphainit;
    s = gkalpha_cpp(beta, alpha, i, k, nbeta, A, Y, ymin, ymax, TCOV, delta, nw);
  }else {
    nu = nuinit;
    nnu = nnuinit;
    s = gkZIP_cpp(beta, nu, i, k, nbeta, nnu, A, Y, TCOV, delta, nw);
  }
  return(s);
}

// ----------------------------------------------------------------------------
// Likelihood multiple
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double likelihoodMult_cpp(arma::vec vparam,
                          List lng, 
                          List lnx,
                          List lnbeta,
                          List ln,
                          List lA,
                          List lY,
                          List lX,
                          List lymin,
                          List lymax, 
                          Nullable<List> lTCOVinit,
                          List lnw,
                          arma::vec vp,
                          arma::mat mk,
                          List lnnu,
                          std::vector<std::string> model){
  List lTCOV(lTCOVinit);  
  if (lTCOVinit.isNotNull()){
    lTCOV= lTCOVinit;
  }
  
  List lparam;
  vec vpcum;
  vpcum = cumsum(vp);
  int itmp = vp.size();
  
  for (int i = 0; i < itmp - 1; ++i){
    lparam.push_back(vparam.subvec(vpcum[i], vpcum[i+1] - 1));
  }

  List ltheta;
  List lbeta;
  Rcpp::List lalpha(vp.size() - 1);
  List ldelta;
  List lnu(vp.size() - 1);
  mat mpsi;
  
  int ilparam = lparam.size();
  for (int nl = 0; nl < ilparam - 1; ++nl){
    vec vtmp = lparam[nl];
    int ilng = lng[nl];
    int ilnx = lnx[nl];
    vec  z = zeros<vec>(1);
    ltheta.push_back(join_cols(z, vtmp.subvec(0, (ilng - 1)*ilnx - 1)));
    vec lnbetatmp = lnbeta[nl];
    lbeta.push_back(vtmp.subvec((ilng - 1)*ilnx, (ilng - 1)*ilnx + sum(lnbetatmp) - 1));
    if (model[nl] == "CNORM"){
      lalpha[nl] = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + ilng - 1);  
    }
    if (model[nl] ==  "ZIP"){
      vec lnnutmp = lnnu[nl];
      lnu[nl] = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + sum(lnnutmp) - 1);
    }
  }
  mpsi = mPsi_cpp(lparam[ilparam - 1], lng);
  
  double a = 0;
  int n = ln[0];
  for (int i = 0; i < n; ++i){
    double s = 0;
    int itmp = mk.n_rows;
    for (int l = 0; l < itmp; ++l){
      double prodprob = 1;
      int itmp2 = mk.n_cols;
      for (int m = 0; m < itmp2; ++m){
        NumericMatrix mTCOV;
        if (lTCOV.size() != 0){
          mTCOV = as<NumericMatrix>(lTCOV[m]);
        }
        
        List lnutmp;
        IntegerVector nnutmp;
        
        if (model[m] == "ZIP"){
          lnutmp = convtolist_cpp(lnu[m], lnnu[m]);
          nnutmp = as<IntegerVector>(lnnu[m]);
        }
        
        NumericVector alpham;
        if (model[m] == "CNORM"){
          alpham = as<NumericVector>(lalpha[m]);
        }
        
        prodprob *= densityChoice_cpp(convtolist_cpp(lbeta[m], lnbeta[m]),
                                      alpham,
                                      i,
                                      mk(l,m),
                                      as<IntegerVector>(lnbeta[m]),
                                      as<NumericMatrix>(lA[m]),
                                      as<NumericMatrix>(lY[m]),
                                      as<double>(lymin[m]),
                                      as<double>(lymax[m]),
                                      mTCOV,
                                      ldelta,
                                      lnw[m],
                                      lnutmp,
                                      nnutmp,
                                      model[m]);
      }   
      s += piikMult_cpp(ltheta, mpsi, i, conv_to<vec>::from(mk.row(l)), lng, lX, mk)*prodprob;
    }
    a += log(s);
  }
  
  return(a);
}

// ----------------------------------------------------------------------------
// LOGIT dif likelihood betak for k and i 
// ----------------------------------------------------------------------------
NumericVector difLbetavkiLOGIT_cpp(List beta,
                                   Nullable<List> delta,
                                   int k,
                                   int i,
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
    // initialize mean
    NumericVector muikt = muikt_cpp(beta[k], nbeta[k], i, period, A, TCOV, delta, nw, k);
    NumericVector tmp1;
    for (int t = 0; t < period; ++t){
      tmp1.push_back(pow(1-1/(1+exp(muikt[t])), Y(i,t))*pow(1/(1+exp(muikt[t])), 1-Y(i,t)));
    }
    double tmp2 = 0;
    for (int t = 0; t < period; ++t){
      NumericVector tmp1minus = tmp1;
      tmp1minus.erase(tmp1minus.begin() + t);
      double ytmp = -1;
      if (Y(i, t) == 1){ ytmp = 1;}
      tmp2 += ytmp*pow(A(i, t), l)/(1+exp(muikt[t]))*(1-1/(1+exp(muikt[t])))*prodvect(tmp1minus);
    }
    betas.push_back(tmp2);
  }
  return(betas); 
}

// ----------------------------------------------------------------------------
// ZIP dif likelihood betak for k and i
// ----------------------------------------------------------------------------
arma::vec difLbetavkiZIP_cpp(List beta,
                                 List nu,
                                 Nullable<List> delta,
                                 int k,
                                 int i,
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
    betas.push_back(difbkl0*py+difbkl*py0);
  }
  return(betas); 
}

// ----------------------------------------------------------------------------
// ZIP dif likelihood nuk for k and i
// ----------------------------------------------------------------------------
arma::vec difLnuvkiZIP_cpp(List beta,
                               List nu,
                               Nullable<List> delta,
                               int k,
                               int i,
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
    
    
    nus.push_back(difbkl0*py+difbkl*py0);
  }
  return(nus); 
}

// ----------------------------------------------------------------------------
// differential betak for k and i
// ----------------------------------------------------------------------------
arma::vec difbetavkialpha_cpp(List beta,
                              NumericVector alpha,
                              Nullable<List> delta,
                              int k,
                              int i,
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
    betas.push_back(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax);
  }
  
  return(betas); 
}

// ----------------------------------------------------------------------------
// differential alphak for k and i
// ----------------------------------------------------------------------------
double difalphavki_cpp(List beta,
                       NumericVector alpha,
                       Nullable<List> delta,
                       int k,
                       int i,
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
  double alphas = difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax;
  return(alphas); 
}

// ----------------------------------------------------------------------------
// dif likelihood thetapsi
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::rowvec difLalphaMult_cpp(arma::vec vparam,
                               List lng, 
                               List lnx,
                               List lnbeta,
                               List ln,
                               List lA,
                               List lY,
                               List lX,
                               List lymin,
                               List lymax, 
                               Nullable<List> lTCOVinit,
                               List lnw,
                               arma::vec vp,
                               arma::mat mk,
                               List lnnu,
                               std::vector<std::string> model){
  List lTCOV(lTCOVinit);  
  if (lTCOVinit.isNotNull()){
    lTCOV= lTCOVinit;
  }
  
  List lparam;
  vec vpcum;
  vpcum = cumsum(vp);
  int itmp =  vp.size();
  
  for (int i = 0; i < itmp - 1; ++i){
    lparam.push_back(vparam.subvec(vpcum[i], vpcum[i+1] - 1));
  }
  
  List ltheta;
  List lbeta;
  List lalpha(vp.size() - 1);
  List ldelta;
  List lnu(vp.size() - 1);
  IntegerVector nuind(vp.size() - 1);
  nuind.fill(-1);  
  IntegerVector cnormind(vp.size() - 1);
  cnormind.fill(-1);
  
  List lthetas;
  List lbetas;
  List lnus(vp.size() - 1);
  List lalphas(vp.size() - 1);
  mat mpsi;
  
  int ilparam = lparam.size();
  for (int nl = 0; nl < ilparam - 1; ++nl){
    vec vtmp = lparam[nl];
    int ilng = lng[nl];
    int ilnx = lnx[nl];
    vec  z = zeros<vec>(1);
    ltheta.push_back(join_cols(z, vtmp.subvec(0, (ilng - 1)*ilnx - 1)));
    vec lnbetatmp = lnbeta[nl];
    lbeta.push_back(vtmp.subvec((ilng - 1)*ilnx, (ilng - 1)*ilnx + sum(lnbetatmp)- 1));
    if (model[nl] == "CNORM"){
      lalpha[nl] = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + ilng - 1);  
      cnormind[nl] = nl;
    }
    if (model[nl] ==  "ZIP"){
      vec lnnutmp = lnnu[nl];
      lnu[nl] = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + sum(lnnutmp) - 1);
      nuind[nl] = nl;
    }
    // initialize the solution list with zero
    vec vtmp1 = vtmp.subvec(0, (ilng - 1)*ilnx - 1);
    vtmp1.zeros();
    lthetas.push_back(vtmp1);
    vtmp1 = vtmp.subvec((ilng - 1)*ilnx, (ilng - 1)*ilnx + sum(lnbetatmp)- 1);
    vtmp1.zeros();
    lbetas.push_back(vtmp1);
    if (model[nl] == "CNORM"){
      vtmp1 = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + ilng - 1);
      vtmp1.zeros();
      lalphas[nl] = vtmp1;
    }
    if (model[nl] ==  "ZIP"){
      vec lnnutmp = lnnu[nl];
      vtmp1 = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + sum(lnnutmp) - 1);
      vtmp1.zeros();
      lnus[nl] = vtmp1;
    }
  }
  mpsi = mPsi_cpp(lparam[ilparam - 1], lng);
  
  vec vpsis = lparam[ilparam - 1];
  vpsis.zeros();
  
  int n = ln[1];

  for (int i = 0; i < n; ++i){
    // common values
    vec vpiikMult(mk.n_rows);
    vec vprodprob;
    mat matgk;
    int itmp = mk.n_rows;
    
    for (int l = 0; l < itmp; ++l){
      vec vtmp(mk.n_cols);
      int itmp2 = mk.n_cols;
      for (int m = 0; m < itmp2; ++m){
        NumericMatrix mTCOV;
        if (lTCOV.size() != 0){
          mTCOV = as<NumericMatrix>(lTCOV[m]);
        } 
        List lnutmp;
        IntegerVector nnutmp;
        
        if (nuind[m] == m){
          lnutmp = convtolist_cpp(lnu[m], lnnu[m]);
          nnutmp = as<IntegerVector>(lnnu[m]);
        }
        
        NumericVector alphaj;
        if (cnormind[m] == m){
          alphaj = as<NumericVector>(lalpha[m]);
        }
        
        vtmp(m) = densityChoice_cpp(convtolist_cpp(lbeta[m], lnbeta[m]),
                                               alphaj,
                                               i,
                                               mk(l,m),
                                               as<IntegerVector>(lnbeta[m]),
                                               as<NumericMatrix>(lA[m]),
                                               as<NumericMatrix>(lY[m]),
                                               as<double>(lymin[m]),
                                               as<double>(lymax[m]),
                                               mTCOV,
                                               ldelta,
                                               lnw[m],
                                               lnutmp,
                                               nnutmp,
                                               model[m]);
      }
      matgk = join_cols(matgk, arma::conv_to<rowvec>::from(vtmp));
      vpiikMult(l) = piikMult_cpp(ltheta, mpsi, i, conv_to<vec>::from(mk.row(l)), lng, lX, mk);
    }
    vprodprob = prod(matgk, 1);
    double deno = as_scalar(conv_to<rowvec>::from(vpiikMult)*vprodprob);
    
    //theta
    for (int j = 0; j < ltheta.size(); ++j){
      int indx = lnx[j];
      int indj = lng[j];
      vec vtmp(indx*(indj-1));
      int indv = 0;
      for (int k = 1; k < indj; ++k){
        uvec indk = find(mk.col(j) == k);
        for (int l = 0; l < indx; ++l){
          double spiiikkjin = sum(vpiikMult(indk));
          double spiikgkjin = as_scalar(conv_to<rowvec>::from(vpiikMult(indk))*vprodprob(indk));
          mat mtmp = lX[j];
          vtmp(indv) = mtmp(i, l)*(spiikgkjin/deno - spiiikkjin);
          indv++; 
        }
      }
      vec vtmptheta = lthetas(j);
      lthetas(j) = vtmptheta + vtmp;
    }
    
    // psi
    vec vtmp(vpsis.size());
    int indv = 0;
    int indng = lng.size();
    for (int j = 0; j < indng - 1; ++j){
      for (int h = j + 1; h < indng; ++h){
        mat mext = extmat_cpp(mpsi, j, h, lng);
        int indc = mext.n_cols - 1;
        int indr = mext.n_rows - 1;
        for (int l  = 0; l < indc; ++l){
          for ( int k = 0; k < indr; ++k){
            uvec indkl = find(mk.col(j) == k+1 && mk.col(h) == l+1);
            double spiiikkjin = sum(vpiikMult(indkl));
            double spiikgkjin = as_scalar(conv_to<rowvec>::from(vpiikMult(indkl))*vprodprob(indkl));
            vtmp(indv) = spiikgkjin/deno - spiiikkjin;
            indv++; 
          }
        }
      }
    }
    vpsis = vpsis + vtmp;
    
    // beta
    for (int j = 0; j < lbeta.size(); ++j){
      vec lnbetaj = lnbeta[j];
      int indv = 0;
      int indj = lng[j];
      vec vtmp(sum(lnbetaj));

      for (int k = 0; k < indj; ++k){
        uvec indk = find(mk.col(j) == k);
        NumericMatrix mTCOV;
        if (lTCOV.size() != 0){
          mTCOV = as<NumericMatrix>(lTCOV[j]);
        }
        List lnutmp;
        IntegerVector nnutmp;

        if (nuind[j] == j){
          lnutmp = convtolist_cpp(lnu[j], lnnu[j]);
          nnutmp = as<IntegerVector>(lnnu[j]);
        }

        NumericVector alphaj;
        if (cnormind[j] == j){
          alphaj = as<NumericVector>(lalpha[j]);
        }

        vec difgk;
        if (model[j] == "CNORM"){
          difgk =  difbetavkialpha_cpp(convtolist_cpp(lbeta[j], lnbeta[j]),
                                           as<NumericVector>(lalpha[j]),
                                           ldelta,
                                           k,
                                           i,
                                           lng[j],
                                           lnx[j],
                                           as<IntegerVector>(lnbeta[j]),
                                           ln[j],
                                           as<NumericMatrix>(lA[j]),
                                           as<NumericMatrix>(lY[j]),
                                           as<NumericMatrix>(lX[j]),
                                           lymin[j],
                                           lymax[j],
                                           as<NumericMatrix>(mTCOV),
                                           lnw[j]);
        }else if (model[j] == "LOGIT"){
          difgk =  difLbetavkiLOGIT_cpp(convtolist_cpp(lbeta[j], lnbeta[j]),
                                           ldelta,
                                           k,
                                           i,
                                           lng[j],
                                           lnx[j],
                                           ln[j],
                                           as<IntegerVector>(lnbeta[j]),
                                           as<NumericMatrix>(lA[j]),
                                           as<NumericMatrix>(lY[j]),
                                           as<NumericMatrix>(lX[j]),
                                           as<NumericMatrix>(mTCOV),
                                           lnw[j]);
        }else if (model[j] == "ZIP"){
          difgk =  difLbetavkiZIP_cpp(convtolist_cpp(lbeta[j], lnbeta[j]),
                                          lnutmp,
                                          ldelta,
                                          k,
                                          i,
                                          lng[j],
                                          lnx[j],
                                          as<IntegerVector>(lnbeta[j]),
                                          nnutmp,
                                          ln[j],
                                          as<NumericMatrix>(lA[j]),
                                          as<NumericMatrix>(lY[j]),
                                          as<NumericMatrix>(lX[j]),
                                          as<NumericMatrix>(mTCOV),
                                          lnw[j]);
        }

        for (int l = 0; l < lnbetaj(k); ++l){
          vec vdif(indk.size());
          int itmp = indk.size();
          for (int ind = 0; ind < itmp; ++ind){
            rowvec vtmp1 = matgk.row(indk(ind));
            double prodtmp = 1;
            for (int jtmp = 0; jtmp < lbeta.size(); ++jtmp){
              if (jtmp != j){
                prodtmp *= vtmp1(jtmp);
              }
            }
            vdif(ind) =  difgk(l)*prodtmp;
          }
          vtmp(indv) = as_scalar(conv_to<rowvec>::from(vpiikMult(indk))*vdif)/deno;
          indv++;
        }
      }
      vec vtmpbeta = lbetas(j);
      lbetas(j) = vtmpbeta + vtmp;
    }

    // same sigma
    for (int j = 0; j < lbeta.size(); ++j){
      if (model[j] == "CNORM"){
      int lngj = lng[j];
      vec vtmp(lngj);
      vec vdif(mk.n_rows);
      int itmp = mk.n_rows;
      for (int ro = 0; ro < itmp; ro++){
        NumericMatrix mTCOV;
        if (lTCOV.size() != 0){
          mTCOV = as<NumericMatrix>(lTCOV[j]);
        }
        List lnutmp;
        IntegerVector nnutmp;

        if (nuind[j] == j){
          lnutmp = convtolist_cpp(lnu[j], lnnu[j]);
          nnutmp = as<IntegerVector>(lnnu[j]);
        }

        NumericVector alphaj;
        if (cnormind[j] == j){
          alphaj = as<NumericVector>(lalpha[j]);
        }

        double difgk =  difalphavki_cpp(convtolist_cpp(lbeta[j], lnbeta[j]),
                                        alphaj,
                                         ldelta,
                                         mk(ro, j),
                                         i,
                                         lng[j],
                                         lnx[j],
                                         as<IntegerVector>(lnbeta[j]),
                                         ln[j],
                                         as<NumericMatrix>(lA[j]),
                                         as<NumericMatrix>(lY[j]),
                                         as<NumericMatrix>(lX[j]),
                                         lymin[j],
                                         lymax[j],
                                         as<NumericMatrix>(mTCOV),
                                         lnw[j]);

        rowvec vtmp1 = matgk.row(ro);
        double prodtmp = 1;
        for (int jtmp = 0; jtmp < lbeta.size(); ++jtmp){
          if (jtmp != j){
            prodtmp *= vtmp1(jtmp);
          }
        }
        vdif(ro) =  difgk*prodtmp;
      }
      for (int ind = 0; ind < lngj; ind++){
        vtmp(ind) = as_scalar(conv_to<rowvec>::from(vpiikMult)*vdif)/deno;
      }
      vec vtmpalpha = lalphas(j);
      lalphas(j) = vtmpalpha + vtmp;
    }
    }

    // nu
    for (int j = 1; j < lbeta.size(); ++j){
      if (model[j] ==  "ZIP"){
        vec vdif;
        int indv = 0;
        int indj = lng[j];
        
        List lnutmp;
        IntegerVector nnutmp;
        if (nuind[j] == j){
          lnutmp = convtolist_cpp(lnu[j], lnnu[j]);
          nnutmp = as<IntegerVector>(lnnu[j]);
        }
        vec vtmp(sum(nnutmp));
        
        for (int k = 0; k < indj; ++k){
          uvec indk = find(mk.col(j) == k);
          
          NumericMatrix mTCOV;
          if (lTCOV.size() != 0){
            mTCOV = as<NumericMatrix>(lTCOV[j]);
          }
          
          NumericVector alphaj;
          if (cnormind[j] == j){
            alphaj = as<NumericVector>(lalpha[j]);
          }
          
          vec difgk =  difLnuvkiZIP_cpp(convtolist_cpp(lbeta[j], lnbeta[j]),
                                        lnutmp,
                                        ldelta,
                                        k,
                                        i,
                                        lng[j],
                                        lnx[j],
                                        as<IntegerVector>(lnbeta[j]),
                                        nnutmp,
                                        ln[j],
                                        as<NumericMatrix>(lA[j]),
                                        as<NumericMatrix>(lY[j]),
                                        as<NumericMatrix>(lX[j]),
                                        as<NumericMatrix>(mTCOV),
                                        lnw[j]);

          for (int l = 0; l < nnutmp(k); ++l){
            vec vdif(indk.size());
            int itmp = indk.size();
            for (int ind = 0; ind < itmp; ++ind){
              rowvec vtmp1 = matgk.row(indk(ind));
              double prodtmp = 1;
              for (int jtmp = 0; jtmp < lnutmp.size(); ++jtmp){
                if (jtmp != j){
                  prodtmp *= vtmp1(jtmp);
                }
              }
              vdif(ind) =  difgk(l)*prodtmp;
            }
            vtmp(indv) = as_scalar(conv_to<rowvec>::from(vpiikMult(indk))*vdif)/deno;
            indv++;
          }
        }
        vec vtmpnu = lnus[j];
        lnus[j] = vtmpnu + vtmp;
      }
    }

  }//i

  vec vparams;
  int lngn = lng.size();
  for (int j = 0; j < lngn; ++j){
    vec th = lthetas[j];
    vec lb = lbetas[j];
    vparams = join_cols(vparams, th, lb);
    if (model[j] == "CNORM"){
      vec la = lalphas[j];
      vparams = join_cols(vparams, la);
    }
    if (model[j] == "ZIP"){
      vec ln = lnus[j];
      vparams = join_cols(vparams, ln);
    }
  }
  vparams = join_cols(vparams, vpsis);
  
  return(conv_to<rowvec>::from(vparams));
}

// ----------------------------------------------------------------------------
// function piik with vector entry
// ----------------------------------------------------------------------------
double piikMultV_cpp(arma::vec vthetapsi,
                    int i,
                    arma::vec vk,
                    List lng,
                    List lX,
                    arma::mat mk){
  List ltheta;
  int indmin = 0;
  int indmax = -1;
  for (int nl = 0; nl < lng.size(); ++nl){
    mat mtmp = lX(nl);
    int lngtmp = lng[nl];
    indmax += mtmp.n_cols*lngtmp;
    ltheta.push_back(vthetapsi.subvec(indmin, indmax));
    indmin = indmax + 1;
  }
  mat mPsi = mPsi_cpp(vthetapsi.subvec(indmin, vthetapsi.size() - 1), lng);
  
  double s = 0;
  int itmp = mk.n_rows;
  for (int l = 0; l < itmp; ++l){
    s += exp(mukMult_cpp(ltheta, mPsi, i, conv_to<vec>::from(mk.row(l)), lng, lX));
  }
  return(exp(mukMult_cpp(ltheta, mPsi, i, vk, lng, lX))/s);
}

// ----------------------------------------------------------------------------
//  Function to find root for theta and psi
// ----------------------------------------------------------------------------
double ThetaPsiPiikMult_cpp(arma::vec vthetapsi,
                            List lng,
                            List lnx,
                            List lnbeta,
                            List ln,
                            List lA,
                            List lY,
                            List lX, 
                            List lymin,
                            List lymax,
                            Nullable<List> lTCOVinit,
                            List lnw,
                            arma::mat taux,
                            arma::mat mk){
  List lTCOV(lTCOVinit);  
  if (lTCOVinit.isNotNull()){
    lTCOV= lTCOVinit;
  }
  
  List ltheta;
  int indmin = 0;
  int indmax = -1;
  for (int nl = 0; nl < lng.size(); ++nl){
    mat mtmp = lX[nl];
    int lngtmp = lng[nl];
    indmax += mtmp.n_cols*(lngtmp - 1);
    int lnx = mtmp.n_cols;
    vec z(lnx);
    z.zeros();
    ltheta.push_back(join_cols(z, vthetapsi.subvec(indmin, indmax)));
    indmin = indmax + 1;
  }
  mat mPsi = mPsi_cpp(vthetapsi.subvec(indmin, vthetapsi.size() - 1), lng);

  double s = 0;
  int n = ln[1];
  for (int i = 0; i < n; ++i){
    vec muk(mk.n_rows);
    int itmp = mk.n_rows;
    for (int ind = 0; ind < itmp; ++ind){
      muk(ind) =  mukMult_cpp(ltheta, mPsi, i, conv_to<vec>::from(mk.row(ind)), lng, lX);
    }
    int itmp2 = mk.n_rows;
    for (int ind = 0; ind < itmp2; ++ind){
      s += taux(ind, i + lA.size())*(muk(ind)-log(sum(exp(muk))));
    }
  }
  
  return(s);
}

// ----------------------------------------------------------------------------
//  diferential of piik by theta and psi
// ----------------------------------------------------------------------------
arma::vec difThetaPsiPiikMult_cpp(arma::vec vthetapsi,
                            List lng,
                            List lnx,
                            List lnbeta,
                            List ln,
                            List lA,
                            List lY,
                            List lX, 
                            List lymin,
                            List lymax,
                            Nullable<List> lTCOVinit,
                            List lnw,
                            arma::mat taux,
                            arma::mat mk){
  List lTCOV(lTCOVinit);  
  if (lTCOVinit.isNotNull()){
    lTCOV= lTCOVinit;
  }
  
  List ltheta;
  int indmin = 0;
  int indmax = -1;
  for (int nl = 0; nl < lng.size(); ++nl){
    mat mtmp = lX[nl];
    int lngtmp = lng[nl];
    indmax += mtmp.n_cols*(lngtmp - 1);
    int lnx = mtmp.n_cols;
    vec z(lnx);
    z.zeros();
    ltheta.push_back(join_cols(z, vthetapsi.subvec(indmin, indmax)));
    indmin = indmax + 1;
  }
  mat mPsi = mPsi_cpp(vthetapsi.subvec(indmin, vthetapsi.size() - 1), lng);
  
  mat mprob = mk;
  int n = ln[1];
  for (int i = 0; i < n; ++i){
    vec vtmp(mk.n_rows);
    int itmp = mk.n_rows;
    for (int l = 0; l < itmp; ++l){
      vtmp(l) = piikMult_cpp(ltheta, mPsi, i, conv_to<vec>::from(mk.row(l)), lng, lX, mk);
    }
    mprob = join_rows(mprob, vtmp);
  }

  // loop for theta
  vec vthetas(1);
  int indv = 0;
  for (int j = 0; j < ltheta.size(); ++j){
    int lngj = lng[j];
    for (int k = 1; k < lngj; ++k){
      int lnxj = lnx[j];
      for (int l = 0; l < lnxj; ++l){
        double stheta = 0;
        for (int i = 0; i < n; ++i){
          int itmp = mk.n_rows;
          for (int ind2 = 0 ; ind2 < itmp; ++ind2){
            if (mk(ind2, j) == k){
              mat mlX = lX[j];
              stheta += mlX(i, l)*(taux(ind2, i + lA.size()) - mprob(ind2, i + lA.size()));
            }
          }
        }
        vthetas(indv) = stheta;
        int indvv = vthetas.size();
        vthetas.resize(indvv + 1);
        indv++;
      }
    }
  } 
  vthetas.resize(indv);
  
  // lopp for psi 
  vec psis(1);
  int indpsi = 0;
  int itmp = lng.size();
  for (int j = 0; j < itmp - 1; ++j){
    for (int h = j + 1; h < itmp; ++h){
      mat mext = extmat_cpp(mPsi, j, h, lng);
      int itmp2 =  mext.n_cols;
      for (int l = 0; l <itmp2 -1; ++l){
        int itmp3 = mext.n_rows;
        for (int k = 0; k < itmp3 - 1; ++k){
          uvec ind2 = find(mk.col(j) == k+1 && mk.col(h) == l+1);
          double stmp = 0;
          int itmp4 = ind2.size();
          for (int ind3 = 0; ind3 < itmp4; ++ind3){
            rowvec vtmp1 = taux.submat(ind2[ind3], lA.size(), ind2[ind3], taux.n_cols - 1);
            rowvec vtmp2 = mprob.submat(ind2[ind3], lA.size(), ind2[ind3], taux.n_cols - 1);
            stmp +=  sum(vtmp1 - vtmp2);
          }
          psis(indpsi) = stmp;
          int indpsit = psis.size();
          psis.resize(indpsit + 1);
          indpsi++;
        }
      }
    }
  }
  psis.resize(indpsi);
  
  vec sol = join_cols(vthetas, psis);
  return(sol);
}

// ----------------------------------------------------------------------------
//  choice of the density in the calculus of the probability for the EM algorithm
// ----------------------------------------------------------------------------
double EMdensityChoice_cpp(List beta,
                           Nullable<NumericVector> sigmainit,
                           int i,
                           int k,
                           IntegerVector nbeta,
                           NumericMatrix A,
                           NumericMatrix Y,
                           double ymin,
                           double ymax, 
                           Nullable<NumericMatrix> TCOV,
                           Nullable<List> delta,
                           int nw,
                           Nullable<List> nuinit,
                           Nullable<IntegerVector> nnuinit,
                           std::string model){
  double s;
  List nu(nuinit);
  IntegerVector nnu(nnuinit);
  NumericVector sigma(sigmainit);
  
  if (model == "LOGIT"){
    s = gkLOGIT_cpp(beta, i, k, nbeta, A, Y, TCOV, delta, nw);
  }else if (model == "CNORM"){
    sigma = sigmainit;
    // inside gkCNORM i = i-1 and k = k-1
    s = gkCNORM_cpp(beta, sigma, i + 1, k + 1, nbeta, A, Y, ymin, ymax, TCOV, delta, nw);
  }else{
    nu = nuinit;
    nnu = nnuinit;
    s = gkZIP_cpp(beta, nu, i, k, nbeta, nnu, A, Y, TCOV, delta, nw);
  }
  return(s);
}

// ----------------------------------------------------------------------------
//  Function rate multiple
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List  ftauxPiikMult_cpp(List lthetainit,
                        arma::mat mPsi,
                        List lbeta,
                        Nullable<List> lsigmainit,
                        List lng,
                        List lnbeta,
                        List ln,
                        List lA,
                        List lY,
                        List lymin,
                        List lymax,
                        Nullable<List> lTCOVinit,
                        List ldelta,
                        List lnw,
                        List lnx,
                        List lX, 
                        arma::mat mk,
                        arma::vec vp,
                        Nullable<List> nuinit,
                        Nullable<List> nnuinit,
                        std::vector<std::string> model){
  List lTCOV(lTCOVinit);  
  if (lTCOVinit.isNotNull()){
    lTCOV= lTCOVinit;
  }
  List lsigma(lsigmainit);  
  if (lsigmainit.isNotNull()){
    lsigma = lsigmainit;
  }
  
  List lnu(nuinit);  
  if (nuinit.isNotNull()){
    lnu= nuinit;
  }
  List lnnu(nnuinit);  
  if (nnuinit.isNotNull()){
    lnnu= nnuinit;
  }
  
  List ltheta;
  int itmp = lng.size();
  for (int nl = 0; nl < itmp; ++nl){
    mat mX = lX[nl];
    int lnx = mX.n_cols;
    vec z(lnx);
    z.zeros();
    vec vtmp = lthetainit[nl];
    ltheta.push_back(join_cols(z, vtmp));  
  }
  
  mat sol = mk;
  mat sol2 = mk;
  int n = ln[0];
  for (int i = 0; i < n; ++i){
    vec vpiikMult(mk.n_rows);
    vec vprod(mk.n_rows);
    int itmp = mk.n_rows;
    for (int l = 0; l < itmp; ++l){
      double prod = 1;
      int itmp2 = mk.n_cols;
      for (int m = 0; m < itmp2; ++m){
        NumericMatrix mTCOV;
        if (lTCOV.size() != 0){
          mTCOV = as<NumericMatrix>(lTCOV[l]);
        } 
        
        List lnutmp;
        IntegerVector nnutmp;
        
        if (model[m] == "ZIP"){
          lnutmp = convtolist_cpp(lnu[m], lnnu[m]);
          nnutmp = as<IntegerVector>(lnnu[m]);
        }
        
        NumericVector sigmam;
        if (model[m] == "CNORM"){
          sigmam = as<NumericVector>(lsigma[m]);
        }
        prod *= EMdensityChoice_cpp(convtolist_cpp(lbeta[m], lnbeta[m]),
                                    sigmam,
                                    i,
                                    mk(l,m),
                                    as<IntegerVector>(lnbeta[m]),
                                    as<NumericMatrix>(lA[m]),
                                    as<NumericMatrix>(lY[m]),
                                    as<double>(lymin[m]),
                                    as<double>(lymax[m]),
                                    mTCOV,
                                    ldelta,
                                    lnw[m],
                                    lnutmp,
                                    nnutmp,
                                    model[m]);
      }
      vprod(l) = prod;
      vpiikMult(l) = piikMult_cpp(ltheta, mPsi, i, conv_to<vec>::from(mk.row(l)), lng, lX, mk);
    }
    sol = join_rows(sol, (vpiikMult % vprod) /as_scalar((conv_to<rowvec>::from(vpiikMult)*vprod)));
    sol2 = join_rows(sol2, vpiikMult);
  }
  
  List lsol;
  lsol.push_back(sol);
  lsol.push_back(sol2);
  
  return(lsol);
}

// ----------------------------------------------------------------------------
//  Likelihood EM
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
double likelihoodMultEM_cpp(List lthetainit,
                            arma::mat mPsi,
                            List lbeta,
                            Nullable<List> lsigmainit,
                            List lng,
                            List lnx,
                            List lnbeta,
                            List ln,
                            List lA,
                            List lY,
                            List lX, 
                            List lymin,
                            List lymax,
                            Nullable<List> lTCOVinit,
                            List ldelta,
                            List lnw,
                            arma::mat mk,
                            arma::mat mprob,
                            arma::vec vp,
                            Nullable<List> nuinit,
                            Nullable<List> nnuinit,
                            std::vector<std::string> model){
  List lTCOV(lTCOVinit);  
  if (lTCOVinit.isNotNull()){
    lTCOV= lTCOVinit;
  }
  List lsigma(lsigmainit);  
  if (lsigmainit.isNotNull()){
    lsigma= lsigmainit;
  }
  List lnu(nuinit);  
  if (nuinit.isNotNull()){
    lnu= nuinit;
  }
  List lnnu(nnuinit);  
  if (nnuinit.isNotNull()){
    lnnu= nnuinit;
  }
  
  List ltheta;
  vec  z = zeros<vec>(1);
  int itmp = lng.size();
  for (int nl = 0; nl < itmp; ++nl){
    vec vtmp = lthetainit[nl];
    ltheta.push_back(join_cols(z, vtmp));  
  }
  
  double a = 0;
  int n = ln[0];
  for (int i = 0; i < n; ++i){
    double s = 0;
    int itmp = mk.n_rows;
    for (int l = 0; l < itmp; ++l){
      double prodprob = 1;
      int itmp2 = mk.n_cols;
      for (int m = 0; m < itmp2; ++m){
        NumericMatrix mTCOV;
        if (lTCOV.size() != 0){
          mTCOV = as<NumericMatrix>(lTCOV[m]);
        } 
        List lnutmp;
        IntegerVector nnutmp;
        if (model[m] == "ZIP"){
          lnutmp = convtolist_cpp(lnu[m], lnnu[m]);
          nnutmp = as<IntegerVector>(lnnu[m]);
        }
        NumericVector sigmam;
        if (model[m] == "CNORM"){
          sigmam = as<NumericVector>(lsigma[m]);
        }
        
        prodprob *= EMdensityChoice_cpp(convtolist_cpp(lbeta[m], lnbeta[m]),
                                    sigmam,
                                    i,
                                    mk(l,m),
                                    as<IntegerVector>(lnbeta[m]),
                                    as<NumericMatrix>(lA[m]),
                                    as<NumericMatrix>(lY[m]),
                                    as<double>(lymin[m]),
                                    as<double>(lymax[m]),
                                    mTCOV,
                                    ldelta,
                                    lnw[m],
                                    lnutmp,
                                    nnutmp,
                                    model[m]);
      }
      s += mprob(l, i + lA.size())*prodprob;
    }
    a += log(s);
  }
  
  return(a);
}

// ----------------------------------------------------------------------------
//  find parameters for the CNORM model
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List EMCNORMparam(arma::vec beta, 
                  int nw, 
                  int ng,
                  int n,
                  int period,
                  int j,
                  arma::vec nbeta,
                  arma::mat Y,
                  List lA,
                  arma::mat taux,
                  arma::vec nbetacum,
                  Nullable<List> lTCOV,
                  Nullable<List> ldelta){
  
  vec newbeta;
  mat A = lA[j];
  double b = 0;
  double stauxk = 0;
  for (int k = 0; k < ng; ++k){
    double staux = 0;
    rowvec a(nbeta[k]);
    a.fill(0);
    mat Ai(nbeta[k], period);
    for (int i = 0; i < n; ++i){
      vec vtmp1 = taux.col(i + lA.size());
      uvec indtmp = find(taux.col(j) == k);
      double sigtaux = sum(vtmp1(indtmp));
      staux += sigtaux;
      for (int t = 0; t < period; ++t){
        for (int kk = 0; kk < nbeta[k]; ++kk){
          Ai(kk, t) = pow(A(i, t), kk);
        }
      }
      rowvec vtmp = Y.row(i);
      a += sigtaux*vtmp*trans(Ai);
      rowvec mtmp = vtmp - trans(beta.subvec(nbetacum[k], nbetacum[k+1]-1))*Ai;
      b += sigtaux*as_scalar(mtmp*trans(mtmp));
    }
    stauxk += staux;
    newbeta = join_cols(newbeta, trans(a*inv(Ai*trans(Ai))/staux));
  }
  vec newsigma(ng);
  for (int i = 0; i < ng; ++i){
    newsigma(i) = sqrt(b / (period * stauxk));
  }
  
  List lsol;
  lsol.push_back(newbeta);
  lsol.push_back(newsigma);
  
  return(lsol);
}
// ----------------------------------------------------------------------------
//  find parameters for the LOGIT model
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List EMLOGITparam(arma::vec beta, 
                  int nw, 
                  int ng,
                  int n,
                  int period,
                  int j,
                  arma::vec nbeta,
                  arma::mat Y,
                  List lA,
                  arma::mat taux,
                  arma::vec nbetacum,
                  Nullable<List> lTCOV,
                  Nullable<List> ldelta){
  
  vec newbeta;
  mat A = lA[j];
  for (int k = 0; k < ng; ++k){
    //double staux = 0;
    vec newbetaIRLS;
    vec betaIRLS = beta.subvec(nbetacum[k], nbetacum[k+1]-1);
    vec precIRLS(betaIRLS.n_elem);
    precIRLS.fill(1);
    int stop = 0;
    
    while(all(abs(precIRLS) > 0.000001) &&  stop < 300){
      stop +=1;
      
      mat Aw(nbeta[k], n*period);
      vec tmpr(n*period);
      vec S(n*period);
      vec tmp2;
      for (int i = 0; i < n; ++i){
        vec vtmp1 = taux.col(i + lA.size());
        uvec indtmp = find(taux.col(j) == k);
        double sigtaux = sum(vtmp1(indtmp));
       // staux += sigtaux;
        vec vtmp(period);
        for (int t = 0; t < period; ++t){
          vec vtmp2(nbeta[k]);
          for (int po = 0; po < nbeta[k]; ++po){
            vtmp2(po) = pow(A(i, t), po);
          }
          double betaAit = as_scalar(trans(betaIRLS)*vtmp2);
          double rhoikt = exp(betaAit)/(1+exp(betaAit));
          for (int kk = 0; kk < nbeta[k]; ++kk){
            Aw(kk, period*i + t) = pow(A(i, t), kk);
          }
          tmpr(period*i + t) = rhoikt*(1-rhoikt);
          S[period*i + t] = betaAit + (Y(i, t) - rhoikt)/(rhoikt*(1-rhoikt));
          vtmp(t) = sigtaux;
        }
        tmp2 = join_cols(tmp2, vtmp);
      }

      mat ZW = diagmat(tmp2 % tmpr);
      mat Q;
      mat R;
      qr(Q, R, trans(Aw*sqrt(ZW)));
      mat Rc = R.submat(0, 0, nbeta[k] - 1, nbeta[k] - 1);
      mat Qc = Q.submat(0, 0, Q.n_rows - 1, nbeta[k] - 1);
      newbetaIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*S);
      
      precIRLS = betaIRLS - newbetaIRLS;
      betaIRLS = newbetaIRLS;
    }
    newbeta = join_cols(newbeta, betaIRLS);
  }
  
  List lsol;
  lsol.push_back(newbeta);
  
  return(lsol);
}
// ----------------------------------------------------------------------------
//  find parameters for the LOGIT model
// ----------------------------------------------------------------------------
// List EMZIPparam(arma::vec beta, 
//                 arma::vec nu,
//                 int nw, 
//                 int ng,
//                 int n,
//                 int period,
//                 int j,
//                 arma::vec nbeta,
//                 arma::vec nnu,
//                 arma::mat Y,
//                 List lA,
//                 arma::mat taux,
//                 arma::vec nbetacum,
//                 arma::vec nnucum,
//                 Nullable<IntegerVector> ndeltacuminit,
//                 Nullable<List> lTCOVinit,
//                 Nullable<List> ldeltainit){
//   
//   List lTCOV(lTCOVinit);  
//   List ldelta(ldeltainit);
//   IntegerVector ndeltacum(ndeltacuminit);
//   if (lTCOVinit.isNotNull()){
//     lTCOV= lTCOVinit;
//     ldelta = ldeltainit;
//     ndeltacum = ndeltacuminit;
//   }
//   
//   NumericVector pi;
//   pi.fill(0);
//   vec newbeta;
//   vec newnu;
//   mat A = lA[j];
//   for (int k = 0; k < ng; ++k){
//     vec newbetaIRLS;
//     vec betaIRLS = beta.subvec(nbetacum[k], nbetacum[k+1]-1);
//     vec newnuIRLS(nnu[k]);
//     vec nuIRLS = nu.subvec(nnucum[k], nnucum[k+1]-1);
//     vec precIRLS(betaIRLS.n_elem + nuIRLS.n_elem);
//     precIRLS.fill(1);
//     int stop = 0;
//     
//     while(all(abs(precIRLS) > 0.000001) &&  stop < 300){
//       stop +=1;
//       mat Aw(nnu[k], n*period);
//       mat Awp(nbeta[k], n*period);
//       vec Sn(n*period);
//       vec Sw(n*period);
//       vec Sp(n*period);
//       vec W(n*period);
//       vec Wp(n*period);
//       vec Z(n*period);
//       
//       for (int i = 0; i < n; ++i){
//         vec vtmp1 = taux.col(i + lA.size());
//         uvec indtmp = find(taux.col(j) == k);
//         double sigtaux = sum(vtmp1(indtmp));
//         for (int t = 0; t < period; ++t){
//           double sikt = fSikt_cpp(pi, NumericVector(beta.begin(), beta.end()), NumericVector(nu.begin(), nu.end()), 
//                                   k, i, t, IntegerVector(nbeta.begin(), nbeta.end()), IntegerVector(nnu.begin(), nnu.end()),
//                                   n, wrap(A), wrap(Y), wrap(lTCOV[j]), wrap(ldelta[j]), nw,
//                                   ndeltacum, period, 
//                                   IntegerVector(nbetacum.begin(), nbetacum.end()),IntegerVector(nnucum.begin(), nnucum.end()));
//           double nuAit = 0;
//           for (int po = 0; po < nnu[k]; ++po){
//             nuAit += pow(A(i, t), po)*nuIRLS[po];
//             Aw(po, period*i + t) = pow(A(i, t), po);
//           }
//           
//           double rhoikt = exp(nuAit)/(1+exp(nuAit));
//           W[period*i + t] = rhoikt*(1-rhoikt);
//           Sn[period*i + t] = nuAit + (sikt - rhoikt)/(rhoikt*(1-rhoikt));
//           
//           double betaAit = 0;
//           for (int po = 0; po < nbeta[k]; ++po){
//             betaAit += pow(A(i, t), po)*betaIRLS[po];
//             Awp(po, period*i + t) = pow(A(i, t), po);
//           }
//           
//           double lambdaikt = exp(betaAit);
//           Sp[period*i + t] = 1-sikt;
//           Sw[period*i + t] = betaAit+Y(i,t)/lambdaikt-1;
//           Wp[period*i + t] = lambdaikt;
//           Z[period*i + t] = sigtaux;
//         }
//       }
//       List lsol2;
//       lsol2.push_back(Sp);
//       
//       return(lsol2);
//       mat Q;
//       mat R;
//       mat ZW = diagmat(Z % W);
//       qr(Q, R, trans(Aw*sqrt(ZW)));
//       
//       mat Rc = R.submat(0, 0, nnu[k] - 1, nnu[k] - 1);
//       mat Qc = Q.submat(0, 0, Q.n_rows - 1, nnu[k] - 1);
//       newnuIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*Sn);
//       
//       mat SpZWp = diagmat(Sp % Z % Wp);
//       qr(Q, R, trans(Awp*sqrt(SpZWp)));
//       Rc = R.submat(0, 0, nbeta[k] - 1, nbeta[k] - 1);
//       Qc = Q.submat(0, 0, Q.n_rows - 1, nbeta[k] - 1);
//       newbetaIRLS = solve(Rc, trans(Qc)*sqrt(SpZWp)*Sw);
//     
//       precIRLS = join_cols(betaIRLS, nuIRLS) - join_cols(newbetaIRLS, newnuIRLS);
//       betaIRLS = newbetaIRLS;
//       nuIRLS = newnuIRLS;
//     }
//     newbeta = join_cols(newbeta, betaIRLS);
//     newnu = join_cols(newnu, nuIRLS);
//   }
//   
//   List lsol;
//   lsol.push_back(newbeta);
//   lsol.push_back(newnu);
//   
//   return(lsol);
// }

// [[Rcpp::export]]
List EMZIPparam(arma::vec beta, 
                arma::vec nu,
                int nw, 
                int ng,
                int n,
                int period,
                int j,
                arma::vec nbeta,
                arma::vec nnu,
                arma::mat Y,
                List lA,
                arma::mat taux,
                arma::vec nbetacum,
                arma::vec nnucum,
                Nullable<IntegerVector> ndeltacum,
                Nullable<List> lTCOVinit,
                Nullable<List> ldeltainit){
  
  
  NumericMatrix  TCOV;
  NumericVector delta;
  
  NumericVector pi;
  pi.fill(0);
  vec newbeta;
  vec newnu;
  mat A = lA[j];
  for (int k = 0; k < ng; ++k){
    vec newbetaIRLS;
    vec betaIRLS = beta.subvec(nbetacum[k], nbetacum[k+1]-1);
    vec newnuIRLS(nnu[k]);
    vec nuIRLS = nu.subvec(nnucum[k], nnucum[k+1]-1);
    vec precIRLS(betaIRLS.n_elem + nuIRLS.n_elem);
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
        vec vtmp1 = taux.col(i + lA.size());
        uvec indtmp = find(taux.col(j) == k);
        double sigtaux = sum(vtmp1(indtmp));
        for (int t = 0; t < period; ++t){
          double sikt = fSikt_cpp(pi, NumericVector(beta.begin(), beta.end()), NumericVector(nu.begin(), nu.end()),
                                  k, i, t, IntegerVector(nbeta.begin(), nbeta.end()), IntegerVector(nnu.begin(), nnu.end()),
                                  n, wrap(A), wrap(Y), TCOV, delta, nw,
                                  ndeltacum, period,
                                  IntegerVector(nbetacum.begin(), nbetacum.end()),IntegerVector(nnucum.begin(), nnucum.end()));
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
          Z[period*i + t] = sigtaux;
        }
      }
      
      mat Q;
      mat R;
      mat ZW = diagmat(Z % W);
      qr(Q, R, trans(Aw*sqrt(ZW)));
      
      mat Rc = R.submat(0, 0, nnu[k] - 1, nnu[k] - 1);
      mat Qc = Q.submat(0, 0, Q.n_rows - 1, nnu[k] - 1);
      newnuIRLS = solve(Rc, trans(Qc)*sqrt(ZW)*Sn);
      
      mat SpZWp = diagmat(Sp % Z % Wp);
      qr(Q, R, trans(Awp*sqrt(SpZWp)));
      Rc = R.submat(0, 0, nbeta[k] - 1, nbeta[k] - 1);
      Qc = Q.submat(0, 0, Q.n_rows - 1, nbeta[k] - 1);
      newbetaIRLS = solve(Rc, trans(Qc)*sqrt(SpZWp)*Sw);
      
      precIRLS = join_cols(betaIRLS, nuIRLS) - join_cols(newbetaIRLS, newnuIRLS);
      betaIRLS = newbetaIRLS;
      nuIRLS = newnuIRLS;
    }
    newbeta = join_cols(newbeta, betaIRLS);
    newnu = join_cols(newnu, nuIRLS);
  }
  
  
  List lsol;
  lsol.push_back(newbeta);
  lsol.push_back(newnu);
  
  return(lsol);
}

// ----------------------------------------------------------------------------
//  EM 
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List EMMult_cpp(List lparam,
                List lng,
                List lnx,
                List lnbeta,
                List lnnu,
                List ln,
                List lA,
                List lY,
                List lX,
                List lymin,
                List lymax,
                Nullable<List> lTCOVinit,
                List lnw,
                arma::mat mk,
                arma::vec vp,
                std::vector<std::string> model,
                int itermax,
                bool EMIRLS){
  
  List lTCOV(lTCOVinit);
  if (lTCOVinit.isNotNull()){
    lTCOV= lTCOVinit;
  }
  double prec = 0.000001;
  int ilparam = lparam.size();
  List ltheta;
  List lbeta;
  List lsigma(ilparam - 1);
  List lnu(ilparam - 1);
  List ldelta;
  mat mpsi;
  vec vthetapsi;
  int nbtheta = 0;
  
  for (int nl = 0; nl < ilparam - 1; ++nl){
    vec vtmp = lparam[nl];
    int ilng = lng[nl];
    int ilnx = lnx[nl];
    //   vec  z = zeros<vec>(1);
    //ltheta.push_back(join_cols(z, vtmp.subvec(0, (ilng - 1)*ilnx - 1)));
    ltheta.push_back(vtmp.subvec(0, (ilng - 1)*ilnx - 1));
    nbtheta += (ilng - 1)*ilnx;
    vec lnbetatmp = lnbeta[nl];
    lbeta.push_back(vtmp.subvec((ilng - 1)*ilnx, (ilng - 1)*ilnx + sum(lnbetatmp) - 1));
    if (model[nl] == "CNORM"){
      lsigma[nl] = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + ilng - 1);  
    }
    if (model[nl] ==  "ZIP"){
      vec lnnutmp = lnnu[nl];
      lnu[nl] = vtmp.subvec((ilng - 1)*ilnx + sum(lnbetatmp), (ilng - 1)*ilnx + sum(lnbetatmp) + sum(lnnutmp) - 1);
    }
    vthetapsi = join_cols(vthetapsi, vtmp.subvec(0, (ilng - 1)*ilnx - 1));
  }
  vec vtmp = lparam[ilparam - 1];
  mpsi = mPsi_cpp(vtmp, lng);
  
  vthetapsi = join_cols(vthetapsi, vtmp);
  
  List ltmp;
  for (int i = 0; i < lng.size(); ++i){
    int itmp = lng[i];
    ltmp.push_back(regspace(0, 1, itmp - 1));
  }
  
  int tour = 1;
  while (tour < itermax){
    // E-step
    List ltmp = ftauxPiikMult_cpp(ltheta,
                                  mpsi,
                                  lbeta,
                                  lsigma,
                                  lng,
                                  lnbeta,
                                  ln,
                                  lA,
                                  lY,
                                  lymin,
                                  lymax,
                                  lTCOVinit,
                                  ldelta,
                                  lnw,
                                  lnx,
                                  lX,
                                  mk,
                                  vp,
                                  lnu,
                                  lnnu,
                                  model);
    mat taux = ltmp[0];
    mat mprob = ltmp[1];
    double like = likelihoodMultEM_cpp(ltheta,
                                       mpsi,
                                       lbeta,
                                       lsigma,
                                       lng,
                                       lnx,
                                       lnbeta,
                                       ln,
                                       lA,
                                       lY,
                                       lX,
                                       lymin,
                                       lymax,
                                       lTCOVinit,
                                       ldelta,
                                       lnw,
                                       mk,
                                       mprob,
                                       vp,
                                       lnu,
                                       lnnu,
                                       model);
    
    Rprintf("iter %3d value ", tour);
    Rprintf("%.6f\n", -like);
    
    List lnewbeta;
    List lnewsigma(ilparam - 1);
    List lnewnu(ilparam - 1);
    
    for (int j = 0; j < lA.size(); ++j){
      mat mA = lA[j];
      int periodj = mA.n_cols;
      
      vec betaj = lbeta[j];
      vec nbetaj = lnbeta[j];
      
      vec  z = zeros<vec>(1);
      vec nbetacumj = join_cols(z, nbetaj);
      nbetacumj = cumsum(nbetacumj);
      
      vec newbeta;
      vec newsigma;
      vec newnu;
      
      if (model[j] == "CNORM"){
        List ltmp = EMCNORMparam(betaj, lnw[j], lng[j], ln[j], periodj, j,
                                 nbetaj, lY[j], lA, taux, nbetacumj, lTCOV, ldelta);
        newbeta = as<arma::vec>(ltmp[0]);
        lnewbeta.push_back(newbeta);
        newsigma = as<arma::vec>(ltmp[1]);
        lnewsigma(j) = newsigma;
      }else if (model[j] == "LOGIT"){
        List ltmp =  EMLOGITparam(betaj, lnw[j], lng[j], ln[j], periodj, j,
                                  nbetaj, lY[j], lA, taux, nbetacumj, lTCOV, ldelta);
        newbeta = as<arma::vec>(ltmp[0]);
        lnewbeta.push_back(newbeta);
      }else if (model[j] == "ZIP"){
        vec nuj = lnu[j];
        vec nnuj = lnnu[j];
        vec nnucumj = join_cols(z, nnuj);
        nnucumj = cumsum(nnucumj);
        IntegerVector ndeltacum;
        List ltmp = EMZIPparam(betaj, nuj, lnw[j], lng[j], ln[j], periodj, j,
                               nbetaj, nnuj, lY[j], lA, taux, nbetacumj, nnucumj, 
                               ndeltacum, lTCOV, ldelta);
        newbeta = as<arma::vec>(ltmp[0]);
        lnewbeta.push_back(newbeta);
        newnu = as<arma::vec>(ltmp[1]);
        lnewnu(j) = newnu;
      }
    }
    
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];
    List tmp = optim(Rcpp::_["par"] = vthetapsi,
                     Rcpp::_["fn"] = Rcpp::InternalFunction(&ThetaPsiPiikMult_cpp),
                     Rcpp::_["gr"] = Rcpp::InternalFunction(&difThetaPsiPiikMult_cpp),
                     Rcpp::_["lng"] = lng,
                     Rcpp::_["lnx"] = lnx,
                     Rcpp::_["lnbeta"] = lnbeta,
                     Rcpp::_["ln"] = ln,
                     Rcpp::_["lA"] = lA,
                     Rcpp::_["lY"] = lY,
                     Rcpp::_["lX"] = lX,
                     Rcpp::_["lymin"] = lymin,
                     Rcpp::_["lymax"] = lymax,
                     Rcpp::_["lTCOVinit"] = lTCOV,
                     Rcpp::_["lnw"] = lnw,
                     Rcpp::_["taux"] = taux,
                     Rcpp::_["mk"] = mk,
                     Rcpp::_["hessian"] =  0,
                     Rcpp::_["method"] = "BFGS",
                     Rcpp::_["control"] = List::create(Named("fnscale")=-1)
    );
    vec newvthetapsi =  Rcpp::as<arma::vec>(tmp[0]);
    
    vec param = vthetapsi;
    vec newparam = newvthetapsi;
    
    int indmin = 0;
    int indmax = -1;
    for (int i = 0; i < lbeta.size(); ++i){
      vec vtmp1 = lbeta[i];
      param = join_cols(param, vtmp1);
      vec vtmp2 = lnewbeta[i];
      newparam = join_cols(newparam, vtmp2);
      //update ltheta;
      mat mtmp = lX[i];
      int lngtmp = lng[i];
      int lnx = mtmp.n_cols;
      indmax += lnx*(lngtmp - 1);
      vec vtmp3 = vthetapsi.subvec(indmin, indmax);
      ltheta[i] = vtmp3;
      indmin = indmax + 1;
    }
    
    for (int i = 0; i < lsigma.size(); ++i){
      if (model[i] == "CNORM"){
        vec vtmp1 = lsigma[i];
        param = join_cols(param, vtmp1);
        vec vtmp2 = lnewsigma[i];
        newparam = join_cols(newparam, vtmp2);
      }
      if (model[i] == "ZIP"){
        vec vtmp1 = lnu[i];
        param = join_cols(param, vtmp1);
        vec vtmp2 = lnewnu[i];
        newparam = join_cols(newparam, vtmp2);
      }
    }
    
    vec vcomp(newparam.size());
    vcomp.fill(prec);
    if (all(abs(newparam-param) < vcomp)){
      tour = itermax + 2;
    }
    
    lbeta = lnewbeta;
    lsigma = lnewsigma;
    lnu = lnewnu;
    vthetapsi = newvthetapsi;
    mpsi = mPsi_cpp(vthetapsi.subvec(nbtheta, vthetapsi.size() - 1), lng);
    tour++;
  }
  
  List sparam;
  int indmin = 0;
  int indmax = -1;
  int ng = lng[0];
  for (int i = 0; i < ng; ++i){
    mat mtmp = lX[i];
    int lngtmp = lng[i];
    int lnx = mtmp.n_cols;
    indmax += lnx*(lngtmp - 1);
    
    vec vtmp1 = vthetapsi.subvec(indmin, indmax);
    vec vtmp2 = lbeta[i];
    vec vtmp3;
    if (model[i] == "CNORM"){
      vtmp3 = as<arma::vec>(lsigma[i]);
    }
    if (model[i] == "ZIP"){
      vtmp3 = as<arma::vec>(lnu[i]);
    }
    sparam.push_back(join_cols(vtmp1, vtmp2, vtmp3));
    indmin = indmax + 1;
  }
  sparam.push_back(vthetapsi.subvec(indmin, vthetapsi.size() - 1));
  
  return(sparam);
}
