#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]
using namespace Rcpp;
using namespace arma;

// factorial function
double facto(double nb){
  double res =1;
  while (nb>1){
    res *= nb;
    nb--;
  }
  return(res);
}

double piik_cpp(NumericVector theta,
                int i,
                int k,
                int ng,
                NumericMatrix X){
  int ntheta = X.ncol();
  
  i = i-1;
  k = k-1;
  
  NumericVector vtmp;
  for (int s = 0; s < ng; ++s){
    double tmp = 0;
    for (int ind = 0; ind < ntheta; ++ind){
      tmp += theta[s*ntheta+ind]*X(i,ind); 
    }
    vtmp.push_back(tmp);
  }
  // to avoid overflow due to exp function
  vtmp = exp(vtmp -  max(vtmp));
  return(vtmp[k]/sum(vtmp));
}

// [[Rcpp::export]]
double piikIntern_cpp(NumericVector theta,
                      int i,
                      int k,
                      int ng,
                      NumericMatrix X){
  int ntheta = X.ncol();
  NumericVector vtmp;
  for (int s = 0; s < ng; ++s){
    double tmp = 0;
    for (int ind = 0; ind < ntheta; ++ind){
      tmp += theta[s*ntheta+ind]*X(i,ind); 
    }
    vtmp.push_back(tmp);
  }
  // to avoid overflow due to exp function
  vtmp = exp(vtmp -  max(vtmp));
  return(vtmp[k]/sum(vtmp));
}
// product of all elements of a vector
double prodvect(NumericVector vec){
  double res = 1;
  for (int i = 0; i < vec.size(); i++){
    res *= vec[i];
  }
  return(res);
}
// ----------------------------------------------------------------------------
// Wit
// ----------------------------------------------------------------------------
double Wit_cpp(Nullable<NumericMatrix> TCOV,
               int period,
               Nullable<List> delta, 
               int nw,
               int i,
               int t,
               int k){
  if (nw == 0){
    return(0);
  }else{
    double a = 0;
    NumericMatrix mTCOV(TCOV.get());
    List deltaL(delta.get());
    NumericVector vtmp = deltaL[k];
    for (int s = 0; s < nw; ++s){
      a += vtmp[s]*mTCOV(i, t + s*period);
    }
    return(a);
  }
}
// compute muikt
// [[Rcpp::export]]
NumericVector muikt_cpp(NumericVector beta,
                        int nbeta,
                        int i,
                        int period,
                        NumericMatrix A, 
                        Nullable<NumericMatrix> TCOV,
                        Nullable<List> delta,
                        int nw,
                        int k){
  NumericVector muikt;
  for (int s = 0; s < period; ++s){
    NumericVector vtmp2;
    for (int po = 0; po < nbeta; ++po){
      vtmp2.push_back(pow(A(i,s), po));
    }
    muikt.push_back(sum(beta*vtmp2)+Wit_cpp(TCOV, period, delta, nw, i, s, k));
  }
  return(muikt);
}
// ----------------------------------------------------------------------------
// Wit for EM algorithm
// ----------------------------------------------------------------------------
double WitEM_cpp(Nullable<NumericMatrix> TCOV,
               int period,
               Nullable<NumericVector> delta, 
               int nw,
               int i,
               int t,
               int k){
  if (nw == 0){
    return(0);
  }else{
    double a = 0;
    NumericVector deltak(delta.get());
    NumericMatrix mTCOV(TCOV.get());
    for (int s = 0; s < nw; ++s){
      a += deltak[s]*mTCOV(i, t + s*period);
    }
    return(a);
  }
}
// subset of a matrix that matches a logical statement
NumericMatrix submat_cpp(NumericMatrix X, LogicalVector condition) { 
  int n=X.nrow(), k=X.ncol();
  NumericMatrix out(n, sum(condition));
  for (int i = 0, j = 0; i < k; i++) {
    if(condition[i]) {
      out(_, j) = X(_, i);
      j = j+1;
    }
  }
  return(out);
}
// [[Rcpp::export]]
double ftheta_cpp(NumericVector theta,
                  NumericMatrix taux,
                  NumericMatrix X,
                  int n,
                  int ng,
                  int period) {
  
  double a = 0;
  int nx = X.cols();
  NumericVector thetatmp (nx);
  NumericVector tmp (ng);
  
  for (int i = 0; i < n; i++){
    for (int k = 0; k < ng; k++){
       for (int s = 0; s < nx; s++){
         thetatmp[s] = theta[k*nx + s];
       }
       tmp[k] = sum(thetatmp*X(i, _));
    }
    for (int k = 0; k < ng; k++){
      a += taux(i, k)*(tmp[k]-log(sum(exp(tmp))));
    }
  }
  return(a);
}
// [[Rcpp::export]]
NumericVector difftheta_cpp(NumericVector theta,
                            NumericMatrix taux,
                            NumericMatrix X,
                            int n,
                            int ng,
                            int period){
  
  int nx = X.cols();
  NumericVector thetas (nx*ng);
  int indt = 0;
  double a;
  NumericVector thetatmp (nx);
  NumericVector tmp (ng);
  
  for (int k = 0; k < ng; k++){
    for (int l = 0; l < nx; l++){
      a = 0;
      for (int i = 0; i < n; i++){
        for (int kt = 0; kt < ng; kt++){
          for (int s = 0; s < nx; s++){
            thetatmp[s] = theta[kt*nx + s];
          }
          tmp[kt] = exp(sum(thetatmp*X(i, _)));
        }  
        a += X(i,l)*(taux(i,k)-tmp[k]/sum(tmp));
      }
      thetas[indt] = a;
      indt += 1;
    }
  }
  return(thetas);
}
// [[Rcpp::export]]
Rcpp::NumericVector thethaIRLS_cpp(Rcpp::NumericVector thetaIRLS,
                                   int n,
                                   int ng,
                                   Rcpp::NumericMatrix X,
                                   arma::mat taux,
                                   int refgr){

  int nx = X.cols();
  //Rcpp::Function piik("piik");
  int stop = 0;
  arma::vec precIRLS(thetaIRLS.length());
  precIRLS.fill(1);

  while(all(abs(precIRLS) > 0.000001) &&  stop < 300){
    stop +=1;

    Rcpp::NumericVector tmpPiik(1);
    Rcpp::NumericVector tmpPiil(1);
    Rcpp::NumericVector tmp4(0);
    Rcpp::IntegerVector ind = Rcpp::seq(0,ng-1);
    arma::mat PIw (n*(ng-1),n*(ng-1));
    arma::mat Xng(n*(ng-1), nx*(ng-1));
    arma::vec tmp2;
    Rcpp::NumericVector tmp5(1);

    int kind = 0;
    ind.erase(refgr-1);
    Xng.zeros();
    PIw.zeros();

    for (auto k = ind.begin(); k != ind.end(); ++k){
      arma::mat PIwtmp(n, (ng-1)*n);
      int lind = 0;
      for (auto l = ind.begin(); l != ind.end(); ++l){
        Rcpp::NumericVector vtmp(nx);
        Rcpp::NumericVector tmp1(0);
        for (int it = 0; it < (ng-1)*nx; ++it){
          vtmp.push_back(thetaIRLS[it]);
        }
        if (*k == *l){
          for (int i = 0; i < n; i++){
            tmpPiik = piik_cpp(vtmp,i+1, *k+1, ng, X);
            tmp1.push_back(tmpPiik[0]*(1-tmpPiik[0]));
            tmp4.push_back(tmpPiik[0]);
          }
        }else{
          for (int i = 0; i < n; ++i){
            tmpPiik = piik_cpp(vtmp, i+1, *k+1, ng, X);
            tmpPiil = piik_cpp(vtmp, i+1, *l+1, ng, X);
            tmp1.push_back(-tmpPiik[0]*tmpPiil[0]);
          }
        }
        PIwtmp.submat(0,lind*n,n-1,(lind+1)*n-1) = arma::diagmat( Rcpp::as<arma::vec>(tmp1));
        lind += 1;
      }
      PIw.submat(kind*n,0,(kind+1)*n-1,(ng-1)*n-1) = PIwtmp;
      Xng.submat(kind*n, kind*nx, (kind+1)*n-1, (kind+1)*nx-1) =  Rcpp::as<arma::mat>(X);
      tmp2 = arma::join_cols(tmp2,taux.col(*k));
      kind += 1;
    }
    arma::colvec Z = tmp2;
    arma::colvec PIm = tmp4;
    arma::vec newthetaIRLS = solve(trans(Xng)*PIw*Xng, trans(Xng)*(PIw*Xng*(Rcpp::as<arma::vec>(thetaIRLS))+Z-PIm));

    precIRLS = Rcpp::as<arma::vec>(thetaIRLS)-newthetaIRLS;
    thetaIRLS = Rcpp::NumericVector(newthetaIRLS.begin(), newthetaIRLS.end());
  }

  return(Rcpp::NumericVector(thetaIRLS.begin(), thetaIRLS.end()));
}
//
//
//
// [[Rcpp::export]]
NumericVector findtheta_cpp(NumericVector theta, 
                           NumericMatrix taux, 
                           NumericMatrix X, 
                           int n, 
                           int ng, 
                           int nx, 
                           int period, 
                           bool EMIRLS, 
                           int refgr){
  NumericVector newtheta(nx);
  if (EMIRLS){
    NumericVector thetaIRLS = theta;
    NumericVector tmp;
    for (int i = (refgr-1)*nx; i < refgr*nx; ++i){
      tmp.push_back(theta[i]);
    }
    tmp = rep(tmp, ng-1);
    thetaIRLS.erase((refgr-1)*nx, refgr*nx);
    NumericVector newthetaIRLS = thethaIRLS_cpp(thetaIRLS-tmp, n, ng,  X, as<mat>(taux), refgr);
    for (int i = 0; i < (ng-1)*nx; ++i){
      newtheta.push_back(newthetaIRLS(i));
    }
  }else{
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];
    List tmp = optim(Rcpp::_["par"] = theta,
                     Rcpp::_["fn"] = Rcpp::InternalFunction(&ftheta_cpp),
                     Rcpp::_["gr"] = Rcpp::InternalFunction(&difftheta_cpp),
                     Rcpp::_["taux"] = taux,
                     Rcpp::_["X"] = X,
                     Rcpp::_["n"] = n,
                     Rcpp::_["ng"] = ng,
                     Rcpp::_["period"] = period,
                     Rcpp::_["hessian"] =  0,
                     Rcpp::_["control"] = List::create(Named("fnscale")=-1)
    );
    newtheta = tmp[0];
  }
  return(newtheta);
}