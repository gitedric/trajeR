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
      Rprintf("%.6f\n", -likelihoodZIP_cpp(NumericVector(vparam.begin(), vparam.end()), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw));
    }
    // E-step
    NumericMatrix zk = ftauxZIP_cpp(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X);
    
    
    NumericVector newbeta;
    NumericVector newnu;
    NumericVector newdelta;
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];
    
    vec vtmp;
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
        NumericVector betadeltak = betak;
        for (int i = 0; i < deltak.size(); i++){
          betadeltak.push_back(deltak[i]);
        }
        
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
        List tmp = optim(Rcpp::_["par"] = betadeltak,
                         Rcpp::_["fn"] = Rcpp::InternalFunction(&QbetadeltakZIP_cpp),
                         Rcpp::_["gr"] = Rcpp::InternalFunction(& difQbetadeltakZIP_cpp),
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
                         Rcpp::_["nw"] = nw,
                         Rcpp::_["hessian"] =  0,
                         Rcpp::_["control"] = List::create(Named("fnscale")=-1)
        );
        vtmpnu = join_cols(vtmpnu, as<arma::vec>(tmpnu[0]));
        NumericVector vtmp;
        vtmp = tmp[0];
        for (int i = 0 ; i < nbeta[k]; ++i){
          newbeta.push_back(vtmp[i]);
        }
        for (int i = 0 ; i < nw; ++i){
          newdelta.push_back(vtmp[i+nbeta[k]-1]);
        }
      }
      newnu = NumericVector(vtmpnu.begin(), vtmpnu.end());
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
