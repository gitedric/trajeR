// ----------------------------------------------------------------------------
  // Q function for delta
// ----------------------------------------------------------------------------
  // [[Rcpp::export]]
// double QdeltaLOGIT_cpp(NumericVector beta,
                          //                       NumericMatrix taux,
                          //                       int n,
                          //                       int ng,
                          //                       IntegerVector nbeta,
                          //                       NumericMatrix A,
                          //                       NumericMatrix Y,
                          //                       Nullable<NumericMatrix> TCOV,
                          //                       Nullable<NumericVector> delta,
                          //                       int nw){
  //   double a = 0;
  //   for (int k = 0; k < ng; ++k){
    //     a += QdeltakLOGIT_cpp(beta, taux, k, n, ng, nbeta, A, Y, TCOV, delta, nw);
    //   }
  //   return(a);
  // }
// ----------------------------------------------------------------------------
  // differential of Q delta
// ----------------------------------------------------------------------------
  // [[Rcpp::export]]
// NumericVector difQdeltaLOGIT_cpp(NumericVector beta,
                                    //                                  NumericMatrix taux,
                                    //                                  int n,
                                    //                                  int ng,
                                    //                                  IntegerVector nbeta,
                                    //                                  NumericMatrix A,
                                    //                                  NumericMatrix Y,
                                    //                                  Nullable<NumericMatrix> TCOV,
                                    //                                  Nullable<NumericVector> delta,
                                    //                                  int nw){
  //   vec res;
  //   for (int k = 0; k < ng; ++k){
    //     NumericVector tmp = difQdeltakLOGIT_cpp(delta, taux, k, n, ng, nbeta, A, Y, TCOV, delta, nw);
    //     res = join_cols(res, as<arma::vec>(tmp));
    //   }
  //   return(Rcpp::NumericVector(res.begin(), res.end()));
  // }