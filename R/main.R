my_env <- new.env(parent = emptyenv())
#################################################################################
# find parameters for a given method
#################################################################################
#' Fitting longitudinal mixture models
#'
#' \code{trajeR} is used to fit longitudinal mixture models. It used 3 types of mixture models : LOGIT, ZIP and censored Normal.
#'
#' Models for trajeR is, by default, a polynomial regression of the time value parameters for each groups. The number fo group is controlled by the integer \code{ng}.
#' We can spcecify the degre of the polynomial shape for each groups by the vector \code{degre}.
#'
#' @param Y Matrix. A matrix containing the variables in the model.
#' @param A Matrix. A matrix containing the time variable data.
#' @param Risk Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
#' with one column  with value 1.
#' @param TCOV Matrix. An optional matrix containing the time covariate that influence the trajectory themselves.
#' By default its value is NULL.
#' @param degre Vector of integer. The degree of every polynomial function.
#' @param degre.nu Vector of integer. The degree of all Poisson part for a ZIP model.
#' @param degre.phi Vector of integer. The degree of  beta parametr for  a BETA model.
#' @param Model String. The model used. The value are LOGIT for a Logit Mixture model,
#'CNORM for a Censored Normal Mixture Model or ZIP for Zero Inflated Poisson Mixture model.
#' @param Method String. Determine the method used for find the parameters of the model.
#' The value are L for the Maximum Likelihood Estimation, EM for Expectation Maximization method
#' with quasi newton method inside, EMIWRLS for Expectation Maximization method with Iterative
#' Weighted Least Square.
#' @param ssigma Logical. By default its value is FALSE. For the CNORM model,
#' indicate if we want the same sigma for all normal density function.
#' @param ymax Real. For the CNORM model, indicate the maximum value of the data. It concern only the model
#' with censored data. By default its value is the maximum value of the data plus 1.
#' @param ymin Real. For the CNORM model, indicate the minimum value of the data. It concern only the model
#' with censored data. By default its value is the maximum value of the data minus 1.
#' @param hessian Logical. Indicate if we want calculate the hessian matrix. Default is FALSE.
#' If the method use is Likelihood, the hessian is calculated by inverting the Information's Fisher Matrix.
#' To avoid numerically singular matrix we find the pseudo inverse matrix by using the \code{ginv} function int he package MASS.
#' If the method is EM or EMIWRLS, the hessian is calculated by using Louis method.
#' @param itermax Integer. Indicate the maximal number of iteration for \code{optim} function or for the EM algorithm.
#' @param paraminit Vector. The vector of initial parameters. By default \code{trajeR} calculate the initial value
#' based of the range or the standard deviation.
#' @param ProbIRLS Logical. Indicate the method to sue in the search of predictor's probability. If TRUE (by default) we use 
#' IRLS method and if FALSE we use optimization method.
#' @param refgr Integer. The number of reference group. By default is 1.
#' @param fct  Function. The definition of the function  f in the definition in nonlinear model.
#' @param diffct Function. The differential of the function f in the nonlinear model.
#' @param nbvar Integer. The number of variable in the nonlinear model.
#' @param  ng.nl Integer. The number of group for a non linear model.
#' @param  nls.lmiter Integer. In the case of non linear model, the maximum number of iterations allowed.
#' @return return an object of class "\code{Trajectory.LOGIT}".
#' The generic accessor functions \code{beta}, \code{delta}, \code{theta}, \code{sd}, \code{tab}, \code{Likelihood},
#' \code{ng}, \code{model} and \code{method} extract various useful features of the value returned by \code{trajeR}.
#'
#' An object of class "\code{Trajectory.LOGIT}" is a list containing at least the following components:
#'
#'\describe{
#' \item{\code{beta}}{a vector of the parameters beta.}
#' \item{\code{delta}}{a vector of the parameter delta. Only if we use time covariate.}
#' \item{\code{theta}}{a vector with the parameter theta if there exist a covariate X that modify
#'   the probability or the probability of group membership.}
#' \item{\code{sd}}{a vector of the standrad deviation of the parameters.}
#' \item{\code{tab}}{a matrix with all the parameters and standard deviation.}
#' \item{\code{Likelihood}}{a real with the Likelihhod obtnaied by the parameters.}
#' \item{\code{ng}}{a integer with the number of group.}
#' \item{\code{model}}{a string with the model used.}
#' \item{\code{method}}{a string with the method used.}
#'}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' load("data/dataNORM01.RData")
#' solL = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
#'               Model="CNORM", Method = "L", ssigma = FALSE, 
#'               hessian = TRUE)
#'               }
trajeR <- function(Y, A, Risk = NULL, TCOV = NULL, degre = NULL, degre.nu = 0, degre.phi = 0,
                   Model, Method = "L", ssigma = FALSE,
                   ymax = max(Y, na.rm = TRUE)+1, ymin = min(Y, na.rm = TRUE)-1, hessian = TRUE,
                   itermax = 100, paraminit = NULL, ProbIRLS = TRUE, refgr = 1,
                   fct = NULL, diffct = NULL, nbvar = NULL, ng.nl = NULL, nls.lmiter = 50){
  ng = ifelse(is.null(ng.nl), length(degre), ng.nl) 
  EMIRLS = ProbIRLS
  X = Risk
  if (is.null(degre)){
    degre = rep(nbvar, ng)
  }else{
    degre = degre + 1
  }
  degre.phi = degre.phi + 1
  n = nrow(Y)
  period = ncol(A)
  if (is.null(X)){
    X = cbind(rep(1, n))
    nx = ncol(X)
    colnames(X) = c("Intercept")
  }else{
    X = cbind(matrix(rep(1, n), ncol = 1), X)
    nx = ncol(X)
    if (is.null(colnames(X))){
      colnames(X) = c("Intercept", paste0("X",2:nx))
    }
  }
  if (is.null(TCOV)){
    nw = 0
    delta = NULL
  }else{
    nw = ncol(TCOV)/period
    delta = lapply(1:ng, function(s){
      rep(0, nw)
    })
  }
  ntheta = nx
  nu = NULL
  if (Model == "CNORM"){
    beta = lapply(1:(ng), function(s){
      c(stats::qnorm((2*(s-1)+1)/(2*ng),mean(Y, na.rm = TRUE),stats::sd(Y, na.rm = TRUE)), rep(0, degre[s]-1))
    })
    sigma = rep(stats::sd(Y, na.rm = TRUE), ng)
    paraff = c(unlist(beta), sigma)
  }else if (Model == "LOGIT"){
    beta = lapply(1:(ng), function(s){
      c(suppressWarnings(log(stats::qnorm((2*(s-1)+1)/(2*ng),mean(Y, na.rm = TRUE),stats::sd(Y, na.rm = TRUE))/(1-stats::qnorm((2*(s-1)+1)/(2*ng),mean(Y, na.rm = TRUE),stats::sd(Y, na.rm = TRUE))))), rep(0, degre[s]-1))
    })
    beta[is.na(beta)]=-5
    paraff = c(unlist(beta))
  }else if (Model =="ZIP"){
    nu = lapply(1:(ng), function(s){
      c(-3, rep(0, degre.nu[s]))
    })
    beta = lapply(1:(ng), function(s){
      c(log(stats::qgamma((2*(s-1)+1)/(2*ng), mean(Y, na.rm = T), mean(Y, na.rm = T)**2/(stats::sd(Y, na.rm = T)**2-mean(Y, na.rm = T)))), rep(0, degre[s]-1))
    })
    paraff = c(unlist(beta), unlist(nu))
  }else  if (Model == "BETA"){
    beta = lapply(1:(ng), function(s){
      c(suppressWarnings(max(log(stats::qnorm((2*(s-1)+1)/(2*ng),mean(Y, na.rm = TRUE),stats::sd(Y, na.rm = TRUE))/(1-stats::qnorm((2*(s-1)+1)/(2*ng),mean(Y, na.rm = TRUE),stats::sd(Y, na.rm = TRUE)))), -5)), rep(0, degre[s]-1))
    })
    beta[is.na(beta)]=-5
    phi = lapply(1:(ng), function(s){
      c(5, rep(0, degre.phi[s]-1))
    })
    paraff = c(unlist(beta), unlist(phi))
  }else if (Model =="POIS"){
    beta = lapply(1:(ng), function(s){
      c(log(stats::qgamma((2*(s-1)+1)/(2*ng), mean(Y, na.rm = T)**2, mean(Y, na.rm = T)**2/(stats::sd(Y, na.rm = T)**2-mean(Y, na.rm = T)))), rep(0, degre[s]-1))
    })
    paraff = c(unlist(beta))
  }else{
    beta = lapply(1:(ng), function(s){
      c(stats::qnorm((2*(s-1)+1)/(2*ng),mean(Y, na.rm = T),stats::sd(Y, na.rm = T)), rep(0, degre[s]-1))
    })
    paraff = c(unlist(beta), rep(stats::sd(Y, na.rm = T), ng))
  }
  nbeta = degre
  theta = rep(0, ng*ntheta)
  X = as.matrix(X)
  pi = sapply(1:ng, function(s){piik(theta = theta, i = 1, k = s, ng = ng, X = X)})
  cat('Starting Values\n')
  if (is.null(paraminit)){
    cat(c(pi, paraff, unlist(delta)))  
  }
  else{
    cat(paraminit)
  }
  cat('\n\n')
  cat("Likelihood\n")
  
  if (Model == "CNORM"){
    res = trajeR.CNORM(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, sigma, delta, pi, Method, ssigma,
                       ymax, ymin, hessian, itermax, paraminit, EMIRLS, refgr)
  }else if (Model == "LOGIT"){
    res = trajeR.LOGIT(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, delta, pi, Method,
                       hessian, itermax, paraminit, EMIRLS, refgr)
  }else if (Model == "ZIP"){
    res = trajeR.ZIP(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, degre.nu, theta, beta, nu, delta, pi, Method,
                     hessian, itermax, paraminit, EMIRLS, refgr)
  }else if (Model == "BETA"){
    res = trajeR.BETA(Y, A, X, TCOV, ng, nx, n, nbeta, degre.phi, nw, ntheta, period, degre, theta, beta, phi, delta, pi, Method,
                      hessian, itermax, paraminit, EMIRLS, refgr)
  }else if (Model == "POIS"){
    res = trajeR.POIS(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, delta, pi, Method,
                      hessian, itermax, paraminit, EMIRLS, refgr)
  }
  else{
    diffctind = 1
    if (is.null(diffct)){
      ffh <- function(beta, t, TCOV){
        return(fct(t, beta, TCOV))
      }
      diffct <- function(t, betak, TCOV){
        return(numDeriv::jacobian(func = ffh,
                                  x = betak, t = t , TCOV = TCOV)
        )
      }
      diffctind = 0
    }
    res = trajeR.NL(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, sigma, pi, Method, ssigma,
                    hessian, itermax, paraminit, EMIRLS, refgr, fct, diffct, diffctind, nls.lmiter)
  }
  return(res)
}
