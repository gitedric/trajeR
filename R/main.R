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
#' @param Y a matrix containing the variables in the model.
#' @param A a matrix containing the time variable data.
#' @param X an optionnal matrix that modifie the probabilty of belong to group. By default its value is a mtrix
#' with one column  with value 1.
#' @param TCOV an optionnal matrix containing the time covariate that influence the trajectory themselve.
#' By default its value is NULL.
#' @param ng integer. The number of group.
#' @param degre vector of integer. The degre of every ploynomial function.
#' @param degre.nu vecor of integer. The degre of all Poisson part for a ZIP model.
#' @param Model the model used. The value are LOGIT for a Logit Mixture model,
#'CNORM for a Censored Normal Mixture Model or ZIP for Zero Inflated Poisson Mixture model.
#' @param Method determine the method used for find the paramters of the model.
#' The value are L for the Maximum Likelihiood Estimation, EM for Expectation Maximization method
#' with quasi newton method inside, EMIWRLS for Expectation Maximization method with Iterative
#' Weighted ??
#' @param ssigma logical. By default its value is FALSE. For the CNORM model,
#' indicate if we want the same sigma for all nomrmal density function.
#' @param ymax real. For the CNORM model, indicate the maximum value of the data. It oncerne only the model
#' with censored data. By default its value is the maximum value of the data plus 1.
#' @param ymin real. For the CNORM model, indicate the minimum value of the data. It oncerne only the model
#' with censored data. By default its value is the maximum value of the data minus 1.
#' @param hessian logical. Indicate if we want calculate the hessian matrix. Default is FALSE.
#' If the method use is Likelihood, the hessian is calculated by inversing the Information's Fisher Matrix.
#' To avoid numerically singular matrix we find the pseudo inverse matrix by using the ginv function int he package MASS.
#' If the method is EM or EMIWRLS, the hessian is calculted by using Louis method.
#' @param itermax integer. Indicate the maximal number of iteration for optim function or for the EM algorithm.
#' @param paraminit vector. The vector of initial parameters. By default trajeR calculate the initial value
#' based of the range or the standrad deviation.
#' @param fct a function. The definition of the function  f in the definition in nonlinear model.
#' @param diffct a function. The differential of the function f in the nonlinear model.
#' @param nbvar integer. The number of variable in the nonlinear model.
#'
#' @return return an object of class "\code{Trajectory.LOGIT}".
#' The generic accessor functions \code{beta}, \code{delta}, \code{theta}, \code{sd}, \code{tab}, \code{Likelihood},
#' \code{ng}, \code{model} and \code{method} extract various useful features of the value returned by \code{trajeR}.
#'
#' An object of class "\code{Trajectory.LOGIT}" is a list containing at least the following components:
#'
#'\describe{
#' \item{\code{beta}}{a vector of the parameters beta.}
#' \item{\code{delta}}{a vector of the parameter delta. Only if we use time covariate.}
#' \item{\code{theta}}{a vector with the paramter theta if there exist a coavriate X that modifie
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
#' load(file = "/data/LOGIT_data01")
#' data=as.matrix(data)
#' solL = trajeR(data[,2:11], data[,12:21], ng = 3, degre=c(0,3,4),
#'               Model="LOGIT", Method = "L", hessian = TRUE,
#'               itermax = 300, paraminit = paraminit)

trajeR <- function(Y, A, Risk = NULL, TCOV = NULL, ng, degre = NULL, degre.nu = 0,
                   Model, Method = "L", ssigma = FALSE,
                   ymax = max(Y)+1, ymin = min(Y)-1, hessian = TRUE,
                   itermax = 100, paraminit = NULL, EMIRLS = TRUE, refgr = 1,
                   fct = NULL, diffct = NULL, nbvar = NULL, nls.lmiter = 50){
  X = Risk
  if (is.null(degre)){
    degre = rep(nbvar, ng)
  }else{
    degre = degre + 1
  }
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
    beta= lapply(1:(ng), function(s){
      c(qnorm((2*(s-1)+1)/(2*ng),mean(Y),sd(Y)), rep(0, degre[s]-1))
    })
  }else if (Model == "LOGIT"){
    beta = lapply(1:(ng), function(s){
      c(suppressWarnings(log(qnorm((2*(s-1)+1)/(2*ng),mean(Y),sd(Y))/(1-qnorm((2*(s-1)+1)/(2*ng),mean(Y),sd(Y))))), rep(0, degre[s]-1))
    })
    beta[is.na(beta)]=-5
  }else if (Model =="ZIP") {
    nu = lapply(1:(ng), function(s){
      c(-3, rep(0, degre.nu[s]))
    })
    ytmp = sort(Y)
    beta= lapply(1:(ng), function(s){
      c(qpois((2*(s-1)+1)/(2*ng),mean(ytmp)), rep(0, degre[s]-1))
    })
  }else{
    beta = lapply(1:(ng), function(s){
      c(qnorm((2*(s-1)+1)/(2*ng),mean(Y),sd(Y)), rep(0, degre[s]-1))
    })
  }
  nbeta = degre
  theta = rep(1, ng*ntheta)
  X = as.matrix(X)
  pi = sapply(1:ng, function(s){piik(theta = theta, i = 1, k = s, ng = ng, X = X)})
  if (is.null(paraminit)){
    cat("Starting Values\n ")
    cat(c(pi, unlist(beta), unlist(nu), rep(sd(Y), ng), unlist(delta)))
    cat('\n\n')
    cat('Likelihood \n')
  }else{
    cat('Starting Values\n')
    cat(paraminit)
    cat('\n\n')
    cat("Likelihood\n")
  }
  if (Model == "CNORM"){
    res = trajeR.CNORM(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, sigma, delta, pi, Method, ssigma,
                       ymax, ymin, hessian, itermax, paraminit, EMIRLS, refgr)
  }else if (Model == "LOGIT"){
    res = trajeR.LOGIT(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, delta, pi, Method,
                       hessian, itermax, paraminit, EMIRLS, refgr)
  }else if (Model == "ZIP"){
    res = trajeR.ZIP(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, degre.nu, theta, beta, nu, delta, pi, Method,
                     hessian, itermax, paraminit, EMIRLS, refgr)
  }else{
    diffctind = 1
    if (is.null(diffct)){
      ffh <- function(beta, t, TCOV){
        return(fct(t, beta, TCOV))
      }
      diffct <- function(t, betak, TCOV){
        return(jacobian(func = ffh,
                        x = betak, t = t , TCOV = TCOV)
        )
      }
      diffctind = 0
    }
    assign("fct", fct, envir=globalenv())
    assign("diffct", diffct, envir=globalenv())
    res = trajeR.NL(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, sigma, pi, Method, ssigma,
                      hessian, itermax, paraminit, EMIRLS, refgr, fct, diffct, diffctind, nls.lmiter)
  }
  return(res)
}
