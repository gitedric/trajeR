#################################################################################
# find parameters for a CNORM model
#################################################################################
#' Title
#'
#' @inheritParams trajeR.LOGIT
#' @param sigma a vector of reals. The values of sigma.
#' @param ssigma logical. A boolean thaht indicates if we had to use same sigma for all groups. By default is FALSE.
#' @param ymin a real. For censored data, it is the minimum value of the data.
#' @param ymax a real. For censored data, it is the maxium value of the data.
#' @export
#'
trajeR.CNORM <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, sigma, delta, pi, Method, ssigma,
                         ymax, ymin, hessian, itermax, paraminit, EMIRLS, refgr){
  nsigma = ng
  if (Method == "L"){
    # initial value for Likelihood's method
    if (is.null(paraminit)){
      sigma = rep(sd(Y),ng)
      paraminit = c(theta, unlist(beta), sigma, unlist(delta))
    }
    paraminitL = c(paraminit[1:(ng*nx+sum(nbeta))], log(paraminit[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]),
                   paraminit[-c(0: (ng*nx+sum(nbeta)+ng))])
    if (ssigma == FALSE){
      # different sigma
      newparam = optim(par = paraminitL, fn = Likelihoodalpha_cpp, gr=difLalpha_cpp,
                       method = "BFGS",
                       hessian = hessian,
                       control = list(fnscale=-1, trace=1, REPORT=1, maxit = itermax),
                       ng = ng, nx = nx, n =n, A = A, Y = Y, X = X, nbeta = nbeta,
                       ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV)
    } else {
      # same sigma
      newparam = optim(par = paraminitL, fn = Likelihoodalpha_cpp, gr=difLalphaunique_cpp,
                       method = "BFGS",
                       hessian = hessian,
                       control = list(fnscale=-1, trace=1, REPORT=1, maxit = itermax),
                       ng = ng, nx = nx, n =n, A = A, Y = Y, X = X, nbeta = nbeta,
                       ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV)
    }
    param = newparam$par
    if (nw!=0){
      param = c(param[(ng*nx+1):(ng*nx+sum(nbeta))],exp(param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]),
                param[(ng*nx+sum(nbeta)+ng+1):(ng*nx+sum(nbeta)+ng+ng*nw)], param[1:(ng*nx)])
    }else{
      param = c(param[(ng*nx+1):(ng*nx+sum(nbeta))],exp(param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]),
                param[1:(ng*nx)])
    }
    if (hessian == TRUE){
      H = newparam$hessian
      Il = ginv(-H)
      SE = sqrt(diag(Il))
      SE = c(SE[-c(1:ng)], SE[1:ng])
    }
  } else if (Method == "EM"){
    # initial value for Likelihood's method
    if (is.null(paraminit)){
      sigma = rep(sd(Y),ng)
      if (nx == 1){
        paraminitEM = c(pi[1:(ng-1)], unlist(beta), sigma, unlist(delta))
      }else{
        paraminitEM = c(theta, unlist(beta), sigma, unlist(delta))
      }
    }else{
      if (nx == 1){paraminitEM = paraminit[-ng]}
    }
    if (max(Y)<=ymax & min(Y)>=ymin){
      if (ssigma == FALSE){
        param = EM_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }else{
        param = EMSigmaunique_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }
    }else{
      if (ssigma == TRUE){
        param = EMcensoredsamesigma(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, delta, nw, itermax, EMIRLS)
      }
      else{
        param = EMCensored_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }
    }
    if (hessian == TRUE){
      SE = IEM(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
        }
    }else{
      SE = NA
    }
    param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
  }
  if (hessian == TRUE){
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = param/SE,
                   Prob = 1-pt(abs(param/SE), n*period-1)+pt(-abs(param/SE), n*period-1)
    )
  }else{
    SE = rep(NA,length(param))
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = SE,
                   Prob =SE
    )
  }
  namedegre = c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta = c()
  for (i in 1:length(nbeta)){
    namebeta = c(namebeta, namedegre[1:nbeta[i]])
  }
  namesigma = paste0("sigma", 1:ng)
  nametheta = rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0){
    namedelta = NULL
  }else{
    namedelta = rep(paste0("TCOV", 1:nw), ng)
  }
  d.names = c(namebeta, namesigma, namedelta, nametheta)
  colnames(d) = c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta = param[1:(sum(nbeta))]
  if (nw == 0){
    delta = NA
    sigma = param[(sum(nbeta)+1):(sum(nbeta)+ng)]
    theta = param[-c(1:(sum(nbeta)+ng))]
  }else{
    delta = param[(sum(nbeta)+ng+1):(sum(nbeta)+ng+nw*ng)]
    sigma = param[(sum(nbeta)+1):(sum(nbeta)+ng)]
    theta = param[-c(1:(sum(nbeta)+ng+nw*ng))]
  }
  method = ifelse(nx != 1, "L", Method)
  res = list(beta = beta,
             sigma = sigma,
             delta = delta,
             theta = theta,
             sd = SE, tab = d, Model = "CNORM",
             groups = ng, Names = d.names, Method = Method, Size = n,
             Likelihood = Likelihood(param = c(theta, beta, sigma, delta), model = "CNORM",
                                     method = method,
                                     ng = ng, nx = nx, n = n,
                                     nbeta = nbeta, nw = nw, A= A, Y = Y, X = X,
                                     TCOV = TCOV, ymin = ymin, ymax = ymax),
             Time = A[1,], degre = degre - 1, min = ymin, max = ymax)
  class(res) = "Trajectory.CNORM"
  return(res)
}
#################################################################################
# find parameters for a LOGIT model
#################################################################################
#' trajeR.LOGIT
#'
#' @param Y a matrix containing the variables in the model.
#' @param A a matrix containing the time variable data.
#' @param X data that modifie the probabilty of belong to group. By default its value is a mtrix
#' with one column  with value 1.
#' @param TCOV a matrix containing the time covariate that influence the trajectory themselve.
#' By default its value is NULL.
#' @param ng integer. The number of group.
#' @param nx integer. The number ov variable X.
#' @param n integer. The number of individuals. It is the number of row of Y too.
#' @param nbeta vector of integer. The number of parameters beta for each group.
#' @param nw integer. The number of time dependant variable. 0 if TCOV is NULL.
#' @param ntheta vector of integer. The number of covariate X plus 1.
#' @param period integer. the number of time period.
#' @param degre vector of integer. The degre of the ploynomial shape.
#' @param theta vector of real. The values of theta parameters.
#' @param beta vector of real. The values of beta parameters.
#' @param delta vector of real. The values of delta parameters.
#' @param pi vector of real. The values of pi parameters.
#' @param Method string. The method used, Likelihhod, Expectation Maximization or
#' Expectation Maximization IRWLS.
#' @param hessian logicial. If TRUE the hessian matrix is calculated. By default its value is FALSE.
#' @param itermax integer. Indicate the maximal number of iteration for optim function or for the EM algorithm.
#' @param paraminit  vector. The vector of initial parameters. By default trajeR calculate the initial value
#' based of the range or the standrad deviation.
#' @param EMIRLS logicial. If TRUE  and if we use predictors for the membership probability we fit the partametrs
#' by using IRLS method int he EM algorithm. By default its value is FALSE.
#' @param refgr integer. Indicate the group reference for the calculus of the membership probability.
#'
#' @return  return a list object.
#' \itemize{
#'   \item beta -  vector of the paramters beta.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the paramter theta if there exist a coavriate X that modifie
#'   the probability or the probability of group membership.
#'   \item sd - vector of the standrad deviation of the parameters.
#'   \item tab - a matrix with all the parameters and standard deviation.
#'   \item Likelihoo -  a real with the Likelihhod obtnaied by the parameters.
#'   \item ng - a integer with the number of group.
#'   \item model - a string with the model used.
#'   \item method -  a string with the method used.
#' }
#' @export
#'
trajeR.LOGIT <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, delta, pi, Method,
                         hessian, itermax, paraminit, EMIRLS, refgr){
  if (Method == "L"){
    # initial value for Likelihood's method
    if (is.null(paraminit)){
      paraminit = c(theta, unlist(beta), unlist(delta))
    }
    newparam = optim(par = paraminit, fn = LikelihoodLogit, gr = difLLogit,
                     method = "BFGS",
                     hessian = hessian,
                     control = list(fnscale=-1, trace = 6, maxit = itermax),
                     ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
                     nw = nw, TCOV = TCOV)
    param = newparam$par
    param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
    if (hessian == TRUE){
      H = newparam$hessian
      Il = ginv(-H)
      SE = sqrt(diag(Il))
      SE = c(SE[-c(1:ng)], SE[1:ng])
    }
  } else if (Method == "EM"){
    # intial value for EM
    if (is.null(paraminit)){
      if (nx == 1){
        paraminitEM = c(pi[1:(ng-1)],unlist(beta), unlist(delta))
      }else{
        paraminitEM = c(theta, unlist(beta), unlist(delta))
      }
    }else{
      if (nx == 1){
        paraminitEM = paraminit[-ng]
      }else{
        paraminitEM = paraminit
      }
    }
    param = EMLogit(paraminitEM, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE = IEML(param, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)], NA)
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:(ng-1))], param[1:(ng-1)], 1-sum(param[1:(ng-1)]))
    }else{
      param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
    }
  } else if (Method == "EMIRLS"){
    # intial value for EM
    if (is.null(paraminit)){
      if (nx == 1){
        paraminitEM = c(pi[1:(ng-1)],unlist(beta), unlist(delta))
      }else{
        paraminitEM = c(theta, unlist(beta), unlist(delta))
      }
    }else{
      if (nx == 1){
        paraminitEM = paraminit[-ng]
      }else{
        paraminitEM = paraminit
      }
    }
    param = EMLogitIRLS(paraminitEM, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE = IEML(param, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:(ng-1))], param[1:(ng-1)], 1-sum(param[1:(ng-1)]))
    }else{
      param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
    }
  }
  if (hessian == TRUE){
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = param/SE,
                   Prob = 1-pt(abs(param/SE), n*period-1)+pt(-abs(param/SE), n*period-1)
    )
  }else{
    SE = rep(NA,length(param))
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = SE,
                   Prob =SE
    )
  }
  namedegre = c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta = c()
  for (i in 1:length(nbeta)){
    namebeta = c(namebeta, namedegre[1:nbeta[i]])
  }
  nametheta = rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0){
    namedelta = NULL
  }else{
    namedelta = rep(paste0("TCOV", 1:nw), ng)
  }
  d.names = c(namebeta, namedelta, nametheta)
  colnames(d) = c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta = param[1:(sum(nbeta))]
  if (nw == 0){
    delta = NA
    theta = param[-c(1:(sum(nbeta)))]
  }else{
    delta = param[(sum(nbeta)+1):(sum(nbeta)+nw*ng)]
    theta = param[-c(1:(sum(nbeta)+nw*ng))]
  }
  method = ifelse(nx != 1, "L", Method)
  res = list(beta = beta,
             delta = delta,
             theta = theta,
             sd = SE, tab = d, Model = "LOGIT",
             groups = ng, Names = d.names, Method = Method, Size = n,
             Likelihood = Likelihood(c(theta, beta, delta), model = "LOGIT", method = method,
                                     ng = ng, nx = nx, n = n, nbeta = nbeta, nw = nw,
                                     A = A, Y = Y, X = X, TCOV = TCOV),
             Time = A[1,], degre = degre - 1)
  class(res) = "Trajectory.LOGIT"
  return(res)
}
#################################################################################
# find parameters for a ZIP model
#################################################################################
#' @export
trajeR.ZIP <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, degre.nu, theta, beta, nu, delta, pi, Method,
                       hessian, itermax, paraminit, EMIRLS, refgr){
  degre.nu = degre.nu + 1
  nnu = degre.nu
#  nu = unlist(c(sapply(1:ng, function(s){
#    c(min(Y)+s*(max(Y)-min(Y))/(ng+1), rep(0, degre.nu[s]-1))
#  })))
  if (Method == "L"){
    # initial value for Likelihood's method
    if (is.null(paraminit)){
      paraminit = c(theta, unlist(beta), unlist(nu), delta)
    }
    newparam = optim(par = paraminit, fn = LikelihoodZIP, gr = difLZIP,
                     method = "BFGS",
                     hessian = hessian,
                     control = list(fnscale=-1, trace=6, maxit = itermax),
                     ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
                     nw = nw, TCOV = TCOV)
    param = newparam$par
    param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
    if (hessian == TRUE){
      H = newparam$hessian
      Il = ginv(-H)
      SE = sqrt(diag(Il))
      SE = c(SE[-c(1:ng)], SE[1:ng])
    }
  } else if (Method == "EM"){
    # intial value for EM
    if (is.null(paraminit)){
      if (nx == 1){
        paraminitEM = c(pi[1:(ng-1)], unlist(beta), unlist(nu), unlist(delta))
      }else{
        paraminitEM = c(theta, unlist(beta), unlist(nu), unlist(delta))
      }
    }else{
      if (nx == 1){
        paraminitEM = paraminit[-ng]
      }else{
        paraminitEM = paraminit
      }
    }
    param = EMZIP(paraminitEM, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE = IEMZIP(param, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, delta, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:(ng-1))], param[1:(ng-1)], 1-sum(param[1:(ng-1)]))
    }else{
      param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
    }
  }else if (Method == "EMIRLS"){
    # intial value for EM
    if (is.null(paraminit)){
      if (nx == 1){
        paraminitEM = c(pi[1:(ng-1)], unlist(beta), unlist(nu), unlist(delta))
      }else{
        paraminitEM = c(theta, unlist(beta), unlist(nu), unlist(delta))
      }
    }else{
      if (nx == 1){
        paraminitEM = paraminit[-ng]
      }else{
        paraminitEM = paraminit
      }
    }
    param = EMZIPIRLS(paraminitEM, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE = IEMZIP(param, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, delta, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:(ng-1))], param[1:(ng-1)], 1-sum(param[1:(ng-1)]))
    }else{
      param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
    }
  }
  if (hessian == TRUE){
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = param/SE,
                   Prob = 1-pt(abs(param/SE), n*period-1)+pt(-abs(param/SE), n*period-1))
  }else{
    SE = rep(NA,length(param))
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = SE,
                   Prob = SE)
  }
  namedegre = c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta = c()
  for (i in 1:length(nbeta)){
    namebeta = c(namebeta, namedegre[1:nbeta[i]])
  }
  nametheta = rep(c("Intercept", colnames(X)[-1]), ng)
  namenu =c()
  for (k in 1:ng){
    namenu = c(namenu, paste0("Nu", k, 1:nnu[k]))
  }
  if (nw == 0){
    namedelta = NULL
  }else{
    namedelta = rep(paste0("TCOV", 1:nw), ng)
  }
  d.names = c(namebeta, namenu, namedelta, nametheta)
  colnames(d) = c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta = param[1:(sum(nbeta))]
  if (nw == 0){
    delta = NA
    nu = param[(sum(nbeta)+1):(sum(nbeta)+sum(nnu))]
    theta = param[-c(1:(sum(nbeta)+sum(nnu)))]
  }else{
    nu = param[(sum(nbeta)+1):(sum(nbeta)+sum(nnu))]
    delta = param[(sum(nbeta)+sum(nnu)+1):(sum(nbeta)+sum(nnu)+nw*ng)]
    theta = param[-c(1:(sum(nbeta)++sum(nnu)+nw*ng))]
  }
  method = ifelse(nx != 1, "L", Method)
  res = list(beta = beta,
             delta = delta,
             theta = theta,
             nu = nu,
             sd = SE, tab = d, Model = "ZIP",
             groups = ng, Names = d.names, Method = Method, Size = n,
             Likelihood = Likelihood(param = c(theta, beta, nu, delta), model = "ZIP", ng =ng, 
                                     nx = nx, n = n, nbeta = nbeta, nw = nw, A = A, Y = Y, X = X, 
                                     TCOV = TCOV, nnu = nnu),
             Time = A[1,], degre = degre - 1, degre.nu = degre.nu - 1)
  class(res) = "Trajectory.ZIP"
  return(res)
}
#################################################################################
# find parameters for Non Linear Model
#################################################################################
#' @export
trajeR.NL <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, sigma, pi, Method, ssigma,
                         hessian, itermax, paraminit, EMIRLS, refgr, fct, diffct, diffctind, nls.lmiter){
  nsigma = ng
  if (Method == "L"){
    # initial value for Likelihood's method
    if (is.null(paraminit)){
      sigma = rep(sd(Y),ng)
      paraminit = c(theta, unlist(beta), sigma)
    }
    paraminitL = c(paraminit[1:(ng*nx+sum(nbeta))], log(paraminit[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]))
    if (ssigma == FALSE){
      # different sigma
      newparam = optim(par = paraminitL, fn = LikelihoodalphaNL, gr = difLalphaNL,
                       method = "BFGS",
                       hessian = hessian,
                       control = list(fnscale=-1, trace=6, maxit = itermax),
                       ng = ng, nx = nx, n =n, A = A, Y = Y, X = X, nbeta = nbeta, diffctind = diffctind,
                       TCOV = TCOV)
    } else {
      # same sigma
      newparam = optim(par = paraminitL, fn = LikelihoodalphaNL, gr = difLalphauniqueNL,
                       method = "BFGS",
                       hessian = hessian,
                       control = list(fnscale=-1, trace=6, maxit = itermax),
                       ng = ng, nx = nx, n =n, A = A, Y = Y, X = X, nbeta = nbeta, diffctind =  diffctind,
                       TCOV = TCOV)
    }
    param = newparam$par
    if (nw!=0){
      param = c(param[(ng*nx+1):(ng*nx+sum(nbeta))],exp(param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]),
                param[(ng*nx+sum(nbeta)+ng+1):(ng*nx+sum(nbeta)+ng+ng*nw)], param[1:(ng*nx)])
    }else{
      param = c(param[(ng*nx+1):(ng*nx+sum(nbeta))],exp(param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]),
                param[1:(ng*nx)])
    }
    if (hessian == TRUE){
      H = newparam$hessian
      Il = ginv(-H)
      SE = sqrt(diag(Il))
      SE = c(SE[-c(1:ng)], SE[1:ng])
    }
  } else if (Method == "EM"){
    # initial value for Likelihood's method
    if (is.null(paraminit)){
      sigma = rep(sd(Y),ng)
      if (nx == 1){
        paraminitEM = c(pi[1:(ng-1)], unlist(beta), sigma)
      }else{
        paraminitEM = c(theta, unlist(beta), sigma)
      }
    }else{
      if (nx == 1){paraminitEM = paraminit[-ng]}
    }
      if (ssigma == FALSE){
        param = EMNL(paraminitEM, ng, nx, nbeta, n, A, Y, X, TCOV, nw, itermax, EMIRLS, fct, diffct, diffctind, nls.lmiter)
      }else{
        param = EMNLSigmaunique(paraminitEM, ng, nx, nbeta, n, A, Y, X, TCOV,  nw, itermax, EMIRLS, fct, diffct, diffctind, nls.lmiter)
      }
    if (hessian == TRUE){
      SE = IEMNL(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr, fct, diffct, diffctind)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)], sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:(ng-1))], param[1:(ng-1)], 1-sum(param[1:(ng-1)]))
    }else{
      param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
    }
  }
  if (hessian == TRUE){
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = param/SE,
                   Prob = 1-pt(abs(param/SE), n*period-1)+pt(-abs(param/SE), n*period-1)
    )
  }else{
    SE = rep(NA,length(param))
    d = data.frame(Estimate =  param,
                   StandardError = SE,
                   TValue = SE,
                   Prob =SE
    )
  }
  namedegre = c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta = c()
  for (i in 1:length(nbeta)){
    namebeta = c(namebeta, namedegre[1:nbeta[i]])
  }
  namesigma = paste0("sigma", 1:ng)
  nametheta = rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0){
    namedelta = NULL
  }else{
    namedelta = rep(paste0("TCOV", 1:nw), ng)
  }
  d.names = c(namebeta, namesigma, namedelta, nametheta)
  colnames(d) = c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta = param[1:(sum(nbeta))]
  if (nw == 0){
    delta = NA
    sigma = param[(sum(nbeta)+1):(sum(nbeta)+ng)]
    theta = param[-c(1:(sum(nbeta)+ng))]
  }else{
    delta = param[(sum(nbeta)+ng+1):(sum(nbeta)+ng+nw*ng)]
    sigma = param[(sum(nbeta)+1):(sum(nbeta)+ng)]
    theta = param[-c(1:(sum(nbeta)+ng+nw*ng))]
  }
    method = ifelse(nx != 1, "L", Method)
  res = list(beta = beta,
             sigma = sigma,
             delta = delta,
             theta = theta,
             sd = SE, tab = d, Model = "CNORM",
             groups = ng, Names = d.names, Method = Method, Size = n,
             Likelihood = Likelihood(param = c(theta, beta, sigma), model = "NL",
                                     ng = ng, nx = nx, n = n,
                                     nbeta = nbeta, nw = nw, A= A, Y = Y, X = X,
                                     TCOV = TCOV, fct),
             Time = A[1,], degre = degre - 1, fct = fct)
  class(res) = "Trajectory.NL"
  return(res)
}
