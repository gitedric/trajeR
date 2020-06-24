#################################################################################
# find parameters for a CNORM model
#################################################################################
#' Internal function to fit CNORM Model
#'
#' @inheritParams trajeR
#' @return  return a object of class Trajectory.CNORM
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item sigma - vector of the parameters sigma.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a coavriate X that modify
#'   the probability or the probability of group membership.
#'   \item sd - vector of the standard deviation of the parameters.
#'   \item tab - a matrix with all the parameters and standard deviation.
#'   \item Model - a string with the model used.
#'   \item groups -  a integer with the number of group.
#'   \item Names - strings with the name of the parameters.
#'   \item Method  -  a string with the method used.
#'   \item Size - a integer with the number of individuals.
#'   \item Likelihood -  a real with the Likelihood obtained by the parameters.
#'   \item Time - a vector with the first row of time values.
#'   \item degre - a vector with the degree of the polynomial shape.
#'   \item min - a real with the minimum value for censored data.
#'   \item max - a real with the maximum value for censored data.
#' }

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
    if (max(Y)<ymax & min(Y)>ymin){
      if (ssigma == FALSE){
        param = EM_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }else{
        param = EMSigmaunique_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }
    }else{
      if (ssigma == TRUE){
        param = EMCensoredSigmaunique_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }
      else{
        param = EMCensored_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }
    }
    if (hessian == TRUE){
      SE = IEM_cpp(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, refgr)
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
#' Internal function to fit LOGIT Model
#'
#' @inheritParams trajeR
#' @return  return a object of class Trajectory.LOGIT
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a coavriate X that modify
#'   the probability or the probability of group membership.
#'   \item sd - vector of the standard deviation of the parameters.
#'   \item tab - a matrix with all the parameters and standard deviation.
#'   \item Model - a string with the model used.
#'   \item groups -  a integer with the number of group.
#'   \item Names - strings with the name of the parameters.
#'   \item Method  -  a string with the method used.
#'   \item Size - a integer with the number of individuals.
#'   \item Likelihood -  a real with the Likelihood obtained by the parameters.
#'   \item Time - a vector with the first row of time values.
#'   \item degre - a vector with the degree of the polynomial shape.
#' }


trajeR.LOGIT <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, delta, pi, Method,
                         hessian, itermax, paraminit, EMIRLS, refgr){
  if (Method == "L"){
    # initial value for Likelihood's method
    if (is.null(paraminit)){
      paraminit = c(theta, unlist(beta), unlist(delta))
    }
    newparam = optim(par = paraminit, fn = likelihoodLOGIT_cpp, gr = difLLOGIT_cpp,
                     method = "BFGS",
                     hessian = hessian,
                     control = list(fnscale=-1, trace=1, REPORT=1, maxit = itermax),
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
    param = EMLOGIT_cpp(paraminitEM, ng, nx, n, nbeta, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE = IEMLOGIT_cpp(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)], sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:ng)], param[1:ng])
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
    param = EMLOGITIRLS_cpp(paraminitEM, ng, nx, n, nbeta, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE =  IEMLOGIT_cpp(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = param = c(param[-c(1:ng)], param[1:ng])
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
#' Internal function to fit ZIP Model
#' @inheritParams trajeR
#' @return  return a object of class Trajectory.ZIP
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a coavriate X that modify
#'   the probability or the probability of group membership.
#'   \item nu - vector of the parameters nu.
#'   \item sd - vector of the standard deviation of the parameters.
#'   \item tab - a matrix with all the parameters and standard deviation.
#'   \item Model - a string with the model used.
#'   \item groups -  a integer with the number of group.
#'   \item Names - strings with the name of the parameters.
#'   \item Method  -  a string with the method used.
#'   \item Size - a integer with the number of individuals.
#'   \item Likelihood -  a real with the Likelihood obtained by the parameters.
#'   \item Time - a vector with the first row of time values.
#'   \item degre - a vector with the degree of the polynomial shape for the Poisson part.
#'   \item degre.nu - a vector with the degree of the polynomial shape for the exceeded zero state.
#' }

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
    newparam = optim(par = paraminit, fn = likelihoodZIP_cpp, gr = difLZIP_cpp,
                     method = "BFGS",
                     hessian = hessian,
                     control = list(fnscale=-1, trace=1, REPORT=1, maxit = itermax),
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
    param = EMZIP_cpp(paraminitEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE = IEMZIP_cpp(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:ng)], param[1:ng])
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
    param = EMZIPIRLS_cpp(paraminitEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE){
      SE = IEMZIP_cpp(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1){
        SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
      }else{
        SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
      }
    }else{
      SE = NA
    }
    if (nx == 1){
      param = c(param[-c(1:ng)], param[1:ng])
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
             Likelihood = Likelihood(param = c(theta, beta, nu, delta), model = "ZIP", 
                                     method =method, ng =ng, 
                                     nx = nx, n = n, nbeta = nbeta, nw = nw, A = A, Y = Y, X = X, 
                                     TCOV = TCOV, nnu = nnu),
             Time = A[1,], degre = degre - 1, degre.nu = degre.nu - 1)
  class(res) = "Trajectory.ZIP"
  return(res)
}
#################################################################################
# find parameters for Non Linear Model
#################################################################################
#' Internal function to fit Non Linear Model
#' @inheritParams trajeR
#' @return  return a object of class Trajectory.NL
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item sigma - vector of the parameters sigma.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a coavriate X that modify
#'   the probability or the probability of group membership.
#'   \item sd - vector of the standard deviation of the parameters.
#'   \item tab - a matrix with all the parameters and standard deviation.
#'   \item Model - a string with the model used.
#'   \item groups -  a integer with the number of group.
#'   \item Names - strings with the name of the parameters.
#'   \item Method  -  a string with the method used.
#'   \item Size - a integer with the number of individuals.
#'   \item Likelihood -  a real with the Likelihood obtained by the parameters.
#'   \item Time - a vector with the first row of time values.
#'   \item degre - a vector with the degree of the polynomial shape.
#'   \item fct - the defintion of the function used int this model.
#' }

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
