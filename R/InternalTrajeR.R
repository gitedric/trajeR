#################################################################################
# find parameters for a CNORM model
#################################################################################
#' Internal function to fit CNORM Model
#'
#' @inheritParams trajeR
#' @param X Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
#' with one column  with value 1.
#' @param ng Integer. The number of groups.
#' @param nx Integer. The number of covariates.
#' @param n Integer. Number of individuals.
#' @param nbeta Vector of integers. Number of beta parameters for each group.
#' @param nw Integer. Number of time dependent covariate.
#' @param ntheta Vector of integers. Number of theta parameters for each group.
#' @param period Integer.
#' @param theta Vector of real. The parameter for calculated the group membership probability.
#' @param beta Vector of real. The beta parameter.
#' @param sigma Vector of real. The sigma parameter.
#' @param delta Vector of real. The delta parameter.
#' @param pi Vector of real. The group membership probability.
#' @param EMIRLS Boolean. True if we use EMIRLS method.
#' @return  return a object of class Trajectory.CNORM
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item sigma - vector of the parameters sigma.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a covariate X that modify
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
                         ymax, ymin, hessian, itermax, paraminit, EMIRLS, refgr) {
  set_tour(1)
  set_storelik(10**100)
  nsigma <- ng
  theta <- theta - theta[1:nx]
  if (Method == "L") {
    theta <- theta[-c(1:nx)]
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      sigma <- rep(stats::sd(Y, na.rm = TRUE), ng)
      paraminit <- c(theta, unlist(beta), sigma, unlist(delta))
    } else {
      paraminit <- paraminit[-c(1:nx)]
    }
    if (ssigma == FALSE) {
      # different sigma
      paraminitL <- c(
        paraminit[1:((ng - 1) * nx + sum(nbeta))], log(paraminit[((ng - 1) * nx + sum(nbeta) + 1):((ng - 1) * nx + sum(nbeta) + ng)]),
        paraminit[-c(1:((ng - 1) * nx + sum(nbeta) + ng))]
      )
      if (!hessian) {
        newparam <- stats::optim(
          par = paraminitL, fn = Likelihoodalpha_cpp, gr = difLalpha_cpp,
          method = "BFGS",
          hessian = hessian,
          control = list(fnscale = -1, trace = 1, REPORT = 1, maxit = itermax),
          ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
          ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV, ssigma = ssigma
        )
      } else {
        newparam <- ucminf::ucminf(
          par = paraminitL,
          fn = LCNORM,
          gr = difLCNORM,
          hessian = 2,
          control = list(stepmax = 10**(-6), maxeval = itermax),
          ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
          ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV, ssigma = ssigma
        )
      }
    } else {
      # same sigma
      paraminitL <- c(
        paraminit[1:((ng - 1) * nx + sum(nbeta))], log(paraminit[(ng - 1) * nx + sum(nbeta) + 1]),
        paraminit[-c(1:((ng - 1) * nx + sum(nbeta) + ng))]
      )
      if (!hessian) {
        newparam <- stats::optim(
          par = paraminitL, fn = Likelihoodalpha_cpp, gr = difLalphaunique_cpp,
          method = "BFGS",
          hessian = hessian,
          control = list(fnscale = -1, trace = 1, REPORT = 1, maxit = itermax),
          ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
          ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV, ssigma = ssigma
        )
      } else {
        newparam <- ucminf::ucminf(
          par = paraminitL,
          fn = LCNORM,
          gr = difLCNORMss,
          hessian = 2,
          control = list(stepmax = 10**(-5), maxeval = itermax),
          ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
          ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV, ssigma = ssigma
        )
      }
    }
    param <- newparam$par
    indtheta <- 1:((ng - 1) * nx)
    indbeta <- ((ng - 1) * nx + 1):((ng - 1) * nx + sum(nbeta))
    if (ssigma) {
      indsigma <- ((ng - 1) * nx + sum(nbeta) + 1)
      if (nw != 0) {
        inddelta <- ((ng - 1) * nx + sum(nbeta) + 2):((ng - 1) * nx + sum(nbeta) + 1 + nw * ng)
      }
    } else {
      indsigma <- ((ng - 1) * nx + sum(nbeta) + 1):((ng - 1) * nx + sum(nbeta) + ng)
      if (nw != 0) {
        inddelta <- ((ng - 1) * nx + sum(nbeta) + ng + 1):((ng - 1) * nx + sum(nbeta) + ng + nw * ng)
      }
    }
    if (ssigma) {
      if (nw != 0) {
        param <- c(param[c(indtheta, indbeta)], rep(exp(param[indsigma]), ng), param[inddelta])
      } else {
        param <- c(param[c(indtheta, indbeta)], rep(exp(param[indsigma]), ng))
      }
    } else {
      sigma <- exp(param[indsigma])
      param[indsigma] <- sigma
    }


    theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
    param <- c(param[-c(1:((ng - 1) * nx))], theta)
    if (hessian == TRUE) {
      invH <- newparam$invhessian
      SE <- sqrt(diag(invH))
      if (ssigma) {
        sdsigma <- rep(exp(newparam$par[indsigma[1]]) * SE[indsigma[1]], ng)
      } else {
        matsigma <- diag(sigma)
        sdsigma <- sqrt(diag(matsigma %*% invH[indsigma, indsigma] %*% matsigma))
      }
      if (nx == 1) {
        sdtmp <- deltaTheta(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        sdbase <- deltaThetaBase(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        if (nw != 0) {
          SE <- c(SE[indbeta], sdsigma, SE[inddelta], sdbase, sdtmp)
        } else {
          SE <- c(SE[indbeta], sdsigma, sqrt(sum(sdtmp**2)), sdtmp)
        }
      } else {
        if (nw != 0) {
          SE <- c(SE[indbeta], sdsigma, SE[inddelta], rep(NA, nx), SE[indtheta])
        } else {
          SE <- c(SE[indbeta], sdsigma, rep(NA, nx), SE[indtheta])
        }
      }
    }
  } else if (Method == "EM") {
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      sigma <- rep(stats::sd(Y), ng)
      if (nx == 1) {
        paraminitEM <- c(pi[1:(ng - 1)], unlist(beta), sigma, unlist(delta))
      } else {
        paraminitEM <- c(theta, unlist(beta), sigma, unlist(delta))
      }
    } else {
      if (nx == 1) {
        paraminitEM <- paraminit[-ng]
      }
    }
    if (max(Y) < ymax & min(Y) > ymin) {
      if (ssigma == FALSE) {
        param <- EM_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      } else {
        param <- EMSigmaunique_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }
    } else {
      if (ssigma == TRUE) {
        param <- EMCensoredSigmaunique_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      } else {
        param <- EMCensored_cpp(paraminitEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, itermax, EMIRLS, refgr)
      }
    }
    if (hessian == TRUE) {
      SE <- IEM_cpp(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, refgr)
      if (nx == 1) {
        if (nw == 0) {
          SE <- c(SE[-c(1:(ng - 1))], SE[1:(ng - 1)], sqrt(sum(SE[1:(ng - 1)]**2)))
        } else {
          SE <- c(
            SE[ng:(ng + sum(nbeta) - 1)], SE[(ng + sum(nbeta) + sum(nw * ng)):(ng + sum(nbeta) + sum(nw * ng) + ng - 1)],
            SE[(ng + sum(nbeta)):(ng + sum(nbeta) + sum(nw * ng) - 1)],
            sqrt(sum(SE[1:(ng - 1)]**2)), SE[1:(ng - 1)]
          )
        }
      } else {
        if (nw == 0) {
          SE <- c(SE[-c(1:((ng - 1) * nx))], rep(0, nx), SE[1:((ng - 1) * nx)])
        } else {
          SE <- c(
            SE[(ng * nx):(ng * nx + sum(nbeta) - 1)], SE[(ng * nx + sum(nbeta) + sum(nw * ng)):(ng * nx + sum(nbeta) + sum(nw * ng) + ng - 1)],
            SE[(ng * nx + sum(nbeta)):(ng * nx + sum(nbeta) + sum(nw * ng) - 1)]
          )
        }
      }
    } else {
      SE <- NA
    }
    param <- c(param[-c(1:(ng * nx))], param[1:(ng * nx)])
  }
  if (hessian == TRUE) {
    if (nx == 1 & Method == "L") {
      paramtmp <- c(param[1:(length(param) - ng * nx)], exp(theta) / sum(exp(theta)))
      prob <- 1 - stats::pt(abs(paramtmp / SE), n * period - 1) + stats::pt(-abs(paramtmp / SE), n * period - 1)
    } else {
      prob <- 1 - stats::pt(abs(param / SE), n * period - 1) + stats::pt(-abs(param / SE), n * period - 1)
    }
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = param / SE,
      Prob = prob
    )
  } else {
    SE <- rep(NA, length(param))
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = SE,
      Prob = SE
    )
  }
  namedegre <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta <- c()
  for (i in 1:length(nbeta)) {
    namebeta <- c(namebeta, namedegre[1:nbeta[i]])
  }
  namesigma <- paste0("sigma", 1:ng)
  nametheta <- rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0) {
    namedelta <- NULL
  } else {
    namedelta <- rep(paste0("TCOV", 1:nw), ng)
  }
  d.names <- c(namebeta, namesigma, namedelta, nametheta)
  colnames(d) <- c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta <- param[1:(sum(nbeta))]
  if (nw == 0) {
    delta <- NA
    sigma <- param[(sum(nbeta) + 1):(sum(nbeta) + ng)]
    theta <- param[-c(1:(sum(nbeta) + ng))]
  } else {
    delta <- param[(sum(nbeta) + ng + 1):(sum(nbeta) + ng + nw * ng)]
    sigma <- param[(sum(nbeta) + 1):(sum(nbeta) + ng)]
    theta <- param[-c(1:(sum(nbeta) + ng + nw * ng))]
  }
  res <- list(
    beta = beta,
    sigma = sigma,
    delta = delta,
    theta = theta,
    sd = SE, tab = d, Model = "CNORM",
    groups = ng, Names = d.names, Method = Method, Size = n,
    Likelihood = Likelihood(
      param = c(theta, beta, sigma, delta), model = "CNORM",
      method = Method,
      ng = ng, nx = nx, n = n,
      nbeta = nbeta, nw = nw, A = A, Y = Y, X = X,
      TCOV = TCOV, ymin = ymin, ymax = ymax
    ),
    Time = A[1, ], degre = degre - 1, min = ymin, max = ymax
  )
  class(res) <- "Trajectory.CNORM"
  return(res)
}
#################################################################################
# find parameters for a LOGIT model
#################################################################################
#' Internal function to fit LOGIT Model
#'
#' @inheritParams trajeR
#' @inheritParams trajeR.CNORM
#' @param X Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
#' with one column  with value 1.
#' @return  return a object of class Trajectory.LOGIT
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a covariate X that modify
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
                         hessian, itermax, paraminit, EMIRLS, refgr) {
  set_tour(1)
  set_storelik(10**100)
  theta <- theta - theta[1:nx]
  if (Method == "L") {
    theta <- theta[-c(1:nx)]
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      paraminit <- c(theta, unlist(beta), unlist(delta))
    } else {
      paraminit <- paraminit[-c(1:nx)]
    }
    if (!hessian) {
      newparam <- stats::optim(
        par = paraminit, fn = likelihoodLOGIT_cpp, gr = difLLOGIT_cpp,
        method = "BFGS",
        hessian = hessian,
        control = list(fnscale = -1, trace = 1, REPORT = 1, maxit = itermax),
        ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
        nw = nw, TCOV = TCOV
      )
      param <- newparam$par
      theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
      param <- c(param[-c(1:((ng - 1) * nx))], theta)
    } else {
      newparam <- ucminf::ucminf(
        par = paraminit,
        fn = LLOGIT,
        gr = difLLOGIT,
        hessian = 2,
        control = list(stepmax = 10**(-5), maxeval = itermax),
        ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
        nw = nw, TCOV = TCOV
      )
      param <- newparam$par
      theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
      param <- c(param[-c(1:((ng - 1) * nx))], theta)
    }
    if (hessian == TRUE) {
      invH <- newparam$invhessian
      SE <- sqrt(diag(invH))
      if (nx == 1) {
        sdtmp <- deltaTheta(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        sdbase <- deltaThetaBase(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        SE <- c(SE[-c(1:((ng - 1) * nx))], sdbase, sdtmp)
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(NA, nx), SE[c(1:((ng - 1) * nx))])
      }
    }
  } else if (Method == "EM") {
    # intial value for EM
    if (is.null(paraminit)) {
      if (nx == 1) {
        paraminitEM <- c(pi[1:(ng - 1)], unlist(beta), unlist(delta))
      } else {
        paraminitEM <- c(theta, unlist(beta), unlist(delta))
      }
    } else {
      if (nx == 1) {
        paraminitEM <- paraminit[-ng]
      } else {
        paraminitEM <- paraminit
      }
    }
    param <- EMLOGIT_cpp(paraminitEM, ng, nx, n, nbeta, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE) {
      SE <- IEMLOGIT_cpp(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1) {
        SE <- c(SE[-c(1:(ng - 1))], SE[1:(ng - 1)], sqrt(sum(SE[1:(ng - 1)]**2)))
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(0, nx), SE[1:((ng - 1) * nx)])
      }
    } else {
      SE <- NA
    }
    if (nx == 1) {
      param <- c(param[-c(1:ng)], param[1:ng])
    } else {
      param <- c(param[-c(1:(ng * nx))], param[1:(ng * nx)])
    }
  } else if (Method == "EMIRLS") {
    # intial value for EM
    if (is.null(paraminit)) {
      if (nx == 1) {
        paraminitEM <- c(pi[1:(ng - 1)], unlist(beta), unlist(delta))
      } else {
        paraminitEM <- c(theta, unlist(beta), unlist(delta))
      }
    } else {
      if (nx == 1) {
        paraminitEM <- paraminit[-ng]
      } else {
        paraminitEM <- paraminit
      }
    }
    param <- EMLOGITIRLS_cpp(paraminitEM, ng, nx, n, nbeta, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE) {
      SE <- IEMLOGIT_cpp(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1) {
        SE <- c(SE[-c(1:(ng - 1))], SE[1:(ng - 1)], sqrt(sum(SE[1:(ng - 1)]**2)))
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(0, nx), SE[1:((ng - 1) * nx)])
      }
    } else {
      SE <- NA
    }
    if (nx == 1) {
      param <- param <- c(param[-c(1:ng)], param[1:ng])
    } else {
      param <- c(param[-c(1:(ng * nx))], param[1:(ng * nx)])
    }
  }
  if (hessian == TRUE) {
    if (nx == 1 & Method == "L") {
      paramtmp <- c(param[1:(length(param) - ng * nx)], exp(theta) / sum(exp(theta)))
      prob <- 1 - stats::pt(abs(paramtmp / SE), n * period - 1) + stats::pt(-abs(paramtmp / SE), n * period - 1)
    } else {
      prob <- 1 - stats::pt(abs(param / SE), n * period - 1) + stats::pt(-abs(param / SE), n * period - 1)
    }
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = param / SE,
      Prob = prob
    )
  } else {
    SE <- rep(NA, length(param))
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = SE,
      Prob = SE
    )
  }
  namedegre <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta <- c()
  for (i in 1:length(nbeta)) {
    namebeta <- c(namebeta, namedegre[1:nbeta[i]])
  }
  nametheta <- rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0) {
    namedelta <- NULL
  } else {
    namedelta <- rep(paste0("TCOV", 1:nw), ng)
  }
  d.names <- c(namebeta, namedelta, nametheta)
  colnames(d) <- c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta <- param[1:(sum(nbeta))]
  if (nw == 0) {
    delta <- NA
    theta <- param[-c(1:(sum(nbeta)))]
  } else {
    delta <- param[(sum(nbeta) + 1):(sum(nbeta) + nw * ng)]
    theta <- param[-c(1:(sum(nbeta) + nw * ng))]
  }
  res <- list(
    beta = beta,
    delta = delta,
    theta = theta,
    sd = SE, tab = d, Model = "LOGIT",
    groups = ng, Names = d.names, Method = Method, Size = n,
    Likelihood = Likelihood(c(theta, beta, delta),
      model = "LOGIT", method = Method,
      ng = ng, nx = nx, n = n, nbeta = nbeta, nw = nw,
      A = A, Y = Y, X = X, TCOV = TCOV
    ),
    Time = A[1, ], degre = degre - 1
  )
  class(res) <- "Trajectory.LOGIT"
  return(res)
}
#################################################################################
# find parameters for a ZIP model
#################################################################################
#' Internal function to fit ZIP Model
#'
#' @inheritParams trajeR
#' @inheritParams trajeR.CNORM
#' @param nu Vector of real. The nu parameter.
#' @param X Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
#' with one column  with value 1.
#' @return  return a object of class Trajectory.ZIP
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a covariate X that modify
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
                       hessian, itermax, paraminit, EMIRLS, refgr) {
  set_tour(1)
  set_storelik(10**100)
  theta <- theta - theta[1:nx]
  degre.nu <- degre.nu + 1
  nnu <- degre.nu
  #  nu = unlist(c(sapply(1:ng, function(s){
  #    c(min(Y)+s*(max(Y)-min(Y))/(ng+1), rep(0, degre.nu[s]-1))
  #  })))
  if (Method == "L") {
    theta <- theta[-c(1:nx)]
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      paraminit <- c(theta, unlist(beta), unlist(nu), delta)
    }
    else {
      paraminit <- paraminit[-c(1:nx)]
    }
    if (!hessian) {
      newparam <- stats::optim(
        par = paraminit, fn = likelihoodZIP_cpp, gr = difLZIP_cpp,
        method = "BFGS",
        hessian = hessian,
        control = list(fnscale = -1, trace = 1, REPORT = 1, maxit = itermax),
        ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
        nw = nw, TCOV = TCOV
      )
      param <- newparam$par
      theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
      param <- c(param[-c(1:((ng - 1) * nx))], theta)
    } else {
      newparam <- ucminf::ucminf(
        par = paraminit,
        fn = LZIP,
        gr = difLZIP,
        hessian = 2,
        control = list(stepmax = 10**(-5), maxeval = itermax),
        ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
        nw = nw, TCOV = TCOV
      )
      param <- newparam$par
      theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
      param <- c(param[-c(1:((ng - 1) * nx))], theta)
    }
    if (hessian == TRUE) {
      invH <- newparam$invhessian
      SE <- sqrt(diag(invH))
      if (nx == 1) {
        sdtmp <- deltaTheta(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        sdbase <- deltaThetaBase(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        SE <- c(SE[-c(1:((ng - 1) * nx))], sdbase, sdtmp)
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(NA, nx), SE[c(1:((ng - 1) * nx))])
      }
    }
  } else if (Method == "EM") {
    # intial value for EM
    if (is.null(paraminit)) {
      if (nx == 1) {
        paraminitEM <- c(pi[1:(ng - 1)], unlist(beta), unlist(nu), unlist(delta))
      } else {
        paraminitEM <- c(rep(0, nx), theta, unlist(beta), unlist(nu), unlist(delta))
      }
    } else {
      if (nx == 1) {
        paraminitEM <- paraminit[-ng]
      } else {
        paraminitEM <- paraminit
      }
    }
    param <- EMZIP_cpp(paraminitEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE) {
      SE <- IEMZIP_cpp(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1) {
        SE <- c(SE[-c(1:(ng - 1))], SE[1:(ng - 1)], sqrt(sum(SE[1:(ng - 1)]**2)))
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(0, nx), SE[1:((ng - 1) * nx)])
      }
    } else {
      SE <- NA
    }
    if (nx == 1) {
      param <- c(param[-c(1:ng)], param[1:ng])
    } else {
      param <- c(param[-c(1:(ng * nx))], param[1:(ng * nx)])
    }
  } else if (Method == "EMIRLS") {
    # intial value for EM
    if (is.null(paraminit)) {
      if (nx == 1) {
        paraminitEM <- c(pi[1:(ng - 1)], unlist(beta), unlist(nu), unlist(delta))
      } else {
        paraminitEM <- c(rep(0, nx), theta, unlist(beta), unlist(nu), unlist(delta))
      }
    } else {
      if (nx == 1) {
        paraminitEM <- paraminit[-ng]
      } else {
        paraminitEM <- paraminit
      }
    }
    param <- EMZIPIRLS_cpp(paraminitEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
    if (hessian == TRUE) {
      SE <- IEMZIP_cpp(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw, refgr)
      if (nx == 1) {
        SE <- c(SE[-c(1:(ng - 1))], SE[1:(ng - 1)], sqrt(sum(SE[1:(ng - 1)]**2)))
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(0, nx), SE[1:((ng - 1) * nx)])
      }
    } else {
      SE <- NA
    }
    if (nx == 1) {
      param <- c(param[-c(1:ng)], param[1:ng])
    } else {
      param <- c(param[-c(1:(ng * nx))], param[1:(ng * nx)])
    }
  }
  if (hessian == TRUE) {
    if (nx == 1 & Method == "L") {
      paramtmp <- c(param[1:(length(param) - ng * nx)], exp(theta) / sum(exp(theta)))
      prob <- 1 - stats::pt(abs(paramtmp / SE), n * period - 1) + stats::pt(-abs(paramtmp / SE), n * period - 1)
    } else {
      prob <- 1 - stats::pt(abs(param / SE), n * period - 1) + stats::pt(-abs(param / SE), n * period - 1)
    }
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = param / SE,
      Prob = prob
    )
  } else {
    SE <- rep(NA, length(param))
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = SE,
      Prob = SE
    )
  }
  namedegre <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta <- c()
  for (i in 1:length(nbeta)) {
    namebeta <- c(namebeta, namedegre[1:nbeta[i]])
  }
  nametheta <- rep(c("Intercept", colnames(X)[-1]), ng)
  namenu <- c()
  for (k in 1:ng) {
    namenu <- c(namenu, paste0("Nu", k, 1:nnu[k]))
  }
  if (nw == 0) {
    namedelta <- NULL
  } else {
    namedelta <- rep(paste0("TCOV", 1:nw), ng)
  }
  d.names <- c(namebeta, namenu, namedelta, nametheta)
  colnames(d) <- c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta <- param[1:(sum(nbeta))]
  if (nw == 0) {
    delta <- NULL
    nu <- param[(sum(nbeta) + 1):(sum(nbeta) + sum(nnu))]
    theta <- param[-c(1:(sum(nbeta) + sum(nnu)))]
  } else {
    nu <- param[(sum(nbeta) + 1):(sum(nbeta) + sum(nnu))]
    delta <- param[(sum(nbeta) + sum(nnu) + 1):(sum(nbeta) + sum(nnu) + nw * ng)]
    theta <- param[-c(1:(sum(nbeta) + +sum(nnu) + nw * ng))]
  }
  method <- ifelse(nx != 1, "L", Method)
  res <- list(
    beta = beta,
    delta = delta,
    theta = theta,
    nu = nu,
    sd = SE, tab = d, Model = "ZIP",
    groups = ng, Names = d.names, Method = Method, Size = n,
    Likelihood = Likelihood(
      param = c(theta, beta, nu, delta), model = "ZIP",
      method = method, ng = ng,
      nx = nx, n = n, nbeta = nbeta, nw = nw, A = A, Y = Y, X = X,
      TCOV = TCOV, nnu = nnu
    ),
    #Time = A[1, ], degre = degre - 1, degre.nu = degre.nu - 1
    Time = c(min(A, na.rm = TRUE), max(A, na.rm = TRUE)), period = period,
    degre = degre - 1, degre.nu = degre.nu - 1
  )
  class(res) <- "Trajectory.ZIP"
  return(res)
}
#################################################################################
# find parameters for a Poisson model
#################################################################################
#' Internal function to fit poisson Model
#'
#' @inheritParams trajeR
#' @inheritParams trajeR.CNORM
#' @param X Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
#' with one column  with value 1.
#' @return  return a object of class Trajectory.Pois
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a covariate X that modify
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
#'   \item degre - a vector with the degree of the polynomial shape for the Poisson part.
#' }

trajeR.POIS <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, delta, pi, Method,
                        hessian, itermax, paraminit, EMIRLS, refgr) {
  set_tour(1)
  set_storelik(10**100)
  theta <- theta - theta[1:nx]
  theta <- theta[-c(1:nx)]
  #  nu = unlist(c(sapply(1:ng, function(s){
  #    c(min(Y)+s*(max(Y)-min(Y))/(ng+1), rep(0, degre.nu[s]-1))
  #  })))
  if (Method == "L") {
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      paraminit <- c(theta, unlist(beta), delta)
    }
    if (!hessian) {
      newparam <- stats::optim(
        par = paraminit, fn = likelihoodPois_cpp, gr = difLPois_cpp,
        method = "BFGS",
        hessian = hessian,
        control = list(fnscale = -1, trace = 1, REPORT = 1, maxit = itermax),
        ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
        nw = nw, TCOV = TCOV
      )
      param <- newparam$par
      theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
      param <- c(param[-c(1:((ng - 1) * nx))], theta)
    }
    # }else{
    #   newparam = ucminf::ucminf(par = paraminit,
    #                             fn = LZIP,
    #                             gr = difLZIP,
    #                             hessian = 2,
    #                             control = list(stepmax=10**(-5), maxeval = itermax),
    #                             ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
    #                             nw = nw, TCOV = TCOV)
    #   param = newparam$par
    #   theta = c(rep(0, nx), param[c(1:((ng-1)*nx))])
    #   param = c(param[-c(1:((ng-1)*nx))], theta)
    # }
    #   if (hessian == TRUE){
    #     invH = newparam$invhessian
    #     SE = sqrt(diag(invH))
    #     if (nx == 1){
    #       sdtmp = deltaTheta(theta, invH[1: (ng-1), 1:(ng - 1)], X, ng - 1)
    #       sdbase = deltaThetaBase(theta, invH[1: (ng-1), 1:(ng - 1)], X, ng - 1)
    #       SE = c(SE[-c(1:((ng-1)*nx))], sdbase, sdtmp)
    #     }else{
    #       SE = c(SE[-c(1:((ng-1)*nx))], rep(NA, nx), SE[c(1:((ng-1)*nx))])
    #     }
    #   }
  }
  # } else if (Method == "EM"){
  #   # intial value for EM
  #   if (is.null(paraminit)){
  #     if (nx == 1){
  #       paraminitEM = c(pi[1:(ng-1)], unlist(beta), unlist(nu), unlist(delta))
  #     }else{
  #       paraminitEM = c(rep(0,nx), theta, unlist(beta), unlist(nu), unlist(delta))
  #     }
  #   }else{
  #     if (nx == 1){
  #       paraminitEM = paraminit[-ng]
  #     }else{
  #       paraminitEM = paraminit
  #     }
  #   }
  #   param = EMZIP_cpp(paraminitEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
  #   if (hessian == TRUE){
  #     SE = IEMZIP_cpp(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw, refgr)
  #     if (nx == 1){
  #       SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
  #     }else{
  #       SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
  #     }
  #   }else{
  #     SE = NA
  #   }
  #   if (nx == 1){
  #     param = c(param[-c(1:ng)], param[1:ng])
  #   }else{
  #     param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
  #   }
  # }else if (Method == "EMIRLS"){
  #   # intial value for EM
  #   if (is.null(paraminit)){
  #     if (nx == 1){
  #       paraminitEM = c(pi[1:(ng-1)], unlist(beta), unlist(nu), unlist(delta))
  #     }else{
  #       paraminitEM = c(rep(0,nx), theta, unlist(beta), unlist(nu), unlist(delta))
  #     }
  #   }else{
  #     if (nx == 1){
  #       paraminitEM = paraminit[-ng]
  #     }else{
  #       paraminitEM = paraminit
  #     }
  #   }
  #   param = EMZIPIRLS_cpp(paraminitEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
  #   if (hessian == TRUE){
  #     SE = IEMZIP_cpp(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw, refgr)
  #     if (nx == 1){
  #       SE = c(SE[-c(1:(ng-1))], SE[1:(ng-1)],  sqrt(sum(SE[1:(ng-1)]**2)))
  #     }else{
  #       SE = c(SE[-c(1:((ng-1)*nx))], rep(0, nx), SE[1:((ng-1)*nx)])
  #     }
  #   }else{
  #     SE = NA
  #   }
  #   if (nx == 1){
  #     param = c(param[-c(1:ng)], param[1:ng])
  #   }else{
  #     param = c(param[-c(1:(ng*nx))], param[1:(ng*nx)])
  #   }
  # }
  if (hessian == TRUE) {
    if (nx == 1 & Method == "L") {
      paramtmp <- c(param[1:(length(param) - ng * nx)], exp(theta) / sum(exp(theta)))
      prob <- 1 - stats::pt(abs(paramtmp / SE), n * period - 1) + stats::pt(-abs(paramtmp / SE), n * period - 1)
    } else {
      prob <- 1 - stats::pt(abs(param / SE), n * period - 1) + stats::pt(-abs(param / SE), n * period - 1)
    }
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = param / SE,
      Prob = prob
    )
  } else {
    SE <- rep(NA, length(param))
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = SE,
      Prob = SE
    )
  }
  namedegre <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta <- c()
  for (i in 1:length(nbeta)) {
    namebeta <- c(namebeta, namedegre[1:nbeta[i]])
  }
  nametheta <- rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0) {
    namedelta <- NULL
  } else {
    namedelta <- rep(paste0("TCOV", 1:nw), ng)
  }
  d.names <- c(namebeta, namedelta, nametheta)
  colnames(d) <- c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta <- param[1:(sum(nbeta))]
  if (nw == 0) {
    delta <- NULL
    theta <- param[-c(1:(sum(nbeta)))]
  } else {
    delta <- param[(sum(nbeta) + 1):(sum(nbeta) + nw * ng)]
    theta <- param[-c(1:(sum(nbeta) + nw * ng))]
  }
  method <- ifelse(nx != 1, "L", Method)
  res <- list(
    beta = beta,
    delta = delta,
    theta = theta,
    sd = SE, tab = d, Model = "POIS",
    groups = ng, Names = d.names, Method = Method, Size = n,
    Likelihood = Likelihood(
      param = c(theta, beta, delta), model = "POIS",
      method = method, ng = ng,
      nx = nx, n = n, nbeta = nbeta, nw = nw, A = A, Y = Y, X = X,
      TCOV = TCOV
    ),
    Time = A[1, ], degre = degre - 1
  )
  class(res) <- "Trajectory.POIS"
  return(res)
}
#################################################################################
# find parameters for Non Linear Model
#################################################################################
#' Internal function to fit Non Linear Model
#'
#' @inheritParams trajeR
#' @inheritParams trajeR.CNORM
#' @param X Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
#' with one column  with value 1.
#' @return  return a object of class Trajectory.NL
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item sigma - vector of the parameters sigma.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a covariate X that modify
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
#'   \item fct - the definition of the function used int this model.
#' }

trajeR.NL <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nw, ntheta, period, degre, theta, beta, sigma, pi, Method, ssigma,
                      hessian, itermax, paraminit, EMIRLS, refgr, fct, diffct, nls.lmiter) {
  nsigma <- ng
  if (Method == "L") {
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      sigma <- rep(stats::sd(Y), ng)
      paraminit <- c(theta, unlist(beta), sigma)
    }
    paraminitL <- c(paraminit[1:(ng * nx + sum(nbeta))], log(paraminit[(ng * nx + sum(nbeta) + 1):(ng * nx + sum(nbeta) + ng)]))
    if (ssigma == FALSE) {
      # different sigma
      newparam <- stats::optim(
        par = paraminitL, fn = LikelihoodalphaNL, gr = difLalphaNL,
        method = "BFGS",
        hessian = hessian,
        control = list(fnscale = -1, trace = 6, maxit = itermax),
        ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
        TCOV = TCOV, fct = fct, diffct = diffct
      )
    } else {
      # same sigma
      newparam <- stats::optim(
        par = paraminitL, fn = LikelihoodalphaNL, gr = difLalphauniqueNL,
        method = "BFGS",
        hessian = hessian,
        control = list(fnscale = -1, trace = 6, maxit = itermax),
        ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
        TCOV = TCOV
      )
    }
    param <- newparam$par
    if (nw != 0) {
      param <- c(
        param[(ng * nx + 1):(ng * nx + sum(nbeta))], exp(param[(ng * nx + sum(nbeta) + 1):(ng * nx + sum(nbeta) + ng)]),
        param[(ng * nx + sum(nbeta) + ng + 1):(ng * nx + sum(nbeta) + ng + ng * nw)], param[1:(ng * nx)]
      )
    } else {
      param <- c(
        param[(ng * nx + 1):(ng * nx + sum(nbeta))], exp(param[(ng * nx + sum(nbeta) + 1):(ng * nx + sum(nbeta) + ng)]),
        param[1:(ng * nx)]
      )
    }
    if (hessian == TRUE) {
      H <- newparam$hessian
      Il <- MASS::ginv(-H)
      SE <- sqrt(diag(Il))
      SE <- c(SE[-c(1:ng)], SE[1:ng])
    }
  } else if (Method == "EM") {
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      sigma <- rep(stats::sd(Y), ng)
      if (nx == 1) {
        paraminitEM <- c(pi[1:(ng - 1)], unlist(beta), sigma)
      } else {
        paraminitEM <- c(theta, unlist(beta), sigma)
      }
    } else {
      if (nx == 1) {
        paraminitEM <- paraminit[-ng]
      }
    }
    if (ssigma == FALSE) {
      param <- EMNL(paraminitEM, ng, nx, nbeta, n, A, Y, X, TCOV, nw, itermax, EMIRLS, fct, diffct, nls.lmiter)
    } else {
      param <- EMNLSigmaunique(paraminitEM, ng, nx, nbeta, n, A, Y, X, TCOV, nw, itermax, EMIRLS, fct, diffct, nls.lmiter)
    }
    if (hessian == TRUE) {
      SE <- IEMNL(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr, fct, diffct)
      if (nx == 1) {
        SE <- c(SE[-c(1:(ng - 1))], SE[1:(ng - 1)], sqrt(sum(SE[1:(ng - 1)]**2)))
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(0, nx), SE[1:((ng - 1) * nx)])
      }
    } else {
      SE <- NA
    }
    if (nx == 1) {
      param <- c(param[-c(1:(ng - 1))], param[1:(ng - 1)], 1 - sum(param[1:(ng - 1)]))
    } else {
      param <- c(param[-c(1:(ng * nx))], param[1:(ng * nx)])
    }
  }
  if (hessian == TRUE) {
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = param / SE,
      Prob = 1 - stats::pt(abs(param / SE), n * period - 1) + stats::pt(-abs(param / SE), n * period - 1)
    )
  } else {
    SE <- rep(NA, length(param))
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = SE,
      Prob = SE
    )
  }
  namedegre <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta <- c()
  for (i in 1:length(nbeta)) {
    namebeta <- c(namebeta, namedegre[1:nbeta[i]])
  }
  namesigma <- paste0("sigma", 1:ng)
  nametheta <- rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0) {
    namedelta <- NULL
  } else {
    namedelta <- rep(paste0("TCOV", 1:nw), ng)
  }
  d.names <- c(namebeta, namesigma, namedelta, nametheta)
  colnames(d) <- c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta <- param[1:(sum(nbeta))]
  if (nw == 0) {
    delta <- NA
    sigma <- param[(sum(nbeta) + 1):(sum(nbeta) + ng)]
    theta <- param[-c(1:(sum(nbeta) + ng))]
  } else {
    delta <- param[(sum(nbeta) + ng + 1):(sum(nbeta) + ng + nw * ng)]
    sigma <- param[(sum(nbeta) + 1):(sum(nbeta) + ng)]
    theta <- param[-c(1:(sum(nbeta) + ng + nw * ng))]
  }
  method <- ifelse(nx != 1, "L", Method)
  res <- list(
    beta = beta,
    sigma = sigma,
    delta = delta,
    theta = theta,
    sd = SE, tab = d, Model = "CNORM",
    groups = ng, Names = d.names, Method = Method, Size = n,
    Likelihood = Likelihood(
      param = c(theta, beta, sigma), model = "NL",
      ng = ng, nx = nx, n = n,
      nbeta = nbeta, nw = nw, A = A, Y = Y, X = X,
      TCOV = TCOV, fct = fct
    ),
    Time = A[1, ], degre = degre - 1, fct = fct
  )
  class(res) <- "Trajectory.NL"
  return(res)
}

#################################################################################
# find parameters for Beta Model
#################################################################################
#' Internal function to fit Beta regression
#'
#' @inheritParams trajeR
#' @inheritParams trajeR.CNORM
#' @param phi Vector of real. The phi parameter.
#' @param nphi Vector of integers. Number of phi parameters for each group.
#' @param X Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
#' with one column  with value 1.
#' @return  return a object of class Trajectory.NL
#' \itemize{
#'   \item beta -  vector of the parameter beta.
#'   \item sigma - vector of the parameters sigma.
#'   \item delta - vector of the parameter delta. Only if we use time covariate.
#'   \item theta - vector with the parameter theta if there exist a covariate X that modify
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
trajeR.BETA <- function(Y, A, X, TCOV, ng, nx, n, nbeta, nphi, nw, ntheta, period, degre, theta, beta, phi, delta, pi, Method,
                        hessian, itermax, paraminit, EMIRLS, refgr) {
  set_tour(1)
  set_storelik(10**100)
  theta <- theta - theta[1:nx]
  if (Method == "L") {
    theta <- theta[-c(1:nx)]
    # initial value for Likelihood's method
    if (is.null(paraminit)) {
      paraminit <- c(theta, unlist(beta), unlist(phi), unlist(delta))
    } else {
      paraminit <- paraminit[-c(1:nx)]
    }
    if (!hessian) {
      newparam <- stats::optim(
        par = paraminit, fn = LikelihoodBETA_cpp, gr = difLBETA_cpp,
        method = "BFGS",
        hessian = hessian,
        control = list(fnscale = -1, trace = 1, REPORT = 1, maxit = itermax),
        ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta, nphi = nphi,
        nw = nw, TCOV = TCOV
      )
      param <- newparam$par
      theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
      param <- c(param[-c(1:((ng - 1) * nx))], theta)
      #param[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))] <- exp(param[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))])
    } else {
      newparam <- ucminf::ucminf(
        par = paraminit,
        fn = LBETA,
        gr = difLBETA,
        hessian = 2,
        control = list(stepmax = 10**(-5), maxeval = itermax),
        ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta, nphi = nphi,
        nw = nw, TCOV = TCOV
      )
      param <- newparam$par
      theta <- c(rep(0, nx), param[c(1:((ng - 1) * nx))])
      param <- c(param[-c(1:((ng - 1) * nx))], theta)
      #param[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))] <- exp(param[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))])
    }
    invH <- NULL
    if (hessian == TRUE) {
      invH <- newparam$invhessian
      SE <- sqrt(diag(invH))
      if (nx == 1) {
        sdtmp <- deltaTheta(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        sdbase <- deltaThetaBase(theta, invH[1:(ng - 1), 1:(ng - 1)], X, ng - 1)
        SE <- c(SE[-c(1:((ng - 1) * nx))], sdbase, sdtmp)
      } else {
        SE <- c(SE[-c(1:((ng - 1) * nx))], rep(NA, nx), SE[c(1:((ng - 1) * nx))])
      }
      #SE[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))] <- SE[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))] * param[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))]
    }
  }
  if (hessian == TRUE) {
    if (nx == 1 & Method == "L") {
      paramtmp <- c(param[1:(length(param) - ng * nx)], exp(theta) / sum(exp(theta)))
      prob <- 1 - stats::pt(abs(paramtmp / SE), n * period - 1) + stats::pt(-abs(paramtmp / SE), n * period - 1)
    } else {
      prob <- 1 - stats::pt(abs(param / SE), n * period - 1) + stats::pt(-abs(param / SE), n * period - 1)
    }
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = param / SE,
      Prob = prob
    )
  } else {
    SE <- rep(NA, length(param))
    d <- data.frame(
      Estimate = param,
      StandardError = SE,
      TValue = SE,
      Prob = SE
    )
  }


  namedegre <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", "Sextic", "Septic", "Octic")
  namebeta <- c()
  for (i in 1:length(nbeta)) {
    namebeta <- c(namebeta, namedegre[1:nbeta[i]])
  }
  namephi <- c()
  for (i in 1:length(nphi)) {
    namephi <- c(namephi, namedegre[1:nphi[i]])
  }
  nametheta <- rep(c("Intercept", colnames(X)[-1]), ng)
  if (nw == 0) {
    namedelta <- NULL
  } else {
    namedelta <- rep(paste0("TCOV", 1:nw), ng)
  }
  d.names <- c(namebeta, namephi, namedelta, nametheta)
  colnames(d) <- c("Estimate", "Std. Error", "T for H0 : Parameter=0", "Prob>|T|")
  beta <- param[1:(sum(nbeta))]
  phi <- param[(sum(nbeta) + 1):(sum(nbeta) + sum(nphi))]
  if (nw == 0) {
    delta <- NA
    theta <- param[-c(1:(sum(nbeta) + sum(nphi)))]
  } else {
    delta <- param[(sum(nbeta) + 1 + sum(nphi)):(sum(nbeta) + nw * ng + sum(nphi))]
    theta <- param[-c(1:(sum(nbeta) + +sum(nphi) + nw * ng))]
  }
  res <- list(
    beta = beta,
    phi = phi,
    delta = delta,
    theta = theta,
    sd = SE, tab = d, Model = "BETA",
    groups = ng, Names = d.names, Method = Method, Size = n,
    Likelihood = Likelihood(c(theta, beta, phi, delta),
      model = "BETA", method = Method,
      ng = ng, nx = nx, n = n, nbeta = nbeta, nphi = nphi, nw = nw,
      A = A, Y = Y, X = X, TCOV = TCOV
    ),
    Time = A[1, ], degre = degre - 1, degre.phi = nphi - 1, invH = invH
  )
  class(res) <- "Trajectory.BETA"
  return(res)
}
