#################################################################################
# pass value inside environment
#################################################################################
get_tour <- function() {
  my_env$tour
}
set_tour <- function(value) {
  old <- my_env$tour
  my_env$tour <- value
  invisible(old)
}
get_storelik <- function() {
  my_env$storelik
}
set_storelik <- function(value) {
  old <- my_env$storelik
  my_env$storelik <- value
  invisible(old)
}
#################################################################################
# Piik
#################################################################################
piik <- function(theta, i, k, ng, X) {
  ntheta <- ncol(X)
  tmp <- exp(sapply(1:ng, function(s) {
    theta[((s - 1) * ntheta + 1):(s * ntheta)] %*% X[i, ]
  }))
  return(tmp[k] / sum(tmp))
}
#################################################################################
# Delta method for theta
#################################################################################
deltaTheta <- function(theta, Ht, X, ng) {
  Dg <- matrix(rep(0, ng**2), ncol <- ng)
  for (k in 1:ng) {
    for (l in 1:ng) {
      if (k == l) {
        Dg[k, l] <- piik(theta, 1, k + 1, ng + 1, X) * (1 - piik(theta, 1, k + 1, ng + 1, X))
      } else {
        Dg[k, l] <- -piik(theta, 1, k + 1, ng + 1, X) * piik(theta, 1, l + 1, ng + 1, X)
      }
    }
  }
  sqrt(diag(Dg %*% Ht %*% t(Dg)))
}

deltaThetaBase <- function(theta, Ht, X, ng) {
  Dg <- c()
  prob <- sapply(1:(ng + 1), function(s) {
    piik(theta, 1, s, ng + 1, X)
  })[-1]
  for (k in 1:ng) {
    Dg <- c(Dg, -prob[k] * (1 - prob[k]) + sum(prob[k] * prob[-k]))
  }
  Dg <- matrix(Dg, nrow = 1)
  sqrt(diag(Dg %*% Ht %*% t(Dg)))
}
#################################################################################
# likelihood
#################################################################################
Likelihood <- function(param, model, method, ng, nx, n, nbeta, nw, A, Y, X, TCOV, ymin = NULL, ymax = NULL, nnu = NULL, fct = NULL, nphi = NULL) {
  Y <- data.matrix(Y)
  A <- data.matrix(A)
  if (model == "CNORM") {
    if (method == "L" | nx != 1) {
      # a = likelihoodCNORM_cpp(param[-c(1:nx)], ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw)
      a <- likelihoodCNORM_cpp(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw)
    } else {
      a <- likelihoodEM_cpp(n, ng, nbeta,
        pi = param[1:(ng)],
        beta = param[(ng + 1):(ng + sum(nbeta))],
        sigma = param[(ng + 1 + sum(nbeta)):(ng + sum(nbeta) + ng)],
        delta = param[(ng + sum(nbeta) + ng + 1):(ng + sum(nbeta) + ng + nw * ng)],
        A, Y, ymin, ymax, TCOV,
        nw
      )
    }
  } else if (model == "LOGIT") {
    if (method == "L" | nx != 1) {
      a <- likelihoodLOGIT_cpp(param[-c(1:nx)], ng, nx, n, nbeta, A, Y, X, TCOV, nw)
    } else {
      a <- likelihoodEMLOGIT_cpp(n, ng, nbeta,
        beta = param[(ng):(ng + sum(nbeta) - 1)],
        pi = c(1 - sum(param[1:(ng - 1)]), param[1:(ng - 1)]),
        A, Y, TCOV,
        delta = param[-c(1:(ng + sum(nbeta) - 1))], nw
      )
    }
  } else if (model == "ZIP") {
    if (method == "L" | nx != 1) {
      a <- likelihoodZIP_cpp(param[-c(1:nx)], ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)
    } else {
      a <- likelihoodEMZIP_cpp(n, ng, nbeta, nnu,
        beta = param[(ng + 1):(ng + sum(nbeta))],
        nu = param[(ng + sum(nbeta) + 1):(ng + sum(nbeta) + sum(nnu))],
        pi = param[1:(ng)],
        A, Y, TCOV,
        delta = param[-c(1:(ng + sum(nbeta) + sum(nnu)))], nw
      )
    }
  } else if (model == "POIS") {
    if (method == "L" | nx != 1) {
      a <- likelihoodPois_cpp(param[-c(1:nx)], ng, nx, nbeta, n, A, Y, X, TCOV, nw)
    } else {
      a <- likelihoodEMZIP_cpp(n, ng, nbeta, nnu,
        beta = param[(ng + 1):(ng + sum(nbeta))],
        nu = param[(ng + sum(nbeta) + 1):(ng + sum(nbeta) + sum(nnu))],
        pi = param[1:(ng)],
        A, Y, TCOV,
        delta = param[-c(1:(ng + sum(nbeta) + sum(nnu)))], nw
      )
    }
  } else if (model == "BETA") {
    if (method == "L" | nx != 1) {
      a <- LikelihoodBETA_cpp(param[-c(1:nx)], ng, nx, nbeta, nphi, n, A, Y, X, TCOV, nw)
    } else {

    }
  } else {
    a <- LikelihoodNL(param, ng, nx, nbeta, n, A, Y, X, TCOV, fct = fct)
  }
  return(a)
}
#################################################################################
# compute the value of Wit for i and t and given k
#################################################################################
Wit <- function(TCOV, period, delta, nw, i, t, k) {
  if (nw == 0) {
    return(0)
  } else {
    return(sum(delta[[k]] * TCOV[i, seq(from = t, to = t + (nw - 1) * period, by = period)]))
  }
}
# compute the value of Wit for i and t and given k for the EM algorithm
WitEM <- function(TCOV, period, delta, nw, i, t, k, ndeltacum) {
  if (nw == 0) {
    return(0)
  } else {
    return(sum(delta[(ndeltacum[k] + 1):(ndeltacum[k + 1])] * TCOV[i, seq(from = t, to = t + (nw - 1) * period, by = period)]))
  }
}
#################################################################################
# Calculate the probability of membership for each data
#################################################################################
#' Membership's probabilities
#'
#' \code{GroupProb} calculate the membership probability of each value of the data.
#'
#' @param Obj Trajectory's object. A trajectory object that is return by \code{trajeR} function.
#' @param Y Matrix. A real matrix. The data.
#' @param A Matrix. A real matrix. The time variable.
#' @param TCOV Matrix. A real matrix. Optional, by default the value is NULL. It contained the time dependent covariate.
#' @param X Matrix. A real matrix. Optional, by default the value is NULL. It contained a covariate that modify the probability membership.
#'
#' @return a real matrix. For each individual i in the data, this matrix contained the membership probability of each group.
#'
#' @export
#'
#' @examples
#' data <- read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data <- as.matrix(data)
#' sol <- trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2, 2), Model = "CNORM", Method = "EM")
#' GroupProb(sol, Y = data[, 2:6], A = data[, 7:11])
GroupProb <- function(Obj, Y, A, TCOV = NULL, X = NULL) {
  Y <- data.matrix(Y)
  A <- data.matrix(A)
  n <- Obj$Size
  ng <- Obj$groups
  nbeta <- Obj$degre + 1
  ymin <- Obj$min
  ymax <- Obj$max
  betatmp <- Obj$beta
  j <- 1
  beta <- list()
  for (i in seq_along(nbeta)) {
    beta[[i]] <- betatmp[j:sum(nbeta[1:i])]
    j <- sum(nbeta[1:i]) + 1
  }
  delta <- Obj$delta
  if (any(is.na(delta))) {
    nw <- 0
    delta <- rep(list(0), ng)
  } else {
    nw <- length(delta) / ng
  }
  theta <- Obj$theta
  if (is.null(X)) {
    X <- cbind(rep(1, n))
  } else {
    X <- cbind(rep(1, n), X)
  }
  res <- c()
  if (Obj$Model == "LOGIT") {
    if (Obj$Method == "L") {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          piik(theta, 1, s, ng, X) * gkLOGIT_cpp(beta, i - 1, s - 1, nbeta, A, Y, TCOV, delta, nw)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    } else {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          theta[s] * gkLOGIT_cpp(beta, i - 1, s - 1, nbeta, A, Y, TCOV, delta, nw)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    }
  } else if (Obj$Model == "CNORM") {
    sigma <- Obj$sigma
    if (Obj$Method == "L") {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          piik(theta, i, s, ng, X) * gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    } else {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          theta[s] * gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    }
  } else if (Obj$Model == "ZIP") {
    nutmp <- Obj$nu
    nnu <- Obj$degre.nu + 1
    j <- 1
    nu <- list()
    for (i in seq_along(nnu)) {
      nu[[i]] <- nutmp[j:sum(nnu[1:i])]
      j <- sum(nnu[1:i]) + 1
    }
    if (Obj$Method == "L") {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          piik(theta, i, s, ng, X) * gkZIP_cpp(beta, nu, i - 1, s - 1, nbeta, nnu, A, Y, TCOV, delta, nw)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    } else {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          theta[s] * gkZIP_cpp(beta, nu, i - 1, s - 1, nbeta, nnu, A, Y, TCOV, delta, nw)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    }
  } else if (Obj$Model == "BETA") {
    phitmp <- Obj$phi
    nphi <- Obj$degre.phi + 1
    j <- 1
    phi <- list()
    for (i in seq_along(nphi)) {
      phi[[i]] <- phitmp[j:sum(nphi[1:i])]
      j <- sum(nphi[1:i]) + 1
    }
    if (Obj$Method == "L") {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          piik(theta, i, s, ng, X) * gkBETA_cpp(beta, phi, i - 1, s - 1, nbeta, nphi, A, Y, TCOV, delta, nw)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    }
  } else {
    if (Obj$Method == "L") {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          piik(theta, i, s, ng, X) * gkNL(beta, sigma, i, s, TCOV, A, Y)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    } else {
      for (i in 1:n) {
        tmp <- sapply(1:ng, function(s) {
          theta[s] * gkNL(beta, sigma, i, s, TCOV, A, Y)
        })
        res <- rbind(res, tmp / sum(tmp))
      }
    }
  }
  colnames(res) <- paste0("Gr", 1:ng)
  return(res)
}
#################################################################################
# Function to find theta in the calculus of the membership probability with predictors
#################################################################################
ftheta <- function(theta, taux, X, n, ng, period) {
  nx <- ncol(X)
  a <- 0
  for (i in 1:n) {
    for (k in 1:ng) {
      tmp <- sapply(1:ng, function(s) {
        theta[((s - 1) * nx + 1):(s * nx)] %*% X[i, ]
      })
      a <- a + taux[i, k] * (tmp[k] - log(sum(exp(tmp))))
    }
  }
  return(a)
}
difftheta <- function(theta, taux, X, n, ng, period) {
  nx <- ncol(X)
  thetas <- c()
  for (k in 1:ng) {
    for (l in 1:nx) {
      a <- 0
      for (i in 1:n) {
        tmp <- exp(sapply(1:ng, function(s) {
          theta[((s - 1) * nx + 1):(s * nx)] %*% X[i, ]
        }))
        a <- a + X[i, l] * (taux[i, k] - tmp[k] / sum(tmp))
      }
      thetas <- c(thetas, a)
    }
  }
  return(thetas)
}
findtheta <- function(theta, taux, X, n, ng, nx, period, EMIRLS, refgr) {
  if (EMIRLS == TRUE) {
    newtheta <- c()
    thetaIRLS <- theta[-c(((refgr - 1) * nx + 1):(nx * refgr))]
    thetaIRLS <- thetaIRLS - theta[c(((refgr - 1) * nx + 1):(nx * refgr))]
    ind <- 1:ng
    ind <- ind[-refgr]
    precIRLS <- 1
    while (any(abs(precIRLS) > 10**(-6))) {
      Xng <- matrix(rep(0, n * (ng - 1) * nx * (ng - 1)), ncol = nx * (ng - 1))
      tmp2 <- c()
      PIw <- c()
      tmp4 <- c()
      kind <- 0
      for (k in ind) {
        kind <- kind + 1
        PIwtmp <- c()
        for (l in ind) {
          tmp1 <- c()
          if (k == l) {
            for (i in 1:n) {
              tmpPiik <- piik(c(rep(0, nx), thetaIRLS), i, k, ng, X)
              tmp1 <- c(tmp1, tmpPiik * (1 - tmpPiik))
              tmp4 <- c(tmp4, tmpPiik)
            }
          } else {
            for (i in 1:n) {
              tmp1 <- c(tmp1, -piik(c(rep(0, nx), thetaIRLS), i, k, ng, X) * piik(c(rep(0, nx), thetaIRLS), i, l, ng, X))
            }
          }
          PIwtmp <- cbind(PIwtmp, diag(tmp1))
        }
        PIw <- rbind(PIw, PIwtmp)
        Xng[((kind - 1) * n + 1):(kind * n), ((kind - 1) * nx + 1):(kind * nx)] <- X
        tmp2 <- c(tmp2, taux[, k])
      }
      rm(PIwtmp)
      Z <- matrix(tmp2, ncol = 1)
      PIm <- matrix(tmp4, ncol = 1)
      newthetaIRLS <- as.vector(solve(t(Xng) %*% PIw %*% Xng, t(Xng) %*% (PIw %*% Xng %*% thetaIRLS + Z - PIm), tol = 10**(-20)))
      precIRLS <- c(thetaIRLS - newthetaIRLS)
      thetaIRLS <- newthetaIRLS
    }
    newtheta <- rep(0, ng * nx)
    newtheta[-c(((refgr - 1) * nx + 1):(nx * refgr))] <- thetaIRLS
  } else {
    newtheta <- stats::optim(
      par = theta, fn = ftheta, gr = difftheta,
      taux = taux, X = X, n = n, ng = ng, period = period,
      control = list(fnscale = -1), hessian = FALSE
    )$par
  }
  return(newtheta)
}


#####################################
# Confidence Intervall
#####################################
#' Compute Confidence Intervals for Predicted Trajectories
#'
#' @description
#' Calculates predicted values and confidence intervals for a given model object
#' at specified time points, based on model. The function computes standard errors and confidence bounds
#' for the predictions.
#'
#' @param Obj An object containing model parameters, including:
#'   \itemize{
#'     \item \code{degre}: Vector of polynomial degrees for each group.
#'     \item \code{varcov}: Variance-covariance matrix of the model coefficients.
#'     \item \code{groups}: Number of groups in the model.
#'     \item \code{Model}: Character string specifying the model type ("CNORM" or "LOGIT").
#'     \item \code{beta}: Vector of model coefficients.
#'   }
#' @param newtime Numeric vector of time points at which to compute predictions
#'   and confidence intervals.
#' @param alpha Numeric value specifying the significance level for the confidence
#'   intervals (default is 0.05, corresponding to 95% confidence intervals).
#'
#' @return A list of class \code{"Trajectory.predict"} containing:
#'   \itemize{
#'     \item \code{newtime}: The input time points.
#'     \item \code{Ypred}: Matrix of predicted values for each group at each time point.
#'     \item \code{SEYpred}: Matrix of standard errors for the predicted values.
#'     \item \code{IC.inf}: Matrix of lower bounds of the confidence intervals.
#'     \item \code{IC.max}: Matrix of upper bounds of the confidence intervals.
#'     \item \code{groups}: Number of groups from the input object.
#'   }
#'
#' @examples
#' \donttest{
#' # Example with a hypothetical model object
#' model <- list(
#'   degre = c(2, 2),
#'   varcov = diag(4),
#'   groups = 2,
#'   Model = "CNORM",
#'   beta = c(1, 0.5, 0.2, 2, 0.3, 0.1)
#' )
#' newtime <- seq(0, 1, by = 0.1)
#' result <- confidenceInt(model, newtime, alpha = 0.05)
#' print(result)
#' }
#'
#' @export
confidenceInt <- function(Obj, newtime, alpha = 0.05) {
  nbeta <- Obj$degre + 1
  nbetacum <- cumsum(c(0, nbeta))

  Ypred <- c()
  SEYpred <- c()
  for (i in 1:Obj$groups) {
    ind <- Obj$groups - 1 + (nbetacum[i] + 1):(nbetacum[i + 1])
    varbeta <- Obj$varcov[ind, ind]
    SEYpred <- cbind(SEYpred, sapply(newtime, function(s) matrix((s)**(0:Obj$degre[i]), nrow = 1) %*% varbeta %*% t(matrix((s)**(0:Obj$degre[i]), nrow = 1))))
    if (Obj$Model == "CNORM") {
      Ypred <- cbind(Ypred, sapply(newtime, function(s) sum(Obj$beta[(nbetacum[i] + 1):(nbetacum[i + 1])] * (s)**(0:Obj$degre[i]))))
    } else if (Obj$Model == "LOGIT") {
      Ypred <- cbind(Ypred, sapply(newtime, function(s) 1 / (1 + exp(-(sum(Obj$beta[(nbetacum[i] + 1):(nbetacum[i + 1])] * (s)**(0:Obj$degre[i])))))))
    }
  }
  if (Obj$Model == "LOGIT") {
    SEYpred <- (Ypred * (1 - Ypred))**2 * SEYpred
  }

  IC.inf <- Ypred - qnorm(1 - alpha / 2) * sqrt(SEYpred)
  IC.max <- Ypred + qnorm(1 - alpha / 2) * sqrt(SEYpred)

  res <- list(
    newtime = newtime,
    Ypred = Ypred,
    SEYpred = SEYpred,
    IC.inf = IC.inf,
    IC.max = IC.max,
    groups = Obj$groups
  )
  class(res) <- "Trajectory.predict"
  return(res)
}

