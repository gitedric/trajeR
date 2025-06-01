# LOGIT
LLOGIT <- function(x, ng, nx, n, A, Y, X, nbeta,
                   nw, TCOV) {
  tour <- get_tour()
  storelik <- get_storelik()
  lik <- -likelihoodLOGIT_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
    nw = nw, TCOV = TCOV
  )
  if (storelik > lik) {
    message(sprintf("iter %3d value ", tour))
    message(sprintf("%.6f\n", lik))
    set_tour(tour + 1)
    set_storelik(lik)
  }
  return(lik)
}
difLLOGIT <- function(x, ng, nx, n, A, Y, X, nbeta,
                      nw, TCOV) {
  return(-difLLOGIT_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
    nw = nw, TCOV = TCOV
  ))
}

# CNORM
LCNORM <- function(x, ng, nx, n, A, Y, X, nbeta,
                   nw, TCOV, ymin, ymax, ssigma) {
  tour <- get_tour()
  storelik <- get_storelik()
  lik <- -Likelihoodalpha_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
    ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV, ssigma = ssigma
  )
  if (storelik > lik) {
    message(sprintf("iter %3d value ", tour))
    message(sprintf("%.6f\n", lik))
    set_tour(tour + 1)
    set_storelik(lik)
  }
  return(lik)
}
difLCNORM <- function(x, ng, nx, n, A, Y, X, nbeta,
                      nw, TCOV, ymin, ymax, ssigma) {
  return(-difLalpha_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
    ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV, ssigma = ssigma
  ))
}
# CNORM same sigma
difLCNORMss <- function(x, ng, nx, n, A, Y, X, nbeta,
                        nw, TCOV, ymin, ymax, ssigma) {
  return(-difLalphaunique_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
    ymin = ymin, ymax = ymax, nw = nw, TCOV = TCOV, ssigma = ssigma
  ))
}

# ZIP
LZIP <- function(x, ng, nx, n, A, Y, X, nbeta, nnu,
                 nw, TCOV) {
  tour <- get_tour()
  storelik <- get_storelik()
  lik <- -likelihoodZIP_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta, nnu = nnu,
    nw = nw, TCOV = TCOV
  )
  if (storelik > lik) {
    message(sprintf("iter %3d value ", tour))
    message(sprintf("%.6f\n", lik))
    set_tour(tour + 1)
    set_storelik(lik)
  }
  return(lik)
}
difLZIP <- function(x, ng, nx, n, A, Y, X, nbeta, nnu,
                    nw, TCOV) {
  return(-difLZIP_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
    nnu = nnu, nw = nw, TCOV = TCOV
  ))
}

# BETA
LBETA <- function(x, ng, nx, n, A, Y, X, nbeta, nphi,
                  nw, TCOV) {
  tour <- get_tour()
  storelik <- get_storelik()
  lik <- -LikelihoodBETA_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta, nphi = nphi,
    nw = nw, TCOV = TCOV
  )
  if (storelik > lik) {
    message(sprintf("iter %3d value ", tour))
    message(sprintf("%.6f\n", lik))
    set_tour(tour + 1)
    set_storelik(lik)
  }
  return(lik)
}
difLBETA <- function(x, ng, nx, n, A, Y, X, nbeta, nphi,
                     nw, TCOV) {
  return(-difLBETA_cpp(x,
    ng = ng, nx = nx, n = n, A = A, Y = Y, X = X, nbeta = nbeta,
    nphi = nphi, nw = nw, TCOV = TCOV
  ))
}
