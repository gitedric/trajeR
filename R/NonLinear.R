#################################################################################
#
# Likelihood
#
#################################################################################
#################################################################################
# Function f(beta, ait)
#################################################################################
#'  Function fait
#'  
#' @param betak Vector of integer.
#' @param i Integer. 
#' @param t Real.
#' @param A Matrix of real.
#' @param TCOV Matrix of real.
#' @param fct Function.
#' @param diffct Function.
#' @return real. Compute the value of the function fct for individual i, time t and group k.
#' @export
fait <- function(betak, i, t, A, TCOV, fct, diffct){
  return(fct(A[i,t], betak, TCOV))
}
#' Differential
#' 
#' @inheritParams fait
#' @return real. Compute the value of the differential function fct for individual i, time t and group k.
#' @export
diffaitbeta <- function(betak, i, t, A, TCOV, fct, diffct){
    return(diffct(A[i,t], betak, TCOV))
}
#################################################################################
# Function g
#################################################################################
gkNL <- function(beta, sigma, i, k, TCOV, A, Y, fct){
  period = ncol(A)
  muikt = sapply(1:period, function(s){
    fait(beta[[k]], i, s, A, TCOV, fct)
    })
  return(max(prod(stats::dnorm((Y[i, ]-muikt)/sigma[k])/sigma[k]), 10**(-300)))
}
#################################################################################
# gk with exponential parameterization
#################################################################################
gkalphaNL <- function(beta, alpha, i, k, TCOV, A, Y, fct){
  period = ncol(A)
  sigma = exp(alpha)
  muikt = sapply(1:period, function(s){
    fait(beta[[k]], i, s, A, TCOV, fct)
    })
  return(max(prod(stats::dnorm((Y[i, ]-muikt)/sigma[k])/sigma[k]), 10**(-300)))
}
#################################################################################
#
# Differential of function by all the parameters
#
#################################################################################
#################################################################################
# dif likelihood theta
#################################################################################
difLthetakalphaNL <- function(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct){
  theta = param[1:(ng*nx)]
  betatmp = param[((ng*nx)+1):(ng*nx+sum(nbeta))]
  alpha = param[-c(1:(ng*nx+sum(nbeta)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  thetas = c()
  for (l in 1:nx){
    a = 0
    for (i in 1:n){
      tmp = exp(sapply(1:ng*nx,function(s){theta[((s-1)*nx+1):(s*nx)]%*%X[i,]}))
      tmp1 = tmp[k]*sum(sapply(1:ng, function(s){tmp[s]*(gkalphaNL(beta, alpha, i, k, TCOV, A, Y, fct) - gkalphaNL(beta, alpha, i, s, TCOV, A, Y, fct))}))
      tmp2 = sum(sapply(1:ng, function(s){tmp[s]*gkalphaNL(beta, alpha, i, s, TCOV, A, Y, fct)}))
      a = a + X[i,l]*tmp1/(sum(tmp)*tmp2)
    }
    thetas = c(thetas, a)
  }
  return(thetas)
}
#################################################################################
# dif likelihood beta with reparameterization alpha
#################################################################################
difLbetakalphaNL <- function(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct, diffct){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[((ng*nx)+1):(ng*nx+sum(nbeta))]
  alpha = param[-c(1:(ng*nx+sum(nbeta)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  betas = c()
  m = rep(0,period)
  for (l in 1:nbeta[k]){
    a = 0
    for (i in 1:n){
      difbkl = 0
      ind = 1:period
      muikt = sapply(1:period, function(s){
        fait(beta[[k]], i, s, A, TCOV, fct)
      })
      m = (Y[i, ]-muikt)*exp(-alpha[k])
      tind = 1
      for (t in ind){
        der = diffaitbeta(beta[[k]], i, t, A, TCOV, fct, diffct)
        difbkl = difbkl + der[l]*exp(-2*alpha[k])*m[t]*stats::dnorm(m[t])*prod(exp(-alpha[k])*stats::dnorm(m[-tind]))
        tind = tind +1
      }
      py = prod(exp(-alpha[k])*stats::dnorm(m[ind]))
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkalphaNL(beta, alpha, i, s, TCOV, A, Y, fct)}))
      a = a + piik(theta, i, k, ng, X)/s*difbkl
    }
    betas = c(betas, a)
  }
  return(betas)
}
#################################################################################
# dif likelihood sigma with reparameterization alpha
#################################################################################
difLsigmakalphaNL <- function(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[((ng*nx)+1):(ng*nx+sum(nbeta))]
  alpha = param[-c(1:(ng*nx+sum(nbeta)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  alphas = c()
  a = 0
  m = rep(0,period)
  for (i in 1:n){
    difbkl = 0
    ind = 1:period
    muikt = sapply(1:period, function(s){
      fait(beta[[k]], i, s, A, TCOV, fct)
    })
    m = (Y[i,]-muikt)*exp(-alpha[k])
    tind = 1
    for (t in ind){
      difbkl = difbkl + (-1+m[t]**2)*exp(-alpha[k])*stats::dnorm(m[t])*prod(exp(-alpha[k])*stats::dnorm(m[-tind]))
      tind = tind + 1
    }
    py = prod(exp(-alpha[k])*stats::dnorm(m))
    s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkalphaNL(beta, alpha, i, s, TCOV, A, Y, fct)}))
    a = a + piik(theta, i, k, ng, X)/s*difbkl
  }
  alphas = c(alphas, a)
  return(alphas)
}
#################################################################################
# dif likelihood sigma with reparameterization alpha same sigma
#################################################################################
difLsigmaalphauniqueNL <- function(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[((ng*nx)+1):(ng*nx+sum(nbeta))]
  alpha = param[-c(1:(ng*nx+sum(nbeta)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  alphas = c()
  a = 0
  m = rep(0,period)
  for (i in 1:n){
    s1 = 0
    s2 = 0
    for (k in 1:ng){
      difbkl = 0
      ind = 1:period
      muikt = sapply(1:period, function(s){
        fait(beta[[k]], i, s, A, TCOV)
      })
        m = (Y[i,]-muikt)*exp(-alpha[k])
        tind = 1
          for (t in ind){
            difbkl = difbkl + (-1+m[t]**2)*exp(-alpha[k])*stats::dnorm(m[t])*prod(exp(-alpha[k])*stats::dnorm(m[-tind]))
            tind = tind + 1
          }
        py = prod(exp(-alpha[k])*stats::dnorm(m))
      s1 = s1 + sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkalphaNL(beta, alpha, i, s, TCOV, A, Y, fct)}))
      s2 = s2 + piik(theta, i, k, ng, X)*difbkl
    }
    a = a + s2/s1
  }
  alphas = rep(a, ng)
  return(alphas)
}
#################################################################################
# dif likelihood sigma with exponential reparameterization
#################################################################################
difLalphaNL <- function(param, ng, nx, nbeta, n, A, Y, X, TCOV, fct, diffct ){
  out =c()
  for (k in 1:ng){
    out = c(out, difLthetakalphaNL(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct))
  }
  for (k in 1:ng){
    out = c(out, difLbetakalphaNL(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct, diffct ))
  }
  for (k in 1:ng){
    out = c(out, difLsigmakalphaNL(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct))
  }
  return(out)
}
#################################################################################
# dif likelihood sigma with exponential reparameterization and unique sigma
#################################################################################
difLalphauniqueNL <- function(param, ng, nx, nbeta, n, A, Y, X, TCOV, fct, diffct ){
  out =c()
  for (k in 1:ng){
    out = c(out, difLthetakalphaNL(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct))
  }
  for (k in 1:ng){
    out = c(out, difLbetakalphaNL(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct ))
  }
  out = c(out, difLsigmaalphauniqueNL(param, k, ng, nx, nbeta, n, A, Y, X, TCOV, fct))
  return(out)
}
#################################################################################
# likelihood sigma with exponential reparameterization
#################################################################################
LikelihoodalphaNL <- function(param, ng, nx, nbeta, n, A, Y, X, TCOV, fct, diffct ){
  theta = param[1:(ng*nx)]
  betatmp = param[((ng*nx)+1):(ng*nx+sum(nbeta))]
  alpha = param[-c(1:(ng*nx+sum(nbeta)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  a = 0
  for (i in 1:n){
    a = a + log(sum(sapply(1:ng, function(s){
      piik(theta, i, s, ng, X)*gkalphaNL(beta, alpha, i, s, TCOV, A, Y, fct)
    }))
    )
  }
  return(a)
}
#################################################################################
# likelihood
#################################################################################
LikelihoodNL <- function(param, ng, nx, nbeta, n, A, Y, X, TCOV, fct){
  theta = param[1:(ng*nx)]
  betatmp = param[((ng*nx)+1):(ng*nx+sum(nbeta))]
  sigma = param[-c(1:(ng*nx+sum(nbeta)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  a = 0
  for (i in 1:n){
    a = a + log(sum(sapply(1:ng, function(s){
      piik(theta, i, s, ng, X)*gkNL(beta, sigma, i, s, TCOV, A, Y, fct)
    }))
    )
  }
  return(a)
}
#################################################################################
#
# EM algorithm
#
#################################################################################
#################################################################################
# Function rate
#################################################################################
ftauxNL <- function(pi, beta, sigma, ng, nbeta, n, nw, nx, A, Y, X, TCOV, fct){
  betatmp = beta
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  tmp2 = c()
  if (nx == 1){
    for (k in 1:ng){
      tmp1 = sapply(1:n, function(s){gkNL(beta, sigma, s, k, TCOV, A, Y, fct)})
      tmp2 = cbind(tmp2, tmp1*pi[k])
    }
  }else{
    for (k in 1:ng){
      tmp1 = sapply(1:n, function(s){gkNL(beta, sigma, s, k, TCOV, A, Y, fct)})
      tmp2 = cbind(tmp2, tmp1*piik(pi, i, k, ng, X))
    }
  }
  tmp3 =c()
  for (k in 1:ng){
    tmp3 = cbind(tmp3, 1/(1+ (rowSums(tmp2)-tmp2[,k])/tmp2[,k]))
  }
  return(tmp3)
}
fctbeta <- function(betak, Y, A, TCOV, zk, period, n, nbetak, fct, diffct){
  tmp1 = c()
  tmp2 = c()
  tmp3 = c()
  for (i in 1:n){
    for (t in 1:period){
      tmp1 = c(tmp1, Y[i,t])
      tmp3 =c(tmp3, fait(betak, i, t, A, TCOV, fct))
    }
    tmp2 = c(tmp2, rep(zk[i], period))
  }
  Fa = matrix(tmp3, ncol = 1)
  Z = diag(tmp2)
  Ynls = matrix(tmp1, ncol = 1)
  return(sqrt(Z)%*%(Ynls-Fa))
}
fctbeta.jac <- function(betak, Y, A, TCOV, zk, period, n, nbetak, fct, diffct){
  tmp1 = c()
  tmp2 = c()
  tmp3 = c()
  for (i in 1:n){
    for (t in 1:period){
      tmp1 = c(tmp1, Y[i,t])
      tmp3 = rbind(tmp3, diffaitbeta(betak, i, t, A, TCOV, fct, diffct))
    }
    tmp2 = c(tmp2, rep(zk[i], period))
  }
  Fa = matrix(tmp3, ncol = nbetak)
  Z = diag(tmp2)
  return(-sqrt(Z)%*%Fa)
}
#################################################################################
# Likelihood calculation for EM algorithm
#################################################################################
likelihoodEMNL <- function(n, ng, nbeta, beta, sigma, pi, A, Y, TCOV, nw, fct){
  likeli = 0
  betatmp = beta
  sigma = sigma
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  for (i in 1:n){
    likeli = likeli + log(sum(sapply(1:ng, function(s){
      pi[s]*gkNL(beta, sigma, i, s, TCOV, A, Y, fct)
    }))
    )
  }
  return(likeli)
}
#################################################################################
# EM different sigma
#################################################################################
EMNL <- function(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, itermax, EMIRLS, fct, diffct , nls.lmiter){
  period = ncol(A)
  if (nx ==1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    sigma = param[-c(1:(ng+sum(nbeta)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = param[-c(1:(ng*nx+sum(nbeta)))]
  }
  nbetacum = cumsum(c(0, nbeta))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      message(paste(-likelihoodEMNL(n, ng, nbeta, beta, sigma, pi, A, Y, TCOV, nw, fct), "\n"))
    }else{
      message(paste(-LikelihoodalphaNL(c(pi, beta, log(sigma)), ng, nx, nbeta, n, A, Y, X, TCOV, fct), "\n"))
    }
    # E-step
    zk = ftauxNL(pi, beta, sigma, ng, nbeta, n, nw, nx, A, Y, X, TCOV, fct)
    newbeta = c()
    newsigma = c()
    # M-step
    for (k in 1:ng){
      betatmp = minpack.lm::nls.lm(par = beta[(nbetacum[k]+1):(nbetacum[k+1])], fn = fctbeta, jac = fctbeta.jac,
                       Y = Y, A = A, TCOV = TCOV, zk = zk[ ,k], 
                       period = period, n =n, nbetak = nbeta[k], fct = fct, diffct = diffct,
                       control = minpack.lm::nls.lm.control(maxiter = nls.lmiter))$par
      newbeta = c(newbeta, betatmp)
      b = 0
      for (i in 1:n){
        Mtmp = c()
        for (t in 1:period){
          Mtmp = c(Mtmp, Y[i,t]-fait(beta[(nbetacum[k]+1):(nbetacum[k+1])], i, t, A, TCOV, fct))
        }
        b = b + zk[i,k]*(t(Mtmp)%*%Mtmp)
      }
      newsigma = c(newsigma, sqrt(b/(period*sum(zk[,k]))))
    }
    if (nx == 1){
      pi = colSums(zk)/n
      stoppi = 0
      refpi = 0
    }else{
      newpi = findtheta(pi, zk, X, n, ng, period, EMIRLS)
      stoppi = (pi-pi[1:ng])-(newpi-newpi[1:ng])
      refpi = rep(0, length(stoppi))
      pi = newpi
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newsigma, stoppi)-c(beta, sigma, refpi))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    sigma = newsigma
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, sigma)
  }else{
    param = c(pi, beta, sigma)
  }
  return(param)
}
#################################################################################
# EM same sigma
#################################################################################
EMNLSigmaunique <- function(param, ng, nx, nbeta, n, A, Y, X, TCOV, nw, itermax, EMIRLS, fct, diffct , nls.lmiter){
  period = ncol(A)
  if (nx ==1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    sigma = param[-c(1:(ng+sum(nbeta)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = param[-c(1:(ng*nx+sum(nbeta)))]
  }
  nbetacum = cumsum(c(0, nbeta))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      message(paste(-likelihoodEMNL(n, ng, nbeta, beta, sigma, pi, A, Y, TCOV, nw, fct), "\n"))
    }else{
      message(paste(-LikelihoodalphaNL(c(pi, beta, log(sigma)), ng, nx, nbeta, n, A, Y, X, TCOV, fct), "\n"))
    }
    # E-step
    zk = ftauxNL(pi, beta, sigma, ng, nbeta, n, nw, nx, A, Y, X, TCOV, fct)
    newbeta = c()
    newsigma = c()
    # M-step
    for (k in 1:ng){
      betatmp = minpack.lm::nls.lm(par = beta[(nbetacum[k]+1):(nbetacum[k+1])], fn = fctbeta, jac = fctbeta.jac,
                       Y =Y, A = A, TCOV = TCOV, zk = zk[ ,k],
                       period = period, n =n, nbetak = nbeta[k], fct = fct, diffct = diffct,
                       control = minpack.lm::nls.lm.control(maxiter = nls.lmiter))$par
      newbeta = c(newbeta, betatmp)
    }
    b = 0
    for (i in 1:n){
      Mtmp = c()
      for (t in 1:period){
        Mtmp = c(Mtmp, Y[i,t]-fait(beta[(nbetacum[k]+1):(nbetacum[k+1])], i, t, A, TCOV, fct))
      }
      b = b + zk[i,k]*(t(Mtmp)%*%Mtmp)
    }
    newsigma = sqrt(b / (period * sum(zk)))
    newsigma = rep(newsigma, ng)
    if (nx == 1){
      pi = colSums(zk)/n
      stoppi = 0
      refpi = 0
    }else{
      newpi = findtheta(pi, zk, X, n, ng, period, EMIRLS)
      stoppi = (pi-pi[1:ng])-(newpi-newpi[1:ng])
      refpi = rep(0, length(stoppi))
      pi = newpi
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newsigma, stoppi)-c(beta, sigma, refpi))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    sigma = newsigma
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, sigma)
  }else{
    param = c(pi, beta, sigma)
  }
  return(param)
}
#################################################################################
# Compute of matrix information
#################################################################################
#################################################################################
# Definition of function
#################################################################################
dif2fait <- function(betak, i, t, A, TCOV, fct, diffct){
    return(numDeriv::jacobian(func = diffct, x = betak, t = A[i,t], TCOV = TCOV))
}
#################################################################################
# matrix differential of beta**2 for t and i
#################################################################################
mbetaNL = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct){
  res = matrix(rep(0,(sum(nbeta)**2)), ncol = sum(nbeta))
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nbeta[k]){
      for (lp in 1:nbeta[k]){
        tmp = c(tmp, taux[i,k]/sigma[k]**2*(dif2fait(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct, diffct)[l, lp]*(Y[i,t]-fait(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct))-diffaitbeta(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct, diffct)[l]*diffaitbeta(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct, diffct)[lp]))
      }
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), (nbetacum[k]+1):(nbetacum[k+1])] = matrix(tmp, ncol = nbeta[k], byrow = TRUE)
  }
  return(res)
}
##################################################################################
# matrix differential of  beta sigma
#################################################################################
mbetasigmaNL = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct){
  res = matrix(rep(0,ng*sum(nbeta)), ncol = ng)
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nbeta[k]){
      tmp = c(tmp, -2*taux[i,k]*diffaitbeta(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct, diffct)[l]*(Y[i,t]-fait(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct))/sigma[k]**3)
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), k] = tmp
  }
  return(res)
}
##################################################################################
# matrix differential of  beta sigma
#################################################################################
mdeltasigmaNL = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct){
  res = matrix(rep(0,ng*nw*ng), ncol = ng)

  return(res)
}
##################################################################################
# matrix differential of  beta delta
#################################################################################
mbetadeltaNL = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct){
  res = matrix(rep(0,ng*sum(nbeta)), ncol = nw)
  for (k in 1:nw){
    tmp = c()
    for (l in 1:nbeta[k]){
      tmp = c(tmp, -2*taux[i,k]*diffaitbeta(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct, diffct)[l]*(Y[i,t]-fait(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct))/sigma[k]**3)
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), k] = tmp
  }
  return(res)
}
##################################################################################
# matrix differential of  sigma**2 for t and i
#################################################################################
msigmaNL = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct){
  tmp = c()
  for (k in 1:ng){
    tmp = c(tmp, -taux[i,k]*(-sigma[k]**2+3*(Y[i,t]-fait(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct))**2)/sigma[k]**4)
  }
  return(diag(tmp))
}
#################################################################################
# matrix differential of delta**2 for t and i
#################################################################################
mdeltaNL = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct){
  res = matrix(rep(0,(sum(nbeta)**2)), ncol = nw*ng)

  return(res)
}
##################################################################################
# definition of function Bikl and Sk
##################################################################################
BiklNL <- function(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct){
  tmp = 0
  for (t in 1:period){
    tmp = tmp + diffaitbeta(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct, diffct)[l]*(Y[i,t]-fait(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct))/sigma[k]**2
  }
  return(tmp)
}
SikNL <- function(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct){
  tmp = 0
  for (t in 1:period){
    tmp = tmp - (sigma[k]**2-(Y[i,t]-fait(beta[(nbetacum[k]+1):nbetacum[k+1]], i, t, A, TCOV, fct))**2)/sigma[k]**3
  }
  return(tmp)
}
DiklNL <- function(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct){
  tmp = 0

  return(tmp)
}
# matrix cov pi betak
covPiBetakNL <- function(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, fct, diffct){
  rcovPiBetak =  matrix(rep(0,(ng-1)*nbeta[k]), ncol = nbeta[k])
  for  (kp in 1:(ng-1)){
    for (l in 1:nbeta[k]){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + BiklNL(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + BiklNL(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      rcovPiBetak[kp,l] =  tmp
    }
  }
  return(rcovPiBetak)
}
# cov betak sigma
covBetaSigmakNL <- function(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct){
  rcovBetaSigmak = matrix(rep(0,nbeta[k]*ng), nrow = nbeta[k])
  for (kp in 1:nbeta[k]){
    for (l in 1:ng){
      tmp = 0
      for (i in 1:n){
        if (k==l){
          tmp = tmp + BiklNL(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*SikNL(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)*taux[i,k]*(1-taux[i,k])
        }
        else{
          tmp = tmp - BiklNL(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*SikNL(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)*taux[i,k]*taux[i,l]
        }
      }
      rcovBetaSigmak[kp,l] = tmp
    }
  }
  return(rcovBetaSigmak)
}
# matrix cov pi deltak
covPiDeltakNL <- function(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, fct, diffct){
  rcovPiDeltak =  matrix(rep(0,(ng-1)*nw), ncol = nw)

  return(rcovPiDeltak)
}
# matrix cov Betak Betal
covBetakBetalNL <- function(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct){
  rcovBetakBetal = matrix(rep(0, (nbeta[k]*nbeta[l])), nrow = nbeta[k])
  if (k==l){
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + BiklNL(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*BiklNL(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*taux[i,k]*(1-taux[i,k])
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp - BiklNL(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*BiklNL(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*taux[i,k]*taux[i,l]
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }
  return(rcovBetakBetal)
}
# matrix cov Betak Deltal
covBetakDeltalNL <- function(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct){
  rcovBetakDeltal = matrix(rep(0, (nbeta[k]*nw)), nrow = nbeta[k])
 
  return(rcovBetakDeltal)
}
# matrix cov Delak Deltal
covDeltakDeltalNL <- function(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct){
  rcovDeltakDeltal = matrix(rep(0, (nw*nw)), nrow = nbeta[k])
  
  return(rcovDeltakDeltal)
}
# cov detak sigma
covDeltaSigmakNL <- function(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct){
  rcovDeltaSigmak = matrix(rep(0,nw*ng), nrow = nw)

  return(rcovDeltaSigmak)
}
##################################################################################
# Main function
##################################################################################
IEMNL <- function(paramEM, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr, fct, diffct){
  nsigma = ng
  period = ncol(A)
  #################################################################################
  # Compute of matrix information
  #################################################################################
  if (nx == 1){
    pi = c(paramEM[1:(ng-1)], 1-sum(paramEM[1:(ng-1)]))
    beta = paramEM[(ng):(ng+sum(nbeta)-1)]
    sigma = paramEM[c((ng+sum(nbeta)):(ng+sum(nbeta)+ng-1))]
    delta = paramEM[-c(1:(ng+sum(nbeta)+ng-1))]
  }else{
    theta = c(paramEM[1:(ng*nx)])
    pi = theta
    # we fit the group ref with 0
    theta = theta - theta[((refgr-1)*nx+1):((refgr-1)*nx+nx)]
    beta = paramEM[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = paramEM[c((ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng))]
    delta = paramEM[-c(1:(ng*nx+sum(nbeta)+ng))]
  }
  taux = ftauxNL(pi, beta, sigma, ng, nbeta, n, nw, nx, A, Y, X, TCOV, fct)
  period = ncol(A)
  nbetacum = cumsum(c(0, nbeta))
  ndelta = rep(nw, ng)
  ndeltacum = cumsum(c(0, ndelta))
  #################################################################################
  # matrix differential of  Pi**2 or theta**2
  #################################################################################
  if (nx == 1){
    mPi = matrix(rep(0,(ng-1)**2), ncol = ng-1)
    for (k in 1:(ng-1)){
      for (l in 1:(ng-1)){
        if (k==l){
          mPi[k,l] = -sum(taux[,k]/pi[k]**2 + taux[,ng]/pi[ng]**2)
        }
        else{
          mPi[k,l] = -sum(taux[,ng]/pi[ng]**2)
        }
      }
    }
  }else{
    mPi = matrix(rep(0,(ng*nx)**2), ncol = ng*nx)
    for (k in 1:ng){
      for (kp in 1:ng){
        tmp1 = c()
        for (l in 1:nx){
          for (lp in 1:nx){
            tmp2 = 0
            if (k==kp){
              for (i in 1:n){
                tmpPiik = piik(theta, i, k, ng, X)
                tmp2 = tmp2 - tmpPiik*(1-tmpPiik)*X[i,l]*X[i,lp]
              }
            }else{
              for (i in 1:n){
                tmp2 = tmp2 + piik(theta, i, k, ng, X)*piik(theta, i, kp, ng, X)*X[i,l]*X[i,lp]
              }
            }
            tmp1 = c(tmp1, tmp2)
          }
        }
        mPi[((k-1)*nx+1):((k-1)*nx+nx), ((kp-1)*nx+1):((kp-1)*nx+nx)] = matrix(tmp1, ncol = nx, byrow =TRUE)
      }
    }
    # we have to delete the part of the reference group
    mPi = mPi[-c(((refgr-1)*nx+1):((refgr-1)*nx+nx)), ]
    mPi = mPi[ ,-c(((refgr-1)*nx+1):((refgr-1)*nx+nx))]
  }
  ##################################################################################
  # matrix -B
  #################################################################################
  B = matrix(rep(0, (sum(nbeta)+ng)**2), ncol =  sum(nbeta)+ng) # we take off the dimension of mPi
  if (nw !=0){
    for (i in 1:n){
      for (t in (1:period)){
        B = B + rbind(cbind(mbetaNL(i,t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw),
                            mbetadeltaNL(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw),
                            mbetasigmaNL(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                      cbind(t(mbetadeltaNL(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                            mdeltaNL(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw),
                            mdeltasigmaNL(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                      cbind(t(mbetasigmaNL(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                            t(mdeltasigmaNL(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                            msigmaNL(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct)))
      }
    }
  }else{
    for (i in 1:n){
      for (t in (1:period)){
        B = B + rbind(cbind(mbetaNL(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct),
                            mbetasigmaNL(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct)),
                      cbind(t(mbetasigmaNL(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct, diffct)),
                            msigmaNL(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, fct)))
      }
    }
  }
  B = cbind(matrix(rep(0, (ng-1)*nx*(sum(nbeta)+ng+nw*ng)), ncol = (ng-1)*nx),
            B)
  B = rbind(cbind(mPi, matrix(rep(0,(sum(nbeta)+ng+nw*ng)*(ng-1)*nx), nrow = (ng-1)*nx)),
            B)
  ##################################################################################
  # matrix cov pi
  ##################################################################################
  if (nx == 1){
    covPi = matrix(rep(0,(ng-1)**2), ncol=ng-1)
    for (k in 1:(ng-1)){
      for (l in 1:(ng-1)){
        tmp = 0
        if (k==l){
          for (i in 1:n){
            tmp = tmp + taux[i,k]*(1-taux[i,k])/pi[k]**2+taux[i,ng]*(1-taux[i,ng])/pi[ng]**2-2*taux[i,k]*taux[i,ng]/(pi[k]*pi[ng])
          }
        }else{
          for (i in 1:n){
            tmp = tmp - taux[i,k]*taux[i,l]/(pi[k]*pi[l])+taux[i,k]*taux[i,ng]/(pi[ng]*pi[k])+taux[i,l]*taux[i,ng]/(pi[l]*pi[ng])+taux[i,ng]*(1-taux[i,ng])/pi[ng]**2
          }
        }
        covPi[k,l] = tmp
      }
    }
  }else{
    # enlever a tuax la bonne colonne
    covPi = matrix(rep(0,((ng-1)*nx)**2), ncol=(ng-1)*nx)
    for (k in 1:(ng-1)){
      for (l in 1:(ng-1)){
        covPitmp = matrix(rep(0, nx**2), ncol =nx)
        for (p in 1:nx){
          for (q in 1:nx){
            tmp = 0
            if (k==l){
              for (i in 1:n){
                tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*X[i,q]
              }
            }else{
              for (i in 1:n){
                tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*X[i,q]
              }
            }
            covPitmp[p, q] = tmp
          }
        }
        covPi[((k-1)*nx+1):((k-1)*nx+nx), ((l-1)*nx+1):((l-1)*nx+nx)] = covPitmp
      }
    }
  }
  ##################################################################################
  # matrix cov pi beta
  ##################################################################################
  if (nx == 1){
    # matrix cov pi beta
    covPiBeta = c()
    for (k in 1:(ng-1)){
      covPiBeta = cbind(covPiBeta, covPiBetakNL(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, fct, diffct))
    }
    rcovPiBetak =  matrix(rep(0,(ng-1)*nbeta[ng]), ncol=nbeta[ng])
    for  (kp in 1:(ng-1)){
      for (l in 1:nbeta[ng]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + BiklNL(i, ng, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
        }
      }
      rcovPiBetak[kp,l] =  tmp
    }
    covPiBeta = cbind(covPiBeta, rcovPiBetak)
  }else{
    covPiBeta = matrix(rep(0,((ng-1)*nx)*sum(nbeta)), ncol=sum(nbeta))
    for (k in 1:(ng-1)){
      for (l in 1:ng){
        covPiBetatmp = matrix(rep(0, nx*nbeta[l]), ncol =nbeta[l])
        for (p in 1:nx){
          for (q in 1:nbeta[l]){
            tmp = 0
            if (k==l){
              for (i in 1:n){
                tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*BiklNL(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)
              }
            }else{
              for (i in 1:n){
                tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*BiklNL(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)
              }
            }
            covPiBetatmp[p, q] = tmp
          }
        }
        covPiBeta[((k-1)*nx+1):((k-1)*nx+nx), (nbetacum[l]+1):(nbetacum[l+1])] = covPiBetatmp
      }
    }
  }
  if (nw !=0){
    ##################################################################################
    # matrix cov pi delta
    ##################################################################################
    # matrix cov pi delta
    if (nx == 1 ){
      covPiDelta = c()
      for (k in 1:(ng-1)){
        covPiDelta = cbind(covPiDelta, covPiDeltakNL(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi))
      }
      rcovPiDeltak =  matrix(rep(0,(ng-1)*nw), ncol=nw)
      for  (kp in 1:(ng-1)){
        for (l in 1:nw){
          tmp = 0
          for (i in 1:n){
            tmp = tmp + DiklNL(i, ng, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
          }
        }
        rcovPiDeltak[kp,l] =  tmp
      }
      covPiDelta = cbind(covPiDelta, rcovPiDeltak)
    }else{
      covPiDelta = matrix(rep(0,((ng-1)*nx)*sum(delta)), ncol=sum(ndelta))
      for (k in 1:(ng-1)){
        for (l in 1:ng){
          covPiDeltatmp = matrix(rep(0, nx*ndelta[l]), ncol =ndelta[l])
          for (p in 1:nx){
            for (q in 1:ndelta[l]){
              tmp = 0
              if (k==l){
                for (i in 1:n){
                  tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*DiklNL(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
                }
              }else{
                for (i in 1:n){
                  tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*DiklNL(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
                }
              }
              covPiDeltatmp[p, q] = tmp
            }
          }
          covPiDelta[((k-1)*nx+1):((k-1)*nx+nx), (ndeltacum[l]+1):(ndeltacum[l+1])] = covPiDeltatmp
        }
      }
    }
    ##################################################################################
    # Matrix cov beta delta
    ##################################################################################
    covBetaDelta = matrix(rep(0,(sum(nbeta)*nw*ng)), nrow = sum(nbeta))
    for (k in 1:ng){
      for (l in 1:ng){
        covBetaDelta[(nbetacum[k]+1):(nbetacum[k+1]),(ndeltacum[l]+1):(ndeltacum[l+1])]  = covBetakDeltalNL(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
      }
    }
    ##################################################################################
    # matrix cov delta
    ##################################################################################
    covDelta = matrix(rep(0,sum(nw*ng)**2), ncol = nw*ng)
    for (k in 1:ng){
      for (l in 1:ng){
        covDelta[(ndeltacum[k]+1):(ndeltacum[k+1]),(ndeltacum[l]+1):(ndeltacum[l+1])] = covDeltakDeltalNL(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
      }
    }
    ##################################################################################
    # Matrix cov delta sigma
    ##################################################################################
    covDeltaSigma = c()
    for (k in 1:ng){
      covDeltaSigma = rbind(covDeltaSigma, covDeltaSigmakNL(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw))
    }
    #### enf if nw==0
  }
  ##################################################################################
  # matrix cov pi sigma
  ##################################################################################
  if (nx == 1){
    covPiSigma = matrix(rep(0,(ng-1)*ng), nrow=ng-1)
    for (k in 1:(ng-1)){
      for (l in 1:ng){
        if ( l!=ng){
          tmp = 0
          if (k==l){
            for (i in 1:n){
              tmp = tmp + taux[i,k]*SikNL(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)*((1-taux[i,k])/pi[k]+taux[i,ng]/pi[ng])
            }
          }
          else{
            for (i in 1:n){
              tmp = tmp + taux[i,l]*SikNL(i,l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)/pi[l]*(-taux[i,k]/pi[k]+taux[i,ng]/pi[ng])
            }
          }
        }else{
          tmp = 0
          for (i in 1:n){
            tmp = tmp + taux[i,ng]*SikNL(i, ng, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)*((1-taux[i,ng])/pi[ng]+taux[i,k]/pi[k])
          }
        }
        covPiSigma[k,l] = tmp
      }
    }
  }else{
    covPiSigma = matrix(rep(0,((ng-1)*nx)*ng), ncol=ng)
    for (k in 1:(ng-1)){
      for (l in 1:ng){
        covPiSigmatmp = matrix(rep(0, nx*ng), ncol = ng)
        for (p in 1:nx){
          tmp = 0
          if (k==l){
            for (i in 1:n){
              tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*SikNL(i,k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)
            }
          }else{
            for (i in 1:n){
              tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*SikNL(i,l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)
            }
          }
          covPiSigmatmp[p, l] = tmp
        }
        covPiSigma[((k-1)*nx+1):((k-1)*nx+nx), 1:ng] = covPiSigmatmp
      }
    }
  }
  ##################################################################################
  # matrix cov beta
  ##################################################################################
  covBeta = matrix(rep(0,sum(nbeta)**2), ncol = sum(nbeta))
  for (k in 1:ng){
    for (l in 1:ng){
      covBeta[(nbetacum[k]+1):(nbetacum[k+1]),(nbetacum[l]+1):(nbetacum[l+1])] = covBetakBetalNL(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct)
    }
  }
  ##################################################################################
  # Matrix cov beta sigma
  ##################################################################################
  covBetaSigma = c()
  for (k in 1:ng){
    covBetaSigma = rbind(covBetaSigma, covBetaSigmakNL(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct, diffct))
  }
  ##################################################################################
  # Matrix cov sigma
  ##################################################################################
  covSigma = matrix(rep(0, ng*ng), ncol =ng)
  for (k in 1:ng){
    for (l in 1:ng){
      tmp = 0
      if (k==l){
        for (i in 1:n){
          tmp = tmp + SikNL(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)**2*taux[i,k]*(1-taux[i,k])
        }
      }else{
        for (i in 1:n){
          tmp = tmp - SikNL(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)*SikNL(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, fct)*taux[i,k]*taux[i,l]
        }
      }
      covSigma[k,l] = tmp
    }
  }
  ##################################################################################
  # Matrix of covariance of score function
  ##################################################################################
  if (nw !=0){
    cov =c()
    cov = cbind(covPi, covPiBeta, covPiDelta, covPiSigma)
    cov = rbind(cov, cbind(t(covPiBeta), covBeta, covBetaDelta, covBetaSigma))
    cov = rbind(cov, cbind(t(covPiDelta), t(covBetaDelta), covDelta, covDeltaSigma))
    cov = rbind(cov, cbind(t(covPiSigma), t(covBetaSigma), t(covDeltaSigma), covSigma))
  }else{
    cov =c()
    cov = cbind(covPi, covPiBeta, covPiSigma)
    cov = rbind(cov, cbind(t(covPiBeta), covBeta, covBetaSigma))
    cov = rbind(cov, cbind(t(covPiSigma), t(covBetaSigma), covSigma))
  }
  ##################################################################################
  # Information matrix of Fisher
  ##################################################################################
  IEM = - B - cov
  return(sqrt(diag(solve(IEM))))
}
