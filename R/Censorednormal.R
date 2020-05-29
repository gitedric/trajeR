#################################################################################
#
# Likelihood
#
#################################################################################
#################################################################################
# Function g
#################################################################################
#' Title
#'
#' @inheritParams trajeR.CNORM
#'
#' @return a real. The porduct of probability for each group.
#' @export
#'
gkCNORM <- function(beta, sigma, i, k, nbeta, A, Y, ymin, ymax, TCOV, delta, nw){
  period = ncol(A)
  muikt = sapply(1:period, function(s){
    sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
  })
  tmp1 =c()
  indmin = (Y[i,]<=ymin)
  if (any(indmin)){
    tmp1 =c(tmp1, pnorm((Y[i, which(indmin)]-muikt[which(indmin)])/sigma[k]))
  }
  indmax = (Y[i,]>=ymax)
  if (any(indmax)){
    tmp1 = c(tmp1, pnorm(-(Y[i, which(indmax)]-muikt[which(indmax)])/sigma[k]))
  }
  ind = !(indmin | indmax)
  if (any(ind)){
    tmp1 = c(tmp1, dnorm((Y[i, which(ind)]-muikt[which(ind)])/sigma[k])/sigma[k])
  }
  return(prod(tmp1))
}
#################################################################################
# gk with exponential parametrization
#################################################################################
gkalpha <- function(beta, alpha, i, k, nbeta, A, Y, ymin, ymax, TCOV, delta, nw){
  period = ncol(A)
  sigma = exp(alpha)
  muikt = sapply(1:period, function(s){
    sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
  })
  tmp1 =c()
  indmin = (Y[i,]<=ymin)
  if (any(indmin)){
    tmp1 =c(tmp1, pnorm((Y[i, which(indmin)]-muikt[which(indmin)])/sigma[k]))
  }
  indmax = (Y[i,]>=ymax)
  if (any(indmax)){
    tmp1 = c(tmp1, pnorm(-(Y[i, which(indmax)]-muikt[which(indmax)])/sigma[k]))
  }
  ind = !(indmin | indmax)
  if (any(ind)){
    tmp1 = c(tmp1, dnorm((Y[i, which(ind)]-muikt[which(ind)])/sigma[k])/sigma[k])
  }
  return(prod(tmp1))
}
#################################################################################
#
# Differential of function by all the parameters
#
#################################################################################

#################################################################################
# dif likelihood theta
#################################################################################
difLthetakalpha <- function(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  alpha = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
  deltatmp = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  thetas = c()
  for (l in 1:nx){
    a = 0
    for (i in 1:n){
      tmp = exp(sapply(1:ng,function(s){theta[((s-1)*nx+1):(s*nx)]%*%X[i,]}))
      tmp1 = tmp[k]*sum(sapply(1:ng, function(s){tmp[s]*(gkalpha(beta, alpha, i, k, nbeta, A, Y, ymin, ymax, TCOV, delta, nw) - gkalpha(beta, alpha, i, s,  nbeta, A, Y, ymin, ymax, TCOV, delta, nw))}))
      tmp2 = sum(sapply(1:ng, function(s){tmp[s]*gkalpha(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)}))
      a = a + X[i,l]*tmp1/(sum(tmp)*tmp2)
    }
    thetas = c(thetas, a)
  }
  return(thetas)
}
#################################################################################
# dif likelihood beta with reparametrization alpha
#################################################################################
difLbetakalpha <- function(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  alpha = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
  deltatmp = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  betas = c()
  m = rep(0,period)
  for (l in 1:nbeta[k]){
    a = 0
    for (i in 1:n){
      difbklmin = 0
      difbklmax = 0
      difbkl = 0
      pymin = 1
      pymax = 1
      py = 1
      indmin = (Y[i,]<=ymin)
      indmax = (Y[i,]>=ymax)
      ind = !(indmin | indmax)
      muikt = sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
      })
      if (any(indmin)){
        m[which(indmin)] = (ymin-muikt[which(indmin)])*exp(-alpha[k])
        mtmp = (ymin-muikt[which(indmin)])*exp(-alpha[k])
        if (sum(indmin)>1){
          tind = 1
          for (t in which(indmin)){
            difbklmin = difbklmin - A[i,t]**(l-1)*exp(-alpha[k])*dnorm(m[t])*prod(pnorm(mtmp[-tind]))
            tind = tind + 1
          }
        } else{
          difbklmin = - A[i,which(indmin)]**(l-1)*exp(-alpha[k])*dnorm(m[which(indmin)])
        }
        pymin = prod(pnorm(m[which(indmin)]))
      }
      if (any(indmax)){
        m[which(indmax)] = (ymax-muikt[which(indmax)])*exp(-alpha[k])
        mtmp =  (ymax-muikt[which(indmax)])*exp(-alpha[k])
        if (sum(indmax)>1){
          tind = 1
          for (t in which(indmax)){
            difbklmax = difbklmax + A[i,t]**(l-1)*exp(-alpha[k])*dnorm(m[t])*prod(pnorm(-mtmp[-tind]))
            tind = tind + 1
          }
        } else{
          difbklmax = A[i,which(indmax)]**(l-1)*exp(-alpha[k])*dnorm(m[which(indmax)])
        }
        pymax = prod(pnorm(-m[which(indmax)]))
      }
      if (any(ind)){
        m[which(ind)] = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
        mtmp = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
        if (sum(ind)>1){
          tind = 1
          for (t in which(ind)){
            difbkl = difbkl + A[i,t]**(l-1)*exp(-2*alpha[k])*m[t]*dnorm(m[t])*prod(exp(-alpha[k])*dnorm(mtmp[-tind]))
            tind = tind +1
          }
        } else{
          difbkl = A[i,which(ind)]**(l-1)*exp(-3*alpha[k])*(Y[i,which(ind)]-muikt[which(ind)])*dnorm(m[which(ind)])
        }
        py = prod(exp(-alpha[k])*dnorm(m[which(ind)]))
      }
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkalpha(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)}))
      a = a + piik(theta, i, k, ng, X)/s*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax)
    }
    betas = c(betas, a)
  }
  return(betas)
}
#################################################################################
# dif likelihood delta with reparametrization alpha
#################################################################################
difLdeltakalpha <- function(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  alpha = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
  deltatmp = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  deltas = c()
  m = rep(0,period)
  for (l in 1:nw){
    a = 0
    for (i in 1:n){
      difbklmin = 0
      difbklmax = 0
      difbkl = 0
      pymin = 1
      pymax = 1
      py = 1
      indmin = (Y[i,]<=ymin)
      indmax = (Y[i,]>=ymax)
      ind = !(indmin | indmax)
      muikt = sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
      })
      if (any(indmin)){
        m[which(indmin)] = (ymin-muikt[which(indmin)])*exp(-alpha[k])
        mtmp = (ymin-muikt[which(indmin)])*exp(-alpha[k])
        if (sum(indmin)>1){
          tind = 1
          for (t in which(indmin)){
            difbklmin = difbklmin - TCOV[i, t + (l-1)*period]*exp(-alpha[k])*dnorm(m[t])*prod(pnorm(mtmp[-tind]))
            tind = tind + 1
          }
        } else{
          difbklmin = - TCOV[i, which(indmin) + (l-1)*period]*exp(-alpha[k])*dnorm(m[which(indmin)])
        }
        pymin = prod(pnorm(m[which(indmin)]))
      }
      if (any(indmax)){
        m[which(indmax)] = (ymax-muikt[which(indmax)])*exp(-alpha[k])
        mtmp =  (ymax-muikt[which(indmax)])*exp(-alpha[k])
        if (sum(indmax)>1){
          tind = 1
          for (t in which(indmax)){
            difbklmax = difbklmax + TCOV[i, t + (l-1)*period]*exp(-alpha[k])*dnorm(m[t])*prod(pnorm(-mtmp[-tind]))
            tind = tind + 1
          }
        } else{
          difbklmax = TCOV[i, which(indmax) + (l-1)*period]*exp(-alpha[k])*dnorm(m[which(indmax)])
        }
        pymax = prod(pnorm(-m[which(indmax)]))
      }
      if (any(ind)){
        m[which(ind)] = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
        mtmp = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
        if (sum(ind)>1){
          tind = 1
          for (t in which(ind)){
            difbkl = difbkl +  TCOV[i, t + (l-1)*period]*exp(-2*alpha[k])*m[t]*dnorm(m[t])*prod(exp(-alpha[k])*dnorm(mtmp[-tind]))
            tind = tind +1
          }
        } else{
          difbkl = TCOV[i, which(ind) + (l-1)*period]*exp(-3*alpha[k])*(Y[i,which(ind)]-muikt[which(ind)])*dnorm(m[which(ind)])
        }
        py = prod(exp(-alpha[k])*dnorm(m[which(ind)]))
      }
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkalpha(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)}))
      a = a + piik(theta, i, k, ng, X)/s*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax)
    }
    deltas = c(deltas, a)
  }
  return(deltas)
}
#################################################################################
# dif likelihood sigma with reparametrization alpha
#################################################################################
difLsigmakalpha <- function(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  alpha = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
  deltatmp = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  alphas = c()
  a = 0
  m = rep(0,period)
  for (i in 1:n){
    difbklmin = 0
    difbklmax = 0
    difbkl = 0
    pymin = 1
    pymax = 1
    py = 1
    indmin = (Y[i,]<=ymin)
    indmax = (Y[i,]>=ymax)
    ind = !(indmin | indmax)
    muikt = sapply(1:period, function(s){
      sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
    })
    if (any(indmin)){
      m[which(indmin)] = (ymin-muikt[which(indmin)])*exp(-alpha[k])
      mtmp = (ymin-muikt[which(indmin)])*exp(-alpha[k])
      if (sum(indmin)>1){
        tind = 1
        for (t in which(indmin)){
          difbklmin = difbklmin -m[t]*dnorm(m[t])*prod(pnorm(mtmp[-tind]))
          tind = tind + 1
        }
      } else{
        difbklmin = -m[which(indmin)]*dnorm(m[which(indmin)])
      }
      pymin = prod(pnorm(m[which(indmin)]))
    }
    if (any(indmax)){
      m[which(indmax)] = (ymax-muikt[which(indmax)])*exp(-alpha[k])
      mtmp = (ymax-muikt[which(indmax)])*exp(-alpha[k])
      if (sum(indmax)>1){
        tind = 1
        for (t in which(indmax)){
          difbklmax = difbklmax + m[t]*dnorm(m[t])*prod(pnorm(-mtmp[-tind]))
          tind = tind + 1
        }
      } else{
        difbklmax =m[which(indmax)]*dnorm(m[which(indmax)])
      }
      pymax = prod(pnorm(-m[which(indmax)]))
    }
    if (any(ind)){
      m[which(ind)] = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
      mtmp = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
      if (sum(ind)>1){
        tind = 1
        for (t in which(ind)){
          difbkl = difbkl + (-1+m[t]**2)*exp(-alpha[k])*dnorm(m[t])*prod(exp(-alpha[k])*dnorm(mtmp[-tind]))
          tind = tind + 1
        }
      } else{
        difbkl = (-1+m[which(ind)]**2)*exp(-alpha[k])*dnorm(m[which(ind)])
      }
      py = prod(exp(-alpha[k])*dnorm(m[which(ind)]))
    }
    s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkalpha(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)}))
    a = a + piik(theta, i, k, ng, X)/s*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax)
  }
  alphas = c(alphas, a)
  return(alphas)
}
#################################################################################
# dif likelihood sigma with reparametrization alpha same sigma
#################################################################################
difLsigmaalphaunique <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  alpha = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
  deltatmp = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  alphas = c()
  a = 0
  m = rep(0,period)
  for (i in 1:n){
    s1 = 0
    s2 = 0
    for (k in 1:ng){
      difbklmin = 0
      difbklmax = 0
      difbkl = 0
      pymin = 1
      pymax = 1
      py = 1
      indmin = (Y[i,]<=ymin)
      indmax = (Y[i,]>=ymax)
      ind = !(indmin | indmax)
      muikt = sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
      })
      if (any(indmin)){
        m[which(indmin)] = (ymin-muikt[which(indmin)])*exp(-alpha[k])
        mtmp = (ymin-muikt[which(indmin)])*exp(-alpha[k])
        if (sum(indmin)>1){
          tind = 1
          for (t in which(indmin)){
            difbklmin = difbklmin -m[t]*dnorm(m[t])*prod(pnorm(mtmp[-tind]))
            tind = tind + 1
          }
        } else{
          difbklmin = -m[which(indmin)]*dnorm(m[which(indmin)])
        }
        pymin = prod(pnorm(m[which(indmin)]))
      }
      if (any(indmax)){
        m[which(indmax)] = (ymax-muikt[which(indmax)])*exp(-alpha[k])
        mtmp = (ymax-muikt[which(indmax)])*exp(-alpha[k])
        if (sum(indmax)>1){
          tind = 1
          for (t in which(indmax)){
            difbklmax = difbklmax + m[t]*dnorm(m[t])*prod(pnorm(-mtmp[-tind]))
            tind = tind + 1
          }
        } else{
          difbklmax =m[which(indmax)]*dnorm(m[which(indmax)])
        }
        pymax = prod(pnorm(-m[which(indmax)]))
      }
      if (any(ind)){
        m[which(ind)] = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
        mtmp = (Y[i,which(ind)]-muikt[which(ind)])*exp(-alpha[k])
        if (sum(ind)>1){
          tind = 1
          for (t in which(ind)){
            difbkl = difbkl + (-1+m[t]**2)*exp(-alpha[k])*dnorm(m[t])*prod(exp(-alpha[k])*dnorm(mtmp[-tind]))
            tind = tind + 1
          }
        } else{
          difbkl = (-1+m[which(ind)]**2)*exp(-alpha[k])*dnorm(m[which(ind)])
        }
        py = prod(exp(-alpha[k])*dnorm(m[which(ind)]))
      }
      s1 = s1 + sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkalpha(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)}))
      s2 = s2 + piik(theta, i, k, ng, X)*(difbklmin*py*pymax+difbklmax*py*pymin+difbkl*pymin*pymax)
    }
    a = a + s2/s1
  }
  alphas = rep(a, ng)
  return(alphas)
}
#################################################################################
# dif likelihood sigma with exponential reparametrization
#################################################################################
difLalpha <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  out =c()
  for ( k in 1:ng){
    out = c(out, difLthetakalpha(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
  }
  for ( k in 1:ng){
    out = c(out, difLbetakalpha(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
  }
  for ( k in 1:ng){
    out = c(out, difLsigmakalpha(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
  }
  if (nw != 0){
    for (k in 1:ng){
      out =c(out, difLdeltakalpha(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
    }
  }
  return(out)
}
#################################################################################
# dif likelihood sigma with exponential reparametrization same sigma
#################################################################################
difLalphaunique <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  out =c()
  for (k in 1:ng){
    out = c(out, difLthetakalpha(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
  }
  for (k in 1:ng){
    out = c(out, difLbetakalpha(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
  }
  out = c(out, difLsigmaalphaunique(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
  if (nw != 0){
    for (k in 1:ng){
      out =c(out, difLdeltakalpha(param, k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw))
    }
  }
  return(out)
}
#################################################################################
# likelihood sigma with exponential reparametrization
#################################################################################
Likelihoodalpha <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  alpha = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
  deltatmp = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  a = 0
  for (i in 1:n){
    a = a + log(sum(sapply(1:ng, function(s){
      piik(theta, i, s, ng, X)*gkalpha(beta, alpha, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)
    }))
    )
  }
  return(a)
}
#################################################################################
# likelihood sigma
#################################################################################
LikelihoodCNORM <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw){
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  sigma = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
  deltatmp = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  a = 0
  for (i in 1:n){
    a = a + log(sum(sapply(1:ng, function(s){
      piik(theta, i, s, ng, X)*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)
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
# Likelihood calculation for EM algorithm
#################################################################################
likelihoodEM <- function(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw){
  likeli = 0
  betatmp = beta
  sigma = sigma
  deltatmp = delta
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  for (i in 1:n){
    likeli = likeli + log(sum(sapply(1:ng, function(s){
      pi[s]*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)
    }))
    )
  }
  return(likeli)
}
#################################################################################
# Function rate
#################################################################################
ftaux <- function(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X){
  betatmp = beta
  sigma = sigma
  deltatmp = delta
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  if (nx == 1){
  tmp2 = c()
    for (k in 1:ng){
      tmp1 = sapply(1:n, function(s){gkCNORM(beta, sigma, s, k, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
      tmp2 = cbind(tmp2, tmp1*pi[k])
    }
  tmp3 =c()
  for (k in 1:ng){
    tmp3 = cbind(tmp3, 1/(1+ (rowSums(tmp2)-tmp2[,k])/tmp2[,k]))
  }
  }else{
    tmp3 = c()
    for (i in 1:n){
      tmp = sapply(1:ng, function(s){
        piik(pi, i, s, ng, X)*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
      tmp3 = rbind(tmp3, tmp/sum(tmp))
    }
  }
  return(tmp3)
}
#################################################################################
# EM not censored
#################################################################################
EM <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, delta, nw, itermax, EMIRLS){
  period = ncol(A)
  nsigma = ng
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    sigma = param[(ng+sum(nbeta)):(ng+sum(nbeta)+ng-1)]
    delta = param[-c(1:(ng+sum(nbeta)+ng-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
    delta = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  }
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEM(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw), "\n"))
    }else{
      cat(paste(-LikelihoodCNORM(c(pi, beta, sigma, delta), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw), "\n"))
    }
    # E-step
    taux = ftaux(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)
    newbeta = c()
    newsigma = c()
    newdelta = c()
    # M-step
    if (nw == 0){
      for (k in 1:ng){
        a = 0
        b = 0
        for (i in 1:n){
          Ai = c()
          for (t in 1:period){
            Ai = cbind(Ai,sapply(0:(nbeta[k]-1), function(s){A[i,t]**(s)}))
          }
          a = a + taux[i,k]*Y[i,]%*%t(Ai)
          Mtmp = Y[i,]-beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai
          b = b + taux[i,k]*(Mtmp%*%t(Mtmp))
        }
        newbeta = c(newbeta, a%*%solve(Ai%*%t(Ai))/sum(taux[,k]))
        newsigma = c(newsigma, sqrt(b/(period*sum(taux[,k]))))
      }
    }else{
      for (k in 1:ng){
        a = 0
        c = 0
        Sw = 0
        b = 0
        for (i in 1:n){
          Ai = c()
          Wi =c()
          for (t in 1:period){
            Ai = cbind(Ai,sapply(0:(nbeta[k]-1), function(s){A[i,t]**(s)}))
            Wi = cbind(Wi, TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
          }
          a = a + taux[i,k]*(Y[i,]%*%t(Ai)-delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi%*%t(Ai))
          c = c + taux[i,k]*(Y[i,]%*%t(Wi)-beta[(nbetacum[k]+1):(nbetacum[k+1])]%*%Ai%*%t(Wi))
          Sw = Sw + taux[i,k]*Wi%*%t(Wi)
          Mtmp = Y[i,]-beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai-delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi
          b = b + taux[i,k]*(Mtmp%*%t(Mtmp))
        }
        newbeta = c(newbeta, a%*%solve(Ai%*%t(Ai))/sum(taux[,k]))
        newdelta =c(newdelta, c%*%solve(Sw))
        newsigma = c(newsigma, sqrt(b/(period*sum(taux[,k]))))
      }
    }
    if (nx == 1){
      pi = colSums(taux)/n
      stoppi = 0
      refpi = 0
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS)
      stoppi = (pi-pi[1:ng])-(newpi-newpi[1:ng])
      refpi = rep(0, length(stoppi))
      pi = newpi
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newsigma, newdelta, stoppi)-c(beta, sigma, delta, refpi))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    sigma = newsigma
    delta = newdelta
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, sigma, delta)
  }else{
    param = c(pi, beta, sigma, delta)
  }
  return(param)
}
#################################################################################
# EM not censored sigma unique
#################################################################################
EMSigmaunique <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, delta, nw, itermax, EMIRLS){
  period = ncol(A)
  nsigma = ng
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    sigma = param[(ng+sum(nbeta)):(ng+sum(nbeta)+ng-1)]
    delta = param[-c(1:(ng+sum(nbeta)+ng-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
    delta = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  }
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour < itermax) {
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEM(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw), "\n"))
    }else{
      cat(paste(-LikelihoodCNORM(c(pi, beta, sigma, delta), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw), "\n"))
    }
    # E-step
    taux = ftaux(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)
    newbeta = c()
    newdelta = c()
    # M-step
    if (nw == 0) {
      for (k in 1:ng) {
        a = 0
        for (i in 1:n) {
          Ai = c()
          for (t in 1:period) {
            Ai = cbind(Ai, sapply(0:(nbeta[k] - 1), function(s) {
              A[i, t] ** (s)
            }))
          }
          a = a + taux[i, k] * Y[i,] %*% t(Ai)
        }
        newbeta = c(newbeta, a %*% solve(Ai %*% t(Ai)) / sum(taux[, k]))
      }
      b = 0
      for (k in 1:ng) {
        for (i in 1:n) {
          Ai = c()
          for (t in 1:period) {
            Ai = cbind(Ai, sapply(0:(nbeta[k] - 1), function(s) {
              A[i, t] ** (s)
            }))
          }
          Mtmp = Y[i,] - beta[(nbetacum[k] + 1):nbetacum[k + 1]] %*% Ai
          b = b + taux[i, k] * (Mtmp %*% t(Mtmp))
        }
      }
      newsigma = sqrt(b / (period * sum(taux)))
      newsigma = rep(newsigma, ng)
    }else{
      for (k in 1:ng) {
        a = 0
        c = 0
        Sw = 0
        for (i in 1:n) {
          Ai = c()
          Wi = c()
          for (t in 1:period) {
            Ai = cbind(Ai, sapply(0:(nbeta[k] - 1), function(s) {A[i, t] ** (s)}))
            Wi = cbind(Wi, TCOV[i, seq(from = t, to = t + (nw - 1) * period, by = period)])
          }
          a = a + taux[i, k] * (Y[i, ] %*% t(Ai) - delta[(ndeltacum[k] + 1):ndeltacum[k + 1]] %*% Wi %*% t(Ai))
          c = c + taux[i, k] * (Y[i, ] %*% t(Wi) - beta[(nbetacum[k] + 1):(nbetacum[k + 1])] %*% Ai %*% t(Wi))
          Sw = Sw + taux[i, k] * Wi %*% t(Wi)
        }
        newbeta = c(newbeta, a %*% solve(Ai %*% t(Ai)) / sum(taux[, k]))
        newdelta = c(newdelta, c %*% solve(Sw))
      }
      b = 0
      for (k in 1:ng) {
        for (i in 1:n) {
          Ai = c()
          Wi = c()
          for (t in 1:period) {
            Ai = cbind(Ai, sapply(0:(nbeta[k] - 1), function(s) {A[i, t]**(s)}))
            Wi = cbind(Wi, TCOV[i, seq(from = t, to = t + (nw - 1) * period, by = period)])
          }
          Mtmp = Y[i,] - newbeta[(nbetacum[k] + 1):nbetacum[k + 1]] %*% Ai - delta[(ndeltacum[k] + 1):ndeltacum[k + 1]] %*% Wi
          b = b + taux[i, k] * (Mtmp %*% t(Mtmp))
        }
      }
      newsigma = sqrt(b / (period * sum(taux)))
      newsigma = rep(newsigma, ng)
    }
    if (nx == 1){
      pi = colSums(taux)/n
      stoppi = 0
      refpi = 0
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx,period, EMIRLS,refgr = 1)
      stoppi = (pi-pi[1:ng])-(newpi-newpi[1:ng])
      refpi = rep(0, length(stoppi))
      pi = newpi
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newsigma, newdelta, stoppi)-c(beta, sigma, delta, refpi))<10**(-6))){
      tour = itermax + 2
    }
    sigma = newsigma
    beta = newbeta
    delta = newdelta
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, sigma, delta)
  }else{
    param = c(pi, beta, sigma, delta)
  }
  return(param)
}
#################################################################################
# EM censored
#################################################################################
EMcensored <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, delta, nw, itermax, EMIRLS){
  period = ncol(A)
  nsigma = ng
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    sigma = param[(ng+sum(nbeta)):(ng+sum(nbeta)+ng-1)]
    delta = param[-c(1:(ng+sum(nbeta)+ng-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
    delta = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  }
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour < itermax){
    if (nx == 1){
      cat(paste(-likelihoodEM(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw), "\n"))
    }else{
      cat(paste(-LikelihoodCNORM(c(pi, beta, sigma, delta), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw), "\n"))
    }
    # E-step
    taux = ftaux(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)
    newbeta = c()
    newsigma = c()
    newdelta = c()
    # M-step
    if (nw == 0){
      for (k in 1:ng){
        a = 0
        b = 0
        for (i in 1:n){
          Ytildei = Y[i,]
          Ytilde2i = Y[i,]**2
          muikt = sapply(1:period, function(s){
            sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,s]**(0:(nbeta[k]-1)))
          })
          alphamax = (ymax-muikt)/sigma[k]
          alphamin = (ymin-muikt)/sigma[k]
          #to avoid numerically division by zero
          alphamax[alphamax>37] = 37
          qmax = dnorm(alphamax)/pnorm(-alphamax)
          qmin = dnorm(alphamin)/pnorm(alphamin)
          Ytildei[which(Y[i,]<=ymin)] = muikt[which(Y[i,]<=ymin)]-sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytildei[which(Y[i,]>=ymax)] = muikt[which(Y[i,]>=ymax)]+sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ytilde2i[which(Y[i,]<=ymin)] = sigma[k]**2*(1-alphamin[which(Y[i,]<=ymin)]*qmin[which(Y[i,]<=ymin)])+muikt[which(Y[i,]<=ymin)]**2-2*muikt[which(Y[i,]<=ymin)]*sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytilde2i[which(Y[i,]>=ymax)] = sigma[k]**2*(1+alphamax[which(Y[i,]>=ymax)]*qmax[which(Y[i,]>=ymax)])+muikt[which(Y[i,]>=ymax)]**2+2*muikt[which(Y[i,]>=ymax)]*sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ai = c()
          for (t in 1:period){
            Ai = cbind(Ai,sapply(0:(nbeta[k]-1), function(s){A[i,t]**(s)}))
          }
          a = a + taux[i,k]*Ytildei%*%t(Ai)
          b = b + taux[i,k]*(t(Ytilde2i)%*%matrix(rep(1, period), ncol=1)-2*Ytildei%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai)+(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai)%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai))
        }
        newbeta = c(newbeta, a%*%solve(Ai%*%t(Ai))/sum(taux[,k]))
        newsigma = c(newsigma, sqrt(b/(period*sum(taux[,k]))))
      }
    }else{
      for (k in 1:ng){
        a = 0
        b = 0
        c = 0
        Sw = 0
        for (i in 1:n){
          Ytildei = Y[i,]
          Ytilde2i = Y[i,]**2
          muikt = sapply(1:period, function(s){
            sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,s]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, s, k, ndeltacum)
          })
          alphamax = (ymax-muikt)/sigma[k]
          alphamin = (ymin-muikt)/sigma[k]
          #to avoid numerically division by zero
          alphamax[alphamax>37] = 37
          qmax = dnorm(alphamax)/pnorm(-alphamax)
          qmin = dnorm(alphamin)/pnorm(alphamin)
          Ytildei[which(Y[i,]<=ymin)] = muikt[which(Y[i,]<=ymin)]-sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytildei[which(Y[i,]>=ymax)] = muikt[which(Y[i,]>=ymax)]+sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ytilde2i[which(Y[i,]<=ymin)] = sigma[k]**2*(1-alphamin[which(Y[i,]<=ymin)]*qmin[which(Y[i,]<=ymin)])+muikt[which(Y[i,]<=ymin)]**2-2*muikt[which(Y[i,]<=ymin)]*sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytilde2i[which(Y[i,]>=ymax)] = sigma[k]**2*(1+alphamax[which(Y[i,]>=ymax)]*qmax[which(Y[i,]>=ymax)])+muikt[which(Y[i,]>=ymax)]**2+2*muikt[which(Y[i,]>=ymax)]*sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ai = c()
          Wi = c()
          for (t in 1:period){
            Ai = cbind(Ai,sapply(0:(nbeta[k]-1), function(s){A[i,t]**(s)}))
            Wi = cbind(Wi, TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
          }
          a = a + taux[i,k]*Ytildei%*%t(Ai)
          c = c + taux[i,k]*(Ytildei%*%t(Wi)-beta[(nbetacum[k]+1):(nbetacum[k+1])]%*%Ai%*%t(Wi))
          Sw = Sw + taux[i,k]*Wi%*%t(Wi)
          b = b + taux[i,k]*(t(Ytilde2i)%*%matrix(rep(1, period), ncol=1)
                             -2*Ytildei%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai+delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi)
                             +(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai+delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi)%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai+delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi))
        }
        newbeta = c(newbeta, a%*%solve(Ai%*%t(Ai))/sum(taux[,k]))
        newdelta =c(newdelta, c/Sw)
        newsigma = c(newsigma, sqrt(b/(period*sum(taux[,k]))))
      }
    }
    if (nx == 1){
      pi = colSums(taux)/n
      stoppi = 0
      refpi = 0
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS)
      stoppi = (pi-pi[1:ng])-(newpi-newpi[1:ng])
      refpi = rep(0, length(stoppi))
      pi = newpi
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newsigma, newdelta, stoppi)-c(beta, sigma, delta, refpi))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    sigma = newsigma
    delta = newdelta
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, sigma, delta)
  }else{
    param = c(pi, beta, sigma, delta)
  }
  return(param)
}
#################################################################################
# EM censored same sigma
#################################################################################
EMcensoredsamesigma <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, delta, nw, itermax, EMIRLS){
  period = ncol(A)
  nsigma = ng
  pi = param[1:(ng - 1)]
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    sigma = param[(ng+sum(nbeta)):(ng+sum(nbeta)+ng-1)]
    delta = param[-c(1:(ng+sum(nbeta)+ng-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng)]
    delta = param[-c(1:(ng*nx+sum(nbeta)+ng))]
  }
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour < itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEM(n, ng, nbeta, beta, sigma, pi, A, Y, ymin, ymax, TCOV, delta, nw), "\n"))
    }else{
      cat(paste(-LikelihoodCNORM(c(pi, beta, sigma, delta), ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw), "\n"))
    }
    # E-step
    taux = ftaux(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)
    newbeta = c()
    newsigma = c()
    newdelta = c()
    # Mstep
    if (nw == 0){
      b = 0
      for (k in 1:ng){
        a = 0
        for (i in 1:n){
          Ytildei = Y[i,]
          Ytilde2i = Y[i,]**2
          muikt = sapply(1:period, function(s){
            sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,s]**(0:(nbeta[k]-1)))
          })
          alphamax = (ymax-muikt)/sigma[k]
          alphamin = (ymin-muikt)/sigma[k]
          #to avoid numerically division by zero
          alphamax[alphamax>37] = 37
          qmax = dnorm(alphamax)/pnorm(-alphamax)
          qmin = dnorm(alphamin)/pnorm(alphamin)
          Ytildei[which(Y[i,]<=ymin)] = muikt[which(Y[i,]<=ymin)]-sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytildei[which(Y[i,]>=ymax)] = muikt[which(Y[i,]>=ymax)]+sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ytilde2i[which(Y[i,]<=ymin)] = sigma[k]**2*(1-alphamin[which(Y[i,]<=ymin)]*qmin[which(Y[i,]<=ymin)])+muikt[which(Y[i,]<=ymin)]**2-2*muikt[which(Y[i,]<=ymin)]*sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytilde2i[which(Y[i,]>=ymax)] = sigma[k]**2*(1+alphamax[which(Y[i,]>=ymax)]*qmax[which(Y[i,]>=ymax)])+muikt[which(Y[i,]>=ymax)]**2+2*muikt[which(Y[i,]>=ymax)]*sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ai = c()
          for (t in 1:period){
            Ai = cbind(Ai,sapply(0:(nbeta[k]-1), function(s){A[i,t]**(s)}))
          }
          a = a + taux[i,k]*Ytildei%*%t(Ai)
          b = b + taux[i,k]*(t(Ytilde2i)%*%matrix(rep(1, period), ncol=1)-2*Ytildei%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai)+(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai)%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai))
        }
        newbeta = c(newbeta, a%*%solve(Ai%*%t(Ai))/sum(taux[,k]))
      }
    }else{
      b = 0
      for (k in 1:ng){
        a = 0
        c = 0
        Sw = 0
        for (i in 1:n){
          Ytildei = Y[i,]
          Ytilde2i = Y[i,]**2
          muikt = sapply(1:period, function(s){
            sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,s]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, s, k, ndeltacum)
          })
          alphamax = (ymax-muikt)/sigma[k]
          alphamin = (ymin-muikt)/sigma[k]
          #to avoid numerically division by zero
          alphamax[alphamax>37] = 37
          qmax = dnorm(alphamax)/pnorm(-alphamax)
          qmin = dnorm(alphamin)/pnorm(alphamin)
          Ytildei[which(Y[i,]<=ymin)] = muikt[which(Y[i,]<=ymin)]-sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytildei[which(Y[i,]>=ymax)] = muikt[which(Y[i,]>=ymax)]+sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ytilde2i[which(Y[i,]<=ymin)] = sigma[k]**2*(1-alphamin[which(Y[i,]<=ymin)]*qmin[which(Y[i,]<=ymin)])+muikt[which(Y[i,]<=ymin)]**2-2*muikt[which(Y[i,]<=ymin)]*sigma[k]*qmin[which(Y[i,]<=ymin)]
          Ytilde2i[which(Y[i,]>=ymax)] = sigma[k]**2*(1+alphamax[which(Y[i,]>=ymax)]*qmax[which(Y[i,]>=ymax)])+muikt[which(Y[i,]>=ymax)]**2+2*muikt[which(Y[i,]>=ymax)]*sigma[k]*qmax[which(Y[i,]>=ymax)]
          Ai = c()
          Wi = c()
          for (t in 1:period){
            Ai = cbind(Ai,sapply(0:(nbeta[k]-1), function(s){A[i,t]**(s)}))
            Wi = cbind(Wi, TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
          }
          a = a + taux[i,k]*(Ytildei%*%t(Ai)-delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi%*%t(Ai))
          c = c + taux[i,k]*(Ytildei%*%t(Wi)-beta[(nbetacum[k]+1):(nbetacum[k+1])]%*%Ai%*%t(Wi))
          Sw = Sw + taux[i,k]*Wi%*%t(Wi)
          b = b + taux[i,k]*(t(Ytilde2i)%*%matrix(rep(1, period), ncol=1)
                             -2*Ytildei%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai+delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi)
                             +(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai+delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi)%*%t(beta[(nbetacum[k]+1):nbetacum[k+1]]%*%Ai+delta[(ndeltacum[k]+1):ndeltacum[k+1]]%*%Wi))
        }
        newbeta = c(newbeta, a%*%solve(Ai%*%t(Ai))/sum(taux[,k]))
        newdelta =c(newdelta, c/Sw)
      }
    }
    newsigma = rep(sqrt(b/(period*n)),nsigma)
    if (nx == 1){
      pi = colSums(taux)/n
      stoppi = 0
      refpi = 0
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS)
      stoppi = (pi-pi[1:ng])-(newpi-newpi[1:ng])
      refpi = rep(0, length(stoppi))
      pi = newpi
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newsigma, newdelta, stoppi)-c(beta, sigma, delta, refpi))<10**(-6))){
      tour = itermax + 2
    }
    sigma = newsigma
    beta = newbeta
    delta = newdelta
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, sigma, delta)
  }else{
    param = c(pi, beta, sigma, delta)
  }
  return(param)
}

#################################################################################
# Compute of matrix information
#################################################################################
#################################################################################
# Definition of function
#################################################################################
# matrix differential of beta**2 for t and i
#################################################################################
mbeta = function(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw){
  res = matrix(rep(0,(sum(nbeta)**2)), ncol = sum(nbeta))
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nbeta[k]){
      for (lp in 1:nbeta[k]){
        tmp = c(tmp, -taux[i,k]*A[i,t]**(l-1)*A[i,t]**(lp-1)/sigma[k]**2)
      }
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), (nbetacum[k]+1):(nbetacum[k+1])] = matrix(tmp, ncol = nbeta[k], byrow = TRUE)
  }
  return(res)
}
#################################################################################
# matrix differential of delta**2 for t and i
#################################################################################
mdelta = function(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw){
  res = matrix(rep(0,((nw*ng)**2)), ncol = nw*ng)
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nw){
      for (lp in 1:nw){
        tmp = c(tmp, -taux[i,k]*TCOV[i, t + (l-1)*period]*TCOV[i, t + (lp-1)*period]/sigma[k]**2)
      }
    }
    res[(ndeltacum[k]+1):(ndeltacum[k+1]), (ndeltacum[k]+1):(ndeltacum[k+1])] = matrix(tmp, ncol = nw, byrow = TRUE)
  }
  return(res)
}
##################################################################################
# matrix differential of  beta delta
#################################################################################
mbetadelta = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw){
  res = matrix(rep(0,ng*nw*sum(nbeta)), ncol = nw*ng)
  for (k in 1:ng){
    tmp = c()
    for (lp in 1:nbeta[k]){
      for (l in 1:nw){
        tmp = c(tmp, -taux[i,k]*TCOV[i, t + (l-1)*period]*A[i,t]**(lp-1)/sigma[k]**2)
      }
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), (ndeltacum[k]+1):(ndeltacum[k+1])] = matrix(tmp, ncol = nw, byrow = TRUE)
  }
  return(res)
}
##################################################################################
# matrix differential of  beta sigma
#################################################################################
mbetasigma = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw){
  res = matrix(rep(0,ng*sum(nbeta)), ncol = ng)
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nbeta[k]){
      muikt = sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,t]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
      tmp = c(tmp, -2*taux[i,k]*A[i,t]**(l-1)*(Y[i,t]-muikt)/sigma[k]**3)
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), k] = tmp
  }
  return(res)
}
##################################################################################
# matrix differential of  delta sigma
#################################################################################
mdeltasigma = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw){
  res = matrix(rep(0,ng*ng*nw), ncol = ng)
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nw){
      muikt = sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,t]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
      tmp = c(tmp, -2*taux[i,k]*TCOV[i, t + (l-1)*period]*(Y[i,t]-muikt)/sigma[k]**3)
    }
    res[(ndeltacum[k]+1):(ndeltacum[k+1]), k] = tmp
  }
  return(res)
}
##################################################################################
# matrix differential of  sigma**2 for t and i
#################################################################################
msigma = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw){
  tmp = c()
  for (k in 1:ng){
    muikt = sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,t]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
    tmp = c(tmp, -taux[i,k]*(-sigma[k]**2+3*(Y[i,t]-muikt)**2)/sigma[k]**4)
  }
  return(diag(tmp))
}
##################################################################################
# defintion of function Bikl and Sk
##################################################################################
Bikl <- function(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  return(sum(A[i,]**(l-1)*(Y[i,]-sapply(1:period, function(s){
    sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,s]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, s, k, ndeltacum)
  }))/sigma[k]**2))
}
Dikl <- function(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  return(sum(TCOV[i, 1:period + (l-1)*period]*(Y[i,]-sapply(1:period, function(s){
    sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,s]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, s, k, ndeltacum)
  }))/sigma[k]**2))
}
Sik <- function(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  return(-sum((sigma[k]**2-(Y[i,]-sapply(1:period, function(s){
    sum(beta[(nbetacum[k]+1):nbetacum[k+1]]*A[i,s]**(0:(nbeta[k]-1)))+WitEM(TCOV, period, delta, nw, i, s, k, ndeltacum)
  }))**2)/sigma[k]**3))
}
# matrix cov pi betak
covPiBetak <- function(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi){
  rcovPiBetak =  matrix(rep(0,(ng-1)*nbeta[k]), ncol = nbeta[k])
  for  (kp in 1:(ng-1)){
    for (l in 1:nbeta[k]){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + Bikl(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + Bikl(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      rcovPiBetak[kp,l] =  tmp
    }
  }
  return(rcovPiBetak)
}
# matrix cov pi deltak
covPiDeltak <- function(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi){
  rcovPiDeltak =  matrix(rep(0,(ng-1)*nw), ncol = nw)
  for  (kp in 1:(ng-1)){
    for (l in 1:nw){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + Dikl(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + Dikl(i, k, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      rcovPiDeltak[kp,l] = tmp
    }
  }
  return(rcovPiDeltak)
}
# cov betak sigma
covBetaSigmak <- function(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  rcovBetaSigmak = matrix(rep(0,nbeta[k]*ng), nrow = nbeta[k])
  for (kp in 1:nbeta[k]){
    for (l in 1:ng){
      tmp = 0
      for (i in 1:n){
        if (k==l){
          tmp = tmp + Bikl(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Sik(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
        else{
          tmp = tmp - Bikl(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Sik(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
      }
      rcovBetaSigmak[kp,l] = tmp
    }
  }
  return(rcovBetaSigmak)
}
# matrix cov Betak Betal
covBetakBetal <- function(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  rcovBetakBetal = matrix(rep(0, (nbeta[k]*nbeta[l])), nrow = nbeta[k])
  if (k==l){
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + Bikl(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Bikl(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp - Bikl(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Bikl(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }
  return(rcovBetakBetal)
}
# matrix cov Betak Deltal
covBetakDeltal <- function(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  rcovBetakDeltal = matrix(rep(0, (nbeta[k]*nw)), nrow = nbeta[k])
  if (k==l){
    for (p in 1:nbeta[k]){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + Bikl(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Dikl(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
        rcovBetakDeltal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nbeta[k]){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          tmp = tmp - Bikl(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Dikl(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
        rcovBetakDeltal[p,q] = tmp
      }
    }
  }
  return(rcovBetakDeltal)
}
# matrix cov Deltak Deltal
covDeltakDeltal <- function(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  rcovDeltakDeltal = matrix(rep(0, (nw)**2), nrow = nw)
  if (k==l){
    for (p in 1:nw){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + Dikl(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Dikl(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
        rcovDeltakDeltal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nw){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          tmp = tmp - Dikl(i, k, p, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Dikl(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
        rcovDeltakDeltal[p,q] = tmp
      }
    }
  }
  return(rcovDeltakDeltal)
}
# cov deltak sigma
covDeltaSigmak <- function(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw){
  rcovDeltaSigmak = matrix(rep(0,nw*ng), nrow = nw)
  for (kp in 1:nw){
    for (l in 1:ng){
      tmp = 0
      if (k==l){
        for (i in 1:n){
          tmp = tmp + Dikl(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Sik(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
      }else{
        for (i in 1:n){
          tmp = tmp - Dikl(i, k, kp, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Sik(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
      }
    }
    rcovDeltaSigmak[kp,l] = tmp
    return(rcovDeltaSigmak)
  }
}
##################################################################################
# Main function
##################################################################################
IEM <- function(paramEM, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw, refgr){
  nsigma = ng
  period = ncol(A)
  #################################################################################
  # Compute of matrix information
  #################################################################################
  if (nx == 1){
    pi = c(paramEM[1:(ng-1)], 1-sum(paramEM[1:(ng-1)]))
    beta = paramEM[(ng):(ng+sum(nbeta)-1)]
    sigma = paramEM[c((ng+sum(nbeta)):(ng+sum(nbeta)+ng-1))]
    delta = paramEM[c(-c(1:(ng+sum(nbeta)+ng-1)))]
  }else{
    theta = c(paramEM[1:(ng*nx)])
    pi = theta
    # we fit the group ref with 0
    theta = theta - theta[((refgr-1)*nx+1):((refgr-1)*nx+nx)]
    beta = paramEM[(ng*nx+1):(ng*nx+sum(nbeta))]
    sigma = paramEM[c((ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+ng))]
    delta = paramEM[c(-c(1:(ng*nx+sum(nbeta)+ng)))]
  }
  taux = ftaux(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)
  period = ncol(A)
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
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
  B = matrix(rep(0, (sum(nbeta)+ng+nw*ng)**2), ncol =  sum(nbeta)+ng+ng*nw) # we take off the dimension of mPi
  if (nw !=0){
    for (i in 1:n){
      for (t in (1:period)){
        B = B + rbind(cbind(mbeta(i,t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw),
                            mbetadelta(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw),
                            mbetasigma(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                      cbind(t(mbetadelta(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                            mdelta(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw),
                            mdeltasigma(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                      cbind(t(mbetasigma(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                            t(mdeltasigma(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                            msigma(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)))
      }
    }
  }else{
    for (i in 1:n){
      for (t in (1:period)){
        B = B + rbind(cbind(mbeta(i,t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw),
                            mbetasigma(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                      cbind(t(mbetasigma(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)),
                            msigma(i,t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)))
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
      covPiBeta = cbind(covPiBeta, covPiBetak(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi))
    }
    rcovPiBetak =  matrix(rep(0,(ng-1)*nbeta[ng]), ncol=nbeta[ng])
    for  (kp in 1:(ng-1)){
      for (l in 1:nbeta[ng]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + Bikl(i, ng, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
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
                tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*Bikl(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
              }
            }else{
              for (i in 1:n){
                tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*Bikl(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
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
        covPiDelta = cbind(covPiDelta, covPiDeltak(k, ng, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi))
      }
      rcovPiDeltak =  matrix(rep(0,(ng-1)*nw), ncol=nw)
      for  (kp in 1:(ng-1)){
        for (l in 1:nw){
          tmp = 0
          for (i in 1:n){
            tmp = tmp + Dikl(i, ng, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
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
                  tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*Dikl(i, k, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
                }
              }else{
                for (i in 1:n){
                  tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*Dikl(i, l, q, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
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
        covBetaDelta[(nbetacum[k]+1):(nbetacum[k+1]),(ndeltacum[l]+1):(ndeltacum[l+1])]  = covBetakDeltal(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
      }
    }
    ##################################################################################
    # matrix cov delta
    ##################################################################################
    covDelta = matrix(rep(0,sum(nw*ng)**2), ncol = nw*ng)
    for (k in 1:ng){
      for (l in 1:ng){
        covDelta[(ndeltacum[k]+1):(ndeltacum[k+1]),(ndeltacum[l]+1):(ndeltacum[l+1])] = covDeltakDeltal(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
      }
    }
    ##################################################################################
    # Matrix cov delta sigma
    ##################################################################################
    covDeltaSigma = c()
    for (k in 1:ng){
      covDeltaSigma = rbind(covDeltaSigma, covDeltaSigmak(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw))
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
              tmp = tmp + taux[i,k]*Sik(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*((1-taux[i,k])/pi[k]+taux[i,ng]/pi[ng])
            }
          }
          else{
            for (i in 1:n){
              tmp = tmp + taux[i,l]*Sik(i,l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)/pi[l]*(-taux[i,k]/pi[k]+taux[i,ng]/pi[ng])
            }
          }
        }else{
          tmp = 0
          for (i in 1:n){
            tmp = tmp + taux[i,ng]*Sik(i, ng, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*((1-taux[i,ng])/pi[ng]+taux[i,k]/pi[k])
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
              tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*Sik(i,k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
            }
          }else{
            for (i in 1:n){
              tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*Sik(i,l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
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
      covBeta[(nbetacum[k]+1):(nbetacum[k+1]),(nbetacum[l]+1):(nbetacum[l+1])] = covBetakBetal(k, l, n, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)
    }
  }

  ##################################################################################
  # Matrix cov beta sigma
  ##################################################################################
  covBetaSigma = c()
  for (k in 1:ng){
    covBetaSigma = rbind(covBetaSigma, covBetaSigmak(k, nbeta, n, ng, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw))
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
          tmp = tmp + Sik(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)**2*taux[i,k]*(1-taux[i,k])
        }
      }else{
        for (i in 1:n){
          tmp = tmp - Sik(i, k, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*Sik(i, l, nbeta, A, Y, period, beta, sigma, taux, nbetacum, TCOV, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
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
