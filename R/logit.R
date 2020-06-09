#################################################################################
#
# Likelihood
#
#################################################################################
#################################################################################
# Function g
#################################################################################
gkLogit <- function(beta, i, k, nbeta, A, Y, TCOV, delta, nw){
  period = ncol(A)
  tmp = sapply(1:period, function(s){
    sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
  })
  tmp1 = (1-1/(1+exp(tmp)))**Y[i,]*(1/(1+exp(tmp)))**(1-Y[i,])
  return(prod(tmp1))
}
#################################################################################
# Differential of likelihood by thetak
#################################################################################
difLthetakLogit <- function(param, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw){
  ntheta = ncol(X)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)))]
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
      tmp = exp(sapply(1:ng,function(s){theta[((s-1)*ntheta+1):(s*ntheta)]%*%X[i,]}))
      tmp1 = tmp[k]*sum(sapply(1:ng, function(s){tmp[s]*(gkLogit(beta, i, k, nbeta, A, Y, TCOV, delta, nw) - gkLogit(beta, i , s, nbeta, A, Y, TCOV, delta, nw))}))
      tmp2 = sum(sapply(1:ng, function(s){tmp[s]*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)}))
      a = a + X[i,l]*tmp1/(sum(tmp)*tmp2)
    }
    thetas = c(thetas, a)
  }
  return(thetas)
}
#################################################################################
# Differential of likelihood by betak
#################################################################################
difLbetakLogit <- function(param, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)))]
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
  for (l in 1:nbeta[k]){
    a = 0
    for (i in 1:n){
      tmp = sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
      })
      tmp1 = (1-1/(1+exp(tmp)))**Y[i,]*(1/(1+exp(tmp)))**(1-Y[i,])
      tmp2 = sum(sapply(1:period, function(s){
        A[i,s]**(l-1)*exp(Y[i,s]*tmp[s])/(1+exp(tmp[s]))**2*(Y[i,s]-(1-Y[i,s])*exp(tmp[s]))*prod(tmp1[-s])
      }))
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)}))
      a = a + piik(theta, i, k, ng, X)*tmp2/s
    }
    betas = c(betas, a)
  }
  return(betas)
}
#################################################################################
# Differential of likelihood by deltak
#################################################################################
difLdeltakLogit <- function(param, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)))]
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
  for (l in 1:nw){
    a = 0
    for (i in 1:n){
      tmp = sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
      })
      tmp1 = (1-1/(1+exp(tmp)))**Y[i,]*(1/(1+exp(tmp)))**(1-Y[i,])
      tmp2 = sum(sapply(1:period, function(s){
        TCOV[i, s + (l-1)*period]*exp(Y[i,s]*tmp[s])/(1+exp(tmp[s]))**2*(Y[i,s]-(1-Y[i,s])*exp(tmp[s]))*prod(tmp1[-s])
      }))
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)}))
      a = a + piik(theta, i, k, ng, X)*tmp2/s
    }
    deltas = c(deltas, a)
  }
  return(deltas)
}
#################################################################################
# Differential of likelihood
#################################################################################
difLLogit <- function(param, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw){
  out =c()
  for (k in 1:ng){
    out = c(out, difLthetakLogit(param, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw))
  }
  for (k in 1:ng){
    out = c(out, difLbetakLogit(param, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw))
  }
  if (nw != 0){
    for (k in 1:ng){
      out = c(out, difLdeltakLogit(param, k, ng, nx, n, nbeta, A, Y, X, TCOV, nw))
    }
  }
  return(out)
}
#################################################################################
# Likelihood
#################################################################################
LikelihoodLogit <- function(param, ng, nx, n, nbeta, A, Y, X, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)))]
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
      piik(theta, i, s, ng, X)*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)
    }))
    )
  }
  return(a)
}


#################################################################################
#
# EM function for LOGIT
#
#################################################################################

#################################################################################
# function rate
#################################################################################
#' Calculate the tau in EM algorithm for a Logit model.
#'
#' @param pi A vector of probability.
#' @param beta A vector.
#' @param ng a interger.
#' @param n a integer.
#' @param nbeta a vector of integer.
#' @param A a matrix of real which is time index.
#' @param Y a matrix of real which is ?.
#' @return A vector wich correspond to tau.
ftauxLogit <-function(pi, beta, ng, n, nbeta, A, Y, TCOV, delta, nw, nx, X){
  j = 1
  betatmp = beta
  beta = list()
  for (i in 1:length(nbeta)){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  deltatmp = delta
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  tmp2 = c()
  if (nx == 1){
    for (k in 1:ng){
      tmp1 = sapply(1:n, function(s){gkLogit(beta, s, k, nbeta, A, Y, TCOV, delta, nw)})
      tmp2 = cbind(tmp2, tmp1*pi[k])
    }
  }else{
    for (k in 1:ng){
      tmp1 = sapply(1:n, function(s){gkLogit(beta, s, k, nbeta, A, Y, TCOV, delta, nw)})
      tmp2 = cbind(tmp2, tmp1*piik(pi, i, k, ng, X))
    }
  }
  tmp3 =c()
  for (k in 1:ng){
    tmp3 = cbind(tmp3, 1/(1+ (rowSums(tmp2)-tmp2[,k])/tmp2[,k]))
  }
  return(tmp3)
}
#################################################################################
# Q function for betak. The goal is to sum all Qbetak after
#################################################################################
Qbetak <- function(beta, taux, k, n, nbeta, A, Y, TCOV, delta, nw){
  period = ncol(A)
  a = 0
  nbetacum = cumsum(c(0, nbeta))
  for (i in 1:n){
    tmp = sapply(1:period, function(s){
      sum(beta*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
    })
    b = sum(Y[i,]*tmp-log(1+exp(tmp)))
    a = a + taux[i,k]*b
  }
  return(a)
}
#################################################################################
# Q function for deltak. The goal is to sum all Qdeltak after
#################################################################################
Qdeltak <- function(beta, taux, k, n, nbeta, A, Y, TCOV, delta, nw){
  period = ncol(A)
  a = 0
  nbetacum = cumsum(c(0, nbeta))
  for (i in 1:n){
    tmp = sapply(1:period, function(s){
      sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
    })
    b = sum(Y[i,]*tmp-log(1+exp(tmp)))
    a = a + taux[i,k]*b
  }
  return(a)
}
#################################################################################
# Differential of betak Q function
#################################################################################
difQbetak <- function(beta, taux, k, n, nbeta, A, Y, TCOV, delta, nw){
  period = ncol(A)
  nbetacum = cumsum(c(0, nbeta))
  betas = c()
  for (l in 1:nbeta[k]){
    a=0
    for (i in 1:n){
      tmp = exp(sapply(1:period, function(s){
        sum(beta*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
      }))
      a = a+ sum(taux[i,k]*A[i,]**(l-1)*(Y[i,]-tmp/(1+tmp)))
    }
    betas = c(betas, a)
  }
  return(betas)
}
#################################################################################
# Differential of deltak Q function
#################################################################################
difQdeltak <- function(delta, taux, k, n, nbeta, A, Y, TCOV, beta, nw){
  period = ncol(A)
  nbetacum = cumsum(c(0, nbeta))
  deltas = c()
  for (l in 1:nw){
    a=0
    for (i in 1:n){
      tmp = exp(sapply(1:period, function(s){
        sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,s]**(0:(nbeta[k]-1))) +  Wit(TCOV, period, delta, nw, i, s, k)
      }))
      a = a+ sum(taux[i,k]*TCOV[i, 1:period + (l-1)*period]*(Y[i,]-tmp/(1+tmp)))
    }
    deltas = c(deltas, a)
  }
  return(deltas)
}
#################################################################################
# Q function for beta
#################################################################################
Qbeta <- function(beta, taux, ng, n, nbeta, A, Y, TCOV, delta, nw){
  nbetacum = cumsum(c(0, nbeta))
  sum(sapply(1:ng, function(s){Qbetak(beta[(nbetacum[s]+1):(nbetacum[s+1])] ,taux, s, n, nbeta, A, Y, TCOV, delta, nw)}))
}
#################################################################################
# Differential of Q beta
#################################################################################
difQbeta <- function(beta, taux, ng, n, nbeta, A, Y, TCOV, delta, nw){
  nbetacum = cumsum(c(0, nbeta))
  tmp =c()
  for (k in 1:ng){
    tmp = c(tmp, difQbetak(beta[(nbetacum[k]+1):nbetacum[k+1]], taux, k, n, nbeta, A, Y, TCOV, delta, nw))
  }
  return(tmp)
}
#################################################################################
# Q function for delta
#################################################################################
Qdelta <- function(delta, taux, ng, n, nbeta, A, Y, TCOV, beta, nw){
  sum(sapply(1:ng, function(s){Qdeltak(beta ,taux, s, n, nbeta, A, Y, TCOV, delta, nw)}))
}
#################################################################################
# Differential of Q delta
#################################################################################
difQdelta <- function(delta, taux, ng, n, nbeta, A, Y, TCOV, beta, nw){
  tmp =c()
  for (k in 1:ng){
    tmp = c(tmp, difQdeltak(beta, taux, k, n, nbeta, A, Y, TCOV, delta, nw))
  }
  return(tmp)
}
#################################################################################
# Likelihood calculation for EM algorithm
#################################################################################
likelihoodEMLOGIT <- function(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw){
  likeli = 0
  betatmp = beta
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
      pi[s]*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)
    }))
    )
  }
  return(likeli)
}
#################################################################################
# EM algorithm
#################################################################################
EMLogit <- function(param, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr){
  period = ncol(A)
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    delta = param[-c(1:(ng+sum(nbeta)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    delta = param[-c(1:(ng*nx+sum(nbeta)))]
    pi[c(((refgr-1)*nx+1):(nx*refgr))] = 0
  }
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEMLOGIT(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw), "\n"))
    }else{
      cat(paste(-LikelihoodLogit(c(pi, beta, delta), ng, nx, n, nbeta, A, Y, X, TCOV, nw), "\n"))
    }
    # E-step
    taux = ftauxLogit(pi, beta, ng, n, nbeta, A, Y, TCOV, delta, nw, nx, X)
    # M-step
    newbeta = c()
    for (k in 1:ng){
      newbeta  = c(newbeta, optim(beta[(nbetacum[k]+1):(nbetacum[k+1])], Qbetak, difQbetak, method = "BFGS", taux=taux, control = list(fnscale=-1),
                                       k=k,n = n, nbeta =  nbeta, A=A, Y = Y, TCOV = TCOV, delta = delta, nw = nw)$par)
    }
    newdelta = c()
    if (nw !=0){
      #newdelta = optim(delta, Qdelta, difQdelta, method = "BFGS", taux=taux, control = list(fnscale=-1),
      #                 ng =ng, n = n, nbeta =  nbeta, A=A, Y = Y, TCOV = TCOV, beta = beta, nw = nw)$par
      for (k in 1:ng){
        newdelta  = c(newdelta, optim(delta[(ndeltacum[k]+1):(ndeltacum[k+1])], Qdeltak, difQdeltak, method = "BFGS", taux=taux, control = list(fnscale=-1),
                                    k=k,n = n, nbeta =  nbeta, A=A, Y = Y, TCOV = TCOV, beta=beta, nw = nw)$par)
      }
    }
    if (nx == 1){
      newpi = colSums(taux)/n
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS, refgr)
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newdelta, newpi)-c(beta, delta, pi))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    delta = newdelta
    pi = newpi
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, delta)
  }else{
    param = c(pi, beta, delta)
  }
  return(param)
}
#################################################################################
# EM algorithm
#################################################################################
EMLogitSame <- function(param, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr){
  period = ncol(A)
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    delta = param[-c(1:(ng+sum(nbeta)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    delta = param[-c(1:(ng*nx+sum(nbeta)))]
    pi[c(((refgr-1)*nx+1):(nx*refgr))] = 0
  }
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEMLOGIT(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw), "\n"))
    }else{
      cat(paste(-LikelihoodLogit(c(pi, beta, delta), ng, nx, n, nbeta, A, Y, X, TCOV, nw), "\n"))
    }
    # E-step
    taux = ftauxLogit(pi, beta, ng, n, nbeta, A, Y, TCOV, delta, nw, nx, X)
    # M-step
    newbeta = c()
    newbeta = constrOptim(theta =beta, f=Qbeta,
                          grad=difQbeta, 
                          method = "BFGS",
                          control =  list(fnscale=-1),
                          ui = rbind(c(0,1,0,-1), c(0,-1,0,1)), ci = c(-0.01,-0.01),
                          taux=taux, ng=ng, n=n, nbeta=nbeta,
                          A=A, Y=Y, TCOV=TCOV, delta=delta, nw=nw,)$par
    newdelta = c()
    if (nw !=0){
      #newdelta = optim(delta, Qdelta, difQdelta, method = "BFGS", taux=taux, control = list(fnscale=-1),
      #                 ng =ng, n = n, nbeta =  nbeta, A=A, Y = Y, TCOV = TCOV, beta = beta, nw = nw)$par
      for (k in 1:ng){
        newdelta  = c(newdelta, optim(delta[(ndeltacum[k]+1):(ndeltacum[k+1])], Qdeltak, difQdeltak, method = "BFGS", taux=taux, control = list(fnscale=-1),
                                      k=k,n = n, nbeta =  nbeta, A=A, Y = Y, TCOV = TCOV, beta=beta, nw = nw)$par)
      }
    }
    if (nx == 1){
      newpi = colSums(taux)/n
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS, refgr)
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newdelta, newpi)-c(beta, delta, pi))<10**(-8))){
      tour = itermax + 2
    }
    beta = newbeta
    delta = newdelta
    pi = newpi
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, delta)
  }else{
    param = c(pi, beta, delta)
  }
  return(param)
}
#################################################################################
# EM algorithm IRLS
#################################################################################
EMLogitIRLS <- function(param, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr){
  period = ncol(A)
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    delta = param[-c(1:(ng+sum(nbeta)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    delta = param[-c(1:(ng*nx+sum(nbeta)))]
    pi[c(((refgr-1)*nx+1):(nx*refgr))] = 0
  }
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEMLOGIT(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw), "\n"))
    }else{
      cat(paste(-LikelihoodLogit(c(pi, beta, delta), ng, nx, n, nbeta, A, Y, X, TCOV, nw), "\n"))
    }
    # E-step
    taux = ftauxLogit(pi, beta, ng, n, nbeta, A, Y, TCOV, delta, nw, nx, X)
    # M-step
    newbeta = c()
    newdelta = c()
    newbetae = c()
    newdeltae = c()
    for (k in 1:ng){
      newbetaIRLS = c()
      betaIRLS = beta[(nbetacum[k]+1):(nbetacum[k+1])]
      newdeltaIRLS = c()
      deltaIRLS = delta[((k-1)*nw+1):(k*nw)]
      precIRLS = 1
      while(all(abs(precIRLS)>10**(-6))){
        tmp = c()
        tmpr = c()
        tmp1 = c()
        tmp2 = c()
        tmpw = c()
        tmp1w = c()
        for (i in 1:n){
          for (t in 1:period){
            tmp = c(tmp, A[i,t]**(0:(nbeta[k]-1)))
            betaAit = sum(betaIRLS*A[i,t]**(0:(nbeta[k]-1)))
            deltaWit = WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
            rhoikt = exp(betaAit + deltaWit)/(1+exp(betaAit + deltaWit))
            tmpr = c(tmpr, rhoikt*(1-rhoikt))
            tmp1 = c(tmp1, betaAit + (Y[i,t]-rhoikt)/(rhoikt*(1-rhoikt)))
            tmpw = c(tmpw, TCOV[i, t + 0:((nw-1)*period)])
            tmp1w = c(tmp1w, deltaWit + (Y[i,t]-rhoikt)/(rhoikt*(1-rhoikt)))
          }
          tmp2 = c(tmp2, rep(taux[i,k], period))
        }
        Aw = matrix(tmp, nrow = nbeta[k], byrow = FALSE)
        W = diag(tmpr)
        Z = diag(tmp2)
        S = matrix(tmp1, ncol =1)
        # QR = qr(Aw%*%(Z*W)%*%t(Aw))
        # Q = qr.Q(QR)
        # R = qr.R(QR)
        # qr.solve(Aw%*%(Z*W)%*%t(Aw),Aw%*%(Z*W)%*%S)
        QR = qr(t(Aw%*%sqrt(Z*W)))
        Q = qr.Q(QR)
        R = qr.R(QR)
        newbetaIRLS = backsolve(R,t(Q)%*%(sqrt(Z*W)%*%S))
        if (nw !=0){
          Ww = matrix(tmpw, nrow = nw, byrow = FALSE)
          Sw = matrix(tmp1w, ncol =1)
          QR = qr(t(Ww%*%sqrt(Z*W)))
          Q = qr.Q(QR)
          R = qr.R(QR)
          newdeltaIRLS = backsolve(R,t(Q)%*%(sqrt(Z*W)%*%Sw))
        }else{
          newdeltaIRLS =c()
        }
        precIRLS = c(betaIRLS-newbetaIRLS, deltaIRLS-newdeltaIRLS)
        betaIRLS = newbetaIRLS
        deltaIRLS = newdeltaIRLS
      }
      newbeta = c(newbeta, betaIRLS)
      newdelta = c(newdelta, deltaIRLS)
    }
    if (nx == 1){
      newpi = colSums(taux)/n
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS, refgr)
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newdelta, newpi)-c(beta, delta, pi))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    delta = newdelta
    pi = newpi
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, delta)
  }else{
    param = c(pi, beta, delta)
  }
  return(param)
}


#################################################################################
# Compute of matrix information
#################################################################################
#################################################################################
# General functions
#################################################################################
# Function exp
#################################################################################
fexp <- function(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw){
  tmp = exp(sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,t]**(0:(nbeta[k]-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum))
  return(tmp/(1+tmp)**2)
}
#################################################################################
# matrix differential of beta**2 for t and i
#################################################################################
mbetaL = function(i, t, ng, nbeta, nbetacum, A, taux, beta, TCOV, period, delta, ndeltacum, nw){
  res = matrix(rep(0,(sum(nbeta))**2), ncol = sum(nbeta))
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nbeta[k]){
      for (lp in 1:nbeta[k]){
        tmp = c(tmp, -taux[i,k]*A[i,t]**(l+lp-2)*fexp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
      }
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), (nbetacum[k]+1):(nbetacum[k+1])] = matrix(tmp, ncol = nbeta[k])
  }
  return(res)
}
#################################################################################
# matrix differential of delta**2 for t and i
#################################################################################
mdeltaL = function(i, t, ng, nbeta, nbetacum, A, taux, beta, TCOV, period, delta, ndeltacum, nw){
  res = matrix(rep(0,(sum(nw*ng))**2), ncol = sum(nw*ng))
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nw){
      for (lp in 1:nw){
        tmp = c(tmp, -taux[i,k]*TCOV[i, t + (l+lp-2)*period]*fexp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
      }
    }
    res[(ndeltacum[k]+1):(ndeltacum[k+1]), (ndeltacum[k]+1):(ndeltacum[k+1])] = matrix(tmp, ncol = nw)
  }
  return(res)
}
##################################################################################
# defintion of function Bikl and Sk
##################################################################################
BLikl <- function(i, k, l, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw){
  tmp = sapply(1:period, function(s){
    exp(sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,s]**(0:(nbeta[k]-1)))+ + WitEM(TCOV, period, delta, nw, i, s, k, ndeltacum))
  })
  return(sum(A[i,]**(l-1)*(Y[i,]-tmp/(1+tmp))))
}
DLikl <- function(i, k, l, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw){
  tmp = sapply(1:period, function(s){
    exp(sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,s]**(0:(nbeta[k]-1)))+ + WitEM(TCOV, period, delta, nw, i, s, k, ndeltacum))
  })
  return(sum(TCOV[i, 1:period + (l-1)*period]*(Y[i,]-tmp/(1+tmp))))
}
# matrix cov Betak Betal
covBetakBetalL <- function(k,l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw){
  rcovBetakBetal = matrix(rep(0, (nbeta[k]*nbeta[l])), nrow = nbeta[k])
  if (k==l){
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + BLikl(i, k, p, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*BLikl(i, k, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp - BLikl(i, k, p, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*BLikl(i, l, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }
  return(rcovBetakBetal)
}
# matrix cov Deltak Deltal
covDeltakDeltalL <- function(k,l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw){
  rcovDeltakDeltal = matrix(rep(0, (nw**2)), nrow = nw)
  if (k==l){
    for (p in 1:nw){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + DLikl(i, k, p, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*DLikl(i, k, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
        rcovDeltakDeltal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nw){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          tmp = tmp - DLikl(i, k, p, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*DLikl(i, l, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
        rcovDeltakDeltal[p,q] = tmp
      }
    }
  }
  return(rcovDeltakDeltal)
}
# matrix cov pi betak
covPiBetakL <- function(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw){
  rcovPiBetak =  matrix(rep(0,(ng-1)*nbeta[k]), ncol=nbeta[k])
  for  (kp in 1:(ng-1)){
    for (l in 1:nbeta[k]){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + BLikl(i, k, l, nbeta, nbetacum, A, Y,  beta, TCOV, period, delta, ndeltacum, nw)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + BLikl(i, k, l, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      rcovPiBetak[kp,l] =  tmp
    }
  }
  return(rcovPiBetak)
}
# matrix cov pi deltak
covPiDeltakL <- function(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw){
  rcovPiDeltak =  matrix(rep(0,(ng-1)*nw), ncol=nw)
  for  (kp in 1:(ng-1)){
    for (l in 1:nw){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + DLikl(i, k, l, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + DLikl(i, k, l, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      rcovPiDeltak[kp,l] =  tmp
    }
  }
  return(rcovPiDeltak)
}
# matrix cov betak deltal
covBetakDeltal <- function(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw){
  rcovBetakDeltal = matrix(rep(0, nbeta[k]*nw), ncol = nw)
  for (p in 1:nbeta[k]){
    for (q in 1 :nw){
      tmp = 0
      if (k == l){
        for (i in 1:n){
          tmp = tmp + BLikl(i, k, p, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*DLikl(i, k, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*(1-taux[i,k])
        }
      }else{
        for (i in 1:n){
          tmp = tmp - BLikl(i, k, p, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*DLikl(i, l, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,k]*taux[i,l]
        }
      }
      rcovBetakDeltal[p,q] = tmp
    }
  }
  return(rcovBetakDeltal)
}
#################################################################################
#################################################################################
IEML <- function(paramEM, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, refgr){
  if (nx == 1){
    pi = c(paramEM[1:(ng-1)], 1-sum(paramEM[1:(ng-1)]))
    beta = paramEM[(ng):(ng+sum(nbeta)-1)]
    delta = paramEM[c(-c(1:(ng+sum(nbeta)-1)))]
  }else{
    theta = c(paramEM[1:(ng*nx)])
    pi = theta
    # we fit the group ref with 0
    theta = theta - theta[((refgr-1)*nx+1):((refgr-1)*nx+nx)]
    beta = paramEM[(ng*nx+1):(ng*nx+sum(nbeta))]
    delta = paramEM[c(-c(1:(ng*nx+sum(nbeta))))]
  }
  taux = ftauxLogit(pi, beta, ng, n, nbeta, A, Y, TCOV, delta, nw, nx, X)
  period = ncol(A)
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  #################################################################################
  # matrix differential of  Pi**2
  #################################################################################
  if (nx == 1){
    mPiL = matrix(rep(0,(ng-1)**2), ncol = ng-1)
    for (k in 1:(ng-1)){
      for (l in 1:(ng-1)){
        if (k==l){
          mPiL[k,l] = -sum(taux[,k]/pi[k]**2 + taux[,ng]/pi[ng]**2)
        }
        else{
          mPiL[k,l] = -sum(taux[,ng]/pi[ng]**2)
        }
      }
    }
  }else{
    mPiL = matrix(rep(0,(ng*nx)**2), ncol = ng*nx)
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
        mPiL[((k-1)*nx+1):((k-1)*nx+nx), ((kp-1)*nx+1):((kp-1)*nx+nx)] = matrix(tmp1, ncol = nx, byrow =TRUE)
      }
    }
    # we have to delete the part of the reference group
    mPiL = mPiL[-c(((refgr-1)*nx+1):((refgr-1)*nx+nx)), ]
    mPiL = mPiL[ ,-c(((refgr-1)*nx+1):((refgr-1)*nx+nx))]
  }
  ##################################################################################
  # matrix -B
  #################################################################################
  B = matrix(rep(0, (sum(nbeta))**2), ncol =  sum(nbeta)) # we take off the dimension of mPi
  tmp = matrix(rep(0, (sum(nw*ng))**2), ncol =  sum(nw*ng)) # we take off the dimension of mPi
  for (i in 1:n){
    for (t in (1:period)){
      B = B + mbetaL(i,t, ng, nbeta, nbetacum, A, taux, beta, TCOV, period, delta, ndeltacum, nw)
      tmp = tmp + mdeltaL(i,t, ng, nbeta, nbetacum, A, taux, beta, TCOV, period, delta, ndeltacum, nw)
    }
  }
  B = cbind(matrix(rep(0, (ng-1)*nx*(sum(nbeta))), ncol = (ng-1)*nx),
            B)
  B = rbind(cbind(mPiL, matrix(rep(0,(sum(nbeta))*(ng-1)*nx), nrow = (ng-1)*nx)),
            B)
  if (nw !=0){
    B = rbind(B, matrix(rep(0,(sum(nbeta)+ng-1)*ng*nw), nrow = ng*nw))
    B = cbind(B, rbind(matrix(rep(0,(sum(nbeta)+ng-1)*ng*nw), ncol = ng*nw),
                       tmp))
  }
  ##################################################################################
  # matrix cov pi
  ##################################################################################
  if (nx == 1){
    covPiL = matrix(rep(0,(ng-1)**2), ncol=ng-1)
    for (k in 1:(ng-1)){
      for (l in 1:(ng-1)){
        tmp = 0
        if (k==l){
          for (i in 1:n){
            tmp = tmp + taux[i,k]*(1-taux[i,k])/pi[k]**2+taux[i,ng]*(1-taux[i,ng])/pi[ng]**2+2*taux[i,k]*taux[i,ng]/(pi[k]*pi[ng])
          }
        }else{
          for (i in 1:n){
            tmp = tmp - taux[i,k]*taux[i,l]/(pi[k]*pi[l])+taux[i,k]*taux[i,ng]/(pi[ng]*pi[k])+taux[i,l]*taux[i,ng]/(pi[l]*pi[ng])+taux[i,ng]*(1-taux[i,ng])/pi[ng]**2
          }
        }
        covPiL[k,l] = tmp
      }
    }
  }else{
    covPiL = matrix(rep(0,((ng-1)*nx)**2), ncol=(ng-1)*nx)
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
                tmp = tmp -  taux[i,k]*taux[i,l]*X[i,p]*X[i,q]
              }
            }
            covPitmp[p, q] = tmp
          }
        }
        covPiL[((k-1)*nx+1):((k-1)*nx+nx), ((l-1)*nx+1):((l-1)*nx+nx)] = covPitmp
      }
    }
  }
  ##################################################################################
  # matrix cov pi beta
  ##################################################################################
  # matrix cov pi beta
  if (nx == 1){
  covPiBetaL = c()
  for (k in 1:(ng-1)){
    covPiBetaL = cbind(covPiBetaL, covPiBetakL(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw))
  }
  for (kp in 1:(ng-1)){
    rcovPiBetakL =  matrix(rep(0,(ng-1)*nbeta[ng]), ncol=nbeta[ng])
    for (l in 1:nbeta[ng]){
      tmp = 0
      for (i in 1:n){
        tmp = tmp + BLikl(i, ng, l, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
      }
    }
    rcovPiBetakL[kp,l] =  tmp
  }
  covPiBetaL = cbind(covPiBetaL, rcovPiBetakL)
  }else{
    covPiBetaL = matrix(rep(0,((ng-1)*nx)*sum(nbeta)), ncol=sum(nbeta))
    for (k in 1:(ng-1)){
      for (l in 1:ng){
        covPiBetatmp = matrix(rep(0, nx*nbeta[l]), ncol =nbeta[l])
        for (p in 1:nx){
          for (q in 1:nbeta[l]){
            tmp = 0
            if (k==l){
              for (i in 1:n){
                tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*BLikl(i, k, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)
              }
            }else{
              for (i in 1:n){
                tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*BLikl(i, l, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)
              }
            }
            covPiBetatmp[p, q] = tmp
          }
        }
        covPiBetaL[((k-1)*nx+1):((k-1)*nx+nx), (nbetacum[l]+1):(nbetacum[l+1])] = covPiBetatmp
      }
    }
  }
  ##################################################################################
  # matrix cov pi delta
  ##################################################################################
  # matrix cov pi delta
  if (nw !=0){
    if (nx == 1 ){
      covPiDeltaL = c()
      for (k in 1:(ng-1)){
        covPiDeltaL = cbind(covPiDeltaL, covPiDeltakL(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw))
      }
      for (kp in 1:(ng-1)){
        rcovPiDeltakL =  matrix(rep(0,(ng-1)*nw), ncol=nw)
        for (l in 1:nw){
          tmp = 0
          for (i in 1:n){
            tmp = tmp + DLikl(i, ng, l, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
          }
        }
        rcovPiDeltakL[kp,l] =  tmp
      }
      covPiDeltaL = cbind(covPiDeltaL, rcovPiDeltakL)
    }else{
      covPiDeltaL = matrix(rep(0,((ng-1)*nx)*sum(delta)), ncol=sum(ndelta))
      for (k in 1:(ng-1)){
        for (l in 1:ng){
          covPiDeltatmp = matrix(rep(0, nx*ndelta[l]), ncol =ndelta[l])
          for (p in 1:nx){
            for (q in 1:ndelta[l]){
              tmp = 0
              if (k==l){
                for (i in 1:n){
                  tmp = tmp + taux[i,k]*(1-taux[i,k])*X[i,p]*DLikl(i, k, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)
                }
              }else{
                for (i in 1:n){
                  tmp = tmp - taux[i,k]*taux[i,l]*X[i,p]*DLikl(i, l, q, nbeta, nbetacum, A, Y, beta, TCOV, period, delta, ndeltacum, nw)
                }
              }
              covPiDeltatmp[p, q] = tmp
            }
          }
          covPiDeltaL[((k-1)*nx+1):((k-1)*nx+nx), (ndeltacum[l]+1):(ndeltacum[l+1])] = covPiDeltatmp
        }
      }
    }
  }
  ##################################################################################
  # matrix cov beta
  ##################################################################################
  # matrix cov Beta
  covBetaL = c()
  for (l in 1:ng){
    covBetakL = c()
    for (k in 1:ng){
      covBetakL= rbind(covBetakL, covBetakBetalL(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw))
    }
    covBetaL = cbind(covBetaL, covBetakL)
  }
  ##################################################################################
  # matrix cov beta delta
  ##################################################################################
  # matrix cov Beta Delta
  if (nw !=0){
    covBetaDeltaL = c()
    for (k in 1:ng){
      covBetaDeltakL = c()
      for (l in 1:ng){
        covBetaDeltakL= rbind(covBetaDeltakL, covBetakDeltal(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw))
      }
      covBetaDeltaL = cbind(covBetaDeltaL, covBetaDeltakL)
    }
  }
  ##################################################################################
  # matrix cov delta
  ##################################################################################
  # matrix cov Delta
  if (nw !=0){
    covDeltaL = c()
    for (l in 1:ng){
      covDeltakL = c()
      for (k in 1:ng){
        covDeltakL= rbind(covDeltakL, covDeltakDeltalL(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw))
      }
      covDeltaL = cbind(covDeltaL, covDeltakL)
    }
  }
  ##################################################################################
  # Matrix of covariance of score function
  ##################################################################################
  cov =c()
  if (nw !=0){
    cov = cbind(covPiL, covPiBetaL, covPiDeltaL)
    cov = rbind(cov, cbind(t(covPiBetaL), covBetaL, covBetaDeltaL))
    cov = rbind(cov, cbind(t(covPiDeltaL), t(covBetaDeltaL), covDeltaL))
  }else{
    cov = cbind(covPiL, covPiBetaL)
    cov = rbind(cov, cbind(t(covPiBetaL), covBetaL))
  }


  ##################################################################################
  # Information matrix of Fisher
  ##################################################################################
  IEM = - B - cov
  return(sqrt(diag(solve(IEM))))
}
