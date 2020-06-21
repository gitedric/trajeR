#################################################################################
#
# Likelihood
#
#################################################################################
#################################################################################
# Function g
#################################################################################
gkZIP <- function(beta, nu, i, k, nbeta, nnu, A, Y, TCOV, delta, nw){
  period = ncol(A)
  tmp = exp(-sapply(1:period, function(s){
    sum(nu[[k]]*A[i,s]**(0:(nnu[k]-1)))
  }))
  rhoikt = 1/(1+tmp)
  lambdaikt = exp(sapply(1:period, function(s){
    sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
  }))
  ind0 = (Y[i,] == 0)
  tmp1 =c()
  if (any(ind0)){
    tmp1 = c(tmp1, rhoikt[which(ind0)]+(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)]))
  }
  ind = !ind0
  if (any(ind)){
    tmp1 = c(tmp1, (1-rhoikt[which(ind)])*lambdaikt[which(ind)]**Y[i,which(ind)]*exp(-lambdaikt[which(ind)])/factorial(Y[i,which(ind)]))
  }
  return(prod(tmp1))
}
#################################################################################
# Differential of likelihood by thetak
#################################################################################
difLthetakZIP <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
  ntheta = nx
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  nu = list()
  for (i in 1:ng){
    nu[[i]] = nutmp[j:sum(nnu[1:i])]
    j = sum(nnu[1:i]) + 1
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
      tmp1 = tmp[k]*sum(sapply(1:ng, function(s){tmp[s]*(gkZIP(beta, nu, i, k, nbeta, nnu, A, Y, TCOV, delta, nw) - gkZIP(beta, nu, i ,s, nbeta, nnu, A, Y, TCOV, delta, nw))}))
      tmp2 = sum(sapply(1:ng, function(s){tmp[s]*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)}))
      a = a + X[i,l]*tmp1/(sum(tmp)*tmp2)
    }
    thetas = c(thetas, a)
  }
  return(thetas)
}
#################################################################################
# Differential of likelihood by betak
#################################################################################
# difLbetakZIP <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
#   period = ncol(A)
#   theta = param[1:(ng*nx)]
#   betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
#   nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
#   deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
#   j = 1
#   beta = list()
#   for (i in 1:ng){
#     beta[[i]] = betatmp[j:sum(nbeta[1:i])]
#     j = sum(nbeta[1:i]) + 1
#   }
#   j = 1
#   nu = list()
#   for (i in 1:ng){
#     nu[[i]] = nutmp[j:sum(nnu[1:i])]
#     j = sum(nnu[1:i]) + 1
#   }
#   ndeltacum = cumsum(c(0, rep(nw, ng)))
#   delta = list()
#   for (i in 1:ng){
#     delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
#   }
#   betas = c()
#   for (l in 1:nbeta[k]){
#     a = 0
#     for (i in 1:n){
#       difbkl0 = 0
#       difbkl = 0
#       py0 = 1
#       py = 1
#       ind0 = (Y[i,] == 0)
#       ind = !ind0
#       tmp = exp(sapply(1:period, function(s){
#         sum(nu[[k]]*A[i,s]**(0:(nnu[k]-1)))
#       }))
#       rhoikt = tmp/(1+tmp)
#       lambdaikt = exp(sapply(1:period, function(s){
#         sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
#       }))
#       if (any(ind0)){
#         if (sum(ind0)>1){
#           tind = which(ind0)
#           for (j in 1:length(tind)){
#             difbkl0 = difbkl0 - A[i,tind][j]**(l-1)*lambdaikt[tind][j]*(1-rhoikt[tind][j])*exp(-lambdaikt[tind][j])*prod(rhoikt[tind][-j]+(1-rhoikt[tind][-j])*exp(-lambdaikt[tind][-j]))
#           }
#         } else{
#           difbkl0 = - A[i,which(ind0)]**(l-1)*lambdaikt[which(ind0)]*(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)])
#         }
#         py0 = prod(rhoikt[which(ind0)]+(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)]))
#       }
#       if (any(ind)){
#         if (sum(ind)>1){
#           tind = which(ind)
#           for (j in 1:length(tind)){
#             difbkl = difbkl + A[i,tind][j]**(l-1)*lambdaikt[tind][j]**(Y[i,tind][j])*(1-rhoikt[tind][j])*exp(-lambdaikt[tind][j])*(Y[i,tind][j]-lambdaikt[tind][j])/factorial(Y[i,tind][j])*prod((1-rhoikt[tind][-j])*lambdaikt[tind][-j]**Y[i,tind][-j]*exp(-lambdaikt[tind][-j])/factorial(Y[i,tind][-j]))
#           }
#         } else{
#           difbkl = (1-rhoikt[which(ind)])*A[i,which(ind)]**(l-1)*lambdaikt[which(ind)]**(Y[i,which(ind)])*exp(-lambdaikt[which(ind)])*(Y[i,which(ind)]-lambdaikt[which(ind)])/factorial(Y[i,which(ind)])
#         }
#         py = prod((1-rhoikt[which(ind)])*lambdaikt[which(ind)]**Y[i,which(ind)]*exp(-lambdaikt[which(ind)])/factorial(Y[i,which(ind)]))
#       }
#       s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)}))
#       a = a + piik(theta, i, k, ng, X)/s*(difbkl0*py+difbkl*py0)
#     }
#     betas = c(betas, a)
#   }
#   return(betas)
# }
difLbetakZIP <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  nu = list()
  for (i in 1:ng){
    nu[[i]] = nutmp[j:sum(nnu[1:i])]
    j = sum(nnu[1:i]) + 1
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
      difbkl0 = 0
      difbkl = 0
      py0 = 1
      py = 1
      ind0 = (Y[i,] == 0)
      ind = !ind0
      tmp = exp(sapply(1:period, function(s){
        sum(nu[[k]]*A[i,s]**(0:(nnu[k]-1)))
      }))
      rhoikt = tmp/(1+tmp)
      lambdaikt = exp(sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
      }))
      if (any(ind0)){
        if (sum(ind0)>1){
          for (t in which(ind0)){
            difbkl0 = difbkl0 - A[i,t]**(l-1)*lambdaikt[t]*(1-rhoikt[t])*exp(-lambdaikt[t])*prod(rhoikt[which(ind0)[which(ind0)!=t]]+(1-rhoikt[which(ind0)[which(ind0)!=t]])*exp(-lambdaikt[which(ind0)[which(ind0)!=t]]))
          }
        } else{
          difbkl0 = - A[i,which(ind0)]**(l-1)*lambdaikt[which(ind0)]*(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)])
        }
        py0 = prod(rhoikt[which(ind0)]+(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)]))
      }
      if (any(ind)){
        if (sum(ind)>1){
          for (t in which(ind)){
            difbkl = difbkl + A[i,t]**(l-1)*lambdaikt[t]**(Y[i,t])*(1-rhoikt[t])*exp(-lambdaikt[t])*(Y[i,t]-lambdaikt[t])/factorial(Y[i,t])*prod((1-rhoikt[which(ind)[which(ind)!=t]])*lambdaikt[which(ind)[which(ind)!=t]]**Y[i,which(ind)[which(ind)!=t]]*exp(-lambdaikt[which(ind)[which(ind)!=t]])/factorial(Y[i,which(ind)[which(ind)!=t]]))
          }
        } else{
          difbkl = (1-rhoikt[which(ind)])*A[i,which(ind)]**(l-1)*lambdaikt[which(ind)]**(Y[i,which(ind)])*exp(-lambdaikt[which(ind)])*(Y[i,which(ind)]-lambdaikt[which(ind)])/factorial(Y[i,which(ind)])
        }
        py = prod((1-rhoikt[which(ind)])*lambdaikt[which(ind)]**Y[i,which(ind)]*exp(-lambdaikt[which(ind)])/factorial(Y[i,which(ind)]))
      }
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)}))
      a = a + piik(theta, i, k, ng, X)/s*(difbkl0*py+difbkl*py0)
    }
    betas = c(betas, a)
  }
  return(betas)
}
#################################################################################
# Differential of likelihood by nuk
#################################################################################
# difLnukZIP <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
#   period = ncol(A)
#   theta = param[1:(ng*nx)]
#   betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
#   nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
#   deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
#   j = 1
#   beta = list()
#   for (i in 1:ng){
#     beta[[i]] = betatmp[j:sum(nbeta[1:i])]
#     j = sum(nbeta[1:i]) + 1
#   }
#   j = 1
#   nu = list()
#   for (i in 1:ng){
#     nu[[i]] = nutmp[j:sum(nnu[1:i])]
#     j = sum(nnu[1:i]) + 1
#   }
#   ndeltacum = cumsum(c(0, rep(nw, ng)))
#   delta = list()
#   for (i in 1:ng){
#     delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
#   }
#   nus = c()
#   for (l in 1:nnu[k]){
#     a = 0
#     for (i in 1:n){
#       difbkl0 = 0
#       difbkl = 0
#       py0 = 1
#       py = 1
#       ind0 = (Y[i,] == 0)
#       ind = !ind0
#       enuikt = exp(sapply(1:period, function(s){
#         sum(nu[[k]]*A[i,s]**(0:(nnu[k]-1)))
#       }))
#       rhoikt = enuikt/(1+enuikt)
#       lambdaikt = exp(sapply(1:period, function(s){
#         sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
#       }))
#       if (any(ind0)){
#         if (sum(ind0)>1){
#           tind = which(ind0)
#           for (j in 1:length(tind)){
#             difbkl0 = difbkl0 + A[i,tind][j]**(l-1)*rhoikt[tind][j]*(1-exp(-lambdaikt[tind][j]))/(1+enuikt[tind][j])*prod(rhoikt[tind][-j]+(1-rhoikt[tind][-j])*exp(-lambdaikt[tind][-j]))
#           }
#         } else{
#           difbkl0 = A[i,which(ind0)]**(l-1)*rhoikt[which(ind0)]*(1-exp(-lambdaikt[which(ind0)]))/(1+enuikt[which(ind0)])
#         }
#         py0 = prod(rhoikt[which(ind0)]+(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)]))
#       }
#       if (any(ind)){
#         if (sum(ind)>1){
#           tind = which(ind)
#           for (j in 1:length(tind)){
#             difbkl = difbkl -A[i,tind][j]**(l-1)*rhoikt[tind][j]*lambdaikt[tind][j]**(Y[i,tind][j])*exp(-lambdaikt[tind][j])/(factorial(Y[i,tind][j])*(1+enuikt[tind][j]))*prod((1-rhoikt[tind][-j])*lambdaikt[tind][-j]**Y[i,tind][-j]*exp(-lambdaikt[tind][-j])/factorial(Y[i,tind][-j]))
#           }
#         } else{
#           difbkl = -A[i,which(ind)]**(l-1)*rhoikt[which(ind)]*lambdaikt[which(ind)]**(Y[i,which(ind)])*exp(-lambdaikt[which(ind)])/(factorial(Y[i,which(ind)])*(1+enuikt[which(ind)]))
#         }
#         py = prod((1-rhoikt[which(ind)])*lambdaikt[which(ind)]**Y[i,which(ind)]*exp(-lambdaikt[which(ind)])/factorial(Y[i,which(ind)]))
#       }
#       s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)}))
#       a = a + piik(theta, i, k, ng, X)/s*(difbkl0*py+difbkl*py0)
#     }
#     nus = c(nus, a)
#   }
#   return(nus)
# }
difLnukZIP <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  nu = list()
  for (i in 1:ng){
    nu[[i]] = nutmp[j:sum(nnu[1:i])]
    j = sum(nnu[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  nus = c()
  for (l in 1:nnu[k]){
    a = 0
    for (i in 1:n){
      difbkl0 = 0
      difbkl = 0
      py0 = 1
      py = 1
      ind0 = (Y[i,] == 0)
      ind = !ind0
      enuikt = exp(sapply(1:period, function(s){
        sum(nu[[k]]*A[i,s]**(0:(nnu[k]-1)))
      }))
      rhoikt = enuikt/(1+enuikt)
      lambdaikt = exp(sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
      }))
      if (any(ind0)){
        if (sum(ind0)>1){
          for (t in which(ind0)){
            difbkl0 = difbkl0 + A[i,t]**(l-1)*rhoikt[t]*(1-exp(-lambdaikt[t]))/(1+enuikt[t])*prod(rhoikt[which(ind0)[which(ind0)!=t]]+(1-rhoikt[which(ind0)[which(ind0)!=t]])*exp(-lambdaikt[which(ind0)[which(ind0)!=t]]))
          }
        } else{
          difbkl0 = A[i,which(ind0)]**(l-1)*rhoikt[which(ind0)]*(1-exp(-lambdaikt[which(ind0)]))/(1+enuikt[which(ind0)])
        }
        py0 = prod(rhoikt[which(ind0)]+(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)]))
      }
      if (any(ind)){
        if (sum(ind)>1){
          for (t in which(ind)){
            difbkl = difbkl -A[i,t]**(l-1)*rhoikt[t]*lambdaikt[t]**(Y[i,t])*exp(-lambdaikt[t])/(factorial(Y[i,t])*(1+enuikt[t]))*prod((1-rhoikt[which(ind)[which(ind)!=t]])*lambdaikt[which(ind)[which(ind)!=t]]**Y[i,which(ind)[which(ind)!=t]]*exp(-lambdaikt[which(ind)[which(ind)!=t]])/factorial(Y[i,which(ind)[which(ind)!=t]]))
          }
        } else{
          difbkl = -A[i,which(ind)]**(l-1)*rhoikt[which(ind)]*lambdaikt[which(ind)]**(Y[i,which(ind)])*exp(-lambdaikt[which(ind)])/(factorial(Y[i,which(ind)])*(1+enuikt[which(ind)]))
        }
        py = prod((1-rhoikt[which(ind)])*lambdaikt[which(ind)]**Y[i,which(ind)]*exp(-lambdaikt[which(ind)])/factorial(Y[i,which(ind)]))
      }
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)}))
      a = a + piik(theta, i, k, ng, X)/s*(difbkl0*py+difbkl*py0)
    }
    nus = c(nus, a)
  }
  return(nus)
}
#################################################################################
# Differential of likelihood by deltak
#################################################################################
# difLdeltakZIP <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
#   period = ncol(A)
#   theta = param[1:(ng*nx)]
#   betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
#   nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
#   deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
#   j = 1
#   beta = list()
#   for (i in 1:ng){
#     beta[[i]] = betatmp[j:sum(nbeta[1:i])]
#     j = sum(nbeta[1:i]) + 1
#   }
#   j = 1
#   nu = list()
#   for (i in 1:ng){
#     nu[[i]] = nutmp[j:sum(nnu[1:i])]
#     j = sum(nnu[1:i]) + 1
#   }
#   ndeltacum = cumsum(c(0, rep(nw, ng)))
#   delta = list()
#   for (i in 1:ng){
#     delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
#   }
#   deltas = c()
#   for (l in 1:nw){
#     a = 0
#     for (i in 1:n){
#       difbkl0 = 0
#       difbkl = 0
#       py0 = 1
#       py = 1
#       ind0 = (Y[i,] == 0)
#       ind = !ind0
#       tmp = exp(sapply(1:period, function(s){
#         sum(nu[[k]]*A[i,s]**(0:(nnu[k]-1)))
#       }))
#       rhoikt = tmp/(1+tmp)
#       lambdaikt = exp(sapply(1:period, function(s){
#         sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
#       }))
#       if (any(ind0)){
#         if (sum(ind0)>1){
#           tind = which(ind0)
#           for (j in 1:length(tind)){
#             difbkl0 = difbkl0 - TCOV[i, tind + (l-1)*period][j]*lambdaikt[tind][j]*(1-rhoikt[tind][j])*exp(-lambdaikt[tind][j])*prod(rhoikt[tind][-j]+(1-rhoikt[tind][-j])*exp(-lambdaikt[tind][-j]))
#           }
#         } else{
#           difbkl0 = - TCOV[i, which(ind0) + (l-1)*period]*lambdaikt[which(ind0)]*(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)])
#         }
#         py0 = prod(rhoikt[which(ind0)]+(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)]))
#       }
#       if (any(ind)){
#         if (sum(ind)>1){
#           tind = which(ind)
#           for (j in 1:length(tind)){
#             difbkl = difbkl + TCOV[i, tind + (l-1)*period][j]*lambdaikt[tind][j]**(Y[i,tind][j])*(1-rhoikt[tind][j])*exp(-lambdaikt[tind][j])*(Y[i,tind][j]-lambdaikt[tind][j])/factorial(Y[i,tind][j])*prod((1-rhoikt[tind][-j])*lambdaikt[tind][-j]**Y[i,tind][-j]*exp(-lambdaikt[tind][-j])/factorial(Y[i,tind][-j]))
#           }
#         } else{
#           difbkl = (1-rhoikt[which(ind)])*TCOV[i, which(ind) + (l-1)*period]*lambdaikt[which(ind)]**(Y[i,which(ind)])*exp(-lambdaikt[which(ind)])*(Y[i,which(ind)]-lambdaikt[which(ind)])/factorial(Y[i,which(ind)])
#         }
#         py = prod((1-rhoikt[which(ind)])*lambdaikt[which(ind)]**Y[i,which(ind)]*exp(-lambdaikt[which(ind)])/factorial(Y[i,which(ind)]))
#       }
#       s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)}))
#       a = a + piik(theta, i, k, ng, X)/s*(difbkl0*py+difbkl*py0)
#     }
#     deltas = c(deltas, a)
#   }
#   return(deltas)
# }
difLdeltakZIP <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
  period = ncol(A)
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  nu = list()
  for (i in 1:ng){
    nu[[i]] = nutmp[j:sum(nnu[1:i])]
    j = sum(nnu[1:i]) + 1
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
      difbkl0 = 0
      difbkl = 0
      py0 = 1
      py = 1
      ind0 = (Y[i,] == 0)
      ind = !ind0
      tmp = exp(sapply(1:period, function(s){
        sum(nu[[k]]*A[i,s]**(0:(nnu[k]-1)))
      }))
      rhoikt = tmp/(1+tmp)
      lambdaikt = exp(sapply(1:period, function(s){
        sum(beta[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
      }))
      if (any(ind0)){
        if (sum(ind0)>1){
          for (t in which(ind0)){
            difbkl0 = difbkl0 - TCOV[i, t + (l-1)*period]*lambdaikt[t]*(1-rhoikt[t])*exp(-lambdaikt[t])*prod(rhoikt[which(ind0)[which(ind0)!=t]]+(1-rhoikt[which(ind0)[which(ind0)!=t]])*exp(-lambdaikt[which(ind0)[which(ind0)!=t]]))
          }
        } else{
          difbkl0 = - TCOV[i, which(ind0) + (l-1)*period]*lambdaikt[which(ind0)]*(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)])
        }
        py0 = prod(rhoikt[which(ind0)]+(1-rhoikt[which(ind0)])*exp(-lambdaikt[which(ind0)]))
      }
      if (any(ind)){
        if (sum(ind)>1){
          for (t in which(ind)){
            difbkl = difbkl + TCOV[i, t + (l-1)*period]*lambdaikt[t]**(Y[i,t])*(1-rhoikt[t])*exp(-lambdaikt[t])*(Y[i,t]-lambdaikt[t])/factorial(Y[i,t])*prod((1-rhoikt[which(ind)[which(ind)!=t]])*lambdaikt[which(ind)[which(ind)!=t]]**Y[i,which(ind)[which(ind)!=t]]*exp(-lambdaikt[which(ind)[which(ind)!=t]])/factorial(Y[i,which(ind)[which(ind)!=t]]))
          }
        } else{
          difbkl = (1-rhoikt[which(ind)])*TCOV[i, which(ind) + (l-1)*period]*lambdaikt[which(ind)]**(Y[i,which(ind)])*exp(-lambdaikt[which(ind)])*(Y[i,which(ind)]-lambdaikt[which(ind)])/factorial(Y[i,which(ind)])
        }
        py = prod((1-rhoikt[which(ind)])*lambdaikt[which(ind)]**Y[i,which(ind)]*exp(-lambdaikt[which(ind)])/factorial(Y[i,which(ind)]))
      }
      s = sum(sapply(1:ng, function(s){piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)}))
      a = a + piik(theta, i, k, ng, X)/s*(difbkl0*py+difbkl*py0)
    }
    deltas = c(deltas, a)
  }
  return(deltas)
}
#################################################################################
# dif likelihood
#################################################################################
difLZIP <- function(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
  out =c()
  for ( k in 1:ng){
    out = c(out, difLthetakZIP(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw))
  }
  for ( k in 1:ng){
    out = c(out, difLbetakZIP(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw))
  }
  for ( k in 1:ng){
    out = c(out, difLnukZIP(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw))
  }
  if (nw != 0){
    for ( k in 1:ng){
      out = c(out, difLdeltakZIP(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw))
    }
  }
  return(out)
}
#################################################################################
# likelihood
#################################################################################
LikelihoodZIP <- function(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
  theta = param[1:(ng*nx)]
  betatmp = param[(ng*nx+1):(ng*nx+sum(nbeta))]
  nutmp = param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
  deltatmp = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  nu = list()
  for (i in 1:ng){
    nu[[i]] = nutmp[j:sum(nnu[1:i])]
    j = sum(nnu[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  a = 0
  for (i in 1:n){
    a = a + log(sum(sapply(1:ng, function(s){
      piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)
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
ftauxZIP <-function(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X){
  j = 1
  betatmp = beta
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  nutmp = nu
  nu = list()
  for (i in 1:ng){
    nu[[i]] = nutmp[j:sum(nnu[1:i])]
    j = sum(nnu[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  deltatmp = delta
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  tmp2 = c()
  if (nx == 1){
    for (k in 1:ng){
      tmp1 = sapply(1:n, function(s){gkZIP(beta, nu, s, k, nbeta, nnu, A, Y, TCOV, delta, nw)})
      tmp2 = cbind(tmp2, tmp1*pi[k])
    }
  }else{
    for (k in 1:ng){
      tmp1 = sapply(1:n, function(s){gkZIP(beta, nu, s, k, nbeta, nnu, A, Y, TCOV, delta, nw)})
      tmp2 = cbind(tmp2, tmp1*piik(pi, i, k, ng, X))
    }
  }
  tmp3 =c()
  for (k in 1:ng){
    tmp3 = cbind(tmp3, 1/(1+(rowSums(tmp2)-tmp2[,k])/tmp2[,k]))
  }
  return(tmp3)
}
fzkSikt <- function(pi, beta, nu, zk, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum){
  if (Y[i,t]>0){
    prob = 0
  }else{
    nuikt = sum(nu[(nnucum[k]+1):(nnucum[k+1])]*A[i,t]**(0:(nnu[k]-1)))
    lambdaikt = exp(
      sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,t]**(0:(nbeta[k]-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum))
    prob = zk[i,k]/(1+exp(-nuikt-lambdaikt))
  }
  return(prob)
}
fSikt <- function(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum){
  if (Y[i,t]>0){
    prob = 0
  }else{
    nuikt = sum(nu[(nnucum[k]+1):(nnucum[k+1])]*A[i,t]**(0:(nnu[k]-1)))
    lambdaikt = exp(sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,t]**(0:(nbeta[k]-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum))
    prob = 1/(1+exp(-nuikt-lambdaikt))
  }
  return(prob)
}
#################################################################################
# Q function for betak, nuk and deltak given a value k.  The goal is to sum all Qbetaknuk after
#################################################################################
QbetakZIP <- function(beta, zk, zkSit, k, nbeta, nnu, n, A, Y, TCOV, delta, nw){
  period = ncol(A)
  a = 0
  for (i in 1:n){
    zik = zk[i,k]
    for (t in 1:period){
      ziksit = zkSit[i,(k-1)*period+t]
      betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
      #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
      a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit))
    }
  }
  return(a)
}
QnukZIP <- function(nu, zk, zkSit, k, nbeta, nnu, n, A, Y){
  period = ncol(A)
  a = 0
  for (i in 1:n){
    zik = zk[i,k]
    for (t in 1:period){
      ziksit = zkSit[i,(k-1)*period+t]
      nuAit = sum(nu*A[i,t]**(0:(nnu-1)))
      a = a + ziksit*nuAit - zik*log(1+exp(nuAit))
    }
  }
  return(a)
}
QdeltakZIP <- function(delta, zk, zkSit, k, nbeta, nnu, n, A, Y, TCOV, beta, nw){
  period = ncol(A)
  a = 0
  for (i in 1:n){
    zik = zk[i,k]
    for (t in 1:period){
      ziksit = zkSit[i,(k-1)*period+t]
      betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + sum(delta*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
      #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
      a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit))
    }
  }
  return(a)
}
#################################################################################
# Differential of betak, nuk and deltak Q function
#################################################################################
difQbetakkZIP <- function(beta, k, zk, zkSit, nbeta, nnu, n, A, Y, TCOV, delta, nw){
  period = ncol(A)
  betas =c()
  for (l in 1:nbeta){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        ziksit = zkSit[i,(k-1)*period+t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
        a = a + A[i,t]**(l-1)*(zik-ziksit)*(Y[i,t]-exp(betaAit))
      }
    }
    betas = c(betas, a)
  }
  return(betas)
}
difQnukkZIP <- function(nu, k, zk, zkSit, nbeta, nnu, n, A, Y){
  period = ncol(A)
  nus =c()
  for (l in 1:nnu){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        ziksit = zkSit[i,(k-1)*period+t]
        nuAit = sum(nu*A[i,t]**(0:(nnu-1)))
        a = a + A[i,t]**(l-1)*(ziksit-zik*exp(nuAit)/(1+exp(nuAit)))
      }
    }
    nus = c(nus, a)
  }
  return(nus)
}
difQdeltakkZIP <- function(delta, k, zk, zkSit, nbeta, nnu, n, A, Y, TCOV, beta, nw){
  period = ncol(A)
  deltas =c()
  for (l in 1:nw){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        ziksit = zkSit[i,(k-1)*period+t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + sum(delta*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
        a = a + TCOV[i, t + (l-1)*period]*(zik-ziksit)*(Y[i,t]-exp(betaAit))
      }
    }
    deltas = c(deltas, a)
  }
  return(deltas)
}
# #################################################################################
# # Q function for beta
# #################################################################################
# QbetaZIP <- function(beta, zk, zkSit, k, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw){
#   sum(sapply(1:ng, function(s){QbetakZIP(beta, zk, zkSit, k, ng, nbeta[s], nnu[s], n, A, Y, TCOV, delta, nw)}))
# }
# #################################################################################
# # Differential of Q beta
# #################################################################################
# difQbetaZIP <- function(beta, k, zk, zkSit, nbeta, nnu, n, A, Y, TCOV, delta, nw, ng){
#   tmp =c()
#   for (k in 1:ng){
#     tmp = c(tmp, difQbetakkZIP(beta, k, zk, zkSit, nbeta[k], nnu[k], n, A, Y, TCOV, delta, nw))
#   }
#   return(tmp)
# }
# #################################################################################
# # Q function for nu
# #################################################################################
# QnuZIP <- function(nu, k, zk, zkSit, nbeta, nnu, n, A, Y){
#   sum(sapply(1:ng, function(s){QnukZIP(nu, zk, zkSit, k, ng, nbeta[s], nnu[s], n, A, Y)}))
# }
# #################################################################################
# # Differential of Q beta
# #################################################################################
# difQnuZIP <- function(nu, k, zk, zkSit, nbeta, nnu, n, A, Y){
#   tmp =c()
#   for (k in 1:ng){
#     tmp = c(tmp, difQnukkZIP(nu, k, zk, zkSit, nbeta[k], nnu[k], n, A, Y))
#   }
#   return(tmp)
# }
#################################################################################
# Likelihood calculation for EM algorithm
#################################################################################
likelihoodEMZIP <- function(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw, nu, nnu){
  likeli = 0
  betatmp = beta
  deltatmp = delta
  j = 1
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  nutmp = nu
  nu = list()
  for (i in 1:ng){
    nu[[i]] = nutmp[j:sum(nnu[1:i])]
    j = sum(nnu[1:i]) + 1
  }
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  delta = list()
  for (i in 1:ng){
    delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
  }
  for (i in 1:n){
    likeli = likeli + log(sum(sapply(1:ng, function(s){
      pi[s]*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)
    }))
    )
  }
  return(likeli)
}

WitEMtmp <- function(TCOV, period, delta, nw, i, t, k){
  if (nw == 0){
    return(0)
  }else{
    return(sum(delta*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)]))
  }
}
difQbetadeltakZIP <- function(betadelta, k, zk, Sikt, nbeta, nnu, n, A, Y, TCOV, nw){
  period = ncol(A)
  beta = betadelta[1:nbeta]
  delta=betadelta[-c(1:nbeta)]
  betadeltas =c()
  for (l in 1:nbeta){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        sikt = Sikt[i,t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
             a = a + A[i,t]**(l-1)*zik*(1-sikt)*(Y[i,t]-exp(betaAit))
      }
    }
    betadeltas = c(betadeltas, a)
  }
  for (l in 1:nw){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        sikt = Sikt[i,t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + sum(delta*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
        a = a + TCOV[i, t + (l-1)*period]*zik*(1-sikt)*(Y[i,t]-exp(betaAit))
      }
    }
    betadeltas = c(betadeltas, a)
  }
  return(betadeltas)
}
QbetadeltakZIP <- function(betadelta, k, zk, Sikt, nbeta, nnu, n, A, Y, TCOV, nw){
  period = ncol(A)
  beta = betadelta[1:nbeta]
  delta=betadelta[-c(1:nbeta)]
  a = 0
  for (i in 1:n){
    zik = zk[i,k]
    for (t in 1:period){
      sikt = Sikt[i,t]
      betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
      #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
      a = a + zik*(1-sikt)*(Y[i,t]*betaAit-exp(betaAit))
    }
  }
  return(a)
}
QnukZIPtmp <- function(nu, zk, Sikt, k, nbeta, nnu, n, A, Y){
  period = ncol(A)
  a = 0
  for (i in 1:n){
    zik = zk[i,k]
    for (t in 1:period){
      sikt = Sikt[i,t]
      nuAit = sum(nu*A[i,t]**(0:(nnu-1)))
      a = a + zik*(sikt*nuAit - log(1+exp(nuAit)))
    }
  }
  return(a)
}
difQnukkZIPtmp <- function(nu, k, zk, Sikt, nbeta, nnu, n, A, Y){
  period = ncol(A)
  nus =c()
  for (l in 1:nnu){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        sikt = Sikt[i,t]
        nuAit = sum(nu*A[i,t]**(0:(nnu-1)))
        a = a + A[i,t]**(l-1)*zik*(sikt-exp(nuAit)/(1+exp(nuAit)))
      }
    }
    nus = c(nus, a)
  }
  return(nus)
}
QbetakZIPtmp <- function(beta, zk, Sikt, k, nbeta, nnu, n, A, Y, TCOV, delta, nw){
  period = ncol(A)
  a = 0
  for (i in 1:n){
    zik = zk[i,k]
    for (t in 1:period){
      sikt = Sikt[i,t]
      betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
      #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
      a = a + zik*(1-sikt)*(Y[i,t]*betaAit-exp(betaAit))
    }
  }
  return(a)
}
difQbetakkZIPtmp <- function(beta, k, zk, Sikt, nbeta, nnu, n, A, Y, TCOV, delta, nw){
  period = ncol(A)
  betas =c()
  for (l in 1:nbeta){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        sikt = Sikt[i,t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
        a = a + A[i,t]**(l-1)*zik*(1-sikt)*(Y[i,t]-exp(betaAit))
      }
    }
    betas = c(betas, a)
  }
  return(betas)
}
QdeltakZIPtmp <- function(delta, zk, Sikt, k, nbeta, nnu, n, A, Y, TCOV, beta, nw){
  period = ncol(A)
  a = 0
  for (i in 1:n){
    zik = zk[i,k]
    for (t in 1:period){
      sikt = Sikt[i,t]
      betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
      #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
      a = a + zik*(1-sikt)*(Y[i,t]*betaAit-exp(betaAit))
    }
  }
  return(a)
}
difQdeltakkZIPtmp <- function(delta, k, zk, Sikt, nbeta, nnu, n, A, Y, TCOV, beta, nw){
  period = ncol(A)
  deltas =c()
  for (l in 1:nw){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        sikt = Sikt[i,t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) +WitEMtmp(TCOV, period, delta, nw, i, t, k)
        a = a + TCOV[i, t + (l-1)*period]*zik*(1-sikt)*(Y[i,t]-exp(betaAit))
      }
    }
    deltas = c(deltas, a)
  }
  return(deltas)
}
#################################################################################
# EM algorithm
#################################################################################
EMZIP <- function(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr){
  period = ncol(A)
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    nu =  param[(ng+sum(nbeta)):(ng+sum(nbeta)+sum(nnu)-1)]
    delta = param[-c(1:(ng+sum(nbeta)+sum(nnu)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    nu =  param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
    delta = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
    pi[c(((refgr-1)*nx+1):(nx*refgr))] = 0
  }
  nbetacum = cumsum(c(0, nbeta))
  nnucum = cumsum(c(0, nnu))
  ndeltacum = cumsum(c(0, rep(nw, ng)))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEMZIP(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw, nu, nnu), "\n"))
    }else{
      cat(paste(-LikelihoodZIP(c(pi, beta, nu, delta), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw), "\n"))
    }
    # E-step
    zk = ftauxZIP(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X)
    zkSit = c()
    for (k in 1:ng){
      tmp3 =c()
      for (i in 1:n){
        tmp4 = c()
        for (t in 1:period){
          tmp4 = c(tmp4, fzkSikt(pi, beta, nu, zk, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum))
        }
        tmp3 = rbind(tmp3, tmp4)
      }
      zkSit = cbind(zkSit, tmp3)
    }
    #M-step
    newbeta = c()
    newnu = c()
    newdelta = c()
    for (k in 1:ng){
      newbeta = c(newbeta, optim(beta[(nbetacum[k]+1):(nbetacum[k+1])], QbetakZIP, difQbetakkZIP, method = "BFGS",
                                 zk = zk, zkSit = zkSit, k = k,
                                 nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
                                 TCOV = TCOV, delta = delta[(ndeltacum[k]+1):(ndeltacum[k+1])], nw = nw, 
                                 control = list(fnscale=-1))$par)
      newnu = c(newnu, optim(nu[(nnucum[k]+1):(nnucum[k+1])], QnukZIP, difQnukkZIP, method = "BFGS",
                             zk = zk, zkSit = zkSit, k = k,
                             nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
                             control = list(fnscale=-1))$par)
      if (nw != 0){
        newdelta = c(newdelta, optim(delta[(ndeltacum[k]+1):(ndeltacum[k+1])], QdeltakZIP, difQdeltakkZIP, method = "BFGS",
                                     zk = zk, zkSit = zkSit, k = k,
                                     nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
                                     TCOV = TCOV, beta = beta[(nbetacum[k]+1):(nbetacum[k+1])], nw = nw,
                                     control = list(fnscale=-1))$par)
      }else{
        newdelta = c()
      }
    }
    if (nx == 1){
      newpi = colSums(zk)/n
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS, refgr)
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newnu, newdelta)-c(beta,nu, delta))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    nu = newnu
    delta = newdelta
    pi = newpi
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, nu, delta)
  }else{
    param = c(pi, beta, nu, delta)
  }
  return(param)
}
#################################################################################
# EM algorithm IRLS
#################################################################################
EMZIPIRLS <- function(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr){
  period = ncol(A)
  if (nx == 1){
    pi = param[1:(ng-1)]
    beta = param[(ng):(ng+sum(nbeta)-1)]
    nu =  param[(ng+sum(nbeta)):(ng+sum(nbeta)+sum(nnu)-1)]
    delta = param[-c(1:(ng+sum(nbeta)+sum(nnu)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    pi = param[1:(ng*nx)]
    beta = param[(ng*nx+1):(ng*nx+sum(nbeta))]
    nu =  param[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
    delta = param[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
    pi[c(((refgr-1)*nx+1):(nx*refgr))] = 0
  }
  nbetacum = cumsum(c(0, nbeta))
  nnucum = cumsum(c(0, nnu))
  tour = 1
  while (tour<itermax){
    ###########################
    # print likelihood for every loop
    ###########################
    if (nx == 1){
      cat(paste(-likelihoodEMZIP(n, ng, nbeta, beta, pi, A, Y, TCOV, delta, nw, nu, nnu), "\n"))
    }else{
      cat(paste(-LikelihoodZIP(c(pi, beta, nu, delta), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw), "\n"))
    }
    # E-step
    zk =  ftauxZIP(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X)
    # M-step
    newbeta = c()
    newnu = c()
    newdelta = c()
    for (k in 1:ng){
      newbetaIRLS = c()
      betaIRLS = beta[(nbetacum[k]+1):nbetacum[k+1]]
      newnuIRLS = c()
      nuIRLS = nu[(nnucum[k]+1):nnucum[k+1]]
      newdeltaIRLS = c()
      deltaIRLS = delta[((k-1)*nw+1):(k*nw)]
      precIRLS = 1
      while(any(abs(precIRLS)>10**(-6))){
        tmp = c()
        tmpr = c()
        tmp1 = c()
        tmp2 = c()
        tmp.p = c()
        tmpr.p =c()
        tmp1.p =c()
        tmp2.p = c()
        tmp.w = c()
        tmp1.w = c()
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            #nuIRLS
            tmp = c(tmp, A[i,t]**(0:(nnu[k]-1)))
            nuAit = sum(nuIRLS*A[i,t]**(0:(nnu[k]-1)))
            rhoikt = exp(nuAit)/(1+exp(nuAit))
            tmpr = c(tmpr, rhoikt*(1-rhoikt))
            tmp1 = c(tmp1, nuAit + (sikt-rhoikt)/(rhoikt*(1-rhoikt)))
            #betaIRLS
            tmp.p = c(tmp.p, A[i,t]**(0:(nbeta[k]-1)))
            betaAit = sum(betaIRLS*A[i,t]**(0:(nbeta[k]-1)))
            deltaWit = WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
            lambdaikt = exp(betaAit + deltaWit)
            tmpr.p = c(tmpr.p, (1-sikt))
            tmp1.p = c(tmp1.p, betaAit+Y[i,t]/lambdaikt-1)
            tmp2.p = c(tmp2.p, lambdaikt)
            #deltaIRLS
            tmp.w = c(tmp.w, TCOV[i, t + 0:((nw-1)*period)])
            tmp1.w = c(tmp1.w, deltaWit+Y[i,t]/lambdaikt-1)
          }
          tmp2 = c(tmp2, rep(zk[i,k], period))
        }
        #nuIRLS
        Aw = matrix(tmp, nrow = nnu[k], byrow = FALSE)
        W = diag(tmpr)
        Z = diag(tmp2)
        S = matrix(tmp1, ncol =1)
        QR = qr(t(Aw%*%sqrt(Z*W)))
        Q = qr.Q(QR)
        R = qr.R(QR)
        newnuIRLS = backsolve(R,t(Q)%*%(sqrt(Z*W)%*%S))
        #betaIRLS
        Awp = matrix(tmp.p, nrow = nbeta[k], byrow = FALSE)
        Wp = diag(tmp2.p)
        Zs = diag(tmp2*tmpr.p)# for Sp%*%Zp in  the defintion 7.6.2
        Sp = matrix(tmp1.p, ncol =1)
        QR = qr(t(Awp%*%sqrt(Zs*Wp)))
        Q = qr.Q(QR)
        R = qr.R(QR)
        newbetaIRLS = backsolve(R,t(Q)%*%(sqrt(Zs*Wp)%*%Sp))
        #deltaIRLS
        if (nw != 0){
          Ww = matrix(tmp.w, nrow = nw, byrow = FALSE)
          Sw = matrix(tmp1.w, ncol =1)
          QR = qr(t(Ww%*%sqrt(Zs*Wp)))
          Q = qr.Q(QR)
          R = qr.R(QR)
          newdeltaIRLS = backsolve(R,t(Q)%*%(sqrt(Zs*Wp)%*%Sw))
        }
        precIRLS = c(betaIRLS-newbetaIRLS, nuIRLS-newnuIRLS, deltaIRLS-newdeltaIRLS)
        betaIRLS = newbetaIRLS
        nuIRLS = newnuIRLS
        deltaIRLS = newdeltaIRLS
      }
      newbeta = c(newbeta, betaIRLS)
      newnu = c(newnu, nuIRLS)
      newdelta = c(newdelta, deltaIRLS)
    }
    if (nx == 1){
      newpi = colSums(zk)/n
    }else{
      newpi = findtheta(pi, taux, X, n, ng, nx, period, EMIRLS, refgr)
    }
    tour = tour + 1
    if (all(abs(c(newbeta, newnu, newdelta)-c(beta,nu,delta))<10**(-6))){
      tour = itermax + 2
    }
    beta = newbeta
    nu = newnu
    delta = newdelta
    pi = newpi
  }
  if (nx == 1){
    param = c(pi[1:(ng-1)], beta, nu, delta)
  }else{
    param = c(pi, beta, nu, delta)
  }
  return(param)
}

#################################################################################
# Compute of matrix information
#################################################################################
#################################################################################
# General functions
#################################################################################
# Function lambda and nu
#################################################################################
lambdaikt <- function(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw){
  return(exp(sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,t]**(0:(nbeta[k]-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)))
}
rhoikt <- function(k, i, t, nnu, nnucum, A, nu){
  tmp = exp(sum(nu[(nnucum[k]+1):(nnucum[k+1])]*A[i,t]**(0:(nnu[k]-1))))
  return(tmp/(1+tmp))
}
#################################################################################
# matrix differential of beta**2 for t and i
#################################################################################
mbetaZIP = function(i, t, ng, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu){
  res = matrix(rep(0,(sum(nbeta))**2), ncol = sum(nbeta))
  for (k in 1:ng){
    tmp = c()
    sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
    for (l in 1:nbeta[k]){
      for (lp in 1:nbeta[k]){
        tmp = c(tmp, -taux[i,k]*(1-sikt)*A[i,t]**(l+lp-2)*lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
      }
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), (nbetacum[k]+1):(nbetacum[k+1])] = matrix(tmp, ncol = nbeta[k])
  }
  return(res)
}
##################################################################################
# matrix differential of  beta delta
#################################################################################
mbetadeltaZIP = function(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu){
  res = matrix(rep(0,ng*nw*sum(nbeta)), ncol = nw*ng)
  for (k in 1:ng){
    tmp = c()
    sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
    for (lp in 1:nbeta[k]){
      for (l in 1:nw){
        tmp = c(tmp, -taux[i,k]*(1-sikt)*TCOV[i, t + (l-1)*period]*A[i,t]**(lp-1)*lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
      }
    }
    res[(nbetacum[k]+1):(nbetacum[k+1]), (ndeltacum[k]+1):(ndeltacum[k+1])] = matrix(tmp, ncol = nw, byrow = TRUE)
  }
  return(res)
}
#################################################################################
# matrix differential of nu**2 for t and i
#################################################################################
mnuZIP = function(i, t, ng, nnu, nnucum, A, taux, nu){
  res = matrix(rep(0,(sum(nnu))**2), ncol = sum(nnu))
  for (k in 1:ng){
    tmp = c()
    for (l in 1:nnu[k]){
      for (lp in 1:nnu[k]){
        tmp = c(tmp, -taux[i,k]*A[i,t]**(l+lp-2)*rhoikt(k, i, t, nnu, nnucum, A, nu)*(1-rhoikt(k, i, t, nnu, nnucum, A, nu)))
      }
    }
    res[(nnucum[k]+1):(nnucum[k+1]), (nnucum[k]+1):(nnucum[k+1])] = matrix(tmp, ncol = nnu[k])
  }
  return(res)
}
#################################################################################
# matrix differential of delta**2 for t and i
#################################################################################
mdeltaZIP = function(i, t, ng, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu){
  res = matrix(rep(0,(sum(nw*ng))**2), ncol = sum(nw*ng))
  for (k in 1:ng){
    tmp = c()
    sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
    for (l in 1:nw){
      for (lp in 1:nw){
        tmp = c(tmp, -taux[i,k]*(1-sikt)*TCOV[i, t + (l+lp-2)*period]*lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
      }
    }
    res[(ndeltacum[k]+1):(ndeltacum[k+1]), (ndeltacum[k]+1):(ndeltacum[k+1])] = matrix(tmp, ncol = nw)
  }
  return(res)
}
##################################################################################
# defintion of function Bikl and Sk
##################################################################################
BPiiklZIP <- function(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum){
  tmp = sapply(1:period, function(s){
    A[i,s]**(l-1)*(Y[i,s]-lambdaikt(k, i, s, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))*(fSikt(pi, beta, nu, k, i, s, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)-1)
  })
  return(sum(tmp))
}
DPiiklZIP <- function(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum){
  tmp = sapply(1:period, function(s){
    TCOV[i, s + (l-1)*period]*(Y[i,s]-lambdaikt(k, i, s, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))*(fSikt(pi, beta, nu, k, i, s, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)-1)
  })
  return(sum(tmp))
}
NiklZIP <- function(pi, beta, nu, k, i, l, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum){
  tmp = sapply(1:period, function(s){
    A[i,s]**(l-1)*(fSikt(pi, beta, nu, k, i, s, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)-rhoikt(k, i, s, nnu, nnucum, A, nu))
  })
  return(sum(tmp))
}
# matrix cov pi betak
covPiBetakZIP <- function(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovPiBetak =  matrix(rep(0,(ng-1)*nbeta[k]), ncol=nbeta[k])
  for  (kp in 1:(ng-1)){
    for (l in 1:nbeta[k]){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + BPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + BPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,k]*(taux[i,kp]/pi[kp]-taux[i,ng]/pi[ng])
        }
      }
      rcovPiBetak[kp,l] =  tmp
    }
  }
  return(rcovPiBetak)
}
# matrix cov pi nuk
covPiNukZIP <- function(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu){
  rcovPiNuk =  matrix(rep(0,(ng-1)*nnu[k]), ncol=nnu[k])
  for  (kp in 1:(ng-1)){
    for (l in 1:nnu[k]){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + NiklZIP(pi, beta, nu, k, i, l, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + NiklZIP(pi, beta, nu, k, i, l, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      rcovPiNuk[kp,l] =  tmp
    }
  }
  return(rcovPiNuk)
}
# matrix cov pi deltak
covPiDeltakZIP <- function(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovPiDeltak =  matrix(rep(0,(ng-1)*nw), ncol=nw)
  for  (kp in 1:(ng-1)){
    for (l in 1:nw){
      tmp = 0
      if (kp == k){
        for (i in 1:n){
          tmp = tmp + DPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      else{
        for (i in 1:n){
          tmp = tmp + DPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
        }
      }
      rcovPiDeltak[kp,l] =  tmp
    }
  }
  return(rcovPiDeltak)
}
# matrix cov Betak Betal
covBetakBetalZIP <- function(k,l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovBetakBetal = matrix(rep(0, (nbeta[k]*nbeta[l])), nrow = nbeta[k])
  if (k==l){
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pikpt = A[i,t]**(p-1)*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            pikqt = A[i,t]**(q-1)*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + pikpt*pikqt*taux[i,k]*((1-taux[i,k])*(1-2*sikt)-sikt*(1-taux[i, k]*sikt))
          }
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nbeta[k]){
      for (q in 1:nbeta[l]){
        tmp = 0
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pikpt = A[i,t]**(p-1)*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            pilqt = A[i,t]**(q-1)*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + pikpt*pilqt*taux[i,k]*taux[i,l]*(sikt-1)*(1-silt)
          }
        }
        rcovBetakBetal[p,q] = tmp
      }
    }
  }
  return(rcovBetakBetal)
}
# matrix cov betak nul
covBetakNulZIP <- function(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovBetakNul = matrix(rep(0, nbeta[k]*nnu[l]), ncol = nnu[l])
  for (p in 1:nbeta[k]){
    for (q in 1 : nnu[l]){
      tmp = 0
      if (k == l){
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pilqt = A[i,t]**(q-1)*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + A[i,t]**(p-1)*pilqt*taux[i,k]*(sikt-1)*(taux[i,k]*sikt+rhoikt(k, i, t, nnu, nnucum, A, nu)*(1-taux[i, k]))
          }
        }
      }else{
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pilpt = A[i,t]**(p-1)*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + A[i,t]**(q-1)*pilpt*taux[i,k]*taux[i,l]*(silt-1)*(sikt-rhoikt(k, i, t, nnu, nnucum, A, nu))
          }
        }
      }
      rcovBetakNul[p,q] = tmp
    }
  }

  return(rcovBetakNul)
}
# matrix cov betak deltal
covBetakDeltaZIP <- function(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovBetakDeltal = matrix(rep(0, nbeta[k]*nw), ncol = nw)
  for (p in 1:nbeta[k]){
    for (q in 1 :nw){
      tmp = 0
      if (k == l){
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pikpt = A[i,t]**(p-1)*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            pikqt =  TCOV[i, t + (q-1)*period]*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + pikpt*pikqt*taux[i,k]*((1-taux[i,k])*(1-2*sikt)-sikt*(1-taux[i, k]*sikt))
          }
        }
      }else{
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pikpt = A[i,t]**(p-1)*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            pilqt =  TCOV[i, t + (q-1)*period]*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + pikpt*pilqt*taux[i,k]*taux[i,l]*(sikt-1)*(1-silt)
          }
        }
      }
      rcovBetakDeltal[p,q] = tmp
    }
  }
  return(rcovBetakDeltal)
}
# matrix cov nuk nul
covNukNulZIP <- function(k,l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovNukNul = matrix(rep(0, (nnu[k]*nnu[l])), nrow = nnu[k])
  if (k==l){
    for (p in 1:nnu[k]){
      for (q in 1:nnu[l]){
        tmp = 0
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            roikt = rhoikt(k, i, t, nnu, nnucum, A, nu)
            tmp = tmp + A[i,t]**(p-1)*A[i,t]**(q-1)*taux[i,k]*(sikt*(1-taux[i,k]*sikt)+roikt**2*(1-taux[i,k])+2*roikt*sikt*(1-taux[i,k]))
          }
        }
        rcovNukNul[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nnu[k]){
      for (q in 1:nnu[l]){
        tmp = 0
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            tmp = tmp +  A[i,t]**(p-1)* A[i,t]**(q-1)*taux[i,k]*taux[i,l]*(-sikt*silt+rhoikt(l, i, t, nnu, nnucum, A, nu)*sikt+rhoikt(k, i, t, nnu, nnucum, A, nu)*silt-rhoikt(k, i, t, nnu, nnucum, A, nu)*rhoikt(l, i, t, nnu, nnucum, A, nu))
          }
        }
        rcovNukNul[p,q] = tmp
      }
    }
  }
  return(rcovNukNul)
}
# matrix cov deltak nul
covDeltakNulZIP <- function(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovDeltakNul = matrix(rep(0, nw*nnu[l]), ncol = nnu[l])
  for (p in 1:nw){
    for (q in 1 : nnu[l]){
      tmp = 0
      if (k == l){
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pilpt = TCOV[i, t + (p-1)*period]*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + A[i,t]**(p-1)*pilpt*taux[i,k]*(sikt-1)*(taux[i,k]*sikt+rhoikt(k, i, t, nnu, nnucum, A, nu)*(1-taux[i, k]))
          }
        }
      }else{
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pilpt =TCOV[i, t + (p-1)*period]*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + A[i,t]**(q-1)*pilpt*taux[i,k]*taux[i,l]*(silt-1)*(sikt-rhoikt(k, i, t, nnu, nnucum, A, nu))
          }
        }
      }
      rcovDeltakNul[p,q] = tmp
    }
  }
  return(rcovDeltakNul)
}
# matrix cov Deltak Deltal
covDeltakDeltalZIP <- function(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
  rcovDeltakDeltal = matrix(rep(0, (nw**2)), nrow = nw)
  if (k==l){
    for (p in 1:nw){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pikpt = TCOV[i, t + (p-1)*period]*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            pikqt =  TCOV[i, t + (q-1)*period]*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + pikpt*pikqt*taux[i,k]*((1-taux[i,k])*(1-2*sikt)-sikt*(1-taux[i, k]*sikt))
          }
        }
        rcovDeltakDeltal[p,q] = tmp
      }
    }
  }else{
    for (p in 1:nw){
      for (q in 1:nw){
        tmp = 0
        for (i in 1:n){
          for (t in 1:period){
            sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
            pikpt = TCOV[i, t + (p-1)*period]*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            pilqt =  TCOV[i, t + (q-1)*period]*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
            tmp = tmp + pikpt*pilqt*taux[i,k]*taux[i,l]*(sikt-1)*(1-silt)
          }
        }
        rcovDeltakDeltal[p,q] = tmp
      }
    }
  }
  return(rcovDeltakDeltal)
}
#################################################################################
#################################################################################
IEMZIP <- function(paramEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, delta, nw, refgr){
  # pi = c(paramEM[1:(ng*nx-1)], 1-sum(paramEM[1:(ng*nx-1)]))
  # beta = paramEM[(ng*nx):(ng*nx+sum(nbeta)-1)]
  # pi = c(paramEM[1:(ng-1)], 1-sum(paramEM[1:(ng-1)]))
  # beta = paramEM[(ng):(ng+sum(nbeta)-1)]
  # delta = paramEM[-c(1:(ng+sum(nbeta)-1))]
  if (nx == 1){
    pi = paramEM[1:(ng-1)]
    beta = paramEM[(ng):(ng+sum(nbeta)-1)]
    nu =  paramEM[(ng+sum(nbeta)):(ng+sum(nbeta)+sum(nnu)-1)]
    delta = paramEM[-c(1:(ng+sum(nbeta)+sum(nnu)-1))]
    pi = c(pi, 1-sum(pi))
  }else{
    theta = paramEM[1:(ng*nx)]
    beta = paramEM[(ng*nx+1):(ng*nx+sum(nbeta))]
    nu =  paramEM[(ng*nx+sum(nbeta)+1):(ng*nx+sum(nbeta)+sum(nnu))]
    delta = paramEM[-(1:(ng*nx+sum(nbeta)+sum(nnu)))]
    theta[c(((refgr-1)*nx+1):(nx*refgr))] = 0
    pi = theta
  }
  taux =  ftauxZIP(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X)
  period = ncol(A)
  nbetacum = cumsum(c(0, nbeta))
  nnucum = cumsum(c(0, nnu))
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
  tmp1 = matrix(rep(0, (sum(nw*ng))**2), ncol =  sum(nw*ng)) # we take off the dimension of mPi
  tmp2 = matrix(rep(0, (sum(nnu))**2), ncol =  sum(nnu)) # we take off the dimension of mPi
  tmp3 = matrix(rep(0, sum(nbeta)*nw*ng), ncol =  nw*ng) # we take off the dimension of mPi
  for (i in 1:n){
    for (t in (1:period)){
      B = B + mbetaZIP(i,t, ng, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu)
      tmp1 = tmp1 + mdeltaZIP(i,t, ng, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu)
      tmp2 = tmp2 + mnuZIP(i, t, ng, nnu, nnucum, A, taux, nu)
      tmp3 = tmp3 +mbetadeltaZIP(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu)
    }
  }
  B = cbind(matrix(rep(0, (ng-1)*nx*(sum(nbeta))), ncol = (ng-1)*nx),
            B,
            matrix(rep(0, (sum(nbeta)*sum(nnu))), ncol = sum(nnu)))
  B = rbind(cbind(mPiL, matrix(rep(0,(sum(nbeta)+sum(nnu))*(ng-1)*nx), nrow = (ng-1)*nx)),
            B)
  B = rbind(B, cbind(matrix(rep(0,(sum(nbeta)+(ng-1)*nx)*sum(nnu)), nrow = sum(nnu)),
                     tmp2))
  if (nw != 0){
    B = rbind(B, matrix(rep(0,(sum(nbeta)+sum(nnu)+ng-1)*ng*nw), nrow = ng*nw))
    B = cbind(B, rbind(matrix(rep(0,(sum(nbeta)+sum(nnu)+ng-1)*ng*nw), ncol = ng*nw),
                       tmp1))
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
            tmp = tmp + taux[i,k]*(1-taux[i,k])/pi[k]**2+taux[i,ng]*(1-taux[i,ng])/pi[ng]**2-2*taux[i,k]*taux[i,ng]/(pi[k]*pi[ng])
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
      covPiBetaL = cbind(covPiBetaL, covPiBetakZIP(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
    }
    for (kp in 1:(ng-1)){
      rcovPiBetakL =  matrix(rep(0,(ng-1)*nbeta[ng]), ncol=nbeta[ng])
      for (l in 1:nbeta[ng]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp + BPiiklZIP(pi, beta, nu, ng, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
        }
        rcovPiBetakL[kp,l] =  tmp
      }
     
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
            if (k == l){
              for (i in 1:n){
                for (t in 1:period){
                  sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
                  pikqt = A[i,t]**(q-1)*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
                  tmp = tmp + X[i,p]*pikqt*taux[i,k]*taux[i,k]*(1-taux[i,k])*(1-sikt)
                }
              }
            }else{
              for (i in 1:n){
                for (t in 1:period){
                  silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
                  pilqt = A[i,t]**(q-1)*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
                  tmp = tmp + X[i,p]*pilqt*taux[i,k]*taux[i,l]*(1+silt)
                }
              }
            }
            covPiBetatmp[p, q] = tmp
          }
        }
        covPiBetaL[((k-1)*nx+1):((k-1)*nx+nx), (nbetacum[l]+1):(nbetacum[l+1])] = covPiBetatmp
      }
    }
    #covPiBetatmp = matrix(rep(0, nx*nbeta[l]), ncol =nbeta[l])
  }
  ##################################################################################
  # matrix cov pi nu
  ##################################################################################
  # matrix cov pi nu
  if (nx == 1){
    covPiNuL = c()
    for (k in 1:(ng-1)){
      covPiNuL = cbind(covPiNuL, covPiNukZIP(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu))
    }
    for (kp in 1:(ng-1)){
      rcovPiNukL =  matrix(rep(0,(ng-1)*nnu[ng]), ncol=nnu[ng])
      for (l in 1:nnu[ng]){
        tmp = 0
        for (i in 1:n){
          tmp = tmp - NiklZIP(pi, beta, nu, ng, i, l, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
        }
        rcovPiNukL[kp,l] =  tmp
      }

    }
    covPiNuL = cbind(covPiNuL, rcovPiNukL)
  }else{
    covPiNuL = matrix(rep(0,((ng-1)*nx)*sum(nnu)), ncol=sum(nnu))
    for (k in 1:(ng-1)){
      for (l in 1:ng){
        covPiNutmp = matrix(rep(0, nx*nnu[l]), ncol =nnu[l])
        for (p in 1:nx){
          for (q in 1:nnu[l]){
            tmp = 0
            if (k==l){
              for (i in 1:n){
                for (t in 1:period){
                sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
                tmp = tmp + X[i,p]*A[i,t]**(q-1)*taux[i,k]*(1-taux[i,k])*(sikt-rhoikt(k, i, t, nnu, nnucum, A, nu))
                }
              }
            }else{
              for (i in 1:n){
                for (t in 1:period){
                silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
                tmp = tmp + X[i,p]*A[i,t]**(q-1)*taux[i,k]*taux[i,l]*(rhoikt(l, i, t, nnu, nnucum, A, nu)-silt)
              }}
            }
            covPiNutmp[p, q] = tmp
          }
        }
        covPiNuL[((k-1)*nx+1):((k-1)*nx+nx), (nnucum[l]+1):(nnucum[l+1])] = covPiNutmp
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
        covPiDeltaL = cbind(covPiDeltaL, covPiDeltakZIP(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
      }
      for (kp in 1:(ng-1)){
        rcovPiDeltakL =  matrix(rep(0,(ng-1)*nw), ncol=nw)
        for (l in 1:nw){
          tmp = 0
          for (i in 1:n){
            tmp = tmp + DPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,ng]*((1-taux[i,ng])/pi[ng]+taux[i,kp]/pi[kp])
          }
        }
        rcovPiDeltakL[kp,l] =  tmp
      }
      covPiDeltaL = cbind(covPiDeltaL, rcovPiDeltakL)
    }else{
      covPiDeltaL = matrix(rep(0,((ng-1)*nx)*sum(ndelta)), ncol=sum(ndelta))
      for (k in 1:(ng-1)){
        for (l in 1:ng){
          covPiDeltatmp = matrix(rep(0, nx*ndelta[l]), ncol =ndelta[l])
          for (p in 1:nx){
            for (q in 1:ndelta[l]){
              tmp = 0
              if (k == l){
                for (i in 1:n){
                  for (t in 1:period){
                    sikt = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
                    pikqt = TCOV[i, t + (q-1)*period]*(Y[i,t]-lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
                    tmp = tmp + X[i,p]*pikqt*taux[i,k]*taux[i,k]*(1-taux[i,k])*(1-sikt)
                  }
                }
              }else{
                for (i in 1:n){
                  for (t in 1:period){
                    silt = fSikt(pi, beta, nu, l, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
                    pilqt = TCOV[i, t + (q-1)*period]*(Y[i,t]-lambdaikt(l, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw))
                    tmp = tmp + X[i,p]*pilqt*taux[i,k]*taux[i,l]*(1+silt)
                  }
                }
              }
              covPiDeltatmp[p, q] = tmp
            }
          }
          covPiDeltaL[((k-1)*nx+1):((k-1)*nx+nx), (ndeltacum[l]+1):(ndeltacum[l+1])] = covPiDeltatmp
        }
      }
      covPiDeltatmp = matrix(rep(0, nx*ndelta[l]), ncol =ndelta[l])
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
      covBetakL= rbind(covBetakL, covBetakBetalZIP(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
    }
    covBetaL = cbind(covBetaL, covBetakL)
  }
  ##################################################################################
  # matrix cov beta nu
  ##################################################################################
  # matrix cov Beta Nu
  covBetaNuL = c()
  for (l in 1:ng){
    covNukL = c()
    for (k in 1:ng){
      covNukL= rbind(covNukL, covBetakNulZIP(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
    }
    covBetaNuL = cbind(covBetaNuL, covNukL)
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
        covBetaDeltakL= rbind(covBetaDeltakL, covBetakDeltaZIP(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
      }
      covBetaDeltaL = cbind(covBetaDeltaL, covBetaDeltakL)
    }
  }
  ##################################################################################
  # matrix cov Nu
  ##################################################################################
  # matrix cov Nu
  covNuL = c()
  for (l in 1:ng){
    covNukL = c()
    for (k in 1:ng){
      covNukL= rbind(covNukL, covNukNulZIP(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
    }
    covNuL = cbind(covNuL, covNukL)
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
        covDeltakL= rbind(covDeltakL, covDeltakDeltalZIP(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
      }
      covDeltaL = cbind(covDeltaL, covDeltakL)
    }
    # matrix cov Delta Nu
    covDeltaNu = c()
    for (l in 1:ng){
      covDeltaNuk = c()
      for (k in 1:ng){
        covDeltaNuk= rbind(covDeltaNuk, covDeltakNulZIP(k, l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu))
      }
      covDeltaNu = cbind(covDeltaNu, covDeltaNuk)
    }
  }
  ##################################################################################
  # Matrix of covariance of score function
  ##################################################################################
  cov =c()
  if (nw !=0){
    cov = cbind(covPiL, covPiBetaL, covPiNuL, covPiDeltaL)
    cov = rbind(cov, cbind(t(covPiBetaL), covBetaL, covBetaNuL, covBetaDeltaL))
    cov = rbind(cov, cbind(t(covPiNuL), t(covBetaNuL), covNuL, t(covDeltaNu)))
    cov = rbind(cov, cbind(t(covPiDeltaL), t(covBetaDeltaL), covDeltaNu, covDeltaL))
  }else{
    cov = cbind(covPiL, covPiBetaL, covPiNuL)
    cov = rbind(cov, cbind(t(covPiBetaL), covBetaL, covBetaNuL))
    cov = rbind(cov, cbind(t(covPiNuL), t(covBetaNuL), covNuL))
  }


  ##################################################################################
  # Information matrix of Fisher
  ##################################################################################
  IEM = - B - cov
  return(sqrt(diag(solve(IEM))))
}

