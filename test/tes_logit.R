library(trajeR)
load(file="/mnt/Travail/These/R/Package_cpp/trajeR/data/dataLOGIT01.RData")
load(file = "/mnt/Travail/These/R/Package/trajeR/tests/LOGIT_data01_sol")


beta = solL$beta
theta=solL$theta
delta=solLTCOV$delta
pi=solEM$theta
ng=3
nbeta = c(0,3,4)+1
Y = data[, 1:10]
A = data[, 11:20]
n=nrow(data)
X = cbind(rep(1, nrow(data)), data[, 21:22])
X= matrix(rep(1, nrow(data)), ncol = 1)
nx= ncol(X)
TCOV = data[, 23:32]
nw=1
TCOV=NULL
nw=0
delta=NULL
itermax=100
hessian = TRUE
EMIRLS=TRUE
refgr=1
ndeltacum=NULL

ndeltacum=c(0,1,2,3)


param=c(pi[-1], beta, delta)

j = 1
betaL = list()
for (i in 1:ng){
  betaL[[i]] = beta[j:sum(nbeta[1:i])]
  j = sum(nbeta[1:i]) + 1
}

taux=ftauxLogit(pi, beta, ng, n, nbeta, A, Y, TCOV, delta, nw, nx, X)


solR = trajeR(Y = data[,1:10], A = data[,11:20],
       ng = 3, degre = c(0,3,4), Risk = data[,c(21,22)],
       Model = "LOGIT", Method = "L", hessian = FALSE,
       itermax = 300)
sol = trajeR(Y = data[,1:10], A = data[,11:20],
              ng = 3, degre = c(0,3,4), 
              Model = "LOGIT", Method = "L", hessian = TRUE,
              itermax = 300)
solEM = trajeR(Y = data[,1:10], A = data[,11:20],
             ng = 3, degre = c(0,3,4), 
             Model = "LOGIT", Method = "EM", hessian = TRUE,
             itermax = 300)
solEMIRLS = trajeR(Y = data[,1:10], A = data[,11:20],
               ng = 3, degre = c(0,3,4), 
               Model = "LOGIT", Method = "EMIRLS", hessian = FALSE,
               itermax = 20)




param=c( 0.3333333, 0.3333333 ,-5, -0.7254226, 0, 0, 0, 1.264463, 0, 0, 0, 0)

mbetaLtmp = function(i, t, ng, nbeta, nbetacum, A, taux, beta, TCOV, period, delta, ndeltacum, nw){
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

i=10
k=1
t=5
fexp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw)

i=i-1
k=k-1
t=t-1
fexp_cpp(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw)


library(trajeR)
LOGIT_datanew02 <- read.csv("/mnt/Travail/These/R/data/LOGIT_datanew02.csv")
d=as.matrix(LOGIT_datanew02)
paraminit = c( 0.5, 0.5, -5,0 ,0, 0.3876904, 0, 0)
s1=trajeR(Y = data[,2:8], A = data[,9:15],
       ng = 2, degre = c(2,2), 
       Model = "LOGIT", Method = "L", hessian = TRUE, paraminit = paraminit)
s2=trajeR(Y = data[,2:8], A = data[,9:15],
       ng = 2, degre = c(2,2), 
       Model = "LOGIT", Method = "EM", hessian = TRUE, paraminit = paraminit)

paramEM=c(s2$theta,s2$beta)
ng=2
nx=1
nbeta=c(2,2)+1
A=d[,9:15]
Y=d[,2:8]
TCOV=NULL
delta=NULL
nw=0
refgr=1
n=nrow(Y)
X = matrix(rep(1, nrow(data)), ncol=1)
nbetacum =c(0,3,6)


IEMLtmp <- function(paramEM, ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, refgr){
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
  print(B)
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

/***R
IEMLtmp(paramEM[-1], ng, nx, n, nbeta, A, Y, X, TCOV, delta, nw, refgr)
IEMLOGIT_cpp(paramEM, ng, nx, nbeta, n, A, Y, X, TCOV, nw, refgr)

i=10
t=2
mbetaLOGIT_cpp(i, t, ng, nbeta, A, beta, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mbetaL(i+1, t+1, ng, nbeta, nbetacum,A,  taux,  beta,TCOV, period, delta, ndeltacum, nw)
*/