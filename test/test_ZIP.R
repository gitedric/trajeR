library(trajeR)
library(MASS)
load(file = "/mnt/Travail/These/R/Package/trajeR/data/ZIP_data01")
data = as.matrix(dataZIP)
load(file = "/mnt/Travail/These/R/Package/trajeR/tests/ZIPdata01_sol")

Y = data[,2:6] 
A = data[,7:11]
ng = 2
nbeta=c(2,2)+1
nnu = c(1,1)+1
refgr = 1
beta = list(solL$beta[1:3],solL$beta[4:6])
nu = list(solL$nu[1:2],solL$nu[3:4])
beta=solL$beta
nu=solL$nu
theta=solL$theta
delta=solLTCOV$delta
pi=solEM$theta
n=nrow(data)
EMRILS = TRUE
X = cbind(rep(1, nrow(data)), data[, 12])
X= matrix(rep(1, nrow(data)), ncol = 1)
nx= ncol(X)

nnucum = c(0,2,4)
nbetacum=c(0,3,6)
TCOV=NULL
nw=0
delta=NULL
itermax=100
hessian = TRUE
EMIRLS=TRUE
refgr=1
ndeltacum=NULL
period=5
param=c(pi, unlist(beta), unlist(nu), unlist(delta))




param=c(solEMTCOV$theta, solEMTCOV$beta, solEMTCOV$nu, solEMTCOV$delta )

taux=ftauxZIP(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X)


paraminit = c(1,1,
              -3.066653,0,0,
              -0.046659,0,0,
              -3,0,
              -3,0)
solL = trajeR(Y = data[,2:6], A = data[,7:11], ng = 2, degre=c(2,2), degre.nu = c(1,1),
              Model="ZIP", Method = "L", hessian = TRUE,
              paraminit = paraminit)

paraminitEM = c(0.5,0.5,
              -3.066653,0,0,
              -0.046659,0,0,
              -3,0,
              -3,0)
paraminitEM = c(0.5,0.5,
                -0.738409,0,0,
                0.489316,0,0,
                -3, 0,
                -3, 0)

EMZIP_cpp(paraminitEM[-1], ng, nx, n, nbeta, nnu, A, Y, X, TCOV=NULL, nw=0, itermax=30, EMIRLS, refgr)
EMZIP(paraminitEM[-1], ng, nx, nbeta, nnu, n, A, Y, X, TCOV=NULL, delta, nw=0, itermax=30, EMIRLS, refgr)

solEM = trajeR(Y = data[,2:6], A = data[,7:11], ng = 2, degre=c(2,2), degre.nu = c(1,1),
              Model="ZIP", Method = "EM", hessian = TRUE, itermax=30,
              paraminit = paraminitEM)
solEMRisk = trajeR(Y = data[,2:6], A = data[,7:11], Risk = data[,12],
                  ng = 2, degre=c(2,2), degre.nu = c(1,1),
                  Model="ZIP", Method = "EM", hessian = TRUE,itermax = 300,
                  paraminit = paraminitR)
solEMIRLSRisk = trajeR(Y = data[,2:6], A = data[,7:11], Risk = data[,12],
                   ng = 2, degre=c(2,2), degre.nu = c(1,1),
                   Model="ZIP", Method = "EMIRLS", hessian = TRUE,itermax = 30,
                   paraminit = paraminitR)

solEMIRLS = trajeR(Y = data[,2:6], A = data[,7:11], ng = 2, degre=c(2,2), degre.nu = c(1,1),
               Model="ZIP", Method = "EMIRLS", hessian = TRUE, 
               paraminit = paraminitEM)

paraminitR = c(1,1,1,1,
              -3.066653,0,0,
              -0.046659,0,0,
              -3,0,
              -3,0)
solLRisk = trajeR(Y = data[,2:6], A = data[,7:11], Risk = data[,12],
                  ng = 2, degre=c(2,2), degre.nu = c(1,1),
                  Model="ZIP", Method = "L", hessian = TRUE,
                  paraminit = paraminit)

paraminit=c(1,1,
            -0.738409,0,0,
            0.489316,0,0,
            -3, 0,
            -3, 0,
            0, 0)
paraminit = c(1,1,
              -3.066653,0,0,
              -0.046659,0,0,
              -3,0,
              -3,0,
              0,0)
solLTCOV = trajeR(Y = data[,2:6], A = data[,7:11], TCOV = data[,13:17],
                  ng = 2, degre=c(2,2), degre.nu = c(1,1),
                  Model="ZIP", Method = "L", hessian = FALSE,
                  paraminit = paraminit)

paraminitR = c(1,1,1,1,
               -3.066653,0,0,
               -0.046659,0,0,
               -3,0,
               -3,0,
               0,0)
solLRiskTCOV = trajeR(Y = data[,2:6], A = data[,7:11], Risk = data[,12], TCOV = data[,13:17],
                  ng = 2, degre=c(2,2), degre.nu = c(1,1),
                  Model="ZIP", Method = "L", hessian = FALSE,
                  paraminit = paraminitR)



TCOV = data[,13:17]
nw=1
itermax=100
a1=optim(par = paraminit, fn = likelihoodZIP_cpp, gr = difLZIP_cpp,
      method = "BFGS",
      hessian = FALSE,
      control = list(fnscale=-1, trace=1, REPORT=1, maxit =itermax),
      ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
      nw = nw, TCOV = TCOV)
a2=optim(par = paraminit, fn = LikelihoodZIP, gr = difLZIP,
      method = "BFGS",
      hessian = FALSE,
      control = list(fnscale=-1, trace=1, REPORT=1, maxit =itermax),
      ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
      nw = nw, TCOV = TCOV)


paraminitEM = c(0.5,0.5,
                -0.738409,0,0,
                0.489316,0,0,
                -3, 0,
                -3, 0,
                0, 0)
paraminitEM = c(0.5,0.5,
              -3.066653,0,0,
              -0.046659,0,0,
              -3,0,
              -3,0,
              0,0)
solEMTCOV = trajeR(Y = data[,2:6], A = data[,7:11], TCOV = data[,13:17],
                   ng = 2, degre=c(2,2), degre.nu = c(1,1),
                   Model="ZIP", Method = "EM", hessian = TRUE,itermax=30,
                   paraminit = paraminitEM)
solEMIRLSTCOV = trajeR(Y = data[,2:6], A = data[,7:11], TCOV = data[,13:17],
                   ng = 2, degre=c(2,2), degre.nu = c(1,1), 
                   Model="ZIP", Method = "EMIRLS", hessian = TRUE,
                   paraminit = paraminitEM)


ndeltacum=c(0,1,2)
delta=c(0.5,1)
TCOV = data[,13:17]
ndeltacum=c(0,1,2)
nw=1
nw=1

QbetadeltakZIP_cpp(c( beta[(nbetacum[k]+1):nbetacum[k+1]],  delta[(ndeltacum[k]+1):ndeltacum[k+1]]),zk,zkSit, k-1,nbeta[k], nnu[k], n, A, Y,TCOV,nw)

QbetakZIP(beta[(nbetacum[k]+1):(nbetacum[k+1])], k, zk, zkSit,  nbeta[k], nnu[k],  n=to, A, Y, TCOV, delta, nw, ndeltacum)
difQbetakkZIP(beta[(nbetacum[k]+1):(nbetacum[k+1])], k, zk, zkSit,  nbeta[k], nnu[k],  n=to, A, Y, TCOV, delta, nw, ndeltacum)
difQbetakZIP_cpp(beta[(nbetacum[k]+1):(nbetacum[k+1])], zk, zkSit, k-1, nbeta[k], nnu[k],  n=to, A, Y, TCOV, delta, nw, ndeltacum)

res=c()
for (i in 1:n){
  for (t in 1:period){
    res=c(res, fzkSikt(pi, beta, nu, taux, k, i, t, ng, nbeta, nnu, n=to, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum) - taux[i,k]*fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n=to, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)    )
  }
}



pi=paraminitEM[1:2]
beta=paraminitEM[3:8]
nu=paraminitEM[9:12]
delta=c(0,0)

to=n

k=1

zk =  ftauxZIP(pi, beta, nu, ng, nbeta, nnu, n=to, A, Y, TCOV, delta, nw, nx, X)



betaIRLS = beta[(nbetacum[k]+1):nbetacum[k+1]]
nuIRLS = nu[(nnucum[k]+1):nnucum[k+1]]
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
tmps=c()
tmpss=c()
tmpsss=c()
Yit=c()
for (i in 1:to){
  for (t in 1:period){
    sikt = fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n=to, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
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
    tmps = c(tmps, sikt)
    tmpss=c(tmpss, rhoikt)
    tmpsss=c(tmpsss, fzkSikt(pi, beta, nu, taux, k, i, t, ng, nbeta, nnu, n=to, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum))
    Yit=c(Yit, Y[i,t])
  }
  tmp2 = c(tmp2, rep(zk[i,k], period))
}
Aw = matrix(tmp, nrow = nnu[k], byrow = FALSE)
Awp = matrix(tmp.p, nrow = nbeta[k], byrow = FALSE)
Z = diag(tmp2)
Awp%*%Z%*%diag(tmpr.p)%*%(matrix(Yit, ncol=1)-matrix(tmp2.p, ncol=1))

Aw%*%Z%*%(matrix(tmps, ncol=1)-matrix(tmpss, ncol=1))


ss=matrix(tmps, ncol=5, byrow = TRUE)
# 
# 
difQbetadeltakZIP_cpp(c(beta[(nbetacum[k]+1):nbetacum[k+1]],delta[(ndeltacum[k]+1):ndeltacum[k+1]]), zk, ss, k-1, nbeta[k], nnu[k],  n=to, A, Y, TCOV, nw)
difQnukZIP_cpp(nu[(nnucum[k]+1):(nnucum[k+1])], zk, ss, k-1, nbeta[k], nnu[k],  n=to, A, Y)
# difQnukkZIPtmp(nu[(nnucum[k]+1):(nnucum[k+1])], k, ss, nbeta[k], nnu[k], n, A, Y)
# difQbetakkZIPtmp(beta[(nbetacum[k]+1):nbetacum[k+1]],k, zk, ss,  nbeta[k], nnu[k],  n=to, A, Y, TCOV,delta[(ndeltacum[k]+1):ndeltacum[k+1]], nw, ndeltacum)
# 
# difQdeltakkZIPtmp(delta[k], k, zk, ss, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):nbetacum[k+1]], nw)
# 



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
optim(nu[(nnucum[k]+1):(nnucum[k+1])], QnukZIP, difQnukkZIP, method = "BFGS",
      zk = zk, zkSit = zkSit, k = k,
      nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
      control = list(fnscale=-1))$par
optim(nu[(nnucum[k]+1):(nnucum[k+1])], QnukZIP_cpp, difQnukZIP_cpp, method = "BFGS",
      zk = zk, Sikt = ss, k = k,
      nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
      control = list(fnscale=-1))$par
Aw = matrix(tmp, nrow = nnu[k], byrow = FALSE)
W = diag(tmpr)
Z = diag(tmp2)
S = matrix(tmp1, ncol =1)
QR = qr(t(Aw%*%sqrt(Z*W)))
Q = qr.Q(QR)
R = qr.R(QR)
backsolve(R,t(Q)%*%(sqrt(Z*W)%*%S))


QnukZIP_cpp(nu[(nnucum[k]+1):(nnucum[k+1])], zk, ss, k-1, nbeta[k], nnu[k],  n=to, A, Y)
difQnukZIP_cpp(nu[(nnucum[k]+1):(nnucum[k+1])], zk, msikt, k-1, nbeta[k], nnu[k],  n=to, A, Y)
difQnukZIP_cpp(nu[(nnucum[k]+1):(nnucum[k+1])], zk, ss, k-1, nbeta[k], nnu[k],  n=to, A, Y)


msikt=matrix(nrow=n, ncol=period)
for (i in 1:n){
  for (t in 1:period){
    msikt[i, t] = fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
  }
} 

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
QnukZIP_cpp(nu[(nnucum[k]+1):(nnucum[k+1])], zk, msikt, k-1, nbeta[k], nnu[k],  n=to, A, Y)
QnukZIP(nu[(nnucum[k]+1):(nnucum[k+1])], zk, zkSit, k, nbeta[k], nnu[k],  n=to, A, Y)

difQdeltakkZIPtmp <- function(delta, k, zk, Sikt, nbeta, nnu, n, A, Y, TCOV, beta, nw){
  period = ncol(A)
  deltas =c()
  for (l in 1:nw){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + sum(delta*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
        a = a + TCOV[i, t + (l-1)*period]*zik*(1-Sikt[i,t])*(Y[i,t]-exp(betaAit))
      }
    }
    deltas = c(deltas, a)
  }
  return(deltas)
}
difQbetakkZIPtmp <- function(beta, k, zk, Sikt, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum){
  period = ncol(A)
  betas =c()
  for (l in 1:nbeta){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
        a = a + A[i,t]**(l-1)*zik*(1-Sikt[i,t])*(Y[i,t]-exp(betaAit))
      }
    }
    betas = c(betas, a)
  }
  return(betas)
}
difQnukkZIPtmp <- function(nu, k, sikt, nbeta, nnu, n, A, Y){
  period = ncol(A)
  nus =c()
  for (l in 1:nnu){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        nuAit = sum(nu*A[i,t]**(0:(nnu-1)))
        #a = a + A[i,t]**(l-1)*zk[i,k]*(sikt[i, (i-1)*period+t]-exp(nuAit)/(1+exp(nuAit)))
        a = a + A[i,t]**(l-1)*zk[i,k]*(sikt[i, t]-exp(nuAit)/(1+exp(nuAit)))
      }
    }
    nus = c(nus, a)
  }
  return(nus)
}

/***R

EMZIP_cpp(paraminitEM[-1], ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax=200, EMIRLS, refgr)

EMZIPIRLS_cpp(paraminitEM[-1], ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax=10, EMIRLS, refgr)
*/

difQbetakkZIPtmp <- function(beta, k, zk, zkSit, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum){
  period = ncol(A)
  betas =c()
  for (l in 1:nbeta){
    a=0
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        ziksit = zkSit[i,(k-1)*period+t]
        print(sum(beta*A[i,t]**(0:(nbeta-1))))
              print(WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum))
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
        a = a + A[i,t]**(l-1)*(zik-ziksit)*(Y[i,t]-exp(betaAit))
      }
    }
    betas = c(betas, a)
  }
  return(betas)
}

delta=c(0,0)
TCOV = data[,13:17]
ndeltacum=c(0,1,2)
nw=1
EMZIP_cpp(paraminitEM[-2], ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
EMZIP(paraminitEM[-1], ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)

/***R
EMZIP_cpp(paraminitEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
*/
  /***R
EMZIP_cpp(paraminitEM[-c(1,13,14)], ng, nx, n, nbeta, nnu, A, Y, X, TCOV=NULL, nw=0, itermax, EMIRLS, refgr)
*/
#####################################
#
# pour EMIRLS
#
#######################################
itermax=100
param=c(0.5,
        -0.738409,0,0,
        0.489316,0,0,
        -3, 0,
        -3, 0,
        0,0)
EMZIPIRLS_cpp(param, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
EMZIPIRLS(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)


sum(abs(ftauxZIP_cpp(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X)-
ftauxZIP(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X))<10**(-6))


nnucum=c(0,2,4)
k=2
i=3
t=2
fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)


res=c()
for (k in 1:ng){
  for (i in 1:n){
    for(t in 1:period){
      res=c(res,abs(fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
      -fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
      )<10**(-6))
    }
  }
}
which(res==0)

fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)

/***R
itermax=2
EMZIPIRLStmp(paraminitEM[-2], ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)

EMZIPIRLS_cpp(paraminitEM[-2], ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
*/
  /***R
paraminitEM = c(0.5,0.5,
                -0.738409,0,0,
                0.489316,0,0,
                -3, 0,
                -3, 0,
                0, 0)
EMZIPIRLS_cpp(paraminitEM[-2], ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)
paraminitEM = c(0.5,0.5,
                -0.738409,0,0,
                0.489316,0,0,
                -3, 0,
                -3, 0)
EMZIPIRLS_cpp(paraminitEM[-2], ng, nx, n, nbeta, nnu, A, Y, X, TCOV=NULL, nw=0, itermax, EMIRLS, refgr)
*/

EMZIPIRLStmp <- function(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr){
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
    param = c(pi, beta, nu, delta)
  
  return(param)
}

















itermax=100
param=c(0.5,
        -0.738409,0,0,
        0.489316,0,0,
        -3, 0,
        -3, 0)


ndeltacum=c(0,1,2)
k=1
difQdeltakZIP_cpp(delta, zk, zkSit, k, nbeta[k],nnu[k],n,A,Y,TCOV,beta[(nbetacum[k]+1):(nbetacum[k+1])],nw,ndeltacum)


s1=0
s2=0
res=c()
for (i in 1:n){
  for (t in 1:period){
    a1=mbetaZIP(i, t, ng, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu)
    s1=a1+s1
    a2=mbetaZIP_cpp(i-1, t-1, ng, nbeta, A, beta,  taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y)
    s2=a2+s2
    
    res=c(res, all(abs(a1-a2)<10**(-6)))
  }
}

which(res==0)

/***R
i=1
t=2

a1=mbetaZIP(i, t, ng, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu)

a2=mbetaZIP_cpp(i-1, t-1, ng, nbeta, A, beta,  taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y)

all(abs(a1-a2)<10**(-6))
*/

  /***R

i=1
t=1
res=c()
for (i in 1:n){
  for (t in 1:period){
a1=mnuZIP(i, t, ng, nnu, nnucum, A, taux, nu)

a2=mnuZIP_cpp(i-1, t-1, ng, nbeta, A, beta,  taux, nbetacum, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n, Y)

res=c(res, all(abs(a1-a2)<10**(-6)))
  }
}
which(res==0)
*/
  /***R
k=2
l=3
  i=1
  
  res=c()
  for (i in 1:n){
    for (k in 1:ng){
      for (l in 1:ng){
        a1=BPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
        a2=BPiiklZIP_cpp(i-1, k-1,l-1, nbeta, A, Y, period, beta,  taux, nbetacum, TCOV, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n)        
        res=c(res, all(abs(a1-a2)<10**(-6)))
        }

    }
    }
  which(res==0)
*/

  
  
  
  /***R
k=2
l=1
i=1
res=c()
for (i in 1:n){
  for (k in 1:ng){
    for (l in 1:ng){
a1=NiklZIP(pi, beta, nu, k, i, l, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
a2=NiklZIP_cpp(i-1, k-1,l-1, nbeta, A, Y, period, beta,  taux, nbetacum, TCOV, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n)
res=c(res, all(abs(a1-a2)<10**(-6)))
    }
    
  }
}
which(res==0)

*/
  
  /***R
k=2
res=c()

  for (k in 1:ng){

covPiNukZIP(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nnu, nu)
covPiNukZIP_cpp(k-1, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)
res=c(res, all(abs(a1-a2)<10**(-6)))
}
which(res==0)
*/  
  /***R
k=1
l=1
res=c()
for (k in 1:ng){
  for (l in 1:ng){
a1=covBetakBetalZIP(k,l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu)
a2=covBetakBetalZIP_cpp(k-1,l-1, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)
res=c(res, all(abs(a1-a2)<10**(-6)))

  
}
}
which(res==0)


*/
  
  
  /***R
k=2
l=1

res=c()
for (k in 1:ng){
  for (l in 1:ng){
a1=covBetakNulZIP(k,l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu)
a2=covBetakNulZIP_cpp(k-1,l-1, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)
res=c(res, all(abs(a1-a2)<10**(-6)))


  }
}
which(res==0)


*/
  
  /***R
k=2
l=2
res=c()
for (k in 1:ng){
  for (l in 1:ng){
a1=covNukNulZIP(k,l, n, nbeta, nbetacum, A, Y, taux, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu)
a2=covNukNulZIP_cpp(k-1,l-1, ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)
res=c(res, all(abs(a1-a2)<10**(-6)))


  }
}
which(res==0)

*/
  
  
  /***R

covPiBetaZIP_cpp(ng, n, nx, nbeta, A, Y, X, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)
*/
  
  /***R

covPiNuZIP_cpp(ng, n, nx, nbeta, A, Y, X, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)
*/  
  

  
  
  
  /***R

covBetaZIP_cpp( ng, n, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, pi, nnucum, nu, nnu)
*/

  
  covPiZIP_cpp(n, ng, nx, pi, X, taux)

/***R
IEMZIP_cpp(param[-2], ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw, refgr)
IEMZIP(paramEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, delta, nw, refgr)
*/
  
  
  paramEM=c(solEM$theta[1], solEM$beta, solEM$nu)
  IEMZIP(paramEM, ng, nx, n, nbeta, nnu, A, Y, X, TCOV=NULL, delta=NULL, nw=0, refgr)



  for (t in 1:5){
    print(rhoikt(k, i, t, nnu, nnucum, A, nu)-rhoikt_cpp(k-1, i-1, t-1, nnu, nnucum, A, nu)    )
    print(fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)-
          fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum))
  }

res=c()
for (i in 1:n){
  res=c(res, BPiiklZIP_cpp(i-1, k-1, l-1, nbeta, A, Y, period, beta, taux, nbetacum, TCOV, delta, ndeltacum, nw, nnucum, nnu, nu, pi, n)
        -BPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum))
}


kp=l

taux[i,kp]*((1-taux[i,kp])/pi[kp]+taux[i,ng-1]/pi[ng-1])

covPiBetakZIPtmp <- function(k, ng, n, nbeta, nbetacum, A, Y, taux, pi, beta, TCOV, period, delta, ndeltacum, nw, nnucum, nu, nnu){
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
          tmp = tmp + BPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)*taux[i,k]*(-taux[i,kp]/pi[kp]+taux[i,ng]/pi[ng])
       #   print(BPiiklZIP(pi, beta, nu, k, l, i, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum))
      
          }
      }
      rcovPiBetak[kp,l] =  tmp
    }
  }
  return(rcovPiBetak)
}






t=5
fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)

zkSit = c()
for (k in 1:ng){
  tmp3 =c()
  for (i in 1:n){
    tmp4 = c()
    for (t in 1:period){
      tmp4 = c(tmp4, fzkSikt(pi, beta, nu, taux, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum))
    }
    tmp3 = rbind(tmp3, tmp4)
  }
  zkSit = cbind(zkSit, tmp3)
}
Sikt = c()
for (k in 1:ng){
  tmp3 =c()
  for (i in 1:n){
    tmp4 = c()
    for (t in 1:period){
      tmp4 = c(tmp4, fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum))
    }
    tmp3 = rbind(tmp3, tmp4)
  }
  Sikt = cbind(Sikt, tmp3)
}

#################################
k=1
i=7
t=2
fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)



k=2
optim(beta[(nbetacum[k]+1):(nbetacum[k+1])], QbetakZIP_cpp, difQbetakZIP_cpp, method = "BFGS",
      zk = zk, zkSit = zkSit, k = k-1,
      nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
      TCOV = TCOV, delta = delta, nw = nw, ndeltacum = ndeltacum,
      control = list(fnscale=-1))$par
optim(beta[(nbetacum[k]+1):(nbetacum[k+1])], QbetakZIP_cpp,  method = "BFGS",
      zk = zk, zkSit = zkSit, k = k-1,
      nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
      TCOV = TCOV, delta = delta, nw = nw, ndeltacum = ndeltacum,
      control = list(fnscale=-1))$par
beta[(nbetacum[k]+1):(nbetacum[k+1])]
nu[(nnucum[k]+1):(nnucum[k+1])]


optim(c(-10,100,-10), QbetakZIP_cpp, difQbetakZIP_cpp, method = "BFGS",
      zk = zk, zkSit = zkSit, k = k-1,
      nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
      TCOV = TCOV, delta = delta, nw = nw, ndeltacum = ndeltacum,
      control = list(fnscale=-1))$par
QbetakZIP_cpp(c(-10,100,-10), zk, zkSit, k-1, nbeta[k], nnu[k],  n, A, Y, TCOV, delta, nw, ndeltacum)

param=c(pi, unlist(beta), unlist(nu), unlist(delta))

QbetakZIP_cpp(beta[(nbetacum[k]+1):(nbetacum[k+1])], zk, zkSit, k-1, nbeta[k], nnu[k],  n, A, Y, TCOV, delta, nw, ndeltacum)
difQbetakZIP_cpp(beta[(nbetacum[k]+1):(nbetacum[k+1])], zk, zkSit, k-1, nbeta[k], nnu[k],  n, A, Y, TCOV, delta, nw, ndeltacum)


QnukZIP_cpp(nu[(nnucum[k]+1):(nnucum[k+1])], zk, zkSit, k-1, nbeta[k], nnu[k],  n, A, Y)
difQnukZIP_cpp(nu[(nnucum[k]+1):(nnucum[k+1])], zk, zkSit, k-1, nbeta[k], nnu[k],  n, A, Y)


EMZIP(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)

paraminit = c(1,1,
              -3.066653,0,0,
              -0.046659,0,0,
              -3,0,
              -3,0)
solL = trajeR(Y = data[,2:6], A = data[,7:11], ng = 2, degre=c(2,2), degre.nu = c(1,1),
              Model="ZIP", Method = "L", hessian = TRUE,
              paraminit = paraminit)

paraminit = c(1,1,1,1,
              -3.066653,0,0,
              -0.046659,0,0,
              -3,0,
              -3,0)
solLRisk = trajeR(Y = data[,2:6], A = data[,7:11], Risk = data[,12],
                  ng = 2, degre=c(2,2), degre.nu = c(1,1),
                  Model="ZIP", Method = "L", hessian = TRUE,
                  paraminit = paraminit)

paraminit=c(1,1,
            -0.738409,0,0,
            0.489316,0,0,
            -3, 0,
            -3, 0,
            0, 0)
solLTCOV = trajeR(Y = data[,2:6], A = data[,7:11], TCOV = data[,13:17],
                  ng = 2, degre=c(2,2), degre.nu = c(1,1),
                  Model="ZIP", Method = "L", hessian = TRUE,
                  paraminit = paraminit)

itermax=100
param=c(0.5,
        -0.738409,0,0,
        0.489316,0,0,
        -3, 0,
        -3, 0)
EMZIP_cpp(param, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)


itermax=100
param=c(0.5,
            -0.738409,0,0,
            0.489316,0,0,
            -3, 0,
            -3, 0)
EMZIPIRLS_cpp(param, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)


EMZIP(param[-2], ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)

r
param=a1

a1faux=a1
a2faux=a2

param=p1
itermax=2
a1=
  EMZIPIRLS(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)


/***R
IEMZIP_cpp(param[-2], ng, nx, nbeta, nnu, n, A, Y, X,TCOV, nw,refgr)
*/

a2=
  EMZIPIRLS_cpp(param, ng, nx, n, nbeta, nnu, A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)


all(abs(a1-a2[-2])<10**(-6))

param=param[-2]
pi = param[1:(ng-1)]
beta = param[(ng):(ng+sum(nbeta)-1)]
nu =  param[(ng+sum(nbeta)):(ng+sum(nbeta)+sum(nnu)-1)]
delta = param[-c(1:(ng+sum(nbeta)+sum(nnu)-1))]
pi = c(pi, 1-sum(pi))

lambdaikt_cpp(k-1, i-1, t-1, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw) 
lambdaikt_cpp(k-1, i-1, t-1, nbeta, nbetacum, A, beta, TCOV=NULL, period, delta=NULL, ndeltacum=NULL, nw=0) 


k=2
i=1
t=5

res=c()
for (i in 1:n)
{
  for (k in 1:ng)
  {
    for (t in 1 :period){
      res=c(res,lambdaikt(k, i, t, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw)-
      lambdaikt_cpp(k-1, i-1, t-1, nbeta, nbetacum, A, beta, TCOV, period, delta, ndeltacum, nw) )     
    }
  }
}
all(abs(res)<10**(-6))
res=c()
for (i in 1:n)
{
  for (k in 1:ng)
  {
    for (t in 1 :period){
res=c(res, rhoikt(k, i, t, nnu, nnucum, A, nu)
-rhoikt_cpp(k-1, i-1, t-1, nnu, nnucum, A, nu))
    }
  }
}
all(abs(res)<10**(-6))



betaL = list(beta[1:3], beta[4:6])
nuL = list(nu[1:2], nu[3:4])
deltaL = list(0.5,1)
delta=unlist(deltaL)
TCOV = data[,13:17]
nw=1
i=2
k=1

res=c()
for (i in 1:n){
  for (k in 1:ng){
    res=c(res,abs(gkZIP(betaL, nuL, i, k, nbeta, nnu, A, Y, TCOV, deltaL, nw)-
                    gkZIP_cpp(betaL,nuL ,i-1, k-1, nbeta, nnu, A, Y, TCOV, deltaL, nw))<10**(-6))
    
  }
}
which(res==0)
exp(muikt_cpp(beta[(nbetacum[k]+1):(nbetacum[k+1])], nbeta[k], i-1, period, A, TCOV, delta, nw, k-1))

nuikt_cpp(nuL[k], nnu[k], i-1,period, A,  k-1)

difLthetakZIP(c(theta, beta, nu, delta), k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)


X = cbind(rep(1, nrow(data)), data[, 12])
X= matrix(rep(1, nrow(data)), ncol = 1)
nx= ncol(X)

k=1
difLthetakZIP(c(theta, beta, nu, delta), k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)
difLthetakZIP_cpp(theta, betaL, nuL, deltaL, k-1, ng, nx, nbeta, nnu, n, A,Y,X,TCOV, nw)

k=1
difLbetakZIPtmp(c(theta, beta, nu, delta), k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)
difLbetakZIP_cpp(theta, betaL, nuL, deltaL, k-1, ng, nx, nbeta, nnu, n, A,Y,X,TCOV, nw)

difLnukZIP(c(theta, beta, nu, delta), k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)
difLnukZIP_cpp(theta, betaL, nuL, deltaL, k-1, ng, nx, nbeta, nnu, n, A,Y,X,TCOV, nw)


difLdeltakZIP(c(theta, beta, nu, delta), k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)
difLdeltakZIP_cpp(theta, betaL, nuL, deltaL, k-1, ng, nx, nbeta, nnu, n, A,Y,X,TCOV, nw)


difLZIP(c(theta, beta, nu, delta), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)-
difLZIP_cpp(c(theta, beta, nu, delta), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)


LikelihoodZIP(c(theta, beta, nu, delta), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)
likelihoodZIP_cpp(c(theta, beta, nu, delta), ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw)

optim(par = paraminit, fn = likelihoodZIP_cpp, gr = difLZIP_cpp,
      method = "BFGS",
      hessian = hessian,
      control = list(fnscale=-1, trace=1, REPORT=1, maxit = itermax),
      ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
      nw = nw, TCOV = TCOV)
optim(par = paraminit, fn = LikelihoodZIP, gr = difLZIP,
      method = "BFGS",
      hessian = hessian,
      control = list(fnscale=-1, trace=1, REPORT=1, maxit = itermax),
      ng = ng, nx = nx, n = n, nnu = nnu, A = A, Y = Y, X = X, nbeta = nbeta,
      nw = nw, TCOV = TCOV)



exp(sapply(1:period, function(s){
  sum(betaL[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, deltaL, nw, i, s, k)
}))


exp(muikt_cpp(betaL[[k]], nbeta[k], i-1, period, A, TCOV, deltaL, nw, k-1))








difLdeltakZIPtmp <- function(param, k, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, nw){
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

exp(sapply(1:period, function(s){
  sum(betaL[[k]]*A[i,s]**(0:(nbeta[k]-1))) + Wit(TCOV, period, delta, nw, i, s, k)
}))
