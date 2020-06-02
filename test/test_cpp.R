library(RcppArmadillo)
library(Rcpp)
library(rbenchmark)
library(MASS)

load("data/dataNORM01.RData")
#######################################
#Valeur initinale
X= cbind(matrix(rep(1, nrow(data)), ncol = 1), data[,11:15])
#X= matrix(rep(1, nrow(data)), ncol = 1)
nx= ncol(X)
n= nrow(data)
nbeta=c(3,3,3)
ng=3
i=12
k=3
ymin=-100000
ymax=1000000
nw=0
A=data[,6:10]
Y=data[,1:5]
itermax=10
EMIRLS=TRUE
TCOV=NULL
nw=0
delta=NULL
period=5
refgr=1
#nw=2
#TCOV=matrix(round(runif(n*10,0,1),2), ncol=10)

#taux = ftaux(solL$theta, solL$beta, solL$sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)

param=c(rep(1,18),0.7631785, 0, 0 ,3.731205 ,0, 0, 6.699231, 0, 0, rep(3.067976,3)) 
taux = ftaux(param[1:18], param[19:27], param[28:30], ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)

######################################"
gkCNORM(beta, sigma, i, k, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)



load("/mnt/Travail/These/R/Package/trajeR/tests/CNORM_data07")
load("/mnt/Travail/These/R/Package/trajeR/tests/CNORM_data07Censored")
data = as.matrix(data)
dataC = as.matrix(dataC)

i=1
k=2
period=10
Y = data[,2:11]
A =data[,12:21]
ng = 3
degre = c(0,3,4)
beta=list(c(1),
          c(1,1,1,1),
          c(1,1,1,1,1))
sigma=rep(1,3)
delta= list(c(1),
            c(1),
            c(1))
TCOV = data[,24:33]
nw=1
nbeta=c(1,4,5)
X= matrix(rep(1, nrow(data)), ncol = 1)
nx= ncol(X)
pi=c(1/3,1/3,1/3)
n=nrow(data)
ymin=-100000
ymax=1000000
# delta=NULL
# TCOV=NULL
# nw=0

k=1
gkCNORM_cpp(beta, sigma,  i,k, nbeta, A, Y, ymin, ymax,  TCOV=NULL, delta=NULL, nw=0)
gkCNORM(beta, sigma, i, k, nbeta, A, Y, ymin, ymax, TCOV=NULL, delta=NULL, nw=0)


Likelihoodalpha_cpp(c(pi, unlist(beta), sigma,unlist(delta)),
                    ng, 
                    nx,
                    nbeta,
                    n,
                    A,
                    Y,
                    X,
                    ymin,
                    ymax, 
                    TCOV,
                    nw)
Likelihoodalpha(c(pi, unlist(beta), sigma,unlist(delta)),
                    ng, 
                    nx,
                    nbeta,
                    n,
                    A,
                    Y,
                    X,
                    ymin,
                    ymax, 
                    TCOV,
                    nw)
a1=difLalphaunique(c(pi, unlist(beta), sigma,unlist(delta)), ng, nx, nbeta, n, A, Y, X, ymin, ymax=30, TCOV, nw)
a2=difLalphaunique_cpp(c(pi, unlist(beta), sigma,unlist(delta)), ng, nx, nbeta, n, A, Y, X, ymin, ymax=30, TCOV, nw)


k=1
difLdeltakalpha(c(pi, unlist(beta), sigma,unlist(delta)), k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw)
difLdeltakalpha_cpp(pi, beta, sigma,delta, k-1, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw)

  


#################################################
sol1 = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
                 Model="CNORM", Method = "L", ssigma = TRUE, 
                 hessian = FALSE)
sol2 = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
                  Model="CNORM", Method = "L", ssigma = FALSE, 
                  hessian = TRUE)
sol1EM = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
                   Model="CNORM", Method = "EM", ssigma = FALSE, 
                   hessian = FALSE)
sol2R = trajeR(data[,1:5], data[,6:10], Risk = data[, 11:15], ng = 3, degre=c(2,2,2), 
              Model="CNORM", Method = "L", ssigma = FALSE, 
              hessian = TRUE)
sol2EMR = trajeR(data[,1:5], data[,6:10], Risk = data[, 11:15], ng = 3, degre=c(2,2,2), 
                Model="CNORM", Method = "EM", ssigma = FALSE, 
                hessian = FALSE)


############################# exempale vignette
load(file = "/mnt/Travail/These/R/Package/trajeR/tests/CNORM_data07_sol")
load(file = "/mnt/Travail/These/R/Package/trajeR/tests/CNORM_data07Censored_sol")
solL = trajeR(Y = data[,2:11], A = data[,12:21],
              ng = 3, degre = c(0,3,4),
              Model = "CNORM", Method = "L",  hessian = TRUE, ssigma = FALSE)
solLs = trajeR(Y = data[,2:11], A = data[,12:21],
              ng = 3, degre = c(0,3,4),
              Model = "CNORM", Method = "L",  hessian = TRUE, ssigma = TRUE)
solEM = trajeR(Y = data[,2:11], A = data[,12:21],
              ng = 3, degre = c(0,3,4),
              Model = "CNORM", Method = "EM",  hessian = TRUE, ssigma = FALSE)
solLRisks = trajeR(Y = data[,2:11], A = data[,12:21], Risk = data[,42:43],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "L", ssigma = TRUE, hessian = FALSE)
solEMRisk = trajeR(Y = data[,2:11], A = data[,12:21], Risk = data[,42:43],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "EM", ssigma = FALSE, hessian = FALSE)
solLRisk = trajeR(Y = data[,2:11], A = data[,12:21], Risk = data[,42:43],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "L", ssigma = FALSE, hessian = FALSE)
solLTCOV2 = trajeR(Y = data[,2:11], A = data[,12:21], TCOV = data[,22:41],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "L", ssigma = FALSE, hessian = FALSE)
solLTCOV2s = trajeR(Y = data[,2:11], A = data[,12:21], TCOV = data[,22:41],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "L", ssigma = TRUE, hessian = FALSE)
solLC = trajeR(Y = dataC[,2:11], A = dataC[,12:21],
               ng = 3, degre = c(0,3,4),
               Model="CNORM", Method = "L", ssigma = FALSE, hessian = FALSE,
               ymin=2, ymax=23)
solEMC = trajeR(Y = dataC[,2:11], A = dataC[,12:21],
               ng = 3, degre = c(0,3,4),
               Model="CNORM", Method = "EM", ssigma = FALSE, hessian = FALSE,
               ymin=2, ymax=23)
solLCs = trajeR(Y = dataC[,2:11], A = dataC[,12:21],
               ng = 3, degre = c(0,3,4),
               Model="CNORM", Method = "L", ssigma = TRUE, hessian = FALSE,
               ymin=2, ymax=23)
solLCTCOV = trajeR(Y = data[,2:11], A = data[,12:21], TCOV = data[,22:31],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "L", ssigma = FALSE, hessian = FALSE,
                   ymin = 2,  ymax = 23)
solLCTCOVs = trajeR(Y = data[,2:11], A = data[,12:21], TCOV = data[,22:31],
                    ng = 3, degre = c(0,3,4),
                    Model = "CNORM", Method = "L", ssigma = TRUE, hessian = FALSE,
                    ymin = 2,  ymax = 23)
solEMCTCOV = trajeR(Y = data[,2:11], A = data[,12:21], TCOV = data[,22:31],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "EM", ssigma = FALSE, hessian = FALSE,
                   ymin = 2,  ymax = 23)
solEMCTCOVs = trajeR(Y = data[,2:11], A = data[,12:21], TCOV = data[,22:31],
                    ng = 3, degre = c(0,3,4),
                    Model = "CNORM", Method = "EM", ssigma = TRUE, hessian = FALSE,
                    ymin = 2,  ymax = 23)
solEMCs = trajeR(Y = dataC[,2:11], A = dataC[,12:21],
                ng = 3, degre = c(0,3,4),
                Model="CNORM", Method = "EM", ssigma = TRUE, hessian = FALSE,
                ymin=2, ymax=23)


library(trajeR)
load("/mnt/Travail/These/R/Package/trajeR/tests/CNORM_data07")
load("/mnt/Travail/These/R/Package/trajeR/tests/CNORM_data07Censored")
data = as.matrix(data)
dataC = as.matrix(dataC)
i=1
k=1
t=1
X= matrix(rep(1, nrow(data)), ncol = 1)
X= cbind(matrix(rep(1, nrow(data)), ncol = 1), data[,42:43])
theta=solLRisk$theta
nx= ncol(X)
Y = data[,2:11]
A = data[,12:21]
TCOV = data[,22:41]
nw=2
delta=solLTCOV2s$delta
TCOV=NULL
nw=0
ng = 3
nbeta = c(0,3,4)+1
n = nrow(data)
ymin = -100000
ymax = 100000
k=1
pi=solEM$theta
beta=solEM$beta
sigma=solEM$sigma
nbetacum=c(0,cumsum(nbeta))
ndeltacum =NULL
delta=NULL
delta=solEMCTCOV$delta
ndeltacum = cumsum(c(0, rep(nw, ng)))
ndeltacum=NULL
pi=theta
taux = ftaux(pi, beta, sigma, ng, nbeta, n, A, Y, ymin, ymax, TCOV, delta, nw, nx, X)
ndeltacum = cumsum(c(0, rep(nw, ng)))
#delta=solLCTCOV$delta

param = c(pi[-1], beta, sigma, delta)
paramEM=param
param = c(pi, beta, sigma, delta)

itermax=3
EMIRLS=TRUE

EMSigmauniquetmp(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, delta, nw, itermax, EMIRLS )
EMSigmaunique_cpp(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV,  nw, itermax, EMIRLS,refgr=1 )

i=1
t=1
mbeta(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mbetaCNORM_cpp(i-1, t-1, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mdelta(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mdeltaCNORM_cpp(i-1, t-1, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mbetadelta(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mbetadeltaCNORM_cpp(i-1, t-1, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mdeltasigma(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mdeltasigmaCNORM_cpp(i-1, t-1, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mdelta(i, t, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
mdeltasigmaCNORM_cpp(i-1, t-1, ng, nbeta, A, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
msigma(i, t, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)
msigmaCNORM_cpp(i-1, t-1, ng, nbeta, A, Y, beta, sigma, taux, nbetacum, TCOV, period, delta, ndeltacum, nw)


EMSigmauniquetmp <- function(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, delta, nw, itermax, EMIRLS){
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
