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
              Model = "CNORM", Method = "EM",  hessian = FALSE, ssigma = FALSE)
solLRisks = trajeR(Y = data[,2:11], A = data[,12:21], Risk = data[,42:43],
                   ng = 3, degre = c(0,3,4),
                   Model = "CNORM", Method = "L", ssigma = TRUE, hessian = FALSE)
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

X= matrix(rep(1, nrow(data)), ncol = 1)
nx= ncol(X)
Y = data[,2:11]
A = data[,12:21]
TCOV = data[,22:31]
nw=1
ng = 3
nbeta = c(0,3,4)+1
n = nrow(data)
ymin = 2
ymax = 23
k=1
j = 1
betatmp=solLCTCOV$beta
beta = list()
for (i in 1:ng){
  beta[[i]] = betatmp[j:sum(nbeta[1:i])]
  j = sum(nbeta[1:i]) + 1
}
ndeltacum = cumsum(c(0, rep(nw, ng)))
deltatmp=solLCTCOV$delta
delta = list()
for (i in 1:ng){
  delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
}

k=1
difLdeltakalpha(c(solLCTCOV$theta, solLCTCOV$beta, solLCTCOV$sigma,solLCTCOV$delta), k, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw)
difLdeltakalpha_cpp(solLCTCOV$theta, beta, solLCTCOV$sigma,delta, k-1, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw)
