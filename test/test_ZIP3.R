library(trajeR)
paraminit = c(0.5,
              1.242196, 0, 0,
              1.769089, 0, 0, -3, 0, -3, 0,0,0
)
param=paraminit
pi = param[1:(ng-1)]
beta = param[(ng):(ng+sum(nbeta)-1)]
nu =  param[(ng+sum(nbeta)):(ng+sum(nbeta)+sum(nnu)-1)]
delta = param[-c(1:(ng+sum(nbeta)+sum(nnu)-1))]
pi = c(pi, 1-sum(pi))
n=nrow(Y)
X=matrix(rep(1,n), ncol=1)
nx=1
ng=2
nbeta=c(3,3)
degre=c(3,3)
Y = data[,ind]
A = data[,ind + 10]
TCOV = data[, ind + 20]
degre.nu = c(1,1)
delta=c(0,0)
nw=1
itermax=100
EMIRLS=FALSE
refgr=1


k=2

 
  siktm = matrix(rep(0, n*period), ncol=period)
  for (i in 1:n){
    for (t in 1:period){
      siktm[i,t] = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
    }
  }
  siktm2 = matrix(rep(0, n*period), ncol=period)
  for (i in 1:n){
    for (t in 1:period){
      siktm2[i,t] = fSikt_cpp(pi, beta, nu, k-1, i-1, t-1, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
    }
  }

  
  
  zk = ftauxZIP(pi, beta, nu, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, nx, X)

  difQbetadeltakZIP(c(beta[(nbetacum[k]+1):(nbetacum[k+1])],delta[(ndeltacum[k]+1):(ndeltacum[k+1])]), k, zk, siktm, nbeta[k], nnu[k], n, A, Y, TCOV, nw)
  difQbetadeltakZIP_cpp(c(beta[(nbetacum[k]+1):(nbetacum[k+1])],delta[(ndeltacum[k]+1):(ndeltacum[k+1])]), zk, sikt, k-1,nbeta[k], nnu[k], n, A, Y, TCOV, nw)  
  
  difQbetakkZIPtmp(beta[(nbetacum[k]+1):(nbetacum[k+1])],k, zk, siktm,  nbeta[k], nnu[k], n, A, Y, TCOV, delta[(ndeltacum[k]+1):(ndeltacum[k+1])], nw)

  QbetakZIPtmp(beta[(nbetacum[k]+1):(nbetacum[k+1])], zk, siktm, k, nbeta[k], nnu[k], n, A, Y, TCOV, delta[(ndeltacum[k]+1):(ndeltacum[k+1])], nw)
  QbetakZIP_cpp(c(beta[(nbetacum[k]+1):(nbetacum[k+1])],delta[(ndeltacum[k]+1):(ndeltacum[k+1])]), k, zk, sikt, nbeta[k], nnu[k], n, A, Y, TCOV, nw, ndeltacum)
  QbetadeltakZIP(c(beta[(nbetacum[k]+1):(nbetacum[k+1])],delta[(ndeltacum[k]+1):(ndeltacum[k+1])]), k, zk, siktm, nbeta[k], nnu[k], n, A, Y, TCOV, nw)
  QdeltakZIPtmp(delta[(ndeltacum[k]+1):(ndeltacum[k+1])], zk, siktm, k, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):(nbetacum[k+1])], nw)
  
  optim(c(beta[(nbetacum[k]+1):(nbetacum[k+1])],delta[(ndeltacum[k]+1):(ndeltacum[k+1])]), 
        QbetadeltakZIP, difQbetadeltakZIP, method = "BFGS",
        zk = zk, Sikt = siktm, k = k,
        nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
        TCOV = TCOV,  nw = nw,
        control = list(fnscale=-1))$par
  
  optim(beta[(nbetacum[k]+1):(nbetacum[k+1])], 
        QbetakZIPtmp, difQbetakkZIPtmp, method = "BFGS",
        zk = zk, Sikt = siktm, k = k,
        nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
        TCOV = TCOV,  nw = nw,  delta = delta[(ndeltacum[k]+1):(ndeltacum[k+1])],
        control = list(fnscale=-1))$par
  
  optim(delta[(ndeltacum[k]+1):(ndeltacum[k+1])], 
        QdeltakZIPtmp, difQdeltakkZIPtmp, method = "BFGS",
        zk = zk, Sikt = siktm, k = k,
        nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
        TCOV = TCOV,  nw = nw, beta=beta[(nbetacum[k]+1):(nbetacum[k+1])],
        control = list(fnscale=-1))$par
  
 optim(nu[(nnucum[k]+1):(nnucum[k+1])], QnukZIPtmp, difQnukkZIPtmp, method = "BFGS",
                         zk = zk, Sikt = siktm, k = k,
                         nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
                         control = list(fnscale=-1))$par
 
 QnukZIPtmp(nu[(nnucum[k]+1):(nnucum[k+1])],zk, sikt, k, nbteta[k], nnu[k],n, A, Y)
 
  trajeR(Y = data[,ind], A = data[,ind + 10],TCOV = data[, ind + 20],
              ng = 2, degre=c(2,2), 
              Model="ZIP", degre.nu = c(1, 1),
              Method = "EM", hessian = FALSE, itermax = 2
  )  
  trajeR(Y = data[,ind], A = data[,ind + 10],TCOV = data[, ind + 20],
         ng = 2, degre=c(2,2), 
         Model="ZIP", degre.nu = c(1, 1),
         Method = "EMIRLS", hessian = FALSE, itermax = 2
  )  

  difQdeltakkZIPtmp(delta[(ndeltacum[k]+1):(ndeltacum[k+1])], k, zk, sikt, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):(nbetacum[k+1])], nw)
  
  difQbetakkZIPtmp(beta[(nbetacum[k]+1):(nbetacum[k+1])], k, zk, sikt, nbeta[k], nnu[k], n, A, Y, TCOV, delta[(ndeltacum[k]+1):(ndeltacum[k+1])], nw)
  
  
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
 max( zkSit[,11:20]-zk[,k]*siktm2)
  
  a1=QdeltakZIPzks( 0.2472527, zk, zkSit, k, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):(nbetacum[k+1])], nw)
  a2=QdeltakZIPtmp( 0.2472527, zk, siktm2, k, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):(nbetacum[k+1])], nw)
  QdeltakZIP_cpp( 0.2472527, zk, siktm2, k-1, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):(nbetacum[k+1])], nw)
  
mat=c()  
  for (i in 1:n){
    for (t in 1:period){
      mat=rbind(mat, cbind(Q1(i,t, 0.2472527, zk, zkSit, k, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):(nbetacum[k+1])], nw),
      Q2(i,t, 0.2472527, zk, siktm2, k, nbeta[k], nnu[k], n, A, Y, TCOV, beta[(nbetacum[k]+1):(nbetacum[k+1])], nw), i, t))
    }
  }
mat[which.max(abs(a1-a2)),]
i=56
t=3
  
  QdeltakZIPzks <- function(delta, zk, zkSit, k, nbeta, nnu, n, A, Y, TCOV, beta, nw){
    period = ncol(A)
    a = c()
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        ziksit = zkSit[i,(k-1)*period+t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + sum(delta*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
        #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
        a = c(a, (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)))
      }
    }
    return(a)
  }
  QdeltakZIPtmp <- function(delta, zk, Sikt, k, nbeta, nnu, n, A, Y, TCOV, beta, nw){
    period = ncol(A)
    a = c()
    for (i in 1:n){
      zik = zk[i,k]
      for (t in 1:period){
        sikt = Sikt[i,t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
        #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
        a = c(a, zik*(1-sikt)*(Y[i,t]*betaAit-exp(betaAit)))
      }
    }
    return(a)
  }
  
  Q1 <- function(i,t,delta, zk, zkSit, k, nbeta, nnu, n, A, Y, TCOV, beta, nw){
    period = ncol(A)
    a = c()

      zik = zk[i,k]
  
        ziksit = zkSit[i,(k-1)*period+t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + sum(delta*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)])
        #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
        a = c(a, (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)))

    return(a)
  }
  Q2 <- function(i,t,delta, zk, Sikt, k, nbeta, nnu, n, A, Y, TCOV, beta, nw){
    period = ncol(A)
    a = c()

      zik = zk[i,k]

        sikt = Sikt[i,t]
        betaAit = sum(beta*A[i,t]**(0:(nbeta-1))) + WitEMtmp(TCOV, period, delta, nw, i, t, k)
        #a = a + (zik-ziksit)*(Y[i,t]*betaAit-exp(betaAit)-log(factorial(Y[i,t])))
        a = c(a, zik*(1-sikt)*(Y[i,t]*betaAit-exp(betaAit)))

    return(a)
  }
  
  fzkSikt(pi, beta, nu, zk, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
  fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
  
  fzkSikt2(pi, beta, nu, zk, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
  fSikt2(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
  
  
  fzkSikt2 <- function(pi, beta, nu, zk, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum){
    if (Y[i,t]>0){
      prob = 0
    }else{
      nuikt = sum(nu[(nnucum[k]+1):(nnucum[k+1])]*A[i,t]**(0:(nnu[k]-1)))
      lambdaikt = exp(
        sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,t]**(0:(nbeta[k]-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum)
        )
      prob = zk[i,k]*exp(nuikt+lambdaikt)/(1+exp(nuikt+lambdaikt))
    }
    return(prob)
  }
  fSikt2 <- function(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum){
    if (Y[i,t]>0){
      prob = 0
    }else{
      nuikt = sum(nu[(nnucum[k]+1):(nnucum[k+1])]*A[i,t]**(0:(nnu[k]-1)))
      lambdaikt = exp(sum(beta[(nbetacum[k]+1):(nbetacum[k+1])]*A[i,t]**(0:(nbeta[k]-1))) + WitEM(TCOV, period, delta, nw, i, t, k, ndeltacum))
      prob = exp(nuikt+lambdaikt)/(1+exp(nuikt+lambdaikt))
    }
    return(prob)
  }
  