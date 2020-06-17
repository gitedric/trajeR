to=n
c(0.5,0.5,
  -0.738409,0,0,
  0.489316,0,0,
  -3, 0,
  -3, 0,
  0, 0)
k=2


pi = c(0.5,0.5)
beta=c(  -0.738409,0,0,
         0.489316,0,0)
nu=c(  -3, 0,
       -3, 0)
delta=c(0,0)
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




ss=matrix(tmps, ncol=5, byrow = TRUE)


difQbetadeltakZIP_cpp(c(beta[(nbetacum[k]+1):nbetacum[k+1]],delta[(ndeltacum[k]+1):ndeltacum[k+1]]), zk, ss, k-1, nbeta[k], nnu[k],  n=to, A, Y, TCOV, nw)



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

optim(c(beta[(nbetacum[k]+1):(nbetacum[k+1])],delta[k]), QbetadeltakZIP_cpp, difQbetadeltakZIP_cpp, method = "BFGS",
      zk = zk, Sikt = ss, k = k-1,
      nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
      TCOV = TCOV,nw = nw,
      control = list(fnscale=-1))$par

Awp = matrix(tmp.p, nrow = nbeta[k], byrow = FALSE)
Wp = diag(tmp2.p)
Zs = diag(tmp2*tmpr.p)# for Sp%*%Zp in  the defintion 7.6.2
Sp = matrix(tmp1.p, ncol =1)
QR = qr(t(Awp%*%sqrt(Zs*Wp)))
Q = qr.Q(QR)
R = qr.R(QR)
backsolve(R,t(Q)%*%(sqrt(Zs*Wp)%*%Sp))
