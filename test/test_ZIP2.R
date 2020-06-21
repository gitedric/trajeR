library(trajeR)

namefile  = "/mnt/Travail/These/R/data/test/ZIP_data_2grTCOV3_"

cpt = 1
name  = paste0(namefile, cpt, ".csv")
data <- read.csv(name, row.names=1)
data = as.matrix(data)
ind = 1:10


paraminit = c(1,1,
              5.65184745935784, 0, 0,
              5.87659981601246, 0, 0,
              -3,0,
              -3,0
              )
sol1=trajeR(Y = data[,ind], A = data[,ind + 10],
           ng = 2, degre=c(2,2), 
           Model="ZIP", degre.nu = c(1, 1),
           Method = "L", hessian = FALSE, itermax = 100
)
sol2=trajeR(Y = data[,ind], A = data[,ind + 10],
           ng = 2, degre=c(2,2), 
           Model="ZIP", degre.nu = c(1, 1),
           Method = "EM", hessian = FALSE, itermax = 100
)
sol3=trajeR(Y = data[,ind], A = data[,ind + 10],
            ng = 2, degre=c(2,2), 
            Model="ZIP", degre.nu = c(1, 1),
            Method = "EMIRLS", hessian = FALSE, itermax = 100
)
sol4=trajeR(Y = data[,ind], A = data[,ind + 10],TCOV = data[, ind + 20],
            ng = 2, degre=c(2,2), 
            Model="ZIP", degre.nu = c(1, 1),
            Method = "L", hessian = FALSE, itermax = 100
)
sol5=trajeR(Y = data[,ind], A = data[,ind + 10],TCOV = data[, ind + 20],
            ng = 2, degre=c(2,2), 
            Model="ZIP", degre.nu = c(1, 1),
            Method = "EM", hessian = FALSE, itermax = 100
)
sol6=trajeR(Y = data[,ind], A = data[,ind + 10],TCOV = data[, ind + 20],
            ng = 2, degre=c(2,2), 
            Model="ZIP", degre.nu = c(1, 1),
            Method = "EMIRLS", hessian = FALSE, itermax = 100
)

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

EMZIP(paraminit, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)
EMZIPtmp(paraminit, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr)
EMZIP_cpp(paraminit, ng, nx, n, nbeta, nnu,  A, Y, X, TCOV, nw, itermax, EMIRLS, refgr)


plot(sol5, Y = data[,ind], A = data[,ind + 10])


EMZIPtmp <- function(param, ng, nx, nbeta, nnu, n, A, Y, X, TCOV, delta, nw, itermax, EMIRLS, refgr){
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
    siktm = matrix(rep(0, n*period), ncol=period)
    for (i in 1:n){
      for (t in 1:period){
        siktm[i,t] = fSikt(pi, beta, nu, k, i, t, ng, nbeta, nnu, n, A, Y, TCOV, delta, nw, ndeltacum, period, nbetacum, nnucum)
      }
    }
    #M-step
    newbeta = c()
    newnu = c()
    newdelta = c()
    newbetadelta = c()
    for (k in 1:ng){
      newnu = c(newnu, optim(nu[(nnucum[k]+1):(nnucum[k+1])], QnukZIP, difQnukkZIP, method = "BFGS",
                             zk = zk, zkSit = zkSit, k = k,
                             nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
                             control = list(fnscale=-1))$par)
      
      newbetadelta = optim(c(beta[(nbetacum[k]+1):(nbetacum[k+1])],delta[(ndeltacum[k]+1):(ndeltacum[k+1])]), 
            QbetadeltakZIP, difQbetadeltakZIP, method = "BFGS",
            zk = zk, Sikt = siktm, k = k,
            nbeta = nbeta[k], nnu = nnu[k], n = n, A = A, Y = Y,
            TCOV = TCOV,  nw = nw, 
            control = list(fnscale=-1))$par
      newbeta= c(newbeta, newbetadelta[1:nbeta[k]])
      newdelta= c(newdelta, newbetadelta[-c(1:nbeta[k])])

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
