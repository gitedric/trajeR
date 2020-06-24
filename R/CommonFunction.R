#################################################################################
# Piik
#################################################################################
piik <- function(theta, i, k, ng, X){
  ntheta = ncol(X)
  tmp = exp(sapply(1:ng,function(s){theta[((s-1)*ntheta+1):(s*ntheta)]%*%X[i,]}))
  return(tmp[k]/sum(tmp))
}
#################################################################################
# likelihood
#################################################################################
Likelihood <- function(param, model, method, ng, nx, n, nbeta, nw, A, Y, X, TCOV, ymin = NULL, ymax = NULL, nnu = NULL, fct = NULL){
  if (model == "CNORM"){
    if (method == "L"){
      a = likelihoodCNORM_cpp(param, ng, nx, nbeta, n, A, Y, X, ymin, ymax, TCOV, nw)
    }else{
      a = likelihoodEM_cpp(n, ng, nbeta, beta=param[(ng+1):(ng+sum(nbeta))],
                           sigma=param[(ng+sum(nbeta)+1):(ng+sum(nbeta)+ng)],
                           pi=param[1:(ng)],
                           A, Y, ymin, ymax, TCOV,
                           delta=param[-c(1:(ng+sum(nbeta)+ng))], nw)
    }
  } else if (model == "LOGIT"){
    if (method == "L"){
      a = likelihoodLOGIT_cpp(param, ng, nx, n, nbeta, A, Y, X, TCOV, nw)
    }else{
      a = likelihoodEMLOGIT_cpp(n, ng, nbeta, 
                              beta=param[(ng+1):(ng+sum(nbeta))],
                              pi=param[1:(ng)],
                              A, Y, TCOV, 
                              delta=param[-c(1:(ng+sum(nbeta)))], nw)
    }
  } else if (model == "ZIP"){
   if (method == "L"){
     a = likelihoodZIP_cpp(param, ng,nx, nbeta, nnu, n, A, Y, X, TCOV, nw)
   }else{
     a = likelihoodEMZIP_cpp(n, ng, nbeta, nnu, 
                             beta=param[(ng+1):(ng+sum(nbeta))],
                             nu=param[(ng+sum(nbeta)+1): (ng+sum(nbeta)+sum(nnu))],
                             pi=param[1:(ng)],
                             A, Y, TCOV,
                             delta=param[-c(1:(ng+sum(nbeta)+sum(nnu)))], nw)
   }
  }else{
    a = LikelihoodNL(param, ng, nx, nbeta, n, A, Y, X, TCOV, fct)
  }
  return(a)
}
#################################################################################
# compute the value of Wit for i and t and given k
#################################################################################
Wit <- function(TCOV, period, delta, nw, i, t, k){
  if (nw == 0){
    return(0)
  }else{
    return(sum(delta[[k]]*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)]))
  }
}
# compute the value of Wit for i and t and given k for the EM algorithm
WitEM <- function(TCOV, period, delta, nw, i, t, k, ndeltacum){
  if (nw == 0){
    return(0)
  }else{
    return(sum(delta[(ndeltacum[k]+1):(ndeltacum[k+1])]*TCOV[i, seq(from = t, to = t+(nw-1)*period, by = period)]))
  }
}
#################################################################################
# Calculate the probaility of memebrship for each data
#################################################################################
#' GroupProb calculate the membership probability of each value of the data.
#'
#' @param Obj a trajectory object that is return by trajeR function.
#' @param Y a real matrix. The data.
#' @param A a real matrix. The time variable.
#' @param TCOV a real matrix. Optionnal, by default the value is NULL. It contained the time dependent covariate.
#' @param X a real matrix. Optionnal, by default the value is NULL. It contained a caovraite that modify the probaility memebrship.
#'
#' @return a real matrix. For each individual i in the data, this matrix contained the membership probability of each group.
#'
#' @export
#'
#' @examples
#' load("data/dataNORM01.RData")
#' solL = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
#'               Model="CNORM", Method = "L", ssigma = FALSE, 
#'               hessian = TRUE)
#' GroupProb(solL, Y=data[,1:5], A=data[,6:10])
GroupProb <- function(Obj, Y, A, TCOV = NULL, X = NULL){
  n = Obj$Size
  ng = Obj$groups
  nbeta = Obj$degre + 1
  ymin = Obj$min
  ymax = Obj$max
  betatmp = Obj$beta
  j = 1
  beta = list()
  for (i in 1:length(nbeta)){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  delta = Obj$delta
  if (any(is.na(delta))){
    nw = 0
    delta = rep(list(0), ng)
  }else{
    nw = length(delta)/ng
  }
  theta = Obj$theta
  if (is.null(X)){
    X = cbind(rep(1, n))
  }
  res = c()
  if (Obj$Model == "LOGIT"){
    if (Obj$Method == "L"){
    for (i in 1:n){
      tmp = sapply(1:ng, function(s){
        piik(theta, 1, s, ng, X)*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)})
      res = rbind(res, tmp/sum(tmp))
      }
    }else{
      for (i in 1:n){
        tmp = sapply(1:ng, function(s){
          theta[s]*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)
        })
        res = rbind(res, tmp/sum(tmp))
      }
    }
  }else if (Obj$Model == "CNORM"){
    sigma = Obj$sigma
    if (Obj$Method == "L"){
      for (i in 1:n){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
        res = rbind(res, tmp/sum(tmp))
      }
    }else{
      for (i in 1:n){
        tmp = sapply(1:ng, function(s){
          theta[s]*gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
      res = rbind(res, tmp/sum(tmp))
      }
    }
  }else if (Obj$Model == "ZIP"){
    nu = Obj$nu
    if (Obj$Method == "L"){
      for (i in 1:n){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, k, nbeta, nnu, A, Y, TCOV, delta, nw)
          })
        res = rbind(res, tmp/sum(tmp))
      }
    }else{
      for (i in 1:n){
        tmp = sapply(1:ng, function(s){
          theta[s]*gkZIP(beta, nu, i, k, nbeta, nnu, A, Y, TCOV, delta, nw)
          })
        res = rbind(res, tmp/sum(tmp))
      }
    }
  }else{
    if (Obj$Method == "L"){
      for (i in 1:n){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkNL(beta, sigma, i, k, TCOV, A, Y)
        })
        res = rbind(res, tmp/sum(tmp))
      }
    }else{
      for (i in 1:n){
        tmp = sapply(1:ng, function(s){
          theta[s]*gkNL(beta, sigma, i, k, TCOV, A, Y)
        })
        res = rbind(res, tmp/sum(tmp))
      }
    }
  }
  colnames(res) = paste0("Gr", 1:ng)
  return(res)
}
#################################################################################
# Function to find theta in the calculus of the membership probability with predictors
#################################################################################
ftheta <- function(theta, taux, X, n, ng, period){
  a = 0
  for (i in 1:n){
    for (k in 1:ng){
      tmp = sapply(1:ng,function(s){theta[((s-1)*nx+1):(s*nx)]%*%X[i,]})
      a = a + taux[i,k]*(tmp[k]-log(sum(exp(tmp))))
    }
  }
  return(a)
}
difftheta <- function(theta, taux, X, n, ng, period){
  nx = ncol(X)
  thetas = c()
  for (k in 1:ng){
    for (l in 1:nx){
      a = 0
      for (i in 1:n){
        tmp = exp(sapply(1:ng,function(s){theta[((s-1)*nx+1):(s*nx)]%*%X[i,]}))
        a = a + X[i,l]*(taux[i,k]-tmp[k]/sum(tmp))
      }
      thetas = c(thetas, a)
    }
  }
  return(thetas)
}
findtheta <- function(theta, taux, X, n, ng, nx, period, EMIRLS, refgr){
  if (EMIRLS == TRUE){
    newtheta = c()
    thetaIRLS = theta[-c(((refgr-1)*nx+1):(nx*refgr))]
   thetaIRLS = thetaIRLS - theta[c(((refgr-1)*nx+1):(nx*refgr))]
    ind = 1:ng
    ind = ind[-refgr]
    precIRLS = 1
    while(any(abs(precIRLS)>10**(-6))){
      Xng = matrix(rep(0,n*(ng-1)*nx*(ng-1)),ncol=nx*(ng-1))
      tmp2 = c()
      PIw = c()
      tmp4 = c()
      kind = 0
      for (k in ind){
        kind = kind + 1
        PIwtmp = c()
        for (l in ind){
          tmp1 = c()
          if (k==l){
            for (i in 1:n){
              tmpPiik = piik(c(rep(0, nx), thetaIRLS), i, k, ng, X)
              tmp1 = c(tmp1, tmpPiik*(1-tmpPiik))
              tmp4 = c(tmp4, tmpPiik)
            }
          }else{
            for (i in 1:n){
              tmp1 = c(tmp1, -piik(c(rep(0, nx), thetaIRLS), i, k, ng, X)*piik(c(rep(0, nx), thetaIRLS), i, l, ng, X))
            }
          }
          PIwtmp = cbind(PIwtmp, diag(tmp1))
        }
        PIw = rbind(PIw, PIwtmp)
        Xng[((kind-1)*n+1):(kind*n), ((kind-1)*nx+1):(kind*nx)] = X
        tmp2 = c(tmp2, taux[,k])
      }
      rm(PIwtmp)
      Z = matrix(tmp2, ncol=1)
      PIm = matrix(tmp4, ncol=1)
      newthetaIRLS = as.vector(solve(t(Xng)%*%PIw%*%Xng, t(Xng)%*%(PIw%*%Xng%*%thetaIRLS+Z-PIm), tol=10**(-20)))
      precIRLS = c(thetaIRLS-newthetaIRLS)
      thetaIRLS= newthetaIRLS
    }
    newtheta = rep(0, ng*nx)
    newtheta[-c(((refgr-1)*nx+1):(nx*refgr))] = thetaIRLS
  }else{
    newtheta = optim(par = theta, fn = ftheta, gr= difftheta,
                     taux=taux, X=X, n=n, ng=ng, period=period,
                     control = list(fnscale=-1), hessian=FALSE)$par
  }
  return(newtheta)
}
