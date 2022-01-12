###################################################################################
###################################################################################
# modification of plot's method for class trajectory
####################################################################################
###################################################################################
#################################################################################
# plot data
#################################################################################
plotTrajCensoredNormalLikelihood <- function(beta, sigma, theta, delta, plotcov, Y, A, X, mean, alpha,
                                             ng, n, Time, col, degre, ymin, ymax, method, ...){
  period = length(Time)
  if (is.null(plotcov)){
    plotcov = matrix(rep(0, period), nrow = 1)
  }
  nbeta = degre + 1
  ntheta = length(theta)/ng
  betatmp = beta
  j = 1
  beta = list()
  for (i in 1:length(nbeta)){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  if (any(is.na(delta))){
    nw = 0
    delta = rep(list(0), ng)
  }else{
    nw = length(delta)/ng
    deltatmp = delta
    delta = list()
    for (i in 1:ng){
      delta[[i]] = deltatmp[((i-1)*nw +1):(i*nw)]
    }
  }
  if (length(col) == 2*ng){
    cols1 = col[1:ng]
    cols2 = col[(ng+1):(2*ng)]
  }else{
    cols1 = grDevices::gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
    cols2 = rep("black",ng)
  }
  pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
  if (is.null(Y) | is.null(A)){
    vec = c()
    for (i in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        sum(beta[[i]]*pas[s]**(0:(nbeta[i]-1)))
      })
      vec = cbind(vec, tmp)
    }
    graphics::matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, ...)
  }else{
    Ygr = cbind(Y, rep(NA, n)) #for store the group membership
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    if (method == "L"){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV=NULL, delta=NULL, nw=0)})
    }else{
      tmp = sapply(1:ng, function(s){
        theta[s]*gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV=NULL, delta=NULL, nw=0)})
    }
    ncol = which.max(tmp/sum(tmp))
    plot(A[1,],Y[1,],type="b", ylim = c(min(Y, na.rm = T),max(Y, na.rm = T)), pch=16, col = cols1[ncol], ...)
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      if (method == "L"){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV=NULL, delta=NULL, nw=0)})
      }else{
        tmp = sapply(1:ng, function(s){
          theta[s]*gkCNORM_cpp(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV=NULL, delta=NULL, nw=0)})
      }
      ncol = which.max(tmp/sum(tmp))
      graphics::lines(A[i,],Y[i,],type="b", pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1)))
      })
      tmp[tmp<ymin] = ymin
      tmp[tmp>ymax] = ymax
      vec = cbind(vec, tmp)
    }
    graphics::matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        graphics::lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        graphics::points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  if (is.matrix(plotcov) == FALSE){
    plotcov = matrix(plotcov, ncol = period*nw, byrow = TRUE)
  }
  if (any(plotcov != 0)){
    pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
    for (k in 1:ng){
      for (row in 1:nrow(plotcov)){
        vec = sapply(1:length(pas), function(s){
          sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))+sum(delta[[k]]*plotcov[row, seq(from = floor(pas[s]), to = floor(pas[s]) + (nw - 1) * period, by = period)]))
        })
        graphics::lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
      }
    }
  }
}
##########################################################################################
# plot LOGIT function
##########################################################################################
plotTrajLOGIT <- function(beta, theta, delta, Y, A, X, ng, n, Time, dec, col, degre, plotcov, mean, alpha){
  period = length(Time)
  Ygr = cbind(Y, rep(NA, n)) #for store the group membership
  if (is.null(plotcov)){
    plotcov = matrix(rep(0, period), nrow = 1)
  }
  nbeta = degre + 1
  ntheta = length(theta)/ng
  betatmp = beta
  j = 1
  beta = list()
  if (any(is.na(delta))){
    nw = 0
    delta = rep(list(0), ng)
  }else{
    nw = length(delta)/ng
  }
  for (i in 1:length(nbeta)){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  if (length(col) == 2*ng){
    cols1 = col[1:ng]
    cols2 = col[(ng+1):(2*ng)]
  }else{
    cols1 = grDevices::gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
    cols2 = rep("black",ng)
  }
  pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
  if (is.null(Y) | is.null(A)){
    vec = c()
    for (i in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[i]]*pas[s]**(0:(nbeta[i]-1))))
      })
      vec = cbind(vec, tmp/(1+tmp))
    }
    graphics::matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = 'Values and predicted trajectories for all groups')
  }else{
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    tmp = sapply(1:ng, function(s){
      piik(theta, 1, s, ng, X)*gkLOGIT_cpp(beta, i-1, s-1, nbeta, A, Y, TCOV=NULL, delta=NULL, nw=0)})
    ncol = which.max(tmp/sum(tmp))
    a = dec * stats::runif(n*period, 0, 1)*2*3.14159
    r = dec * stats::runif(n*period, 0, 0.01)
    decx = r * cos(a)
    decy = r * sin(a)
    d = max(decx)
    dy = max(decy)
    xlim = c(min(A)-d, max(A)+d)
    plot(A[1,]+decx[1:period],Y[1,]+decy[1:period],type="b",
         ylim = c(-0.2, 1.2), xlim = xlim,
         pch=16, col = cols1[ncol],
         xlab = 'Time', ylab = 'Value', yaxt="n", main = 'Values and predicted trajectories for all groups')
    graphics::axis(side=2, at=seq(0,1,0.25), labels = seq(0,1,0.25))
    graphics::rect(xleft = xlim[1]-d, ybottom = 1-dy, xright = xlim[2]+d, ytop = 1+dy,
         col = "gray91", border = NA)
    graphics::rect(xleft = xlim[1]-d, ybottom = -dy, xright = xlim[2]+d, ytop = dy,
         col = "gray91", border = NA)
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkLOGIT_cpp(beta, i-1, s-1, nbeta, A, Y, TCOV=NULL, delta=NULL, nw=0)})
      ncol = which.max(tmp/sum(tmp))
      graphics::lines(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i],type="b", pch=16, col = grDevices::adjustcolor(cols1[ncol], alpha.f = alpha))
      graphics::points(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i], pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))))
      })
      vec = cbind(vec, tmp/(1+tmp))
    }
    graphics::matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        graphics::lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        graphics::points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  if (is.matrix(plotcov) == FALSE){
    plotcov = matrix(plotcov, ncol = period)
  }
  if (any(plotcov != 0)){
    pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
    for (k in 1:ng){
      for (row in 1:nrow(plotcov)){
        vec = sapply(1:length(pas), function(s){
          exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))+sum(delta[[k]]*plotcov[row,floor(pas[s])])))
        })
        vec = vec/(1+vec)
        graphics::lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
      }
    }
  }
}
##########################################################################################
# plot POIS function
##########################################################################################
plotTrajPOIS <- function(beta, theta, delta, Y, A, X, TCOV, ng, n, Time, dec, col, degre, plotcov, mean, alpha, method){
  period = length(Time)
  Ygr = cbind(Y, rep(NA, n)) #for store the group membership
  if (is.null(plotcov)){
    plotcov = matrix(rep(0, period), nrow = 1)
  }
  nbeta = degre + 1
  ntheta = length(theta)/ng
  betatmp = beta
  j = 1
  beta = list()
  if (any(is.na(delta))){
    nw = 0
    delta = rep(list(0), ng)
  }else{
    nw = length(delta)/ng
  }
  for (i in 1:length(nbeta)){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  if (length(col) == 2*ng){
    cols1 = col[1:ng]
    cols2 = col[(ng+1):(2*ng)]
  }else{
    cols1 = grDevices::gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
    cols2 = rep("black",ng)
  }
  pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
  if (is.null(Y) | is.null(A)){
    vec = c()
    for (i in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[i]]*pas[s]**(0:(nbeta[i]-1))))
      })
      vec = cbind(vec, tmp)
    }
    graphics::matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = 'Values and predicted trajectories for all groups')
  }else{
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    tmp = sapply(1:ng, function(s){
      piik(theta, 1, s, ng, X)*gkPois_cpp(beta, i-1, s-1, nbeta, A, Y, TCOV=NULL, delta=NULL, nw=0)})
    ncol = which.max(tmp/sum(tmp))
    a = dec * stats::runif(n*period, 0, 1)*2*3.14159
    r = dec * stats::runif(n*period, 0, 0.01)
    decx = r * cos(a)
    decy = r * sin(a)
    d = max(decx)
    dy = max(decy)
    xlim = c(min(A), max(A))
    ylim = c(min(Y), max(Y))
    plot(A[1,]+decx[1:period],Y[1,]+decy[1:period],type="b",
         ylim = c(min(Y, na.rm = T),max(Y, na.rm = T)),
         pch=16, col = cols1[ncol],
         xlab = 'Time', ylab = 'Value', yaxt="n", main = 'Values and predicted trajectories for all groups')
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      if (method == "L"){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkPois_cpp(beta, i-1, s-1, nbeta, A, Y, TCOV, delta, nw)})
      }else{
        tmp = sapply(1:ng, function(s){
          theta[s]*gkPois_cpp(beta, i-1, s-1, nbeta, A, Y, TCOV, delta, nw)})
      }
      ncol = which.max(tmp/sum(tmp))
      graphics::lines(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i],type="b", pch=16, col = grDevices::adjustcolor(cols1[ncol], alpha.f = alpha))
      graphics::points(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i], pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))))
      })
      vec = cbind(vec, tmp)
    }
    graphics::matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        graphics::lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        graphics::points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  if (is.matrix(plotcov) == FALSE){
    plotcov = matrix(plotcov, ncol = period)
  }
  if (any(plotcov != 0)){
    pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
    for (k in 1:ng){
      for (row in 1:nrow(plotcov)){
        vec = sapply(1:length(pas), function(s){
          exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))+sum(delta[[k]]*plotcov[row,floor(pas[s])])))
        })
        graphics::lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
      }
    }
  }
}
##########################################################################################
# plot ZIP function
##########################################################################################
plotTrajZIP <- function(beta, nu, theta, delta, Y, A, X, TCOV, ng, n, Time, dec, col, degre, plotcov, mean, alpha, method, degre.nu){
  period = length(Time)
  if (is.null(plotcov)){
    plotcov = matrix(rep(0, period), nrow = 1)
  }
  nbeta = degre + 1
  nnu = degre.nu + 1
  ntheta = length(theta)/ng
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
  if (any(is.na(delta))){
    nw = 0
    delta = rep(list(0), ng)
  }else{
    nw = length(delta)/ng
    deltatmp =delta
    ndeltacum = cumsum(c(0, rep(nw, ng)))
    delta = list()
    for (i in 1:ng){
      delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
    }
  }
  if (length(col) == 2*ng){
    cols1 = col[1:ng]
    cols2 = col[(ng+1):(2*ng)]
  }else{
    cols1 = grDevices::gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
    cols2 = rep("black",ng)
  }
  pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
  if (is.null(Y) | is.null(A)){
    vec = c()
    for (i in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[i]]*pas[s]**(0:(nbeta[i]-1))))
      })
      vec = cbind(vec, tmp/(1+tmp))
    }
    graphics::matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = 'Values and predicted trajectories for all groups')
  }else{
    Ygr = cbind(Y, rep(NA, n)) #for store the group membership
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    if (method == "L"){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkZIP_cpp(beta, nu, 0, s-1, nbeta, nnu, A, Y, TCOV, delta, nw)})
    }else{
      tmp = sapply(1:ng, function(s){
        theta[s]*gkZIP_cpp(beta, nu, 0, s-1, nbeta, nnu, A, Y, TCOV, delta, nw)})
    }
    ncol = which.max(tmp/sum(tmp))
    a = dec * stats::runif(n*period, 0, 1)*2*3.14159
    r = dec * stats::runif(n*period, 0, 0.01)
    decx = r * cos(a)
    decy = r * sin(a)
    d = max(decx)
    dy = max(decy)
    xlim = c(min(A), max(A))
    ylim = c(min(Y), max(Y))
    plot(A[1,]+decx[1:period],Y[1,]+decy[1:period],type="b",
         ylim = c(min(Y, na.rm = T),max(Y, na.rm = T)),
         pch=16, col = cols1[ncol],
         xlab = 'Time', ylab = 'Value', yaxt="n", main = 'Values and predicted trajectories for all groups')
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      if (method == "L"){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkZIP_cpp(beta, nu, i-1, s-1, nbeta, nnu, A, Y, TCOV, delta, nw)})
      }else{
        tmp = sapply(1:ng, function(s){
          theta[s]*gkZIP_cpp(beta, nu, i-1, s-1, nbeta, nnu, A, Y, TCOV, delta, nw)})
      }
      ncol = which.max(tmp/sum(tmp))
      graphics::lines(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i],type="b", pch=16, col = grDevices::adjustcolor(cols1[ncol], alpha.f = alpha))
      graphics::points(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i], pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))))
      })
      tmp2 = sapply(1:length(pas), function(s){
        exp(sum(nu[[k]]*pas[s]**(0:(nnu[k]-1))))/(1+exp(sum(nu[[k]]*pas[s]**(0:(nnu[k]-1)))))
      })
      # plot average of the ZIP distirbution
      vec = cbind(vec, (1-tmp2)*tmp)
      #vec = cbind(vec, tmp2 +(1-tmp2)*tmp)
    }
    graphics::matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        graphics::lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        graphics::points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  if (is.matrix(plotcov) == FALSE){
    plotcov = matrix(plotcov, ncol = period, byrow=TRUE)
  }
  if (any(plotcov != 0)){
    pas = 1:period
    for (k in 1:ng){
      for (row in 1:nrow(plotcov)){
        tmp = sapply(1:length(pas), function(s){
          exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))) +sum(delta[[k]]*plotcov[row,s]))
        })
        tmp2 = sapply(1:length(pas), function(s){
          exp(sum(nu[[k]]*pas[s]**(0:(nnu[k]-1))))/(1+exp(sum(nu[[k]]*pas[s]**(0:(nnu[k]-1)))))
        })
        # plot average of the ZIP distirbution
        vec = (1-tmp2)*tmp
        graphics::lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
      }
    }
  }
}
#################################################################################
plotTrajNL <- function(beta, sigma, theta, plotcov, Y, A, X, mean, alpha,
                                             ng, n, Time, col, degre, ymin, ymax, method, fct, main, TCOV){
  period = length(Time)
  if (is.null(plotcov)){
    plotcov = matrix(rep(0, period), nrow = 1)
  }
  nbeta = degre + 1
  ntheta = length(theta)/ng
  betatmp = beta
  j = 1
  beta = list()
  nw =0
  for (i in 1:length(nbeta)){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  if (length(col) == 2*ng){
    cols1 = col[1:ng]
    cols2 = col[(ng+1):(2*ng)]
  }else{
    cols1 = grDevices::gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
    cols2 = rep("black",ng)
  }
  pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
  if (is.null(Y) | is.null(A)){
    vec = c()
    for (i in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        fct(pas[s], beta[[i]], TCOV)
      })
      vec = cbind(vec, tmp)
    }
    graphics::matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = main)
  }else{
    Ygr = cbind(Y, rep(NA, n)) #for store the group membership
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    if (method == "L"){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkNL(beta, sigma, i, s, TCOV, A, Y, fct)})
    }else{
      tmp = sapply(1:ng, function(s){
        theta[s]*gkNL(beta, sigma, i, s, TCOV, A, Y, fct)})
    }
    ncol = which.max(tmp/sum(tmp))
    plot(A[1,],Y[1,],type="b", ylim = c(min(Y),max(Y)), pch=16, col = cols1[ncol],
         xlab = 'Time', ylab = 'Value', main = main)
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      if (method == "L"){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkNL(beta, sigma, i, s, TCOV, A, Y, fct)})
      }else{
        tmp = sapply(1:ng, function(s){
          theta[s]*gkNL(beta, sigma, i, s, TCOV, A, Y, fct)})
      }
      ncol = which.max(tmp/sum(tmp))
      graphics::lines(A[i,],Y[i,],type="b", pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        fct(pas[s], beta[[k]], TCOV)
      })
      vec = cbind(vec, tmp)
    }
    graphics::matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        graphics::lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        graphics::points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  # if (is.matrix(plotcov) == FALSE){
  #   plotcov = matrix(plotcov, ncol = period)
  # }
  # if (any(plotcov != 0)){
  #   pas = 1:period
  #   for (k in 1:ng){
  #     for (row in 1:nrow(plotcov)){
  #       vec = sapply(1:length(pas), function(s){
  #         sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))+sum(delta[[k]]*plotcov[row,s]))
  #       })
  #       lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
  #     }
  #   }
  # }
}
##########################################################################################
# plot BETA function
##########################################################################################
plotTrajBETA <- function(beta, phi, theta, delta, Y, A, X, TCOV, ng, n, Time, dec, col, degre, plotcov, mean, alpha, method, degre.phi){
  period = length(Time)
  if (is.null(plotcov)){
    plotcov = matrix(rep(0, period), nrow = 1)
  }
  nbeta = degre + 1
  nphi = degre.phi + 1
  ntheta = length(theta)/ng
  j = 1
  betatmp = beta
  beta = list()
  for (i in 1:ng){
    beta[[i]] = betatmp[j:sum(nbeta[1:i])]
    j = sum(nbeta[1:i]) + 1
  }
  j = 1
  phitmp = phi
  phi = list()
  for (i in 1:ng){
    phi[[i]] = phitmp[j:sum(nphi[1:i])]
    j = sum(nphi[1:i]) + 1
  }
  if (any(is.na(delta))){
    nw = 0
    delta = rep(list(0), ng)
  }else{
    nw = length(delta)/ng
    deltatmp =delta
    ndeltacum = cumsum(c(0, rep(nw, ng)))
    delta = list()
    for (i in 1:ng){
      delta[[i]] = deltatmp[(ndeltacum[i]+1):(ndeltacum[i+1])]
    }
  }
  if (length(col) == 2*ng){
    cols1 = col[1:ng]
    cols2 = col[(ng+1):(2*ng)]
  }else{
    cols1 = grDevices::gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
    cols2 = rep("black",ng)
  }
  pas = seq(Time[1],Time[period],(Time[period]-Time[1])/100)
  if (is.null(Y) | is.null(A)){
    vec = c()
    for (i in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[i]]*pas[s]**(0:(nbeta[i]-1))))
      })
      vec = cbind(vec, tmp/(1+tmp))
    }
    graphics::matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    #graphics::matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = 'Values and predicted trajectories for all groups')
  }else{
    Ygr = cbind(Y, rep(NA, n)) #for store the group membership), max(Y))
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    if (method == "L"){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkBETA_cpp(beta, phi, 0, s-1, nbeta, nphi, A, Y, TCOV, delta, nw)})
    }else{
      tmp = sapply(1:ng, function(s){
        theta[s]*gkBETA_cpp(beta, phi, 0, s-1, nbeta, nphi, A, Y, TCOV, delta, nw)})
    }
    ncol = which.max(tmp/sum(tmp))
    xlim = c(min(A), max(A))
    ylim = c(min(Y, na.rm = TRUE), max(Y, na.rm = TRUE))
    plot(A[1,], Y[1,],type="b",
         ylim = c(min(Y, na.rm = T),max(Y, na.rm = T)),
         pch=16, col = cols1[ncol],
         xlab = 'Time', ylab = 'Value', main = 'Values and predicted trajectories for all groups')
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      if (method == "L"){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkBETA_cpp(beta, phi, i-1, s-1, nbeta, nphi, A, Y, TCOV, delta, nw)})
      }else{
        tmp = sapply(1:ng, function(s){
          theta[s]*gkBETA_cpp(beta, phi, i-1, s-1, nbeta, nphi, A, Y, TCOV, delta, nw)})
      }
      ncol = which.max(tmp/sum(tmp))
      graphics::lines(A[i,], Y[i,],type="b", pch=16, col = grDevices::adjustcolor(cols1[ncol], alpha.f = alpha))
      graphics::points(A[i,],Y[i,], pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))))
      })
      # plot average of the ZIP distirbution
      vec = cbind(vec, tmp/(1+tmp))
      #vec = cbind(vec, tmp2 +(1-tmp2)*tmp)
    }
    graphics::matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        graphics::lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        graphics::points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  if (is.matrix(plotcov) == FALSE){
    plotcov = matrix(plotcov, ncol = period, byrow=TRUE)
  }
  if (any(plotcov != 0)){
    pas = 1:period
    for (k in 1:ng){
      for (row in 1:nrow(plotcov)){
        tmp = sapply(1:length(pas), function(s){
          exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))) +sum(delta[[k]]*plotcov[row,s]))
        })
        # plot average of the ZIP distirbution
        vec = tmp/(1+tmp)
        graphics::lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
      }
    }
  }
}
###################################################################################
# Generic function to plot the trajectory
####################################################################################
#' plot trajectory
#'
#' @param Obj an object of class "\code{Trajectory}".
#' @param ... optional parameters
#'
#' @return a graphic.
#' @export
#'
#' @examples
#' data = read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data = as.matrix(data)
#' sol = trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2,2), Model = "CNORM", Method = "EM")
#' plotrajeR(sol)

plotrajeR <- function(Obj, ...){
  UseMethod("plotrajeR", Obj)
}

###################################################################################
# modification of plot's method for class trajectory.CNORM
####################################################################################
#' plot CNORM trajectory
#'
#' @param Obj an object of class "\code{Trajectory.CNORM}".
#' @param plotcov an optionnal vector or matrix with the same length as the time period. Default value is NULL.
#' @param col an optionnal vector. The vecotr of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default valme is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optionnal real. Indiciate the alpha channel of the points color.
#' @param ... optional parameters
#'
#' @return a graphic.
#' @export
#'
plotrajeR.Trajectory.CNORM <- function(Obj,  plotcov = NULL, col = "black", Y = NULL, A = NULL, Risk = NULL, mean = FALSE, alpha = 1,
                                   ...){
  plotTrajCensoredNormalLikelihood(beta = Obj$beta, sigma = Obj$sigma,  theta = Obj$theta, delta = Obj$delta, plotcov =plotcov,
                                   Y = Y, A = A, X = Risk, mean = mean, alpha = alpha,
                                   ng = Obj$groups, n = Obj$Size, Time = Obj$Time, col = col, degre = Obj$degre,
                                   ymin = Obj$min, ymax = Obj$max, method = Obj$Method,
                                   ...)
}
###################################################################################
# modification of plot's method for class trajectory.LOGIT
####################################################################################
#' plot LOGIT trajectory
#'
#' @param Obj an object of class "\code{Trajectory.LOGIT}".
#' @param plotcov an optional vector or matrix with the same length as the time period. Default value is NULL.
#' @param dec an optional real. It precise the shift to draw the data points.
#' @param col an optional vector. The vector of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default value is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optional real. Indicate the alpha channel of the points color.
#' @param ... optional parameters
#'
#' @return a graphic.
#' @export
#'

plotrajeR.Trajectory.LOGIT <- function(Obj, plotcov = NULL, dec = 1, col = "black", Y = NULL, A = NULL, Risk = NULL, mean = FALSE, alpha = 1, ...){
  plotTrajLOGIT(beta = Obj$beta, theta = Obj$theta, delta = Obj$delta, plotcov = plotcov,
                Y = Y, A = A, X = Risk, mean = mean, alpha = alpha,
                ng = Obj$groups, n = Obj$Size, Time = Obj$Time, dec = dec, col = col, degre = Obj$degre)
}
###################################################################################
# modification of plot's method for class trajectory.POIS
####################################################################################
#' plot POIS trajectory
#'
#' @param Obj an object of class "\code{Trajectory.POIS}".
#' @param plotcov an optional vector or matrix with the same length as the time period. Default value is NULL.
#' @param dec an optional real. It precise the shift to draw the data points.
#' @param col an optional vector. The vector of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default value is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optional real. Indicate the alpha channel of the points color.
#' @param ... optional parameters
#'
#' @return a graphic.
#' @export
#'

plotrajeR.Trajectory.POIS <- function(Obj, plotcov = NULL, dec = 0, col = "black", Y = NULL, A = NULL, Risk = NULL,  TCOV = NULL, mean = FALSE, alpha = 1, ...){
  plotTrajPOIS(beta = Obj$beta, theta = Obj$theta, delta = Obj$delta, plotcov = plotcov,
                Y = Y, A = A, X = Risk, TCOV = TCOV, mean = mean, alpha = alpha, method = Obj$Method,
                ng = Obj$groups, n = Obj$Size, Time = Obj$Time, dec = dec, col = col, degre = Obj$degre)
}
###################################################################################
# modification of plot's method for class trajectory.ZIP
####################################################################################
#' plot ZIP trajectory
#'
#' @param Obj an object of class "\code{Trajectory.LOGIT}".
#' @param plotcov an optional vector or matrix with the same length as the time period. Default value is NULL.
#' @param dec an optional real. It precise the shift to draw the data points.
#' @param col an optional vector. The vector of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default value is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optional real. Indicate the alpha channel of the points color.
#' @param ... optional parameters
#'
#' @return a graphic.
#' @export
#'
plotrajeR.Trajectory.ZIP <- function(Obj, plotcov = NULL, dec = 1, col = "black", Y = NULL, A = NULL, Risk = NULL, TCOV = NULL, mean = FALSE, alpha = 1, ...){
  plotTrajZIP(beta = Obj$beta, nu = Obj$nu, theta = Obj$theta, delta = Obj$delta, plotcov = plotcov,
              Y = Y, A = A, X = Risk, TCOV = TCOV, mean = mean, alpha = alpha,
              ng = Obj$groups, n = Obj$Size, Time = Obj$Time, dec = dec, col = col,
              degre = Obj$degre, method = Obj$Method, degre.nu = Obj$degre.nu)
}
###################################################################################
# modification of plot's method for class trajectory.NL
####################################################################################
#' plot Non Linear trajectory
#'
#' @param Obj an object of class "\code{Trajectory.LOGIT}".
#' @param plotcov an optional vector or matrix with the same length as the time period. Default value is NULL.
#' @param col an optional vector. The vector of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default value is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optional real. Indicate the alpha channel of the points color.
#' @param ... optional parameters
#'
#' @return a graphic.
#' @export
#'
plotrajeR.Trajectory.NL <- function(Obj,  plotcov = NULL, col = "black", Y = NULL, A = NULL, Risk = NULL, mean = FALSE, alpha = 1,  TCOV = NULL, ...){
  plotTrajNL(beta = Obj$beta, sigma = Obj$sigma,  theta = Obj$theta, plotcov =plotcov,
                                   Y = Y, A = A, X = Risk, mean = mean, alpha = alpha,
                                   ng = Obj$groups, n = Obj$Size, Time = Obj$Time, col = col, degre = Obj$degre,
                                   ymin = Obj$min, ymax = Obj$max, method = Obj$Method, fct = Obj$fct,
                                   TCOV = TCOV)
}
###################################################################################
# modification of plot's method for class trajectory.BETA
####################################################################################
#' plot BETA trajectory
#'
#' @param Obj an object of class "\code{Trajectory.LOGIT}".
#' @param plotcov an optional vector or matrix with the same length as the time period. Default value is NULL.
#' @param col an optional vector. The vector of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default value is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optional real. Indicate the alpha channel of the points color.
#' @param ... optional parameters
#'
#' @return a graphic.
#' @export
#'
plotrajeR.Trajectory.BETA <- function(Obj, plotcov = NULL, col = "black", Y = NULL, A = NULL, Risk = NULL, TCOV = NULL, mean = FALSE, alpha = 1, ...){
  plotTrajBETA(beta = Obj$beta, phi = Obj$phi, theta = Obj$theta, delta = Obj$delta, plotcov = plotcov,
              Y = Y, A = A, X = Risk, TCOV = TCOV, mean = mean, alpha = alpha,
              ng = Obj$groups, n = Obj$Size, Time = Obj$Time, dec = 1, col = col,
              degre = Obj$degre, method = Obj$Method, degre.phi = Obj$degre.phi)
}