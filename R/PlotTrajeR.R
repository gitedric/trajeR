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
    cols1 = gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
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
    matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, ...)
  }else{
    Ygr = cbind(Y, rep(NA, n)) #for store the group membership
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    if (method == "L"){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
    }else{
      tmp = sapply(1:ng, function(s){
        theta[s]*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
    }
    ncol = which.max(tmp/sum(tmp))
    plot(A[1,],Y[1,],type="b", ylim = c(min(Y),max(Y)), pch=16, col = cols1[ncol], ...)
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      if (method == "L"){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
      }else{
        tmp = sapply(1:ng, function(s){
          theta[s]*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
      }
      # tmp = sapply(1:ng, function(s){
      #   piik(theta, i, s, ng, X)*gkCNORM(beta, sigma, i, s, nbeta, A, Y, ymin, ymax, TCOV, delta, nw)})
      ncol = which.max(tmp/sum(tmp))
      lines(A[i,],Y[i,],type="b", pch=16, col = cols1[ncol])
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
      #lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 4)
    }
    matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  if (is.matrix(plotcov) == FALSE){
    plotcov = matrix(plotcov, ncol = period)
  }
  if (any(plotcov != 0)){
    pas = 1:period
    for (k in 1:ng){
      for (row in 1:nrow(plotcov)){
        vec = sapply(1:length(pas), function(s){
          sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))+sum(delta[[k]]*plotcov[row,s]))
        })
        lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
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
    cols1 = gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
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
    matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = 'Values and predicted trajectories for all groups')
  }else{
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    tmp = sapply(1:ng, function(s){
      piik(theta, 1, s, ng, X)*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)})
    ncol = which.max(tmp/sum(tmp))
    a = dec * runif(n*period, 0, 1)*2*3.14159
    r = dec * runif(n*period, 0, 0.01)
    decx = r * cos(a)
    decy = r * sin(a)
    d = max(decx)
    dy = max(decy)
    xlim = c(min(A)-d, max(A)+d)
    plot(A[1,]+decx[1:period],Y[1,]+decy[1:period],type="b",
         ylim = c(-0.2, 1.2), xlim = xlim,
         pch=16, col = cols1[ncol],
         xlab = 'Time', ylab = 'Value', yaxt="n", main = 'Values and predicted trajectories for all groups')
    axis(side=2, at=seq(0,1,0.25), labels = seq(0,1,0.25))
    rect(xleft = xlim[1]-d, ybottom = 1-dy, xright = xlim[2]+d, ytop = 1+dy,
         col = "gray91", border = NA)
    rect(xleft = xlim[1]-d, ybottom = -dy, xright = xlim[2]+d, ytop = dy,
         col = "gray91", border = NA)
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkLogit(beta, i, s, nbeta, A, Y, TCOV, delta, nw)})
      ncol = which.max(tmp/sum(tmp))
      lines(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i],type="b", pch=16, col = adjustcolor(cols1[ncol], alpha.f = alpha))
      points(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i], pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1))))
      })
      vec = cbind(vec, tmp/(1+tmp))
    }
    matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
      }
    }
  }
  if (is.matrix(plotcov) == FALSE){
    plotcov = matrix(plotcov, ncol = period)
  }
  if (any(plotcov != 0)){
    pas = Time
    for (k in 1:ng){
      for (row in 1:nrow(plotcov)){
        tmp = sapply(1:length(pas), function(s){
          exp(sum(beta[[k]]*pas[s]**(0:(nbeta[k]-1)))+sum(delta[[k]]*plotcov[row,s]))
        })
        vec = tmp/(1+tmp)
        lines(pas, vec, type="l", pch=16,col = "white", lwd = 2)
        lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
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
    cols1 = gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
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
    matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = 'Values and predicted trajectories for all groups')
  }else{
    Ygr = cbind(Y, rep(NA, n)) #for store the group membership
    if (is.null(X)){
      X = cbind(rep(1, n))
    }
    if (method == "L"){
      tmp = sapply(1:ng, function(s){
        piik(theta, i, s, ng, X)*gkZIP(beta, nu, 1, s, nbeta, nnu, A, Y, TCOV, delta, nw)})
    }else{
      tmp = sapply(1:ng, function(s){
        theta[s]*gkZIP(beta, nu, 1, s, nbeta, nnu, A, Y, TCOV, delta, nw)})
    }
    ncol = which.max(tmp/sum(tmp))
    a = dec * runif(n*period, 0, 1)*2*3.14159
    r = dec * runif(n*period, 0, 0.01)
    decx = r * cos(a)
    decy = r * sin(a)
    d = max(decx)
    dy = max(decy)
    xlim = c(min(A), max(A))
    ylim = c(min(Y), max(Y))
    plot(A[1,]+decx[1:period],Y[1,]+decy[1:period],type="b",
         ylim = ylim, xlim = xlim,
         pch=16, col = cols1[ncol],
         xlab = 'Time', ylab = 'Value', yaxt="n", main = 'Values and predicted trajectories for all groups')
    Ygr[1,period+1] = ncol
    for (i in 2:n){
      if (method == "L"){
        tmp = sapply(1:ng, function(s){
          piik(theta, i, s, ng, X)*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)})
      }else{
        tmp = sapply(1:ng, function(s){
          theta[s]*gkZIP(beta, nu, i, s, nbeta, nnu, A, Y, TCOV, delta, nw)})
      }
      ncol = which.max(tmp/sum(tmp))
      lines(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i],type="b", pch=16, col = adjustcolor(cols1[ncol], alpha.f = alpha))
      points(A[i,]+decx[((i-1)*period+1):(i*period)],Y[i,]+decy[i], pch=16, col = cols1[ncol])
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
      vec = cbind(vec, tmp2 +(1-tmp2)*tmp)
    }
    matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
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
        vec = tmp2 +(1-tmp2)*tmp
        lines(pas, vec, type="l", pch=16,col = cols2[k], lwd = 2, lty=2)
      }
    }
  }
}
#################################################################################
plotTrajNL <- function(beta, sigma, theta, plotcov, Y, A, X, mean, alpha,
                                             ng, n, Time, col, degre, ymin, ymax, method, fct, main){
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
    cols1 = gray.colors(ng, start =0.3, end = 0.7, gamma = 2.2, alpha = 0.5)
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
    matplot(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1, xlab = 'Time', ylab = 'Value', main = main)
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
      lines(A[i,],Y[i,],type="b", pch=16, col = cols1[ncol])
      Ygr[i,period+1] = ncol
    }
    vec = c()
    for (k in 1:ng){
      tmp = sapply(1:length(pas), function(s){
        fct(pas[s], beta[[k]], TCOV)
      })
      vec = cbind(vec, tmp)
    }
    matlines(pas, vec, type="l", pch=16, col = cols2, lwd = 4, lty = 1)
    for (k in 1:ng){
      if (mean == TRUE){
        lines(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k], lty =2)
        points(A[1,], colMeans(Ygr[Ygr[,period+1] == k,-(period+1)]), col = cols2[k])
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
###################################################################################
# modification of plot's method for class trajectory.CNORM
####################################################################################
#' plot CNORM trajectory
#'
#'#' @param Obj an object of class "\code{Trajectory.LOGIT}".
#' @param plotcov an optionnal vector or matrix with the same length as the time period. Default value is NULL.
#' @param dec an optionnal real. It precise the shift to draw the data points.
#' @param col an optionnal vector. The vecotr of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default valme is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optionnal real. Indiciate the alpha channel of the points color.
#' @param ...
#'
#' @return a graphic.
#' @export
#'
#' @examples
#' plot(solL)
plot.Trajectory.CNORM <- function(Obj,  plotcov = NULL, col = "black", Y = NULL, A = NULL, X = NULL, mean = FALSE, alpha = 1,
                                   ...){
  plotTrajCensoredNormalLikelihood(beta = Obj$beta, sigma = Obj$sigma,  theta = Obj$theta, delta = Obj$delta, plotcov =plotcov,
                                   Y = Y, A = A, X = X, mean = mean, alpha = alpha,
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
#' @param plotcov an optionnal vector or matrix with the same length as the time period. Default value is NULL.
#' @param dec an optionnal real. It precise the shift to draw the data points.
#' @param col an optionnal vector. The vecotr of colors. It must contain a color for each trajectory and each points of groups.
#' Its length is the double of the number of group. Default valme is a grayscale.
#' @inheritParams trajeR
#' @param mean an optional logicial. Indicate if the mean of ech group and time value must be draw.
#' @param alpha on optionnal real. Indiciate the alpha channel of the points color.
#' @param ...
#'
#' @return a graphic.
#' @export
#'
#' @examples
#' plot(solL)
plot.Trajectory.LOGIT <- function(Obj, plotcov = NULL, dec = 1, col = "black", Y = NULL, A = NULL, X = NULL, mean = FALSE, alpha = 1, ...){
  plotTrajLOGIT(beta = Obj$beta, theta = Obj$theta, delta = Obj$delta, plotcov = plotcov,
                Y = Y, A = A, X = X, mean = mean, alpha = alpha,
                ng = Obj$groups, n = Obj$Size, Time = Obj$Time, dec = dec, col = col, degre = Obj$degre)
}
###################################################################################
# modification of plot's method for class trajectory.ZIP
####################################################################################
#' plot ZIP trajectory
#'
#' @export
#'
plot.Trajectory.ZIP <- function(Obj, plotcov = NULL, dec = 1, col = "black", Y = NULL, A = NULL, X = NULL, TCOV = NULL, mean = FALSE, alpha = 1, ...){
  plotTrajZIP(beta = Obj$beta, nu = Obj$nu, theta = Obj$theta, delta = Obj$delta, plotcov = plotcov,
              Y = Y, A = A, X = X, TCOV = TCOV, mean = mean, alpha = alpha,
              ng = Obj$groups, n = Obj$Size, Time = Obj$Time, dec = dec, col = col,
              degre = Obj$degre, method = Obj$Method, degre.nu = Obj$degre.nu)
}
###################################################################################
# modification of plot's method for class trajectory.NL
####################################################################################
#' plot Non Linear trajectory
#'
#' @export
#'
plot.Trajectory.NL <- function(Obj,  plotcov = NULL, col = "black", Y = NULL, A = NULL, X = NULL, mean = FALSE, alpha = 1, main = 'Values and predicted trajectories for all groups', ...){
  plotTrajNL(beta = Obj$beta, sigma = Obj$sigma,  theta = Obj$theta, plotcov =plotcov,
                                   Y = Y, A = A, X = X, mean = mean, alpha = alpha,
                                   ng = Obj$groups, n = Obj$Size, Time = Obj$Time, col = col, degre = Obj$degre,
                                   ymin = Obj$min, ymax = Obj$max, method = Obj$Method, fct = Obj$fct,
                                   main)
}
