###################################################################################
# modification of print's method for class trajectory
####################################################################################
####################################################################################
####################################################################################
# general functions
####################################################################################
Row <- function(x, n, pad = 1) {
  foo <- function(i, x, n) {
    fmt <- paste0("%", n[i], "s")
    sprintf(fmt, as.character(x[i]))
  }
  rowc <- sapply(seq_along(x), foo, x = x, n = n)
  paste0(" ", paste(paste0(rep(" ", pad), rowc, rep(" ", pad)),
                    collapse = " "),
         " ")
}
SepLine1 <- function(n, pad = 1) {
  tmp <- lapply(n, function(x, pad) paste(rep("-", x + (2* pad)),
                                          collapse = ""),
                pad = pad)
  paste0("-", paste(tmp, collapse = "-"), "-")
}
SepLine2 <- function(n, pad = 1) {
  tmp <- lapply(n, function(x, pad) paste(rep(" ", x + (2* pad)),
                                          collapse = ""),
                pad = pad)
  paste0(" ", paste(tmp, collapse = " "), " ")
}
###################################################################################
# modification of print's method for class trajectory.CNORM
####################################################################################
#' Print CNORM
#'
#' Print mehtod for an object of class "\code{Trajectory.CNORM}".
#'
#' @param Obj an object of class "\code{Trajectory.CNORM}".
#' @param ...
#'
#' @return  The print of Obj.
#' @export
#'
#' @examples
#' load(file='data/CNORM_data01')
#' solL
print.Trajectory.CNORM <- function(Obj, ...){
  # definiton of different sizes
  n= Obj$Size
  ng = Obj$groups
  nbeta = Obj$degre + 1
  ntheta = length(Obj$theta)/ng
  if (any(is.na(Obj$delta))){
    ndelta = 0
  }else{
    ndelta = length(Obj$delta)/ng
  }
  nbetatmp = c(0, nbeta)
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(ndelta, ng)))
  indbeta = 1:sum(nbeta)
  indsigma = (sum(nbeta)+1):(sum(nbeta)+ng)
  if (any(is.na(Obj$delta))){
    inddelta = 0
  }else{
    inddelta = (sum(nbeta)+ng+1):(sum(nbeta)+ng+ndelta*ng)
  }
  indtheta = (sum(nbeta) +ng+ ndelta*ng + 1):(sum(nbeta) +ng+ ndelta*ng + ntheta*ng)
  #for the changement of group in the final table
  indcum=c(0, nbetacum[-1]+(1:ng)*ndelta)
  # adding number of groups
  group = rep("", sum(nbeta)+ng+ndelta*ng + ntheta*ng)
  #number for beta
  for (i in 1:ng){
    group[nbetacum[i]+1] = as.character(i)
  }
  # number for sigma
  for (i in 1:ng){
    group[nbetacum[ng+1]+i] = as.character(i)
  }
  # number for theta
  for (i in 1:(ng)){
    group[indcum[ng+1]+ng+1+(i-1)*ntheta] = as.character(i)
  }
  # create a data frame for print
  dfprint = cbind(group, Parameter = Obj$Names, Obj$tab)
  dfprint[,2] = as.character(dfprint[,2])
  dfprint[,1] = as.character(dfprint[,1])
  if (ntheta == 1){
    if (Obj$Method == "L"){
      dfprint[indtheta, 3] = exp(dfprint[indtheta, 3])/sum(exp(dfprint[indtheta, 3]))
    }
    dfprint[indtheta, 2] = paste0("pi", 1:ng)
  }else{
    # we store the value of the baseline
    dfprint[indtheta,3] = dfprint[indtheta,3] - dfprint[indtheta[1:ntheta],3]
    dfprint = dfprint[-indtheta[2:ntheta],]
    # we change the name of group 1 to baseline
    dfprint[indcum[ng+1]+ng+1,2] = "Baseline"
  }
  # we reorganize the data frame by group ordering
  tmp =c()
  for (i in 1:ng){
    tmp = rbind(tmp, dfprint[(nbetacum[i]+1):(nbetacum[i+1]),])
    if (ndelta !=0){
      tmp = rbind(tmp, dfprint[(sum(nbeta)+ng+ndeltacum[i]+1):(sum(nbeta)+ng+ndeltacum[i+1]),])
    }
  }
  for (i in 1:ng){
    tmp = rbind(tmp, dfprint[sum(nbeta) + i,])
  }
  if (ntheta == 1){
    tmp = rbind(tmp, dfprint[indtheta,])
  }else{
    tmp = rbind(tmp, dfprint[indtheta[1:(length(indtheta)-ntheta+1)],])
  }
  dfprint = tmp
  dfprint[,3] = round(as.numeric(dfprint[,3]), 5)
  dfprint[,5] = round(as.numeric(dfprint[,5]), 5)
  dfprint[,6] = round(as.numeric(dfprint[,6]), 5)
  dfprint[,4] = round(as.numeric(dfprint[,4]), 5)
  fObj  = format(dfprint)
  strings = apply( dfprint,2,  function(x) unlist(dfprint))[1,]
  widths <- nchar(strings)
  names = c("Group" ,"Parameter", "Estimate", "Std. Error", "T for H0:", "Prob>|T|")
  widths <- pmax(nchar(strings), nchar(names))
  csum <- sum(widths + 1) - 1
  cat(paste("Call TrajeR with"), Obj$groups, "groups and a", paste(Obj$degre, collapse=","),
      "degrees of polynomial shape of trajectory.\n")
  cat("Model : Censored Normal\n")
  if (Obj$Method == "L"){
    cat("Method : Likelihood \n \n")
  } else if (Obj$Method == "EM"){
    cat("Method : Expectation-maximization \n \n")
  }else{
    cat("Method : Expectation-maximization with IWRLS\n \n")
  }
  esp = 1
  sep1 <- SepLine1(widths, pad = esp)
  sep2 <- SepLine2(widths, pad = esp)
  namestmp = colnames(dfprint)
  namestmp[5] = "T for H0:"
  writeLines(Row(namestmp, widths, esp))
  writeLines(Row(c("","","","","param.=0",""), widths, esp))
  writeLines(sep1)
  # write beta and delta
  for (i in 1:(ng-1)){
    for (j in 1:(nbeta[i]+ndelta)){
      writeLines(Row(dfprint[indcum[i]+j , ], widths, esp))
    }
    writeLines(sep2)
  }
  for (j in 1:(nbeta[ng]+ndelta)){
    writeLines(Row(dfprint[indcum[ng]+j, ], widths, esp))
  }
  writeLines(sep1)
  # write sigma
  for (i in 1:(ng-1)){
    writeLines(Row(dfprint[sum(nbeta)+ndelta*ng + i, ], widths, esp))
  }
  writeLines(Row(dfprint[sum(nbeta)+ndelta*ng + ng, ], widths, esp))
  writeLines(sep1)
  # write theta or pi
  if (ntheta == 1){
    for (i in 1:ng){
      if (ndelta !=0){
        writeLines(Row(dfprint[sum(nbeta)+ng+ndelta*ng+i, ], widths, esp))
      }else{
        writeLines(Row(dfprint[sum(nbeta)+ng+ndelta*ng+i, ], widths, esp))
      }
    }
  }else{
    writeLines(Row(dfprint[indcum[ng+1]+ng+1, ], widths, esp))
    writeLines(sep2)
    if (ng>2){
        for (j in 1:((ng-2)*ntheta)){
          writeLines(Row(dfprint[indcum[ng+1]+ng+1+j, ], widths, esp))
        }
        writeLines(sep2)
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+ng+(ng-2)*ntheta+1+j, ], widths, esp))
      }
    }else{
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+ng+1+j, ], widths, esp))
      }
    }
  }
  writeLines(sep1)
  cat("Likelihood :", Obj$Likelihood)
}
###################################################################################
# modification of print's method for class trajectory.ZIP
####################################################################################
#' Print ZIP
#'
#' Print method for an object of class "\code{Trajectory.ZIP}".
#' 
#' @export
#' 
#' @return  The print of Obj.
print.Trajectory.ZIP <- function(Obj, ...){
  # definiton of different sizes
  n= Obj$Size
  ng = Obj$groups
  nbeta = Obj$degre + 1
  ntheta = length(Obj$theta)/ng
  nnu = Obj$degre.nu + 1
  if (any(is.na(Obj$delta))){
    ndelta = 0
  }else{
    ndelta = length(Obj$delta)/ng
  }
  nbetatmp = c(0, nbeta)
  nbetacum = cumsum(c(0, nbeta))
  nnucum = cumsum(c(0, nnu))
  ndeltacum = cumsum(c(0, rep(ndelta, ng)))
  indbeta = 1:sum(nbeta)
  indnu = (sum(nbeta)+1):(sum(nbeta)+sum(nnu))
  if (any(is.na(Obj$delta))){
    inddelta = 0
  }else{
    inddelta = (sum(nbeta)+sum(nnu)+1):(sum(nbeta)+sum(nnu)+ndelta*ng)
  }
  indtheta = (sum(nbeta) + sum(nnu) + ndelta*ng + 1):(sum(nbeta) + sum(nnu) + ndelta*ng + ntheta*ng)
  #for the changement of group in the final table
  indcum=c(0, nbetacum[-1]+nnucum[-1]+(1:ng)*ndelta)
  # adding number of groups
  group = rep("", sum(nbeta)+sum(nnu)+ndelta*ng + ntheta*ng)
  for (i in 1:ng){
    group[nbetacum[i]+1] = as.character(i)
    #group[nbetacum[ng+1]+nnucum[i]+1] = as.character(i)
    group[indcum[ng+1]+1+(i-1)*ntheta] = as.character(i)
  }
  # create a data frame for print
  dfprint = cbind(group, Parameter = Obj$Names, Obj$tab)
  dfprint[,2] = as.character(dfprint[,2])
  dfprint[,1] = as.character(dfprint[,1])
  if (ntheta == 1){
    if (Obj$Method == "L"){
      dfprint[indtheta, 3] = exp(dfprint[indtheta, 3])/sum(exp(dfprint[indtheta, 3]))
    }
    dfprint[indtheta, 2] = paste0("pi", 1:ng)
  }else{
    # we store the value of the baseline
    dfprint[indtheta,3] = dfprint[indtheta,3] - dfprint[indtheta[1:ntheta],3]
    dfprint = dfprint[-indtheta[2:ntheta],]
    # we change the name of group 1 to baseline
    dfprint[indcum[ng+1]+1,2] = "Baseline"
  }
  # we reorganize the data frame by group ordering
  tmp =c()
  for (i in 1:ng){
    tmp = rbind(tmp, dfprint[(nbetacum[i]+1):(nbetacum[i+1]),])
    tmp = rbind(tmp,  dfprint[(nbetacum[ng+1]+nnucum[i]+1):(nbetacum[ng+1]+nnucum[i+1]),])
    if (ndelta !=0){
      tmp = rbind(tmp, dfprint[(sum(nbeta)+sum(nnu)+ndeltacum[i]+1):(sum(nbeta)+sum(nnu)+ndeltacum[i+1]),])
    }
  }
  if (ntheta == 1){
    tmp = rbind(tmp, dfprint[indtheta,])
  }else{
    tmp = rbind(tmp, dfprint[indtheta[1:(length(indtheta)-ntheta+1)],])
  }
  dfprint = tmp
  dfprint[,3] = round(as.numeric(dfprint[,3]), 5)
  dfprint[,5] = round(as.numeric(dfprint[,5]), 5)
  dfprint[,6] = round(as.numeric(dfprint[,6]), 5)
  dfprint[,4] = round(as.numeric(dfprint[,4]), 5)
  fObj  = format(dfprint)
  strings = apply( dfprint,2,  function(x) unlist(dfprint))[1,]
  widths <- nchar(strings)
  names = c("Group" ,"Parameter", "Estimate", "Std. Error", "T for H0:", "Prob>|T|")
  widths <- pmax(nchar(strings), nchar(names))
  csum <- sum(widths + 1) - 1
  cat(paste("Call TrajeR with"), Obj$groups, "groups and a", paste(Obj$degre, collapse=","),
      "degrees of polynomial shape of trajectory.\n")
  cat("Model : Zero Inflated Poisson\n")
  if (Obj$Method == "L"){
    cat("Method : Likelihood \n \n")
  } else if (Obj$Method == "EM"){
    cat("Method : Expectation-maximization \n \n")
  }else{
    cat("Method : Expectation-maximization with IWRLS\n \n")
  }
  esp = 1
  sep1 <- SepLine1(widths, pad = esp)
  sep2 <- SepLine2(widths, pad = esp)
  namestmp = colnames(dfprint)
  namestmp[5] = "T for H0:"
  writeLines(Row(namestmp, widths, esp))
  writeLines(Row(c("","","","","param.=0",""), widths, esp))
  writeLines(sep1)
  for (i in 1:(ng-1)){
    for (j in 1:(nbeta[i]+nnu[i]+ndelta)){
      writeLines(Row(dfprint[indcum[i]+j , ], widths, esp))
    }
    writeLines(sep2)
  }
  for (j in 1:(nbeta[ng]+nnu[ng]+ndelta)){
    writeLines(Row(dfprint[indcum[ng]+j, ], widths, esp))
  }
  writeLines(sep1)
  if (ntheta == 1){
    for (i in 1:ng){
      if (ndelta !=0){
        writeLines(Row(dfprint[sum(nbeta)+sum(nnu)+ndelta*ng+i, ], widths, esp))
      }else{
        writeLines(Row(dfprint[sum(nbeta)+sum(nnu)+ndelta*ng+i, ], widths, esp))
      }
    }
  }else{
    writeLines(Row(dfprint[indcum[ng+1]+1, ], widths, esp))
    writeLines(sep2)
    if (ng>2){
      for (i in 2:(ng-1)){
        for (j in 1:((ng-2)*ntheta)){
          writeLines(Row(dfprint[indcum[ng+1]+1+j, ], widths, esp))
        }
        writeLines(sep2)
      }
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+(ng-2)*ntheta+1+j, ], widths, esp))
      }
    }else{
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+1+j, ], widths, esp))
      }
    }
  }
  writeLines(sep1)
  cat("Likelihood :", Obj$Likelihood)
}
###################################################################################
# modification of print's method for class trajectory.LOGIT
####################################################################################
#' Print LOGIT
#'
#' Print mehtod for an object of class "\code{Trajectory.LOGIT}".
#'
#' @param Obj an object of class "\code{Trajectory.LOGIT}".
#' @param ...
#'
#' @return  The print of Obj.
#' @export
#'
#' @examples
#' load(file='data/LOGIT_data01')
#' solL
print.Trajectory.LOGIT <- function(Obj, ...){
  # definiton of different sizes
  n= Obj$Size
  ng = Obj$groups
  nbeta = Obj$degre + 1
  ntheta = length(Obj$theta)/ng
  if (any(is.na(Obj$delta))){
    ndelta = 0
  }else{
    ndelta = length(Obj$delta)/ng
  }
  nbetatmp = c(0, nbeta)
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(ndelta, ng)))
  indbeta = 1:sum(nbeta)
  if (any(is.na(Obj$delta))){
    inddelta = 0
  }else{
    inddelta = (sum(nbeta)+1):(sum(nbeta)+ndelta*ng)
  }
  indtheta = (sum(nbeta) + ndelta*ng + 1):(sum(nbeta) + ndelta*ng + ntheta*ng)
  #for the changement of group in the final table
  indcum=c(0, nbetacum[-1]+(1:ng)*ndelta)
  # adding number of groups
  group = rep("", sum(nbeta)+ndelta*ng + ntheta*ng)
  for (i in 1:ng){
    group[nbetacum[i]+1] = as.character(i)
  }
  for (i in 1:(ng)){
    group[indcum[ng+1]+1+(i-1)*ntheta] = as.character(i)
  }
  # create a data frame for print
  dfprint = cbind(group, Parameter = Obj$Names, Obj$tab)
  dfprint[,2] = as.character(dfprint[,2])
  dfprint[,1] = as.character(dfprint[,1])
  if (ntheta == 1){
    if (Obj$Method == "L"){
      dfprint[indtheta, 3] = exp(dfprint[indtheta, 3])/sum(exp(dfprint[indtheta, 3]))
    }
    dfprint[indtheta, 2] = paste0("pi", 1:ng)
  }else{
    # we store the value of the baseline
    dfprint[indtheta,3] = dfprint[indtheta,3] - dfprint[indtheta[1:ntheta],3]
    dfprint = dfprint[-indtheta[2:ntheta],]
    # we change the name of group 1 to baseline
    dfprint[indcum[ng+1]+1,2] = "Baseline"
  }
  # we reorganize the data frame by group ordering
  tmp =c()
  for (i in 1:ng){
    tmp = rbind(tmp, dfprint[(nbetacum[i]+1):(nbetacum[i+1]),])
    if (ndelta !=0){
      tmp = rbind(tmp, dfprint[(sum(nbeta)+ndeltacum[i]+1):(sum(nbeta)+ndeltacum[i+1]),])
    }
  }
  if (ntheta == 1){
    tmp = rbind(tmp, dfprint[indtheta,])
  }else{
    tmp = rbind(tmp, dfprint[indtheta[1:(length(indtheta)-ntheta+1)],])
  }
  dfprint = tmp
  dfprint[,3] = round(as.numeric(dfprint[,3]), 5)
  dfprint[,5] = round(as.numeric(dfprint[,5]), 5)
  dfprint[,6] = round(as.numeric(dfprint[,6]), 5)
  dfprint[,4] = round(as.numeric(dfprint[,4]), 5)
  fObj  = format(dfprint)
  strings = apply( dfprint,2,  function(x) unlist(dfprint))[1,]
  widths <- nchar(strings)
  names = c("Group" ,"Parameter", "Estimate", "Std. Error", "T for H0:", "Prob>|T|")
  widths <- pmax(nchar(strings), nchar(names))
  csum <- sum(widths + 1) - 1
  cat(paste("Call TrajeR with"), Obj$groups, "groups and a", paste(Obj$degre, collapse=","),
      "degrees of polynomial shape of trajectory.\n")
  cat("Model : Logit\n")
  if (Obj$Method == "L"){
    cat("Method : Likelihood \n \n")
  } else if (Obj$Method == "EM"){
    cat("Method : Expectation-maximization \n \n")
  }else{
    cat("Method : Expectation-maximization with IRLS\n \n")
  }
  esp = 1
  sep1 <- SepLine1(widths, pad = esp)
  sep2 <- SepLine2(widths, pad = esp)
  namestmp = colnames(dfprint)
  namestmp[5] = "T for H0:"
  writeLines(Row(namestmp, widths, esp))
  writeLines(Row(c("","","","","param.=0",""), widths, esp))
  writeLines(sep1)
  for (i in 1:(ng-1)){
    for (j in 1:(nbeta[i]+ndelta)){
      writeLines(Row(dfprint[indcum[i]+j , ], widths, esp))
    }
    writeLines(sep2)
  }
  for (j in 1:(nbeta[ng]+ndelta)){
    writeLines(Row(dfprint[indcum[ng]+j, ], widths, esp))
  }
  writeLines(sep1)
  if (ntheta == 1){
    for (i in 1:ng){
      if (ndelta !=0){
        writeLines(Row(dfprint[sum(nbeta)+ndelta*ng+i, ], widths, esp))
      }else{
        writeLines(Row(dfprint[sum(nbeta)+ndelta*ng+i, ], widths, esp))
      }
    }
  }else{
    writeLines(Row(dfprint[indcum[ng+1]+1, ], widths, esp))
    writeLines(sep2)
    if (ng>2){
      for (i in 2:(ng-1)){
        for (j in 1:((ng-2)*ntheta)){
          writeLines(Row(dfprint[indcum[ng+1]+1+j, ], widths, esp))
        }
        writeLines(sep2)
      }
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+(ng-2)*ntheta+1+j, ], widths, esp))
      }
    }else{
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+1+j, ], widths, esp))
      }
    }
  }
  writeLines(sep1)
  cat("Likelihood :", Obj$Likelihood)
}
###################################################################################
# modification of print's method for class trajectory.ZIP
####################################################################################
#' print NL trajectory
#'
#' @export
#'
print.Trajectory.NL <- function(Obj, ...){
  # definiton of different sizes
  n= Obj$Size
  ng = Obj$groups
  nbeta = Obj$degre + 1
  ntheta = length(Obj$theta)/ng
  if (any(is.na(Obj$delta))){
    ndelta = 0
  }else{
    ndelta = length(Obj$delta)/ng
  }
  nbetatmp = c(0, nbeta)
  nbetacum = cumsum(c(0, nbeta))
  ndeltacum = cumsum(c(0, rep(ndelta, ng)))
  indbeta = 1:sum(nbeta)
  indsigma = (sum(nbeta)+1):(sum(nbeta)+ng)
  if (any(is.na(Obj$delta))){
    inddelta = 0
  }else{
    inddelta = (sum(nbeta)+ng+1):(sum(nbeta)+ng+ndelta*ng)
  }
  indtheta = (sum(nbeta) +ng+ ndelta*ng + 1):(sum(nbeta) +ng+ ndelta*ng + ntheta*ng)
  #for the changement of group in the final table
  indcum=c(0, nbetacum[-1]+(1:ng)*ndelta)
  # adding number of groups
  group = rep("", sum(nbeta)+ng+ndelta*ng + ntheta*ng)
  #number for beta
  for (i in 1:ng){
    group[nbetacum[i]+1] = as.character(i)
  }
  # number for sigma
  for (i in 1:ng){
    group[nbetacum[ng+1]+i] = as.character(i)
  }
  # number for theta
  for (i in 1:(ng)){
    group[indcum[ng+1]+ng+1+(i-1)*ntheta] = as.character(i)
  }
  # create a data frame for print
  dfprint = cbind(group, Parameter = Obj$Names, Obj$tab)
  dfprint[,2] = as.character(dfprint[,2])
  dfprint[,1] = as.character(dfprint[,1])
  if (ntheta == 1){
    if (Obj$Method == "L"){
      dfprint[indtheta, 3] = exp(dfprint[indtheta, 3])/sum(exp(dfprint[indtheta, 3]))
    }
    dfprint[indtheta, 2] = paste0("pi", 1:ng)
  }else{
    # we store the value of the baseline
    dfprint[indtheta,3] = dfprint[indtheta,3] - dfprint[indtheta[1:ntheta],3]
    dfprint = dfprint[-indtheta[2:ntheta],]
    # we change the name of group 1 to baseline
    dfprint[indcum[ng+1]+ng+1,2] = "Baseline"
  }
  # we reorganize the data frame by group ordering
  tmp =c()
  for (i in 1:ng){
    tmp = rbind(tmp, dfprint[(nbetacum[i]+1):(nbetacum[i+1]),])
    if (ndelta !=0){
      tmp = rbind(tmp, dfprint[(sum(nbeta)+ng+ndeltacum[i]+1):(sum(nbeta)+ng+ndeltacum[i+1]),])
    }
  }
  for (i in 1:ng){
    tmp = rbind(tmp, dfprint[sum(nbeta) + i,])
  }
  if (ntheta == 1){
    tmp = rbind(tmp, dfprint[indtheta,])
  }else{
    tmp = rbind(tmp, dfprint[indtheta[1:(length(indtheta)-ntheta+1)],])
  }
  dfprint = tmp
  dfprint[,3] = round(as.numeric(dfprint[,3]), 5)
  dfprint[,5] = round(as.numeric(dfprint[,5]), 5)
  dfprint[,6] = round(as.numeric(dfprint[,6]), 5)
  dfprint[,4] = round(as.numeric(dfprint[,4]), 5)
  fObj  = format(dfprint)
  strings = apply( dfprint,2,  function(x) unlist(dfprint))[1,]
  widths <- nchar(strings)
  names = c("Group" ,"Parameter", "Estimate", "Std. Error", "T for H0:", "Prob>|T|")
  widths <- pmax(nchar(strings), nchar(names))
  csum <- sum(widths + 1) - 1
  cat(paste("Call TrajeR with"), Obj$groups, "groups and a", paste(Obj$degre, collapse=","),
      "degrees of polynomial shape of trajectory.\n")
  cat("Model : Non Linear\n")
  if (Obj$Method == "L"){
    cat("Method : Likelihood \n \n")
  } else if (Obj$Method == "EM"){
    cat("Method : Expectation-maximization \n \n")
  }else{
    cat("Method : Expectation-maximization with IWRLS\n \n")
  }
  esp = 1
  sep1 <- SepLine1(widths, pad = esp)
  sep2 <- SepLine2(widths, pad = esp)
  namestmp = colnames(dfprint)
  namestmp[5] = "T for H0:"
  writeLines(Row(namestmp, widths, esp))
  writeLines(Row(c("","","","","param.=0",""), widths, esp))
  writeLines(sep1)
  # write beta and delt
  for (i in 1:(ng-1)){
    for (j in 1:(nbeta[i]+ndelta)){
      writeLines(Row(dfprint[indcum[i]+j , ], widths, esp))
    }
    writeLines(sep2)
  }
  for (j in 1:(nbeta[ng]+ndelta)){
    writeLines(Row(dfprint[indcum[ng]+j, ], widths, esp))
  }
  writeLines(sep1)
  # write sigma
  for (i in 1:(ng-1)){
    writeLines(Row(dfprint[sum(nbeta)+ndelta*ng + i, ], widths, esp))
  }
  writeLines(Row(dfprint[sum(nbeta)+ndelta*ng + ng, ], widths, esp))
  writeLines(sep1)
  # write theta or pi
  if (ntheta == 1){
    for (i in 1:ng){
      if (ndelta !=0){
        writeLines(Row(dfprint[sum(nbeta)+ng+ndelta*ng+i, ], widths, esp))
      }else{
        writeLines(Row(dfprint[sum(nbeta)+ng+ndelta*ng+i, ], widths, esp))
      }
    }
  }else{
    writeLines(Row(dfprint[indcum[ng+1]+ng+1, ], widths, esp))
    writeLines(sep2)
    if (ng>2){
      for (i in 2:(ng-1)){
        for (j in 1:((ng-2)*ntheta)){
          writeLines(Row(dfprint[indcum[ng+1]+ng+1+j, ], widths, esp))
        }
        writeLines(sep2)
      }
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+ng+(ng-2)*ntheta+1+j, ], widths, esp))
      }
    }else{
      for (j in 1:ntheta){
        writeLines(Row(dfprint[indcum[ng+1]+ng+1+j, ], widths, esp))
      }
    }
  }
  writeLines(sep1)
  cat("Likelihood :", Obj$Likelihood)
}
