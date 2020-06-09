#' Average Posterior Porbability
#'
#' Calculate the Average Posterior Porbability. Average Posterior Probability
#' (AvePP) the average posterior probability of membership for each group for
#' those individuals that were assigned to.
#'
#' @param sol An object of type trajectory.
#' @param Y A matrix containing the variables in the model.
#' @param A A matrix containing the time variable data.
#' @param X An optionnal matrix that modifie the probabilty of belong to group.
#'   By default its value is a mtrix with one column  with value 1.
#'
#' @return A vector of real. The average posterior probability.
#' @export
#'
#' @examples
#' load("data/dataNORM01.RData")
#' sol = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), Model="CNORM",
#' Method = "L", ssigma = FALSE, hessian = FALSE)
#' AvePP(sol, Y = data[, 1;5], A = data[, 6:10])
AvePP <- function(sol, Y, A, X = NULL){
  Xt = cbind(matrix(rep(1, sol$Size), ncol = 1), X)
  if ((length(sol$theta) == sol$groups) & (sol$Method == "L")){
    sol$theta = exp(sol$theta)/sum(exp(sol$theta))
    sol$Method = "EM"
  }
  prob = GroupProb(sol, Y, A, X = Xt)
  gr = sapply(1:nrow(prob), function(s){which.max(prob[s,])})
  res = c()
  for (i in 1:sol$groups){
    res = c(res, mean(apply(prob[gr == i,], 1, max)))
  }
  return(res)
}

#'  Odds of Correct Classification
#'  
#'  Calculate Odds of Correct Classification. The Odds of Correct Classification for
#'  group k (OCCj) is the ratio between the odds of a correct
#'  classification into group j on the basis of the posterior probability rule and
#'  the odds of correct assignment based on random assignments with the probability
#'  of assignment to group j is done with pik, the probability estimate by
#'  the model.
#'  
#' @inheritParams AvePP
#'
#' @return A vector of real. The Odds of Correct Classification.
#' @export
#'
#' @examples
#' load("data/dataNORM01.RData")
#' sol = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), Model="CNORM",
#' Method = "L", ssigma = FALSE, hessian = FALSE)
#' OCC(sol, Y = data[, 1;5], A = data[, 6:10])
OCC <- function(sol, Y, A, X = NULL){
  tmp = AvePP(sol, Y, A, X)
  return(tmp/(1-tmp)/(sol$theta/(1-sol$theta)))
}

#' Assignment proportion
#'
#' Calculate the proportion of individuals in a given group. That is the ratio of
#' the number of individuals in one group and all the individuals.
#' @inheritParams AvePP
#'
#' @return A vector of real. The  proportion.
#' @export
#'
#' @examples
#' load("data/dataNORM01.RData")
#' sol = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), Model="CNORM",
#' Method = "L", ssigma = FALSE, hessian = FALSE)
#' propAssign(sol, Y = data[, 1;5], A = data[, 6:10])
propAssign <- function(sol, Y, A, X = NULL){
  Xt = cbind(matrix(rep(1, sol$Size), ncol = 1), X)
  prob = GroupProb(sol, Y, A, X = Xt)
  gr = sapply(1:nrow(prob), function(s){which.max(prob[s,])})
  return(table(gr)/sol$Size)
}

#' Confidence interval
#'
#' Calculate the confidence interval of the probabilities with a bootstrap
#' method. We have to specify the number of the repetition of bootstrap and the
#' degree of confidence.
#'
#' @inheritParams AvePP
#' @param nb An integer. The numbers of repetition in the bootstrap method.
#' @param alpha A number. The degree of confidence of the interval.
#' @return A vector of real. The two bounds of the confidence interval given a
#'   degree of confidence.
#' @export
#'
#' @examples
#' load("data/dataNORM01.RData")
#' sol = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), Model="CNORM",
#' Method = "L", ssigma = FALSE, hessian = FALSE)
#' ConfIntT(sol, Y = data[, 1;5], A = data[, 6:10])
ConfIntT <- function(sol, Y, A, X = NULL, nb = 10000, alpha = 0.98){
  Xt = cbind(matrix(rep(1, sol$Size), ncol = 1), X)
  theta = sol$tab[(length(c(sol$beta,sol$delta))+sol$groups):(length(c(sol$beta,sol$delta, sol$theta))+sol$groups-1),1]
  sdthet = sol$tab[(length(c(sol$beta,sol$delta))+sol$groups):(length(c(sol$beta,sol$delta, sol$theta))+sol$groups-1),2]
  boottheta = sapply(1:sol$groups, function(s){rnorm(nb, theta[s], sdthet[s])})
  prob = exp(boottheta)/rowSums(exp(boottheta))
  sapply(1:sol$groups, function(s){quantile(prob[,s], probs = c((1-alpha)/2,1-(1-alpha)/2))})
}

#' Adequacy
#'
#' Calculate the summary of the five method above.
#'
#' @inheritParams AvePP
#' @param nb An integer. The numbers of repetition in the bootstrap method.
#' @param alpha A number. The degree of confidence of the interval.
#' @return A table of real. A table with 5 rows: the estimate probabilities, the
#'  two bounds of the confidence interval, the proportion of assignment, the
#'  Average Posterior Probability and the Odds of Correct Classification.
#' @export
#'
#' @examples
#' load("data/dataNORM01.RData")
#' sol = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), Model="CNORM",
#' Method = "L", ssigma = FALSE, hessian = FALSE)
#' adequacy (sol, Y = data[, 1;5], A = data[, 6:10])
adequacy <- function(sol, Y, A, X = NULL, nb = 10000, alpha = 0.98){
  tab = rbind(exp(sol$theta)/sum(exp(sol$theta)),
              ConfIntT(sol, Y, A, X = NULL, nb = 10000, alpha = 0.98),
              propAssign(sol, Y, A, X = NULL),
              AvePP(solSC,  Y = data[,1:5], A = data[,6:10]),
              OCC(solSC, Y=data[,1:5], A=data[,6:10])
              )
  rownames(tab) = c("Prob. est.", "CI inf.", "CI sup.", "Prop.", "AvePP", "OCC")
  return(tab)
}
