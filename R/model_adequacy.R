#' Average Posterior Probability
#'
#' Calculate the Average Posterior Probability. Average Posterior Probability
#' (AvePP) is the average posterior probability of membership for each group for
#' those individuals that were assigned to.
#'
#' @param sol Trajectory's object. An object of type Trajectory.
#' @param Y Matrix. A matrix containing the variables in the model.
#' @param A Matrix. A matrix containing the time variable data.
#' @param X Matrix. An optional matrix that modifies the probability of belong to group.
#'   By default its value is a one column matrix with value 1.
#'
#' @return A vector of reals. The average posterior probability.
#' @export
#'
#' @examples
#' data <- read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data <- as.matrix(data)
#' sol <- trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2, 2), Model = "CNORM", Method = "EM")
#' AvePP(sol, Y = data[, 2:6], A = data[, 7:11])
AvePP <- function(sol, Y, A, X = NULL) {
  # Xt = cbind(matrix(rep(1, sol$Size), ncol = 1), X)
  Xt <- X
  #
  # A VOIR
  #
  # if ((length(sol$theta) == sol$groups) & (sol$Method == "L")){
  #   sol$theta = exp(sol$theta)/sum(exp(sol$theta))
  #   sol$Method = "EM"
  # }
  prob <- GroupProb(sol, Y, A, X = Xt)
  gr <- as.vector(sapply(1:nrow(prob), function(s) {
    which.max(prob[s, ])
  }))
  res <- c()
  for (i in 1:sol$groups) {
    res <- c(res, mean(apply(prob[gr == i, ], 1, max)))
  }
  return(res)
}

#' Odds of Correct Classification
#'
#' Calculate Odds of Correct Classification. The Odds of Correct Classification for
#' group k (OCCj) is the ratio between the odds of a correct
#' classification into group j on the basis of the posterior probability rule and
#' the odds of correct assignment based on random assignments with the probability
#' of assignment to group j is the probability estimate by
#' the model.
#'
#' @inheritParams AvePP
#'
#' @return A vector of reals. The Odds of Correct Classification.
#' @export
#'
#' @examples
#' data <- read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data <- as.matrix(data)
#' sol <- trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2, 2), Model = "CNORM", Method = "EM")
#' OCC(sol, Y = data[, 2:6], A = data[, 7:11])
OCC <- function(sol, Y, A) {
  tmp <- AvePP(sol, Y, A, X = NULL)
  if (sol$Method == "L") {
    prob <- exp(sol$theta) / sum(exp(sol$theta))
  } else {
    prob <- sol$theta
  }
  return(tmp / (1 - tmp) / (prob / (1 - prob)))
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
#' data <- read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data <- as.matrix(data)
#' sol <- trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2, 2), Model = "CNORM", Method = "EM")
#' propAssign(sol, Y = data[, 2:6], A = data[, 7:11])
propAssign <- function(sol, Y, A) {
  # Xt = matrix(rep(1, sol$Size), ncol = 1)
  # prob = GroupProb(sol, Y, A, X = Xt)
  prob <- GroupProb(sol, Y, A, X = NULL)
  gr <- as.vector(sapply(1:nrow(prob), function(s) {
    which.max(prob[s, ])
  }))
  tab <- matrix(sapply(1:sol$groups, function(s) {
    sum(gr == s)
  }), nrow = 1)
  colnames(tab) <- 1:sol$groups
  return(tab / sol$Size)
}

#' Confidence interval
#'
#' Calculate the confidence interval of the probabilities with bootstrap
#' method. We have to specify the number of the repetitions of bootstrap and the
#' degree of confidence.
#'
#' @inheritParams AvePP
#' @param nb An integer. The number of repetitions in the bootstrap method.
#' @param alpha A number. The degree of confidence of the interval.
#' @return A vector of reals. The two bounds of the confidence interval given a
#'   degree of confidence.
#' @export
#'
#' @examples
#' data <- read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data <- as.matrix(data)
#' sol <- trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2, 2), Model = "CNORM", Method = "EM")
#' ConfIntT(sol, Y = data[, 2:6], A = data[, 7:11])
ConfIntT <- function(sol, Y, A, nb = 10000, alpha = 0.98) {
  Xt <- cbind(matrix(rep(1, sol$Size), ncol = 1))
  vtmp <- c(sol$beta, sol$delta, sol$phi, sol$nu, sol$sigma)
  vtmp <- vtmp[!is.na(vtmp)]
  indmin <- length(vtmp) + 1
  indmax <- indmin + length(c(sol$theta)) - 1
  theta <- sol$tab[indmin : indmax, 1]
  sdthet <- sol$tab[indmin : indmax, 2]
  boottheta <- sapply(1:sol$groups, function(s) {
    stats::rnorm(nb, theta[s], sdthet[s])
  })
  prob <- boottheta
  if (sol$Method == "L") {
    prob <- exp(boottheta) / rowSums(exp(boottheta))
  }
  sapply(1:sol$groups, function(s) {
    stats::quantile(prob[, s], probs = c((1 - alpha) / 2, 1 - (1 - alpha) / 2))
  })
}

#' Adequacy of the model
#'
#' Calculate the summary of the five methods : assignment proportion, average posterior probability, confidence interval, odds of Correct Classification.
#'
#' @inheritParams AvePP
#' @param nb Integer. The numbers of repetitions in the bootstrap method.
#' @param alpha  Real. The degree of confidence of the interval.
#' @return A table of reals. A table with 5 rows: the estimate probabilities, the
#'  two bounds of the confidence interval, the proportion of assignment, the
#'  Average Posterior Probability and the Odds of Correct Classification.
#' @export
#'
#' @examples
#' data <- read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data <- as.matrix(data)
#' sol <- trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2, 2), Model = "CNORM", Method = "EM")
#' adequacy(sol, Y = data[, 2:6], A = data[, 7:11])
adequacy <- function(sol, Y, A, nb = 10000, alpha = 0.98) {
  tab <- rbind(
    exp(sol$theta) / sum(exp(sol$theta)),
    ConfIntT(sol, Y, A, nb = 10000, alpha = 0.98),
    propAssign(sol, Y, A),
    AvePP(sol, Y, A),
    OCC(sol, Y, A)
  )
  rownames(tab) <- c("Prob. est.", "CI inf.", "CI sup.", "Prop.", "AvePP", "OCC")
  return(tab)
}
