#' Profiles of each group
#'
#' \code{GroupProfiles} calculate the profil of a group regarding covariate. It is a
#' cross tabulation of individual level trajectory group assignments with
#' individual level characteristic that might be associated with trajectory
#' group membership.
#'
#' @param sol A object of type trajectory.
#' @param Y A matrix containing the variables in the model.
#' @param A A matrix containing the time variable data.
#' @param X An optionnal matrix that modifie the probabilty of belong to group.
#'   By default its value is a mtrix with one column  with value 1.
#'
#' @return A table of real.
#' @export
#'
#' @examples
#' load("data/dataNORM01.RData")
#' sol = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), Model="CNORM",
#' Method = "L", ssigma = FALSE, hessian = FALSE)
#' GroupProfles(sol, Y = data[,1:5], A = data[, 6:10])
GroupProfiles <- function(sol, Y, A, X){
  Xt = cbind(matrix(rep(1, sol$Size), ncol = 1), X)
  prob = GroupProb(sol, Y = Y, A = A, X = Xt)
  gr = sapply(1:nrow(prob), function(s){which.max(prob[s,])})
  tab = c()
  for (j in 1:ncol(X)){
    tabl = c()
    for (i in 1:sol$groups){
      tabl = c(tabl, mean(X[gr == i, j]))
    }
    tab = rbind(tab, tabl)
  }
  colnames(tab) = paste0("Gr ", 1:sol$groups)
  rownames(tab) = colnames(X)
  return(tab)
}
