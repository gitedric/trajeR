#' Profiles of each group
#'
#' \code{GroupProfiles} calculate the profile of a group regarding covariate. It is a
#' cross tabulation of individual level trajectory group assignments with
#' individual level characteristic that might be associated with trajectory
#' group membership.
#'
#' @param sol Trajectory's object. A object of type trajectory.
#' @param Y Matirx. A matrix containing the variables in the model.
#' @param A Matrix. A matrix containing the time variable data.
#' @param X Matrix. An optional matrix that modify the probability of belong to group.
#'   By default its value is a matrix with one column  with value 1.
#'
#' @return A table of real.
#' @export
#'
#' @examples
#' \dontrun{
#' load("data/dataNORM01.RData")
#' solL = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
#'               Model="CNORM", Method = "L", ssigma = FALSE, 
#'               hessian = TRUE)
#' GroupProfles(sol, Y = data[,1:5], A = data[, 6:10])
#' }
GroupProfiles <- function(sol, Y, A, X){
 # Xt = cbind(matrix(rep(1, sol$Size), ncol = 1), X)
  Xt = X
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
