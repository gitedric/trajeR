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
#' data = read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data = as.matrix(data)
#' sol = trajeR(Y = data[, 2:6], A = data[, 7:11], Risk = data[,12], 
#'              degre = c(2,2), Model = "CNORM", Method = "L")
#' GroupProfiles(sol, Y = data[, 2:6], A = data[, 7:11], X = data[,12])
GroupProfiles <- function(sol, Y, A, X){
  prob = GroupProb(sol, Y = Y, A = A, X = X)
  gr = sapply(1:nrow(prob), function(s){which.max(prob[s,])})
  if (!is.matrix((X))){
    X = matrix(X)
  }
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
