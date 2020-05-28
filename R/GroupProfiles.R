#' GroupProfiles
#'
#' Calculate the group profil regarding covariate. It is a cross tabulation
#'  of individual level trajectory group assignments with individual
#'  level characteristic that might be associated with trajectory group membership.
#'
#' @param sol a object of type trajectory.
#' @param Y a matrix containing the variables in the model.
#' @param A a matrix containing the time variable data.
#' @param X an optionnal matrix that modifie the probabilty of belong to group. By default its value is a mtrix
#' with one column  with value 1.
#'
#' @return a table of real.
#' @export
#'
#' @examples
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
