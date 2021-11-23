# -----------------------------------
# compute the BIC function to a model
# -----------------------------------
#' BIC function to an trajectory object
#'
#' Calculate the BIC value to an trajectory object.
#'
#' @param sol Trajectory's object. An object of type trajectory.
#' @return A real. 
#' @export
#' @examples
#' \dontrun{
#' load("data/dataNORM01.RData")
#' solL = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
#'               Model="CNORM", Method = "L", ssigma = FALSE, 
#'               hessian = TRUE)
#' trajeBIC(solL)
#' }
trajeRBIC <- function(sol){
  -2*sol$Likelihood + log(sol$Size)*(nrow(sol$tab) - 1)
}
# -----------------------------------
# compute the AIC function to a model
# -----------------------------------
#' AIC function to an trajectory object
#'
#' Calculate the AIC value to an trajectory object.
#'
#' @param sol Trajectory's object. An object of type trajectory.
#' @return A real. 
#' @export
#' @examples
#' \dontrun{
#' load("data/dataNORM01.RData")
#' solL = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
#'               Model="CNORM", Method = "L", ssigma = FALSE, 
#'               hessian = TRUE)
#' trajeAIC(solL)
#' }
trajeRAIC <- function(sol){
  -2*sol$Likelihood + 2*(nrow(sol$tab) - 1)
}
# -----------------------------------
# compute the SH function to a model
# -----------------------------------
#' SH function to an trajectory object
#'
#' Calculate the Slope Heuristic value to a list of trajectory objects.
#'
#' @param l List. A list of objects of type trajectory.
#' @return A vector of real. 
#' @export
#' @examples
#' \dontrun{
#' load("data/dataNORM01.RData")
#' solL = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
#'               Model="CNORM", Method = "L", ssigma = FALSE, 
#'               hessian = TRUE)
#' trajeSH(solL)
#' }
trajeRSH <- function(l){
  data = data.frame(model = l[[1]]$groups, pen = (nrow(l[[1]]$tab) - 1)/l[[1]]$Size, complexity = nrow(l[[1]]$tab) - 1,
  contrast = -l[[1]]$Likelihood)
  for (i in 2:length(l)){
    data = rbind(data, c(l[[i]]$groups, (nrow(l[[i]]$tab) - 1)/l[[i]]$Size, nrow(l[[i]]$tab) - 1,
                 -l[[i]]$Likelihood))
  }
  return(capushe::DDSE(data))
}
