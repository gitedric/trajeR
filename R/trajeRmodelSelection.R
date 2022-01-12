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
#' data = read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data = as.matrix(data)
#' sol = trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2,2), Model = "CNORM", Method = "EM")
#' trajeRBIC(sol)
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
#' data = read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data = as.matrix(data)
#' sol = trajeR(Y = data[, 2:6], A = data[, 7:11], degre = c(2,2), Model = "CNORM", Method = "EM")
#' trajeRAIC(sol)
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
#' data = read.csv(system.file("extdata", "CNORM2gr.csv", package = "trajeR"))
#' data = as.matrix(data)
#' degre = list(c(2,2), c(1,1), c(1,2), c(2,1), c(0,0), c(0,1), c(1,0), c(0,0), c(0,2), c(2,0))
#' sol = list()
#' for (i in 1:10){
#'   sol[[i]] = trajeR(Y = data[, 2:6], A = data[, 7:11], 
#'                     degre = degre[[i]], Model = "CNORM", Method = "EM")
#'   }
#' trajeRSH(sol)
trajeRSH <- function(l){
  data = data.frame(model = l[[1]]$groups, pen = (nrow(l[[1]]$tab) - 1)/l[[1]]$Size, complexity = nrow(l[[1]]$tab) - 1,
  contrast = -l[[1]]$Likelihood)
  for (i in 2:length(l)){
    data = rbind(data, c(l[[i]]$groups, (nrow(l[[i]]$tab) - 1)/l[[i]]$Size, nrow(l[[i]]$tab) - 1,
                 -l[[i]]$Likelihood))
  }
  return(capushe::DDSE(data))
}
