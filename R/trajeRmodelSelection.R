# -----------------------------------
# compute the BIC function to a model
# -----------------------------------
trajeRBIC <- function(sol){
  -2*sol$Likelihood + log(sol$Size)*(nrow(sol$tab) - 1)
}
# -----------------------------------
# compute the AIC function to a model
# -----------------------------------
trajeRAIC <- function(sol){
  -2*sol$Likelihood + 2*(nrow(sol$tab) - 1)
}
# -----------------------------------
# compute the SH function to a model
# -----------------------------------
trajeRSH <- function(l){
  data = data.frame(model = l[[1]]$groups, pen = (nrow(l[[1]]$tab) - 1)/l[[1]]$Size, complexity = nrow(l[[1]]$tab) - 1,
  contrast = -l[[1]]$Likelihood)
  for (i in 2:length(l)){
    data = rbind(data, c(l[[i]]$groups, (nrow(l[[i]]$tab) - 1)/l[[i]]$Size, nrow(l[[i]]$tab) - 1,
                 -l[[i]]$Likelihood))
  }
  return(DDSE(data))
}
