% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajeRmodelSelection.R
\name{trajeRBIC}
\alias{trajeRBIC}
\title{BIC function to an trajectory object}
\usage{
trajeRBIC(sol)
}
\arguments{
\item{sol}{An object of type trajectory.}
}
\value{
A real.
}
\description{
Calculate the BIC value to an trajectory object.
}
\examples{
load("data/dataNORM01.RData")
solL = trajeR(data[,1:5], data[,6:10], ng = 3, degre=c(2,2,2), 
              Model="CNORM", Method = "L", ssigma = FALSE, 
              hessian = TRUE)
trajeBIC(solL)
}
