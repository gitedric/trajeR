% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InternalTrajeR.R
\name{trajeR.CNORM}
\alias{trajeR.CNORM}
\title{Internal function to fit CNORM Model}
\usage{
trajeR.CNORM(
  Y,
  A,
  X,
  TCOV,
  ng,
  nx,
  n,
  nbeta,
  nw,
  ntheta,
  period,
  degre,
  theta,
  beta,
  sigma,
  delta,
  pi,
  Method,
  ssigma,
  ymax,
  ymin,
  hessian,
  itermax,
  paraminit,
  EMIRLS,
  refgr
)
}
\arguments{
\item{Y}{a matrix containing the variables in the model.}

\item{A}{a matrix containing the time variable data.}

\item{TCOV}{an optionnal matrix containing the time covariate that influence the trajectory themselve.
By default its value is NULL.}

\item{ng}{integer. The number of group.}

\item{degre}{vector of integer. The degre of every ploynomial function.}

\item{Method}{determine the method used for find the paramters of the model.
The value are L for the Maximum Likelihiood Estimation, EM for Expectation Maximization method
with quasi newton method inside, EMIWRLS for Expectation Maximization method with Iterative
Weighted Least Square.}

\item{ssigma}{logical. By default its value is FALSE. For the CNORM model,
indicate if we want the same sigma for all nomrmal density function.}

\item{ymax}{real. For the CNORM model, indicate the maximum value of the data. It oncerne only the model
with censored data. By default its value is the maximum value of the data plus 1.}

\item{ymin}{real. For the CNORM model, indicate the minimum value of the data. It oncerne only the model
with censored data. By default its value is the maximum value of the data minus 1.}

\item{hessian}{logical. Indicate if we want calculate the hessian matrix. Default is FALSE.
If the method use is Likelihood, the hessian is calculated by inversing the Information's Fisher Matrix.
To avoid numerically singular matrix we find the pseudo inverse matrix by using the ginv function int he package MASS.
If the method is EM or EMIWRLS, the hessian is calculted by using Louis method.}

\item{itermax}{integer. Indicate the maximal number of iteration for optim function or for the EM algorithm.}

\item{paraminit}{vector. The vector of initial parameters. By default trajeR calculate the initial value
based of the range or the standrad deviation.}

\item{refgr}{integer. The number of reference group. By default is 1.}
}
\value{
return a object of class Trajectory.CNORM
\itemize{
  \item beta -  vector of the parameter beta.
  \item sigma - vector of the parameters sigma.
  \item delta - vector of the parameter delta. Only if we use time covariate.
  \item theta - vector with the parameter theta if there exist a coavriate X that modify
  the probability or the probability of group membership.
  \item sd - vector of the standard deviation of the parameters.
  \item tab - a matrix with all the parameters and standard deviation.
  \item Model - a string with the model used.
  \item groups -  a integer with the number of group.
  \item Names - strings with the name of the parameters.
  \item Method  -  a string with the method used.
  \item Size - a integer with the number of individuals.
  \item Likelihood -  a real with the Likelihood obtained by the parameters.
  \item Time - a vector with the first row of time values.
  \item degre - a vector with the degree of the polynomial shape.
  \item min - a real with the minimum value for censored data.
  \item max - a real with the maximum value for censored data.
}
}
\description{
Internal function to fit CNORM Model
}
