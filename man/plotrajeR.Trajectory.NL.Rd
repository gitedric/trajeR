% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotTrajeR.R
\name{plotrajeR.Trajectory.NL}
\alias{plotrajeR.Trajectory.NL}
\title{plot Non Linear trajectory}
\usage{
\method{plotrajeR}{Trajectory.NL}(
  Obj,
  plotcov = NULL,
  col = "black",
  Y = NULL,
  A = NULL,
  Risk = NULL,
  mean = FALSE,
  alpha = 1,
  TCOV = NULL,
  ...
)
}
\arguments{
\item{Obj}{an object of class "\code{Trajectory.LOGIT}".}

\item{plotcov}{an optional vector or matrix with the same length as the time period. Default value is NULL.}

\item{col}{an optional vector. The vector of colors. It must contain a color for each trajectory and each points of groups.
Its length is the double of the number of group. Default value is a grayscale.}

\item{Y}{Matrix. A matrix containing the variables in the model.}

\item{A}{Matrix. A matrix containing the time variable data.}

\item{Risk}{Matrix. An optional matrix that modify the probability of belong to group. By default its value is a matrix
with one column  with value 1.}

\item{mean}{an optional logical. Indicate if the mean of ech group and time value must be draw.}

\item{alpha}{on optional real. Indicate the alpha channel of the points color.}

\item{TCOV}{Matrix. An optional matrix containing the time covariate that influence the trajectory themselves.
By default its value is NULL.}

\item{...}{optional parameters}
}
\value{
a graphic.
}
\description{
plot Non Linear trajectory
}
