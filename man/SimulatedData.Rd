\name{SimulatedData}
\alias{SimulatedData}

\title{Simulated data with irregularly spaced time points and different number of within-subject observations.}
\description{This data set contains simulated repeated measurements data on \eqn{n=30} subjects, each with a different number of time points and irregularly spaced time points. There are a total of \eqn{N=238} observations and \eqn{p=100} time-varying covariates.  
 }
\usage{data(SimulatedData)}
\format{
A dataframe with the following columns. 
\describe{
    \item{\code{y}:}{\eqn{N=238} responses for all \eqn{n=30} subjects.}
    \item{\code{t}:}{\eqn{N=238} time points for all \eqn{n=30} subjects.}
    \item{\code{id}:}{identifier for subject (1-30).}
    \item{\code{X1-X100}:}{columns 4-103 are the 100 covariates for all \eqn{N} observations.}
}
}

\references{ 
Bai, R., Boland, M. R., and Chen, Y. (2020). "Fast algorithms and theory for high-dimensional Bayesian varying coefficient models." \emph{arXiv preprint arXiv:1907.06477}.
}