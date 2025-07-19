\name{SimulatedData}
\alias{SimulatedData}

\title{Simulated data for illustration}
\description{This is a simulated dataset for illustration. It contains a total of \eqn{N=436} observations at irregularly spaced time points for \eqn{n=50} subjects. There are \eqn{p=100} covariates. 
 }

\usage{data(SimulatedData)}

\details{This simulated dataset contains \eqn{N=436} observations for \eqn{n=50} subjects, with \eqn{p=100} covariates. The first column \code{y} gives the response variables, the second column \code{t} gives the observation times, the third column \code{id} gives the unique IDs for each of the 50 subjects, and columns 4-103 (\code{x1}, ..., \code{x100}) give the covariate values.

This synthetic dataset is a slight modification from Experiment 2 in Section 5.1 of Bai et al. (2023). We use \eqn{p=100} for illustration, instead of \eqn{p=500} as in the paper.
}

\references{
Bai, R., Boland, M. R., and Chen, Y. (2023). "Scalable high-dimensional Bayesian varying coefficient models with unknown within-subject covariance." \emph{Journal of Machine Learning Research}, \bold{24}:1-49.
}