\name{robustified_NVC_SSL}
\alias{robustified_NVC_SSL}

\title{Robustified Nonparametric Varying Coefficient Spike-and-Slab Lasso (robustified NVC-SSL)}

\description{This function implements the robustified NVC-SSL model of Bai et al. (2020) for repeated measures data. The robustified NVC-SSL model aims to estimate the functions of time, \eqn{\beta_1(t), \ldots, \beta_p(t)} in the model,
\deqn{y_i(t_{ij}) = \sum_{k=1}^{p} x_{ik}(t_{ij}) \beta_k(t_{ij}) + \varepsilon_{ij},}
where the subscript \eqn{ij} refers to the \eqn{j}th observation for the \eqn{i}th subject.

The robustified NVC-SSL model can be used to consistently estimate the functions \eqn{\beta_k(t), k=1, \ldots, p}, when the within-subject covariance matrix is possibly misspecified. The data does \emph{not} need to be regularly spaced. In particular, we have \eqn{n} subjects with \eqn{n_i, i=1, \ldots, n} within-study observations each, for a total of \eqn{N=\sum_{i=1}^{n} n_i} observations. 
}

\usage{robustified_NVC_SSL(y, t, id, X, df=8, 
                           working.cov=c("AR1","CS","ind"), frac.power=0.5, 
                           lambda0=c(seq(from=5,to=30,by=5),seq(from=40,to=100,by=10)), 
                           lambda1=1, a=1, b=ncol(X), tol = 1e-6, print.iteration = TRUE)}
\arguments{
  \item{y}{\eqn{N \times 1} vector of all observed responses, \eqn{(y_{11}, \ldots, y_{1 n_1}, \ldots, y_{n1}, \ldots, y_{n n_n})'}.}
  \item{t}{\eqn{N \times 1} vector of all time points, \eqn{(t_{11}, \ldots, t_{1 n_1}, \ldots, t_{n1}, \ldots, t_{n n_n})'}.}
  \item{id}{\eqn{N \times 1} vector of subject identifiers (\eqn{1} through \eqn{n}).}
  \item{X}{\eqn{N \times p} design matrix of possibly time-varying covariates. The \eqn{k}th column of \eqn{X} is \eqn{(x_{1k}(t_{11}), \ldots, x_{1k}(t_{1 n_1}), \ldots, x_{nk}(t_{n1}), \ldots, x_{nk}(t_{n n_n}))'}.}
  \item{df}{number of basis functions to use in each basis expansion. Default is \code{df=8}, but the user may specify degrees of freedom as any integer greater than or equal to 3. }
  \item{working.cov}{working within-subject covariance structure with which to fit the robustified NVC-SSL model. The options are: \code{"AR1"} for first-order autoregressive, \code{"CS"} for compound symmetry, and \code{"ind"} for independent errors (i.e. no within-subject temporal correlations). The default is \code{"AR1"}.}
  \item{lambda0}{ladder of spike hyperparameters in the spike-and-slab group lasso (SSGL) prior for dynamic posterior exploration. The user may specify either a scalar or a vector. The default  is \code{lambda0=c(5,10,15,20,25,30,40,50,60,70,80,90,100)}. }
  \item{lambda1}{slab hyperparameter in the SSGL prior, which is held fixed. The default is \code{lambda1=1}.}
  \item{a}{shape hyperparameter for the \eqn{B(a,b)} prior on the mixing proportion \eqn{\theta}. Default is \code{a=1}.}
  \item{b}{shape hyperparameter for the \eqn{B(a,b)} prior on the mixing proportion \eqn{\theta}. Default is \code{b=ncol(X)}.}
  \item{frac.power}{The fractional power \eqn{\xi \in (0,1)} to use in order to fit the robustified NVC-SSL model. \code{frac.power} must be strictly between 0 and 1.}
\item{tol}{convergence threshold for the EM algorithm. Default is \code{tol=1e-6}.}
\item{print.iteration}{flag for printing the current value of \code{lambda0} in the ladder of spike hyperparameters. Default is \code{TRUE}.}
}

\value{The function returns a list containing the following components:
\item{t.ordered}{all \eqn{N} time points in order from smallest to largest. Needed for plotting.}
\item{beta.hat}{\eqn{N \times p} matrix of the function estimates. The \eqn{k}th column is the function estimate \eqn{\beta_k(t_{ij})} evaluated at the \eqn{N} time observations \eqn{t_{ij}, i=1,...,n, j=1,...,n_i}.}
\item{gamma.hat}{\eqn{df \times p} estimated basis coefficients \eqn{\gamma}. Needed for prediction.}
\item{intercept}{Estimate of an intercept (or grand mean). Needed for prediction.}
\item{classifications}{\eqn{p \times 1} vector of binary variables for the \eqn{p} covariates. "1" indicates that the covariate was selected, and "0" indicates that it was not selected.}
\item{AICc}{AIC with correction for small sample size. This value can be used for tuning either the fractional power \code{frac.power} or the degrees of freedom \code{df} or both. The user should pick the model which obtains the lowest \eqn{AIC_c}.}
}

\references{ 
Bai, R., Boland, M. R., and Chen, Y. (2020). "Fast algorithms and theory for high-dimensional Bayesian varying coefficient models." \emph{arXiv preprint arXiv:1907.06477}.
}

\examples{

#####################################
# Example on the synthetic data set #
#####################################
data(SimulatedData)
attach(SimulatedData)
y = SimulatedData$y
t = SimulatedData$t
id = SimulatedData$id
X = SimulatedData[,4:103]

## Training set. First 20 subjects
training.indices = which(id <= 20)
y.train = y[training.indices]
t.train = t[training.indices]
id.train = id[training.indices]
X.train = X[training.indices,]

## Test set. Last 10 subjects
test.indices = which(id > 20)
y.test = y[test.indices]
t.test = t[test.indices]
id.test = id[test.indices]
X.test = X[test.indices, ]

## Set seed
set.seed(123)

## Fit robustified NVC-SSL model with 8 degrees of freedom
## (so we are estimating 800 unknown basis coefficients),
## working AR(1) covariance structure, and fractional power 0.6

robust.nvcssl.mod = robustified_NVC_SSL(y=y.train, t=t.train, 
                                       id=id.train, X=X.train, df=8,
                                       working.cov="AR1",
                                       frac.power=0.6)

## AIC with correction for sample size. Can use to tune df or frac.power
robust.nvcssl.mod$AICc

## The first 6 functionals were classified as nonzero
robust.nvcssl.mod$classifications 

## Make plot of the 6 active functionals

oldpar <- par(no.readonly = TRUE)

t.ordered = robust.nvcssl.mod$t.ordered
beta.hat = robust.nvcssl.mod$beta.hat

par(mfrow=c(3,2), mar=c(4,4,4,4))
plot(t.ordered, beta.hat[,1], lwd=3, type='l', col='red', 
     xlab="T", ylim=c(-10,10), ylab=expression(beta[1]))
plot(t.ordered, beta.hat[,2], lwd=3, type='l', col='red', 
     xlab="T", ylim=c(-8,8), ylab=expression(beta[2]))
plot(t.ordered, beta.hat[,3], lwd=3, type='l', col='red', 
     xlab="T", ylim=c(-7,4), ylab=expression(beta[3]))
plot(t.ordered, beta.hat[,4], lwd=3, type='l', col='red', 
     xlab="T", ylim=c(-6,4), ylab=expression(beta[4]))
plot(t.ordered, beta.hat[,5], lwd=3, type='l', col='red', 
     xlab="T", ylim=c(0,15), ylab=expression(beta[5]))
plot(t.ordered, beta.hat[,6], lwd=3, type='l', col='red', 
     xlab="T", ylim=c(-5,0), ylab=expression(beta[6]))

par(oldpar)

## Prediction on the test set
y.preds = predict_NVC(robust.nvcssl.mod, t.new=t.test, id.new=id.test, X.new=X.test)

## MSPE for test set
robust.nvcssl.mod.mspe = mean((y.test-y.preds[,2])^2)
robust.nvcssl.mod.mspe
}