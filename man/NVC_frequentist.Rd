\name{NVC_frequentist}
\alias{NVC_frequentist}

\title{Nonparametric Varying Coefficient Model Fitting with Frequentist Penalties}

\description{This function fits frequentist sparse nonparametric varying coefficients (NVC) models for repeated measures data using basis expansions with regularization on groups of basis coefficients. NVC models aim to estimate the functions of time, \eqn{\beta_1(t), \ldots, \beta_p(t)}, in the model,
\deqn{y_i(t_{ij}) = \sum_{k=1}^{p} x_{ik}(t_{ij}) \beta_k(t_{ij}) + \varepsilon_{ij},}
where the subscript \eqn{ij} refers to the \eqn{j}th observation for the \eqn{i}th subject.

Unlike the Bayesian NVC-SSL model, this function does \emph{not} model the within-subject temporal correlations. The data does \emph{not} need to be regularly spaced. In particular, we have \eqn{n} subjects with \eqn{n_i, i=1, \ldots, n} within-study observations each, for a total of \eqn{N=\sum_{i=1}^{n} n_i} observations. 
}

\usage{NVC_frequentist(y, t, X, df=8, penalty=c("gLasso","gSCAD","gMCP"), lambda = NULL)}

\arguments{
  \item{y}{\eqn{N \times 1} vector of all observed responses, \eqn{(y_{11}, \ldots, y_{1 n_1}, \ldots, y_{n1}, \ldots, y_{n n_n})'}.}
  \item{t}{\eqn{N \times 1} vector of all time points, \eqn{(t_{11}, \ldots, t_{1 n_1}, \ldots, t_{n1}, \ldots, t_{n n_n})'}.}
  \item{X}{\eqn{N \times p} design matrix of possibly time-varying covariates. The \eqn{k}th column of \eqn{X} is \eqn{(x_{1k}(t_{11}), \ldots, x_{1k}(t_{1 n_1}), \ldots, x_{nk}(t_{n1}), \ldots, x_{nk}(t_{n n_n}))'}.}
  \item{df}{number of basis functions to use in each basis expansion. Default is \code{df=8}, but the user may specify degrees of freedom as any integer greater than or equal to 3. }
  \item{penalty}{group regularization method to use on the groups of basis coefficients. \code{"gLasso"} is the group lasso penalty, \code{"gSCAD"} is the group SCAD penalty, and \code{"gMCP"} is the group minimax concave penalty (MCP).}
  \item{lambda}{Optional grid for tuning the regularization parameter \eqn{\lambda} from cross-validation. If the user does not provide this, the program chooses the grid automatically. }
}

\value{The function returns a list containing the following components:
\item{t.ordered}{all \eqn{N} time points in order from smallest to largest. Needed for plotting.}
\item{beta.hat}{\eqn{N \times p} matrix of the function estimates. The \eqn{k}th column is the function estimate \eqn{\beta_k(t_{ij})} evaluated at the \eqn{N} time observations \eqn{t_{ij}, i=1,...,n, j=1,...,n_i}.}
\item{gamma.hat}{\eqn{df \times p} estimated basis coefficients \eqn{\gamma}. Needed for prediction.}
\item{intercept}{Estimate of an intercept (or grand mean). Needed for prediction.}
\item{classifications}{\eqn{p \times 1} vector of binary variables for the \eqn{p} covariates. "1" indicates that the covariate was selected, and "0" indicates that it was not selected.}
\item{AICc}{AIC with correction for small sample size. This value can be used for tuning the degrees of freedom \code{df}. The user should pick the model which obtains the lowest \eqn{AIC_c}.}
}

\references{ 
Bai, R., Boland, M. R., and Chen, Y. (2020). "Fast algorithms and theory for high-dimensional Bayesian varying coefficient models." \emph{arXiv preprint arXiv:1907.06477}.

Huang, J., Breheny, P., and Ma, S. (2012). "A selective review of group selection in high-dimensional models." \emph{Statistical Science}, \bold{27}: 481-499.
}

\examples{

#####################################
# Example on the synthetic data set #
#####################################
data(SimulatedData)
attach(SimulatedData)
y = SimulatedData$y
t = SimulatedData$t
X = SimulatedData[,4:103]

## Set seed
set.seed(123)

## Fit frequentist NVC model with group MCP penalty

nvc.gMCP.mod = NVC_frequentist(y=y, t=t, X=X, penalty="gMCP")

## AIC with correction for sample size. Can use this to tune df
nvc.gMCP.mod$AICc

## The first 6 functionals were classified as nonzero
nvc.gMCP.mod$classifications

## Make plot of the 6 active functionals

oldpar <- par(no.readonly = TRUE)

t.ordered = nvc.gMCP.mod$t.ordered
beta.hat = nvc.gMCP.mod$beta.hat

par(mfrow=c(3,2), mar=c(4,4,4,4))
plot(t.ordered, beta.hat[,1], lwd=3, type='l', col='purple', 
     xlab="T", ylim=c(-10,10), ylab=expression(beta[1]))
plot(t.ordered, beta.hat[,2], lwd=3, type='l', col='purple', 
     xlab="T", ylim=c(-8,8), ylab=expression(beta[2]))
plot(t.ordered, beta.hat[,3], lwd=3, type='l', col='purple', 
     xlab="T", ylim=c(-7,4), ylab=expression(beta[3]))
plot(t.ordered, beta.hat[,4], lwd=3, type='l', col='purple', 
     xlab="T", ylim=c(-6,4), ylab=expression(beta[4]))
plot(t.ordered, beta.hat[,5], lwd=3, type='l', col='purple', 
     xlab="T", ylim=c(0,15), ylab=expression(beta[5]))
plot(t.ordered, beta.hat[,6], lwd=3, type='l', col='purple', 
     xlab="T", ylim=c(-5,0), ylab=expression(beta[6]))

par(oldpar)

\donttest{
#############################
# Example on the yeast cell #
# cycle data set            #
#############################

data(yeast)
attach(yeast)
y = SimulatedData$y
t = SimulatedData$t
id = SimulatedData$id
X = SimulatedData[,4:103]

y = yeast$mRNA
t = yeast$Time
id = yeast$Gene
X = yeast[,c(4:99)]

## Set seed
set.seed(12345)

## Fit frequentist NVC model with group SCAD penalty
nvc.gSCAD.mod = NVC_frequentist(y, t, X, penalty="gSCAD")

## TF's classified as significant
nvc.gSCAD.mod$classifications

## Plot the first TF's functional

oldpar <- par(no.readonly = TRUE)

t.ordered = nvc.gSCAD.mod$t.ordered
beta1.hat = nvc.gSCAD.mod$beta.hat[,1]

plot(t.ordered, beta1.hat, type='l', lwd=5, xlab="Time", 
     ylab=expression(beta[1]), main=colnames(X)[1])
     
par(oldpar)
}
}