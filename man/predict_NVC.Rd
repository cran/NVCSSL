\name{predict_NVC}
\alias{predict_NVC}

\title{Prediction with Nonparametric Varying Coefficient Models}

\description{This function takes a model fit by a varying coefficient model,
\deqn{y_i(t_{ij}) = \sum_{k=1}^{p} x_{ik}(t_{ij}) \beta_k(t_{ij}) + \varepsilon_{ij},}
and predicts the values \eqn{\hat{y}(t_{new})} on a new data set. This function can be used with NVC models fit using the \code{NVC_SSL}, \code{robustified_NVC_SSL}, and \code{NVC_frequentist} functions.
}

\usage{predict_NVC(NVC.mod, t.new, id.new, X.new) }

\arguments{
  \item{NVC.mod}{a model output by \code{NVC_SSL}, \code{robustified_NVC_SSL}, or \code{NVC_frequentist}.}
  \item{t.new}{a vector of new time points.}
  \item{id.new}{a vector of identifiers for the new subjects.}
  \item{X.new}{a design matrix for the new observations.}
}

\value{The function returns a matrix with two columns. The first column \code{id} contains the identifiers for the new subjects. The second column \code{y.hat} contains the predictions \eqn{\hat{y}} for all the new observations.}

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

## Fit NVC model with group lasso penalty
nvc.gLasso.mod = NVC_frequentist(y=y.train, t=t.train, X=X.train, penalty="gLasso")

## Prediction on test set
y.preds = predict_NVC(nvc.gLasso.mod, t.new=t.test, id.new=id.test, X.new=X.test)

## MSPE
nvc.gLasso.mspe = mean((y.test-y.preds[,2])^2)
nvc.gLasso.mspe
}