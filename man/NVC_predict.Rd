\name{NVC_predict}
\alias{NVC_predict}

\title{Prediction for nonparametric varying coefficient (NVC) models}

\description{This is a function to predict the responses \eqn{y(t_{new})} for new subjects
at new time points \eqn{t_{new}} with new covariates \eqn{X_{new}}. The function accepts an estimated NVC model that was fit using either the \code{NVC_SSL} or \code{NVC_frequentist} functions and returns the predicted \eqn{y(t)}'s. This function can be used for either out-of-sample predictions or for in-sample predictions if the "new" subjects are the same as the ones used to obtain the fitted NVC model.}


\usage{
NVC_predict(NVC_mod, t_new, id_new, X_new) 
}

\arguments{
  \item{NVC_mod}{an object with a fitted NVC model returned by the \code{NVC_SSL} or \code{NVC_frequentist} function}
  \item{t_new}{vector of new observation times}
    \item{id_new}{vector of new labels, where a label corresponds to one of the new subjects}
  \item{X_new}{new design matrix with columns \eqn{[X_1, \ldots, X_p]} where the \eqn{k}th column corresponds to the \eqn{k}th covariate. \code{X_new} must have the \eqn{p} columns, i.e. the same number of varying coefficients estimated by \code{NVC_mod}.}
}

\value{The function returns a list containing the following components:
    \item{id}{vector of each \eqn{i}th subject's label}
  \item{time}{vector of each \eqn{j}th observation time for each \eqn{i}th subject}
  \item{y_pred}{vector of predicted responses corresponding to each \eqn{j}th observation time for each \eqn{i}th subject}
}

\references{
Bai, R., Boland, M. R., and Chen, Y. (2023). "Scalable high-dimensional Bayesian varying coefficient models with unknown within-subject covariance." \emph{Journal of Machine Learning Research}, \bold{24}:1-49.
}

\examples{
## Load simulated data
data(SimulatedData)
attach(SimulatedData)
y = SimulatedData$y
t = SimulatedData$t
id = SimulatedData$id
X = SimulatedData[,4:103]

## Fit frequentist penalized NVC model with the group lasso penalty. 
## No need to specify an 'id' argument when using NVC_frequentist() function.

NVC_gLASSO_mod = NVC_frequentist(y=y, t=t, X=X, penalty="gLASSO")

## Make in-sample predictions. Here, we DO need to specify 'id' argument

NVC_gLASSO_predictions = NVC_predict(NVC_gLASSO_mod, t_new=t, id_new=id, X_new=X)

## Subjects
NVC_gLASSO_predictions$id

## Observation times
NVC_gLASSO_predictions$time

## Predicted responses
NVC_gLASSO_predictions$y_pred


\donttest{
## Fit NVC-SSL model to the data instead. Here, we do need to specify id

NVC_SSL_mod = NVC_SSL(y=y, t=t, id=id, X=X)
NVC_SSL_predictions = NVC_predict(NVC_SSL_mod, t_new = t, id_new=id, X_new=X)

## Subjects
NVC_SSL_predictions$id

## Observation times
NVC_SSL_predictions$time

## Predicted responses
NVC_SSL_predictions$y_pred
}
}