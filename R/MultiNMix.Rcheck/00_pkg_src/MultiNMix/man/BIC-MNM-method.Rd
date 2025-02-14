% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MNM_class_methods.R
\name{BIC,MNM-method}
\alias{BIC,MNM-method}
\title{BIC Method for "MNM" Class}
\usage{
\S4method{BIC}{MNM}(object)
}
\arguments{
\item{object}{An object of class \code{"MNM"}.}
}
\value{
A numeric value representing the BIC of the fitted model.
}
\description{
Computes the Bayesian Information Criterion (BIC) for an object of class \code{"MNM"}.
BIC is a metric used for model comparison, accounting for goodness of fit and the size of the dataset.

The formula for BIC is:
\deqn{BIC = -2 \cdot \log L + \log(n) \cdot k}
where:
\itemize{
\item \eqn{\log L} is the log-likelihood of the model.
\item \eqn{n} is the number of observations in the dataset.
\item \eqn{k} is the number of parameters in the model.
}
}
