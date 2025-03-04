% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_envir.R
\name{clean_envir}
\alias{clean_envir}
\title{Clean Environment of Temporary Objects}
\usage{
clean_envir(env_rm = .GlobalEnv)
}
\arguments{
\item{env_rm}{An environment from which to remove objects. Defaults to the global environment (\code{.GlobalEnv}). This may be set to a different environment for more controlled cleanup.}
}
\value{
The function does not return any values.
}
\description{
This function removes temporary objects from a specified environment. It is designed for use in contexts where intermediate objects are created and should be cleaned up after execution to prevent clutter or unintended behavior.
}
\details{
This is an auxiliary function that removes temporary objects created by Nimble from the environment.
}
