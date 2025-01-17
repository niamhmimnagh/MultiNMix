#' @title Multi-Species N-Mixture (MNM) Model Class
#'
#' @description
#' The `"MNM"` class represents a multi-species N-mixture model object fitted using Nimble.
#' It contains key model outputs, including parameter estimates, input data, predictions,
#' log-likelihood, convergence diagnostics, and model summaries.
#'
#' @slot summary A `data.frame` summarizing the model output, including posterior means, standard deviations, and convergence diagnostics (e.g., Rhat values).
#' @slot n_parameters A numeric value indicating the number of parameters estimated in the model.
#' @slot estimates Mean estimates for all monitored parameters.
#' @slot fitted_Y An array containing the fitted values (predicted responses) for the model.
#' @slot data An array containing the input data used to fit the model.
#' @slot logLik A numeric value representing the log-likelihood of the model.
#' @slot n_converged A numeric value indicating the number of parameters with Rhat < 1.1, showing convergence.
#' @slot plot A list containing traceplots and density plots for all monitored variables.
#'
#' @importFrom methods new
#' @export
setClass("MNM",
         slots = list(
           summary="data.frame",
           n_parameters = "numeric",
           fitted_Y="array",
           estimates = "list",
           data = "array",
           logLik = "numeric",
           n_converged = "numeric",
           plot="list"
         ))

#' @title Log-Likelihood Method for "MNM" Class
#'
#' @description
#' Retrieves the log-likelihood value from an object of class `"MNM"`.
#'
#' @param object An object of class `"MNM"`.
#' @return A numeric value representing the log-likelihood of the fitted model.
#' @export
setMethod("logLik", "MNM", function(object) {
  if (!is.null(object@logLik)) {
    return(object@logLik)
  } else {
    stop("Log-likelihood not available in Nimble output.")
  }
})

#' @title AIC Method for "MNM" Class
#'
#' @description
#' Computes the Akaike Information Criterion (AIC) for an object of class `"MNM"`.
#' AIC is a metric used for model comparison, balancing goodness of fit and model complexity.
#'
#' The formula for AIC is:
#' \deqn{AIC = -2 \cdot \log L + 2 \cdot k}
#' where:
#' \itemize{
#'   \item \eqn{\log L} is the log-likelihood of the model.
#'   \item \eqn{k} is the number of parameters in the model.
#' }
#'
#' @param object An object of class `"MNM"`.
#' @return A numeric value representing the AIC of the fitted model.
#' @export
setMethod("AIC", "MNM", function(object) {
  logLik_val <- object@logLik
  k <- object@n_parameters
  -2 * logLik_val + 2 * k
})

#' @title BIC Method for "MNM" Class
#'
#' @description
#' Computes the Bayesian Information Criterion (BIC) for an object of class `"MNM"`.
#' BIC is a metric used for model comparison, accounting for goodness of fit and the size of the dataset.
#'
#' The formula for BIC is:
#' \deqn{BIC = -2 \cdot \log L + \log(n) \cdot k}
#' where:
#' \itemize{
#'   \item \eqn{\log L} is the log-likelihood of the model.
#'   \item \eqn{n} is the number of observations in the dataset.
#'   \item \eqn{k} is the number of parameters in the model.
#' }
#'
#' @param object An object of class `"MNM"`.
#' @return A numeric value representing the BIC of the fitted model.
#' @export
setMethod("BIC", "MNM", function(object) {
  logLik_val <- object@logLik
  k <- object@n_parameters
  n <- nrow(object@data)
  -2 * logLik_val + log(n) * k
})

#' Convergence Check Generic Function
#'
#' A generic function for checking the convergence of an `"MNM"` object.
#' @param object An object of class `"MNM"`.
#' @return A numeric value indicating the number of parameters that meet the convergence criterion (Rhat < 1.1).
#' @export
setGeneric("check_convergence", function(object) standardGeneric("check_convergence"))

#' @title Convergence Check for "MNM" Class
#'
#' @description
#' Checks the convergence status of an object of class `"MNM"` based on convergence diagnostics (e.g., Rhat values).
#' Returns the number of parameters that meet the convergence criterion (Rhat < 1.1).
#'
#' @param object An object of class `"MNM"`.
#' @return A numeric value indicating the number of converged parameters.
#' @export
setMethod("check_convergence", "MNM", function(object) {
  object@n_converged
})

#' @title Predict Fitted Values for "MNM" Class
#'
#' @description
#' Returns the predicted (fitted) values, denoted as \eqn{Y}, from an object of class `"MNM"`.
#'
#' @param object An object of class `"MNM"`.
#' @return An array containing the predicted (fitted) values from the model.
#' @export
setGeneric("predictY", function(object) standardGeneric("predictY"))

#' @describeIn predictY Predict fitted values for the "MNM" class.
setMethod("predictY", "MNM", function(object) {
  object@fitted_Y
})

#' @title Density Plot Method for "MNM" Class
#'
#' @description
#' Generates the density plot for a specified parameter stored in the MNM object.
#'
#' @param x An object of class `"MNM"`.
#' @param param The name of the parameter to plot (e.g., `"N[8, 1]"`).
#' @param ... Additional arguments (not used).
#' @return No return value; displays the density plot for the specified parameter.
#' @examples
#' # Calling the density function
#' # density(y, "N[10, 1]")
#'
#' @export
setMethod("density", "MNM", function(x, param, ...) {
  if (!param %in% names(x@plot)) {
    stop(paste("Parameter", param, "not found in the MNM object."))
  }

  if (!"density" %in% names(x@plot[[param]])) {
    stop(paste("Density plot not available for parameter", param))
  }

  # Call the density plot function for the specified parameter
  x@plot[[param]]$density()
})

#' @title Trace Plot Generic and Method for "MNM" Class
#'
#' @description
#' The `tracePlot` function generates trace plots for parameters in objects. This
#' documentation includes the generic and the method for the `"MNM"` class.
#'
#' @param x An object of class `"MNM"`.
#' @param param The name of the parameter to plot (e.g., `"N[8, 1]"`).
#' @param ... Additional arguments (not used for `"MNM"`).
#' @return No return value; displays a trace plot for the specified parameter.
#' @examples
#'
#' # Assuming `y` is an object of class "MNM" with plots stored
#' # tracePlot(y, "N[8, 1]")   Generates a trace plot for parameter N[8, 1]
#'
#' @export
setGeneric("tracePlot", function(x, param, ...) standardGeneric("tracePlot"))

#' @describeIn tracePlot Method for objects of class `"MNM"`.
setMethod("tracePlot", "MNM", function(x, param, ...) {
  if (!param %in% names(x@plot)) {
    stop(paste("Parameter", param, "not found in the MNM object."))
  }

  if (!"trace" %in% names(x@plot[[param]])) {
    stop(paste("Trace plot not available for parameter", param))
  }

  # Call the trace plot function for the specified parameter
  x@plot[[param]]$trace()
})
