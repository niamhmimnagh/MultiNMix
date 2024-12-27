#' @title List of Valid Nimble Distributions and Their Arguments
#'
#' @description
#' A named list where each key is a valid Nimble distribution, and the value is a vector of the required arguments.
#'
#' @format
#' A named list of length 19.
#'
#' @details
#' This object is used internally to validate user-specified prior distributions and their arguments in Nimble models.
#'
# List of valid Nimble distributions and their required arguments
distribution_args <- list(
  dnorm = c("mean", "sd"),
  dexp = c("rate"),
  dgamma = c("shape", "rate"),
  dbeta = c("shape1", "shape2"),
  dunif = c("min", "max"),
  dlnorm = c("meanlog", "sdlog"),
  dbern = c("prob"),
  dpois = c("lambda"),
  dbinom = c("size", "prob"),
  dcat = c("prob"),
  dmnorm = c("mean", "cov"),
  dwish = c("scale", "df"),
  dchisq = c("df"),
  dinvgamma = c("shape", "scale"),
  dt = c("mu", "sigma", "df"),
  dweib = c("shape", "scale"),
  ddirch = c("alpha"),
  dmulti = c("size", "prob"),
  dmvt = c("mean", "cov", "df")
)


#' @title Validate User-Specified Prior for a Single Prior
#'
#' @description Validates a single user-defined prior to ensure it specifies a valid distribution supported by Nimble and has the correct parameters for the distribution.
#'
#' @param prior_name A string representing the name of the prior (e.g., `prior_mean`).
#' @param prior_string A string specifying the prior distribution (e.g., `'dnorm(0, 0.001)'`).
#'
#' @return Returns `TRUE` if the prior is valid. Throws an error if the prior is invalid.
#'
#' @examples
#' validate_prior("prior_mean", "dnorm(0, 0.001)")
#'
#' @export
# Function to validate a single prior string
validate_prior <- function(prior_name, prior_string) {
  # Extract the distribution name
  match <- regexpr("^[a-zA-Z]+", prior_string)
  dist_name <- regmatches(prior_string, match)

  # Check if the distribution name is valid
  if (!dist_name %in% names(distribution_args)) {
    stop(sprintf(
      "Unsupported distribution: '%s' in '%s'. For a full list of supported distributions, see the help page.",
      dist_name, prior_name
    ))
  }

  # Extract the arguments inside parentheses
  args_match <- regexpr("\\((.*)\\)", prior_string)
  args_string <- regmatches(prior_string, args_match)
  args_string <- sub("^\\((.*)\\)$", "\\1", args_string)  # Remove outer parentheses

  # Split arguments while respecting brackets (square brackets, parentheses)
  args <- strsplit(args_string, ",(?![^\\[\\(]*[\\]\\)])", perl = TRUE)[[1]]
  args <- trimws(args)  # Remove whitespace

  # Get the expected arguments for the distribution
  expected_args <- distribution_args[[dist_name]]

  # Check if the correct number of arguments are provided
  if (length(args) != length(expected_args)) {
    stop(sprintf(
      "Incorrect number of arguments for distribution '%s' in '%s'.\nExpected: %s\nProvided: %s",
      dist_name, prior_name, paste(expected_args, collapse = ", "), paste(args, collapse = ", ")
    ))
  }

  # Success
  return(TRUE)
}

#' @title Validate User-Specified Priors for Nimble Models
#'
#' @description Validates a set of user-defined priors to ensure that they specify valid distributions supported by Nimble and have the correct parameters for each distribution.
#'
#' @details The function checks the following:
#'   - Whether the specified distribution is supported by Nimble.
#'   - Whether the correct number of parameters is provided for the distribution.
#'
#' If any prior is invalid, the function throws an informative error with details about the issue.
#'
#' @param priors A named list of priors specified as strings, where the name is the prior name (e.g., `prior_mean`) and the value is the prior distribution (e.g., `'dnorm(0, 0.001)'`).
#'
#' @return Returns `TRUE` if all priors are valid. Throws an error if any prior is invalid.
#'
#' @seealso
#' - Nimble documentation for a full list of supported distributions.
#'
#' @references
#' - NIMBLE Development Team (2021). NIMBLE: An R Package for Programming with BUGS Models. <https://r-nimble.org/>
#'
#' @export
#'
#'
#'

# Function to validate all priors
validate_all_priors <- function(priors) {
  lapply(names(priors), function(name) {
    validate_prior(name, priors[[name]])
  })
}
