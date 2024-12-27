#' @title Fit a Multi-Species N-Mixture (MNM) Model in Nimble
#'
#' @description Fits a Multi-Species N-Mixture (MNM) model to observed count data using Nimble, with options to include autocorrelation (AR) and/or hurdle components for advanced ecological modeling.
#' @details This function implements the Bayesian MNM model to estimate latent species abundances and inter-species correlations based on observed count data. Extensions include options for incorporating autocorrelation (AR) to account for temporal dependencies and a hurdle model to handle zero inflation in the data. The input data and covariates must conform to specific dimensional requirements as described below.
#'
#' The MNM model (Mimnagh et al., 2022) builds upon Royle's (2004) N-mixture model by allowing simultaneous modeling of multiple species, enabling inference about inter-species relationships and correlations.
#'
#' @param Y Array of observed count data:
#'   - Dimensions for standard MNM and hurdle models: \eqn{R \times T \times S}
#'   - Dimensions for MNM with AR or both AR and hurdle components: \eqn{R \times T \times S \times K}
#'   - `R`: Number of sites
#'   - `T`: Number of replicates
#'   - `S`: Number of species
#'   - `K`: Number of time periods (required if `AR = TRUE`)
#' @param AR Logical. Indicates whether to include an autocorrelation component in the model. Defaults to `FALSE`.
#' @param Hurdle Logical. Indicates whether to include a hurdle component in the model. Defaults to `FALSE`.
#' @param ... Additional arguments for prior distributions. Supported priors include:
#'   - `prior_detection_probability`, `prior_precision`, `prior_mean`, `prior_mean_AR`, `prior_sd_AR`, `prior_hurdle`.
#'   - Supported distributions include: `dnorm`, `dexp`, `dgamma`, `dbeta`, `dunif`, `dlnorm`, `dbern`, `dpois`, `dbinom`, `dcat`, `dmnorm`, `dwish`, `dchisq`, `dinvgamma`, `dt`, `dweib`, `ddirch`, `dmulti`, `dmvt`. Refer to the [Nimble documentation](https://r-nimble.org) for details.
#' @param Xp An array of covariates affecting detection probability, with dimensions (R, S, P1), where:
#'   - `R`: Number of sites
#'   - `S`: Number of species
#'   - `P1`: Number of detection-related covariates
#'   See examples for implementation details.
#' @param Xn An array of covariates affecting abundance, with dimensions (R, S, P2), where:
#'   - `R`: Number of sites
#'   - `S`: Number of species
#'   - `P2`: Number of abundance-related covariates
#'
#' @returns An MNM object that contains the following components:
#'   - summary: Nimble model summary (mean, standard deviation, standard error, quantiles, effective sample size and Rhat value for all monitored values)
#'   - n_parameters: Number of parameters in the model (for use in calculating information criteria)
#'   - data: Observed abundances
#'   - fitted_Y: Predicted values for Y. Posterior predictive checks can be performed by comparing fitted_Y with the observed data.
#'   - logLik: Log-likelihood of the observed data (Y) given the model parameters.
#'   - n_converged: Number of parameters with successful convergence (Rhat < 1.1).
#'
#' @references
#' - Royle, J. A. (2004). N-mixture models for estimating population size from spatially replicated counts. Biometrics, 60(1), 108-115.
#' - Mimnagh, N., Parnell, A., Prado, E., & Moral, R. D. A. (2022). Bayesian multi-species N-mixture models for unmarked animal communities. Environmental and Ecological Statistics, 29(4), 755-778.
#' @seealso
#' - `simulateData`: For generating example datasets compatible with this function.
#' - `MNM`: For details on creation of covariate arrays Xp and Xn.

#' @note
#' Ensure that the dimensions of `Y`, `Xp`, and `Xn` match the requirements specified above. Mismatched dimensions will result in errors during model fitting.#'

#' @examples
#' # Example 1: Fit a standard MNM model
#' Y <- array(data = rpois(60, lambda = 5), dim = c(3, 5, 4))  # Simulated counts
#' Xp <- array(data = rnorm(60), dim = c(3, 4, 2))  # Detection covariates
#' Xn <- array(data = rnorm(60), dim = c(3, 4, 2))  # Abundance covariates
#'
#' \dontrun{model <- MNM_fit(Y = Y, AR = FALSE, Hurdle = FALSE, Xp = Xp, Xn = Xn)}
#'
#' #' # Example 2: Fit an MNM model with AR-1 component
#' Y <- array(data = rpois(180, lambda = 5), dim = c(3, 5, 4, 3))  # Simulated counts
#' Xp <- array(data = rnorm(180), dim = c(3, 4, 3, 2))  # Detection covariates
#' Xn <- array(data = rnorm(180), dim = c(3, 4, 3, 2))  # Abundance covariates
#'
#' \dontrun{model <- MNM_fit(Y = Y, AR = TRUE, Hurdle = FALSE, Xp = Xp, Xn = Xn)}
#'
#' # Example 3: Fit an MNM model with user-specified prior distributions
#' Y <- array(data = rpois(60, lambda = 5), dim = c(3, 5, 4))  # Simulated counts
#' Xp <- array(data = rnorm(60), dim = c(3, 4, 2))  # Detection covariates
#' Xn <- array(data = rnorm(60), dim = c(3, 4, 2))  # Abundance covariates
#'
#' \dontrun{model <- MNM_fit(Y = Y, AR = FALSE, Hurdle = TRUE, Xp = Xp, Xn = Xn,
#'                           prior_detection_probability="dnorm(0.01,0.01)")}
#' # Access traceplots and density plots
#' \dontrun{tracePlot(y, "N[10, 1]")}
#' \dontrun{density(y, "N[10, 1]")}
#'
#' @export
MNM_fit<-function(Y=NULL, AR=FALSE, Hurdle=FALSE, Xp=NULL, Xn=NULL, ...){

  if(AR==TRUE & Hurdle==TRUE){
    if(!is.array(Y)|!length(dim(Y))==4){
      stop("Observations Y must be an array with dimensions R,T,S,K")
    }
    mnm_object<-MNM_Hurdle_AR(Y, Xp=Xp, Xn=Xn,...)

  }
  else if(AR==TRUE & Hurdle==FALSE){
    if(!is.array(Y)|!length(dim(Y))==4){
      stop("Observations Y must be an array with dimensions R,T,S,K")
    }
    mnm_object<-MNM_AR(Y, Xp=Xp, Xn=Xn, ...)

  }
  else if(AR==FALSE & Hurdle==TRUE){
    if(!is.array(Y)|!length(dim(Y))==3){
      stop("Observations Y must be an array with dimensions R,T,S")
    }
    mnm_object<-MNM_Hurdle(Y, Xp=Xp, Xn=Xn, ...)

  }
  else{
    if(!is.array(Y)|!length(dim(Y))==3){
      stop("Observations Y must be an array with dimensions R,T,S")
    }

    mnm_object <- MNM(Y, Xp=Xp, Xn=Xn, ...)  # Call MNM to get the result
    }

  # Run cleanup after MNM completes
  clean_envir(env_rm = .GlobalEnv)
  # Return the MNM object
  return(mnm_object)
}



