#' @title Fit a Multi-Species N-Mixture Model with Hurdle Component using Nimble
#' @description This function fits a multi-species N-mixture (MNM) model incorporating a Hurdle component to handle zero-inflated data, allowing for robust estimation of abundance and detection probabilities across multiple species and sites.
#'
#' @details
#' This function uses the Nimble framework to fit a Hurdle model, which combines a truncated Poisson distribution for non-zero counts with a separate process for modeling zero counts. The model is particularly suitable for ecological data with excess zeros, such as species occurrence data.
#'
#' The model supports covariates influencing both abundance and detection probabilities, and outputs posterior distributions for model parameters, derived quantities, and predicted values. Convergence diagnostics and posterior predictive checks can also be performed using the returned results.
#'
#' @param Y Array of observed counts, with dimensions (R, T, S, K), where:
#'   - R: Number of sites.
#'   - T: Number of repeated counts (replicates).
#'   - S: Number of species.
#' @param iterations Integer. Number of iterations to be used in the JAGS model. Defaults to 60,000.
#' @param burnin Integer. Number of iterations to be discarded as burn-in. Defaults to 20,000.
#' @param thin Integer. Thinning interval for the MCMC chains. Defaults to 10.
#' @param ... Additional arguments passed for prior distribution specification. Supported distributions include dnorm, dexp, dgamma, dbeta, dunif, dlnorm, dbern, dpois, dbinom, dcat, dmnorm, dwish, dchisq, dinvgamma, dt, dweib, ddirch, dmulti, dmvt. Default prior distributions are:
#' \itemize{
#'  \item prior_detection_probability: prior distribution for the detection probability intercept (`gamma`). Default is `'dnorm(0, 0.001)'`.
#'  \item prior_precision: prior distribution for the precision matrix for the species-level random effect. Default is `'dwish(Omega[1:S,1:S], df)'`.
#'  \item prior_mean: prior distribution for the mean of the species-level random effect (`mu`). Default is `'dnorm(0,0.001)'`.
#' \item prior_hurdle: prior distribution for `theta`, the probability of structural zero in hurdle models. Default is `'dbeta(1,1)'`.
#' \item prior_mean_AR: prior distribution for the mean of the autoregressive random effect (`phi`). Default is `'dnorm(0,0.001)'`.
#' \item prior_sd_AR: prior distribution for the standard deviation of the autoregressive random effect (`phi`). Default is `'dexp(1)'`.
#'}
#'See Nimble (r-nimble.org) documentation for distribution details.
#' @param Xp Array of detection covariates with dimensions `(R, S, P1)`, where:
#'   - `R`: Number of sites.
#'   - `S`: Number of species.
#'   - `P1`: Number of detection probability covariates.
#' @param Xn Array of abundance covariates with dimensions `(R, S, P2)`, where:
#'   - `R`: Number of sites.
#'   - `S`: Number of species.
#'   - `P2`: Number of abundance covariates.
#' @param verbose Control the level of output displayed during function execution. Default is TRUE.
#' @returns An MNM object that contains the following components:
#'   - summary: Nimble model summary statistics (mean, standard deviation, standard error, quantiles, effective sample size and Rhat value for all monitored values)
#'   - n_parameters: Number of parameters in the model (for use in calculating information criteria).
#'   - data: Observed abundances.
#'   - fitted_Y: Predicted values for Y. Posterior predictive checks can be performed by comparing fitted_Y with the observed data.
#'   - logLik: Log-likelihood of the observed data (Y) given the model parameters.
#'   - n_converged: Number of parameters with successful convergence (Rhat < 1.1).
#'   - plot: traceplots and density plots for all monitored variables.

#' @references
#' - Royle, J. A. (2004). N-mixture models for estimating population size from spatially replicated counts. Biometrics, 60(1), 108-115.
#' - Mimnagh, N., Parnell, A., Prado, E., & Moral, R. D. A. (2022). Bayesian multi-species N-mixture models for unmarked animal communities. Environmental and Ecological Statistics, 29(4), 755-778.
#' @seealso
#' - `simulateData`: For generating example datasets compatible with this function.
#' - `MNM`: For details on creating covariate arrays Xp and Xn.
#'
#' @note
#' Ensure that the dimensions of `Y`, `Xp`, and `Xn` match the requirements specified above. Mismatched dimensions will result in errors during model fitting.
#' @examples
#' # Example 1:
#' Y <- array(rpois(100, lambda = 5), dim = c(10, 5, 2))
#' Xp <- array(runif(100), dim = c(10, 2, 5))
#' Xn <- array(runif(100), dim = c(10, 2, 3))
#'
#' \donttest{model <- MNM_Hurdle(Y = Y, Xp = Xp, Xn = Xn)}
#' # nimble creates auxiliary functions that may be removed after model run is complete
#' # str_objects <- ls(pattern = "^str")
#' # rm(list = str_objects)
#' # Accessing results
#' \donttest{print(model@summary)}
#'
#'#' data(birds)
#'
#' # Example 2: North American Breeding Bird Data
#' # Data must first be reformatted to an array of dimension (R,T,S,K)
#' R <- 15
#' T <- 10
#' S <- 10
#' K <- 4
#' # Ensure data is ordered consistently
#' birds <- birds[order(birds$Route, birds$Year, birds$English_Common_Name), ]
#'
#' # Create a 4D array with proper dimension
#' Y <- array(NA, dim = c(R, T, S, K))
#'
#' # Map route, species, and year to indices
#' route_idx <- as.numeric(factor(birds$Route))
#' species_idx <- as.numeric(factor(birds$English_Common_Name))
#' year_idx <- as.numeric(factor(birds$Year))
#'
#' # Populate the array
#' stop_data <- as.matrix(birds[, grep("^Stop", colnames(birds))])
#'
#' for (i in seq_len(nrow(birds))) {
#'   Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
#'   }
#'
#'   # Assign dimnames
#'   dimnames(Y) <- list(
#'     Route = sort(unique(birds$Route)),
#'       Stop = paste0("Stop", 1:T),
#'         Species = sort(unique(birds$English_Common_Name)),
#'           Year = sort(unique(birds$Year))
#'           )
#'
#' # Selecting only 5 bird species and 1 year for analysis:
#' Y<-Y[,,1:5,1]
#'
#' \donttest{model<-MNM_fit(Y=Y, AR=FALSE, Hurdle=TRUE, iterations=5000, burnin=1000)}
#' @import abind
#' @export
#'
#'
#'
#'
MNM_Hurdle<-function(Y=NULL,iterations=60000, burnin=20000, thin=10,  Xp=NULL, Xn=NULL, verbose=TRUE, ...){
  if(is.null(Y)){
    stop("Error: No data entered. Please provide Y: an array of dimension (R,T,S).")
  }

  if(!all(is.numeric(as.vector(Y)))){
    stop("Error: Non-numeric elements present in Y. Please confirm that Y contains only observational counts.")
  }

  R=dim(Y)[1]
  T=dim(Y)[2]
  S=dim(Y)[3]

  print(paste0("Observations entered correspond to ", R, " sites, ", T, " sampling occasions, and ", S, " species. If this appears to be incorrect, please re-format your data into an array of dimension (R,T,S)."))

   if(burnin>iterations){
    stop("The number of iterations discarded as burn-in cannot be greater than the total number of iterations.")
  }

  if(iterations<5000){
    print(paste0("Warning: Using too few iterations may result in a model that fails to converge."))
  }


  # Capture additional arguments
  additional_args <- list(...)

    # Retrieve model code
    code <- do.call(MNM_control, c(list(model = "Hurdle", Xp = Xp, Xn = Xn), additional_args))
    nimble::nimbleOptions(showCompilerOutput = FALSE, verboseErrors = FALSE, verbose=FALSE, clearNimbleFunctionsAfterCompiling=TRUE, clearCompiled=TRUE)
    model_code <- eval(parse(text = code))

    if(verbose==TRUE){
      print("Building model ...")
    }

    # Dynamically construct the data list and initial values list, and ensure correct dimensionality of Xp and Xn, if present
    data_list <- list(Y = Y, Omega = base::diag(S))
    inits_list <- list(
      theta = 0.5,
      lambda = matrix(10, R, S),
      gamma = rep(0, S),
      a = matrix(0, R, S),
      precision = base::diag(S) * (S+1),
      z = apply(Y, c(1, 3), function(z) ifelse(any(z > 0), 1, 0)),
      count = apply(Y, c(1,3), max)+1)
    constants_list<-list(S = S, R = R, T = T, df = S+1)


    if(!is.null(Xp) & length(dim(Xp))==length(dim(Y))){
      dim_Xp<-dim(Xp)[length(dim(Y))]
      data_list$Xp <- Xp
      inits_list$beta_p <- matrix(0, nrow = dim_Xp, ncol = S)
      constants_list$dim_Xp<-dim_Xp
    } else if(!is.null(Xp) & length(dim(Xp!=length(dim(Y))))){
      stop("Model building stopped: Xp should be an array of dimension (R,S,P) where P is the number of probability-level covariates")
    } else{
      dim_Xp<-0
    }

    if(!is.null(Xn) & length(dim(Xn))==length(dim(Y))){
      dim_Xn<-dim(Xn)[length(dim(Y))]
      data_list$Xn <- Xn  # Include Xn only if it is not NULL
      inits_list$beta_n <- matrix(0, nrow = dim_Xn, ncol = S)
      constants_list$dim_Xn<-dim_Xn
    } else if(!is.null(Xn) & length(dim(Xn!=length(dim(Y))))){
      stop("Model building stopped: Xn should be an array of dimension (R,S,P) where P is the number of abundance-level covariates")
    } else{
      dim_Xn<-0
    }


    # Define the Nimble model
    nimbleModel <- nimble::nimbleModel(code=model_code,
                                       data=data_list,
                                       constants=constants_list,
                                       inits=inits_list,
                                       check=FALSE,
                                       calculate=FALSE)

    # Configure and build the MCMC in the environment
    if(verbose==TRUE){
        print("Building MCMC object ... ")
    }
    mcmcConf <- nimble::configureMCMC(nimbleModel)

    # Dynamically building the list of parameters to monitor
    mcmcConf$addMonitors(c("mu", "correlation", "covariance","theta", "N", "sigma", "probability", "Y_pred", "a","gamma", "precision" ))
    if (dim_Xp > 0) {
      mcmcConf$addMonitors(c("beta_p"))
    }
    if (dim_Xn > 0) {
      mcmcConf$addMonitors(c("beta_n"))
    }
    mcmc <- nimble::buildMCMC(mcmcConf)

    # Compile the model and MCMC objects:
    if(verbose==TRUE){
          print("Compiling model ... ")
    }
    compiledModel <- nimble::compileNimble(nimbleModel)
    if(verbose==TRUE){
         print("Compiling MCMC object ... ")
    }
    compiledMCMC <- nimble::compileNimble(mcmc, project=nimbleModel)


    # Run the MCMC
    if(verbose==TRUE){
        print("Running MCMC object ... ")
    }
    results <- nimble::runMCMC(compiledMCMC,
                               niter=iterations,
                               nburnin=burnin,
                               thin=thin,
                               nchains=4)




  # extract parameter estimates for monitored parameters
  all_samples <- do.call(rbind, results)
  extract_parameter <- function(results, param_name, dim = NULL) {
    all_samples <- do.call(rbind, results)
    param_columns <- grep(paste0("^", param_name), colnames(all_samples))

    if (length(param_columns) == 0) {
      warning(paste("No columns found for parameter:", param_name))
      return(NULL)
    }

    # Ensure the input to colMeans has at least two dimensions
    param_data <- all_samples[, param_columns, drop = FALSE]
    param_means <- colMeans(param_data)

    if (!is.null(dim)) {
      return(array(param_means, dim = dim))
    }
    return(param_means)
  }


  monitored_params <- unique(sub("\\[.*", "", colnames(all_samples)))
  param_means_list <- lapply(monitored_params, function(param) {
    extract_parameter(results, param)
  })

  names(param_means_list) <- monitored_params

  # Convert results to mcmc.list and extract into 3D array for monitoring
  mcmc_list <- coda::mcmc.list(lapply(results, coda::as.mcmc))
  mcmc_array <- as.array(coda::mcmc.list(lapply(results, as.mcmc)))
  mcmc_array <- aperm(mcmc_array, c(1, 3, 2))  # Rearrange dimensions

  # Run the monitor function to get rhat values
  monitor_results <- rstan::monitor(mcmc_array, warmup=0, print=FALSE)
  monitor_results <- as.data.frame(monitor_results)
  monitor_results <- monitor_results[, 1:10]

  # extract parameters for log likelihood calculation
  N<-array(round(param_means_list$N), dim=c(R,S))
  prob<-array(param_means_list$probability, dim=c(R,S))
  total_log_likelihood<-0

  calculate_log_likelihood <- function() {
    for (s in 1:S) {
      for (i in 1:R) {
        for (t in 1:T) {
          total_log_likelihood <- total_log_likelihood +
            stats::dbinom(Y[i, t, s], size = round(N[i,s]), prob = prob[i,s], log = TRUE)
        }
      }
    }
    return(total_log_likelihood)
  }
  log_likelihood <- as.numeric(calculate_log_likelihood())
  generate_plot_data <- function(results, param_names) {
    # Convert results to mcmc.list
    mcmc_list <- coda::mcmc.list(lapply(results, coda::as.mcmc))

    # Create a list to store samples
    plot_data <- list()

    for (param in param_names) {
      # Extract samples for the parameter
      param_samples <- mcmc_list[[1]][, param]

      # Save samples in the list
      plot_data[[param]] <- param_samples
    }

    return(plot_data)
  }

  # Generate plot data
  param_names <- colnames(mcmc_list[[1]])
  plot_data <- generate_plot_data(results, param_names)

  generate_plot_functions <- function(plot_data) {
    plot_functions <- list()

    for (param in names(plot_data)) {
      # Extract samples for the parameter
      param_samples <- plot_data[[param]]

      # Create a closure to lock in the current parameter and its samples
      trace_plot <- local({
        samples <- param_samples
        param_name <- param
        function() {
          plot(
            samples, type = "l",
            xlab = "Iteration", ylab = "Value",
            main = paste("Trace Plot for", param_name)
          )
        }
      })

      density_plot <- local({
        samples <- param_samples
        param_name <- param
        function() {
          plot(
            stats::density(samples),
            xlab = "Value", ylab = "Density",
            main = paste("Density Plot for", param_name)
          )
        }
      })

      # Store the functions
      plot_functions[[param]] <- list(trace = trace_plot, density = density_plot)
    }

    return(plot_functions)
  }



  # Generate plot functions
  plot_functions <- generate_plot_functions(plot_data)



  mnm_object <- new("MNM",
                    summary = as.data.frame(monitor_results),
                    n_parameters=2*S+1+(R*S)+((S*(S+1))/2)+(S*dim_Xp)+(S*dim_Xn),
                    estimates=param_means_list,
                    data = Y,
                    plot=plot_functions,
                    fitted_Y = array(round(param_means_list$Y_pred), dim = c(R, T, S)),
                    logLik = log_likelihood,
                    n_converged = sum(monitor_results[, "Rhat"] < 1.1) )
  return(mnm_object)
}
