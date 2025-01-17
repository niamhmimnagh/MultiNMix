#' @title Generate Nimble Code for Multi-Species N-Mixture (MNM) Models
#'
#' @description This function generates Nimble code for fitting various Multi-Species N-Mixture (MNM) models. These models can include standard MNM, hurdle components, autoregressive (AR) components, or both hurdle and AR components.
#'
#' @details The generated Nimble code can be used to fit Bayesian hierarchical MNM models. Users can specify prior distributions for key model parameters and provide covariate arrays that influence detection probability and abundance. The function validates prior specifications and adapts model code based on the selected model type. Supported models include:
#'   - **MNM**: Standard MNM model
#'   - **Hurdle**: MNM model with a hurdle component
#'   - **AR**: MNM model with an autoregressive component
#'   - **HurdleAR**: MNM model with both hurdle and AR components
#'
#'   This is an internal function, needed for implementing MNM models using the MNM_fit function. While the MNM_fit function allows the user to define prior distributions and linear covariates on detection probability and abundance, if the user wishes to implement more complex models in Nimble (via unsupported prior distributions or non-linear covariate effects), the Nimble model code may be extracted as in examples below, and modified by the user.
#'
#' @param model Character. Specifies the type of MNM model to generate. Options are:
#'   - `"MNM"`: Standard MNM model (default)
#'   - `"Hurdle"`: Hurdle MNM model
#'   - `"AR"`: Autoregressive MNM model
#'   - `"HurdleAR"`: Hurdle MNM model with autoregression
#' @param prior_detection_probability Character. Prior distribution for the detection probability intercept (`gamma`). Default is `'dnorm(0, 0.001)'`.
#' @param prior_precision Character. Prior distribution for the precision matrix for the species-level random effect. Default is `'dwish(Omega[1:S,1:S], df)'`.
#' @param prior_mean Character. Prior distribution for the mean of the species-level random effect (`mu`). Default is `'dnorm(0,0.001)'`.
#' @param prior_hurdle Character. Prior distribution for `theta`, the probability of structural zero in hurdle models. Default is `'dbeta(1,1)'`.
#' @param prior_mean_AR Character. Prior distribution for the mean of the autoregressive random effect (`phi`). Default is `'dnorm(0,0.001)'`.
#' @param prior_sd_AR Character. Prior distribution for the standard deviation of the autoregressive random effect (`phi`). Default is `'dexp(1)'`.
#' @param Xp Array. Covariates influencing detection probability. Dimensions depend on the model:
#'   - \eqn{R \times S \times P} for models without AR
#'   - \eqn{R \times S \times K \times P} for AR models
#'   where:
#'   - \eqn{R}: Number of sites
#'   - \eqn{S}: Number of species
#'   - \eqn{K}: Number of time points
#'   - \eqn{P}: Number of covariates
#' @param Xn Array. Covariates influencing abundance. Dimensions depend on the model:
#'   - \eqn{R \times S \times P} for models without AR
#'   - \eqn{R \times S \times K \times P} for AR models
#' @return Character. Nimble code for the specified MNM model, which can be used for further analysis or fitting.
#'
#' @examples
#' # Example
#' # In order to implement scenarios involving nonlinear covariate
#' # effects or complex interaction terms,
#' # the user is invited to extract and modify the MNM Nimble code:
#' model<-MNM_control(model="MNM")
#' cat(model)
#'
#' @import abind
#' @import nimble
#' @export
MNM_control <- function(model = "MNM",
                        prior_detection_probability = 'dnorm(0, 0.001)',
                        prior_precision = 'dwish(Omega[1:S,1:S], df)',
                        prior_mean = 'dnorm(0,0.001)',
                        prior_mean_AR='dnorm(0,0.001)',
                        prior_sd_AR='dexp(1)',
                        prior_hurdle='dbeta(1,1)',
                        Xp = NULL,
                        Xn = NULL) {

  prior_list<-list(prior_detection_probability = prior_detection_probability,
                   prior_precision = prior_precision,
                   prior_mean = prior_mean,
                   prior_mean_AR = prior_mean_AR,
                   prior_sd_AR = prior_sd_AR,
                   prior_hurdle=prior_hurdle)
  validate_all_priors(prior_list)

  # Initialize dimensions
  if (model == "MNM" | model == "Hurdle") {
    dim_Xp <- if (!is.null(Xp)) dim(Xp)[3] else 0
    dim_Xn <- if (!is.null(Xn)) dim(Xn)[3] else 0
  } else if (model == "AR" | model == "HurdleAR") {
    dim_Xp <- if (!is.null(Xp)) dim(Xp)[4] else 0
    dim_Xn <- if (!is.null(Xn)) dim(Xn)[4] else 0
  }



  prob_covariates <- if (dim_Xp > 1) {
    'logit(probability[i,s]) <- gamma[s] + sum(beta_p[1:dim_Xp, s] * Xp[i,s,1:dim_Xp])'
  } else if (dim_Xp == 1) {
    'logit(probability[i,s]) <- gamma[s] + beta_p[1, s] * Xp[i,s,1]'
  } else {
    'logit(probability[i,s]) <- gamma[s]'
  }


  prob_beta <- if (dim_Xp > 0) {
    'for (s in 1:S) {
      for (p in 1:dim_Xp) {
        beta_p[p, s] ~ dnorm(0, 0.01)
      }
    }'
  } else {
    ''
  }



  if(model=="MNM"){
    prob_covariates <- if (dim_Xp > 1) {
      'logit(probability[i,s]) <- gamma[s] + sum(beta_p[1:dim_Xp, s] * Xp[i,s,1:dim_Xp])'
    } else if (dim_Xp == 1) {
      'logit(probability[i,s]) <- gamma[s] + beta_p[1, s] * Xp[i,s,1]'
    } else {
      'logit(probability[i,s]) <- gamma[s]'
    }


    prob_beta <- if (dim_Xp > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xp) {
        beta_p[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }

    abund_covariates <- if (dim_Xn > 1) {
      'log(lambda[i,s]) <- a[i,s] + sum(beta_n[1:dim_Xn, s] * Xn[i,s,1:dim_Xn])'
    } else if (dim_Xn == 1) {
      'log(lambda[i,s]) <- a[i,s] + beta_n[1, s] * Xn[i,s,1]'
    } else {
      'log(lambda[i,s]) <- a[i,s]'
    }


    abund_beta <- if (dim_Xn > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xn) {
        beta_n[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }

  # Generate model code
    code <- paste0('nimbleCode({
      for (s in 1:S) {
        # Priors for detection probability coefficients


        for (i in 1:R) {
          N[i, s] ~ dpois(lambda[i, s])
        ', abund_covariates, '
        ', prob_covariates, '
          for (t in 1:T) {
            Y[i, t, s] ~ dbin(probability[i,s], N[i, s])
            Y_pred[i, t, s] ~ dbin(probability[i,s], N[i, s])
          }
        }
      }

      # Random effects a with mean vector mu and variance-covariance matrix
      for (i in 1:R) {
        a[i, 1:S] ~ dmnorm(mu[1:S], precision[1:S, 1:S])
      }

      # Wishart prior on precision
      precision[1:S, 1:S] ~ ', prior_precision, '
      covariance[1:S, 1:S] <- inverse(precision[1:S, 1:S])

      # Correlations and standard deviations
      for (s in 1:S) {
        sigma[s] <- sqrt(covariance[s, s])
        for (s1 in 1:S) {
          correlation[s, s1] <- covariance[s, s1] / sqrt(covariance[s, s] * covariance[s1, s1])
        }
      }

      # Priors for beta_p and beta_n (species-level parameters)
      # Priors
      ', prob_beta, '

      for (s in 1:S) {
        gamma[s] ~ ', prior_detection_probability, '
        mu[s] ~ ', prior_mean, '
      }
     ', abund_beta, '
      }
  )'
    )}




  else if(model=="Hurdle"){
    prob_covariates <- if (dim_Xp > 1) {
      'logit(probability[i,s]) <- gamma[s] + sum(beta_p[1:dim_Xp, s] * Xp[i,s,1:dim_Xp])'
    } else if (dim_Xp == 1) {
      'logit(probability[i,s]) <- gamma[s] + beta_p[1, s] * Xp[i,s,1]'
    } else {
      'logit(probability[i,s]) <- gamma[s]'
    }


    prob_beta <- if (dim_Xp > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xp) {
        beta_p[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }

    abund_covariates <- if (dim_Xn > 1) {
      'log(lambda[i,s]) <- a[i,s] + sum(beta_n[1:dim_Xn, s] * Xn[i,s,1:dim_Xn])'
    } else if (dim_Xn == 1) {
      'log(lambda[i,s]) <- a[i,s] + beta_n[1, s] * Xn[i,s,1]'
    } else {
      'log(lambda[i,s]) <- a[i,s]'
    }


    abund_beta <- if (dim_Xn > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xn) {
        beta_n[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }
  code<-paste0('
    nimbleCode({
  # Likelihood
  for(s in 1:S){
    for (i in 1:R) {
        # Zero-altered abundance model
        ', abund_covariates, '
        ', prob_covariates, '
        z[i,s] ~ dbern(theta)               # Binary indicator for structural zeros
        count[i,s] ~ T(dpois(lambda[i,s]),1,)             # Poisson abundance for ZAP
        N[i,s] <- z[i,s] * count[i,s]        # Apply zero-alteration

      for(t in 1:T){
        Y[i,t,s] ~ dbin(probability[i,s],N[i,s])
        Y_pred[i,t,s]~ dbin(probability[i,s],N[i,s])

      }
    }
  }

  # Random effects a with mean vector mu, and variance-covariance matrix cov
      for (i in 1:R) {
        a[i, 1:S] ~ dmnorm(mu[1:S], precision[1:S, 1:S])
      }


  # Wishart prior on precision with df=S+1 and diagonal matrix Omega
  precision[1:S,1:S] ~', prior_precision,'
  covariance[1:S,1:S]<-inverse(precision[1:S, 1:S])

    # Correlations and standard deviations
      for (s in 1:S) {
        sigma[s] <- sqrt(covariance[s, s])
        for (s1 in 1:S) {
          correlation[s, s1] <- covariance[s, s1] / sqrt(covariance[s, s] * covariance[s1, s1])
        }
      }



  theta~', prior_hurdle,'
  ', prob_beta, '
  for(s in 1:S){
    gamma[s] ~',prior_detection_probability, '
    mu[s] ~', prior_mean,'
  }
    ', abund_beta, '
  })
')}

  else if(model=="AR"){
    prob_covariates <- if (dim_Xp > 1) {
      'logit(probability[i,s, k]) <- gamma[s] + sum(beta_p[1:dim_Xp, s] * Xp[i,s,k,1:dim_Xp])'
    } else if (dim_Xp == 1) {
      'logit(probability[i,s,k]) <- gamma[s] + beta_p[1, s] * Xp[i,s,k,1]'
    } else {
      'logit(probability[i,s,k]) <- gamma[s]'
    }


    prob_beta <- if (dim_Xp > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xp) {
        beta_p[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }

    abund_covariates_k <- if (dim_Xn > 1) {
      'log(lambda[i,s,k]) <- a[i,s] + sum(beta_n[1:dim_Xn, s] * Xn[i,s,k,1:dim_Xn])+ phi[s]*log(N[i,s,k-1]+1)'
    } else if (dim_Xn == 1) {
      'log(lambda[i,s,k]) <- a[i,s] + beta_n[1, s] * Xn[i,s,k,1]+ phi[s]*log(N[i,s,k-1]+1)'
    } else {
      'log(lambda[i,s,k]) <- a[i,s]+ phi[s]*log(N[i,s,k-1]+1)'
    }

    abund_covariates_1 <- if (dim_Xn > 1) {
      'log(lambda[i,s,1]) <- a[i,s] + sum(beta_n[1:dim_Xn, s] * Xn[i,s,1,1:dim_Xn])'
    } else if (dim_Xn == 1) {
      'log(lambda[i,s,1]) <- a[i,s] + beta_n[1, s] * Xn[i,s,1,1]'
    } else {
      'log(lambda[i,s,1]) <- a[i,s]'
    }

    abund_beta <- if (dim_Xn > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xn) {
        beta_n[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }


    code<-paste0('
nimbleCode({
  for(s in 1:S){
    for (i in 1:R) {
      for(k in 1:K){
          ', prob_covariates, '

      }
      ', abund_covariates_1, '
      N[i,s,1]~dpois(lambda[i,s,1])

      for(k in 2:K){
        ', abund_covariates_k, '
         N[i,s,k] ~ dpois(lambda[i,s,k])
      }
    }
  }

        # Loop over time points
        for(i in 1:R){
          for(s in 1:S){
            for(t in 1:T){
              for(k in 1:K){
                Y[i,t,s,k] ~ dbin(probability[i,s,k], N[i,s,k])
                Y_pred[i,t,s,k]~ dbin(probability[i,s,k],N[i,s,k])
              }
            }
        }
     }

  # Random effects a with mean vector mu, and variance-covariance matrix cov
  for(i in 1:R){
    a[i,1:S] ~ dmnorm(mu[1:S], precision[1:S,1:S])
  }


  # Wishart prior on precision with df=S+1 and diagonal matrix Omega
  precision[1:S,1:S] ~', prior_precision,'
  covariance[1:S,1:S]<-inverse(precision[1:S,1:S])

      # Correlations and standard deviations
      for (s in 1:S) {
        sigma[s] <- sqrt(covariance[s, s])
        for (s1 in 1:S) {
          correlation[s, s1] <- covariance[s, s1] / sqrt(covariance[s, s] * covariance[s1, s1])
        }
      }


  # Priors
  ', prob_beta, '
  ', abund_beta, '
  for(s in 1:S){
    gamma[s] ~',prior_detection_probability, '
    mu[s] ~', prior_mean,'
    phi[s] ~ dnorm(muPhi,tauPhi)
  }

  muPhi~',prior_mean_AR,'
  sdPhi~',prior_sd_AR,'
  tauPhi<-pow(sdPhi, -2)
})')}


  else if(model=="HurdleAR"){
    prob_covariates <- if (dim_Xp > 1) {
      'logit(probability[i,s, k]) <- gamma[s] + sum(beta_p[1:dim_Xp, s] * Xp[i,s,k,1:dim_Xp])'
    } else if (dim_Xp == 1) {
      'logit(probability[i,s,k]) <- gamma[s] + beta_p[1, s] * Xp[i,s,k,1]'
    } else {
      'logit(probability[i,s,k]) <- gamma[s]'
    }


    prob_beta <- if (dim_Xp > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xp) {
        beta_p[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }

    abund_covariates_k <- if (dim_Xn > 1) {
      'log(lambda[i,s,k]) <- a[i,s] + sum(beta_n[1:dim_Xn, s] * Xn[i,s,k,1:dim_Xn])+ phi[s]*log(N[i,s,k-1]+1)'
    } else if (dim_Xn == 1) {
      'log(lambda[i,s,k]) <- a[i,s] + beta_n[1, s] * Xn[i,s,k,1]+ phi[s]*log(N[i,s,k-1]+1)'
    } else {
      'log(lambda[i,s,k]) <- a[i,s]+ phi[s]*log(N[i,s,k-1]+1)'
    }

    abund_covariates_1 <- if (dim_Xn > 1) {
      'log(lambda[i,s,1]) <- a[i,s] + sum(beta_n[1:dim_Xn, s] * Xn[i,s,1,1:dim_Xn])'
    } else if (dim_Xn == 1) {
      'log(lambda[i,s,1]) <- a[i,s] + beta_n[1, s] * Xn[i,s,1,1]'
    } else {
      'log(lambda[i,s,1]) <- a[i,s]'
    }

    abund_beta <- if (dim_Xn > 0) {
      'for (s in 1:S) {
      for (p in 1:dim_Xn) {
        beta_n[p, s] ~ dnorm(0, 0.01)
      }
    }'
    } else {
      ''
    }

    code<-paste0('
nimbleCode({
  for(s in 1:S){
    for (i in 1:R) {
      for(k in 1:K){
          ', prob_covariates, '
      }
      ', abund_covariates_1, '
        z[i,s,1] ~ dbern(theta)               # Binary indicator for structural zeros
        count[i,s,1] ~ T(dpois(lambda[i,s,1]),1,)             # Poisson abundance for ZAP
        N[i,s,1] <- z[i,s,1] * count[i,s,1]        # Apply zero-alteration


      for(k in 2:K){
      ', abund_covariates_k, '
        z[i,s,k] ~ dbern(theta)               # Binary indicator for structural zeros
        count[i,s,k] ~ T(dpois(lambda[i,s,k]),1,)             # Poisson abundance for ZAP
        N[i,s,k] <- z[i,s,k] * count[i,s,k]        # Apply zero-alteration
      }
    }
  }

        # Loop over time points
        for(i in 1:R){
          for(s in 1:S){
            for(t in 1:T){
              for(k in 1:K){
                Y[i,t,s,k] ~ dbin(probability[i,s,k], N[i,s,k])
                Y_pred[i,t,s,k]~ dbin(probability[i,s,k],N[i,s,k])

              }
            }
        }
     }

  # Random effects a with mean vector mu, and variance-covariance matrix cov
  for(i in 1:R){
    a[i,1:S] ~ dmnorm(mu[1:S], precision[1:S,1:S])
  }


  # Wishart prior on precision with df=S+1 and diagonal matrix Omega
  precision[1:S,1:S] ~', prior_precision,'
  covariance[1:S,1:S]<-inverse(precision[1:S,1:S])


      # Correlations and standard deviations
      for (s in 1:S) {
        sigma[s] <- sqrt(covariance[s, s])
        for (s1 in 1:S) {
          correlation[s, s1] <- covariance[s, s1] / sqrt(covariance[s, s] * covariance[s1, s1])
        }
      }

  theta~', prior_hurdle,'
  ', prob_beta, '
  ', abund_beta, '
  # Priors
  for(s in 1:S){
    gamma[s] ~',prior_detection_probability, '
    mu[s] ~', prior_mean,'
    phi[s] ~ dnorm(muPhi,tauPhi)
  }

  muPhi~',prior_mean_AR,'
  sdPhi~',prior_sd_AR,'
  tauPhi<-pow(sdPhi, -2)
})')
  }

else{
  print("Unrecognised model specification. Possible model specifications include 'MNM', 'Hurdle', 'AR', 'HurdleAR'")
}
  return(code)
}

