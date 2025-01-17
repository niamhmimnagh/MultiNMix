#' @title Simulate Data for Multi-Species N-Mixture Models
#'
#' @description Simulates multi-species correlated abundance data for various Multi-Species N-Mixture (MNM) model types, including standard MNM, Hurdle, AR (autoregressive), and HurdleAR models.
#' @details This function generates abundance data for multi-species N-mixture models under different configurations:
#'   - **MNM**: Standard multi-species N-mixture model.
#'   - **Hurdle**: Includes a hurdle component to model zero-inflated data.
#'   - **AR**: Includes an autoregressive (AR) component for temporal dependencies.
#'   - **HurdleAR**: Combines hurdle and AR components for zero-inflation and temporal dependencies.
#' The output includes observed and true abundances, detection probabilities, latent variables, and covariance information for the random effects.
#'
#' @param model Character. Specifies the model type. Options are `"MNM"`, `"Hurdle"`, `"AR"`, `"HurdleAR"`. Default is `"MNM"`.
#' @param R Integer. Number of sites. Default is `10`.
#' @param T Integer. Number of replicates. Default is `5`.
#' @param S Integer. Number of species. Default is `2`.
#' @param K Integer. Number of time points (used for AR models). Default is `4`.
#' @param prob Character. Specifies the range of detection probabilities:
#'   - `"small"`: Detection probabilities < 0.4.
#'   - `"large"`: Detection probabilities > 0.5.
#'   - `"all"`: Detection probabilities between 0.01 and 0.99 (default).
#' @param abundance Character. Specifies the abundance size:
#'   - `"small"`: Latent species abundance between 0 and 50.
#'   - `"large"`: Latent species abundance between 0 and 700.
#'   Default is `"small"`.
#' @param theta Numeric. Probability of zero-inflation (used for hurdle models). Default is `0.5`.
#'
#' @return A list containing:
#'   - **Y**: Array of observed abundances.
#'   - **N**: Array of true abundances.
#'   - **p**: Array of detection probabilities.
#'   - **Sigma**: Covariance matrix for the multivariate normal variable `a`.
#'   - **mu**: Mean vector for the multivariate normal variable `a`.
#'   - **lambda**: Latent abundance rate parameter.
#'   - **correlation**: Correlation matrix derived from `Sigma`.
#'   - **R, T, S, K**: Number of sites, sampling occasions, species, and time points.
#'   - Additional elements depending on the model type:
#'     - **phi**: Autoregression parameter (AR and HurdleAR models).
#'     - **muPhi**: Mean of the autoregressive parameter (AR and HurdleAR models).
#'     - **varPhi**: Variance of the autoregressive parameter (AR and HurdleAR models).
#'     - **zeros**: Matrix of zero-indicators for hurdle models.
#'     - **theta**: Zero-inflation parameter for hurdle models.
#' @examples
#' # Simulate data for a standard MNM model
#' data <- simulateData(model = "MNM", R = 10, S = 3, T = 5, prob = "all",
#' abundance = "small")
#'
#' # Simulate data for a hurdle model
#' data <- simulateData(model = "Hurdle", R = 10, S = 3, T = 5, prob = "large",
#' abundance = "large", theta = 0.3)
#'
#' # Simulate data for an autoregressive model
#' data <- simulateData(model = "AR", R = 10, S = 2, T = 5, K = 4, prob = "small",
#' abundance = "small")
#'
#' # Simulate data for a hurdle autoregressive model
#' data <- simulateData(model = "HurdleAR", R = 10, S = 3, T = 5, K = 4, prob = "all",
#' abundance = "large", theta = 0.5)
#'
#' @seealso
#' - `simulateData_MNM`: Helper function for simulating standard MNM data.
#' - `simulateData_Hurdle`: Helper function for simulating hurdle MNM data.
#' - `simulateData_AR`: Helper function for simulating AR MNM data.
#' - `simulateData_Hurdle_AR`: Helper function for simulating hurdle AR MNM data.
#'
#' @export
#'
simulateData<-function(model="MNM", R=10, S=2, T=5, prob="all", abundance="small", K=4, theta=0.5){
  if(model=="MNM"){
    simulateData_MNM(S, R, T, prob, abundance)
  }
  else if (model=="Hurdle"){
    simulateData_Hurdle(S,R,T, prob, abundance, theta)
  }
  else if (model=="AR"){
    simulateData_AR(S,R,T,K, prob, abundance)
  }
  else if (model=="HurdleAR"){
    simulateData_Hurdle_AR(S,R,T,K, prob, abundance, theta)
  }
}

simulateData_MNM<-function(S,R,T, prob, abundance){
  N<-lambda<-matrix(ncol=S, nrow=R)
  Y <-  array(NA, dim = c(R,T,S), dimnames = list(NULL, NULL, paste("species", 1:S)))

  # Probability of detection
  if(prob=="small") p <- array(stats::runif(R*S, 0.01, 0.4), dim=c(R,S))
  else if (prob=="large") p <- array(stats::runif(R*S, 0.5, 0.99), dim=c(R,S))
  else if (prob=="all") p<-array(stats::runif(R*S, 0.01, 0.99), dim=c(R,S))

  # Mean and Covariance matrix for MVN variable a
  Sigma<-clusterGeneration::genPositiveDefMat(S, rangeVar=c(0.2, 1), covMethod="unifcorrmat")[["Sigma"]]
  correlation<-stats::cov2cor(Sigma)
  mu<-rep(ifelse(abundance=="small", 2, 4),S)

  a<-mvtnorm::rmvnorm(R, mean=mu, sigma=Sigma)
  lambda<-exp(a)

  ## observed and true abundances
  for(i in 1:R){
    for(s in 1:S){
      N[i,s]<-stats::rpois(1, lambda[i,s])
      for(t in 1:T){
        Y[i,t,s]<-stats::rbinom(1, N[i,s], p[i,s])
      }
    }
  }

  ylist<-list("Y"=Y,
              "N"=N,
              "R"=R,
              "T"=T,
              "S"=S,
              "p"=p,
              "Sigma"=Sigma,
              "lambda"=lambda,
              "a"=a,
              "mu"=mu,
              "correlation"=correlation)
  return(ylist)
}


simulateData_Hurdle <- function(S,R,T, prob, abundance, theta){
  N<- matrix(ncol=S, nrow=R)
  Y<-array(NA, dim = c(R,T,S), dimnames = list(NULL, NULL, paste("species", 1:S)))


  zeroInflation <- matrix(stats::rbinom(R*S, size=1, prob=1-theta), ncol=S, nrow=R)


  if(prob=="small") p <- array(stats::runif(R*S, 0.01, 0.4), dim=c(R,S))
  else if (prob=="large") p <- array(stats::runif(R*S, 0.5, 0.99), dim=c(R,S))
  else if (prob=="all") p<-array(stats::runif(R*S, 0.01, 0.99), dim=c(R,S))


  mu<-rep(ifelse(abundance=="small", 2, 4),S)



  Sigma <- clusterGeneration::genPositiveDefMat(S, rangeVar=c(0.2, 1), covMethod="unifcorrmat")[["Sigma"]]
  correlation <- stats::cov2cor(Sigma)


  a <- mvtnorm::rmvnorm(R, mean=mu, sigma=Sigma)
  lambda <- exp(a)

  ## observed and true abundances
  for(i in 1:R){
    for(s in 1:S){
      N[i,s] <- ifelse(zeroInflation[i,s]==0,
                       0,
                       extraDistr::rtpois(n=1, lambda=lambda[i,s],a=0))


      for(t in 1:T){
        Y[i,t,s] <- stats::rbinom(1, N[i,s], p[i,s])
      }
    }
  }

  ylist<-list("Y"=Y,
              "N"=N,
              "R"=R,
              "T"=T,
              "S"=S,
              "p"=p,
              "Sigma"=Sigma,
              "lambda"=lambda,
              "a"=a,
              "mu"=mu,
              "correlation"=correlation,
              "zeros"=zeroInflation,
              "theta"=theta)
  return(ylist)
}



simulateData_AR<-function(S,R,T,K, prob, abundance){
  Y<-array(NA, dim=c(R,T,S,K), dimnames=list(NULL, NULL, paste("species", 1:S), NULL))
  N<-lambda<-array(NA, dim=c(R,S,K), dimnames=list(NULL, paste("species", 1:S), NULL))
  gamma<-vector(length=S)


  # Probability of detection
  if(prob=="small") p <- array(stats::runif(R*S, 0.01, 0.4), dim=c(R,S))
  else if (prob=="large") p <- array(stats::runif(R*S, 0.5, 0.99), dim=c(R,S))
  else if (prob=="all") p<-array(stats::runif(R*S, 0.01, 0.99), dim=c(R,S))

  # Normally dist. autocorrelation coefficient
  muPhi<-0
  sigmaPhi<-0.2
  phi<-stats::rnorm(S, muPhi, sigmaPhi)


  # MVN random effect
  Omega<-diag(1, nrow=S, ncol=S) # Scale matrix for wishart distribution
  Sigma<-clusterGeneration::genPositiveDefMat(S, rangeVar=c(0.2, 1), covMethod="unifcorrmat")[["Sigma"]]
  correlation<-stats::cov2cor(Sigma)
  mu<-rep(ifelse(abundance=="small", 2, 4),S)

  # Standard Deviations
  sigma<-vector(length=S)
  for(i in 1:S){
    sigma[i]<-sqrt(Sigma[i,i])
  }


  a <- mvtnorm::rmvnorm(R, mean=mu, sigma=Sigma)


  for(i in 1:R){
    for(s in 1:S){
      lambda[i,s,1]<-exp(a[i,s])
      N[i,s,1]<-stats::rpois(1, lambda[i,s,1])

      for(k in 2:K){
        lambda[i,s,k]<-exp(a[i,s]+phi[s]*log(N[i,s,k-1]+1))
        N[i,s,k]<-stats::rpois(1, lambda[i,s,k])
      }
    }
  }

  for(i in 1:R){
    for(t in 1:T){
      for(s in 1:S){
        for(k in 1:K){
          Y[i,t,s,k]<-stats::rbinom(1, size=N[i,s,k], prob=p[i,s])
        }
      }
    }
  }

  ylist<-list("Y"=Y,
              "N"=N,
              "R"=R,
              "T"=T,
              "S"=S,
              "K"=K,
              "p"=p,
              "Sigma"=Sigma,
              "lambda"=lambda,
              "a"=a,
              "mu"=mu,
              "correlation"=correlation,
              "phi"=phi,
              "muPhi"=muPhi,
              "varPhi"=sigmaPhi)
  return(ylist)
}



simulateData_Hurdle_AR<-function(S,R,T,K, prob, abundance, theta){
  Y<-array(NA, dim=c(R,T,S,K), dimnames=list(NULL, NULL, paste("species", 1:S), NULL))
  N<-lambda<-array(NA, dim=c(R,S,K), dimnames=list(NULL, paste("species", 1:S), NULL))


  # Zero-counts
  zeroInflation <- array(stats::rbinom(R*S*K, size=1, prob=1-theta), dim=c(R,S,K))

  # Probability of Detection
  if(prob=="small") p <- array(stats::runif(R*S, 0.01, 0.4), dim=c(R,S))
  else if (prob=="large") p <- array(stats::runif(R*S, 0.5, 0.99), dim=c(R,S))
  else if (prob=="all") p<-array(stats::runif(R*S, 0.01, 0.99), dim=c(R,S))


  # Uniform dist. autocorrelation coefficient
  muPhi<-0
  sigmaPhi<-0.2
  phi<-stats::rnorm(S, muPhi, sigmaPhi)


  # MVN random effect
  Omega<-diag(1, nrow=S, ncol=S) # Scale matrix for wishart distribution

  Sigma<-clusterGeneration::genPositiveDefMat(S, rangeVar=c(0.2, 1), covMethod="unifcorrmat")[["Sigma"]]
  correlation<-stats::cov2cor(Sigma)
  mu<-rep(ifelse(abundance=="small", 2, 4),S)

  # Standard Deviations
  sigma<-vector(length=S)
  for(i in 1:S){
    sigma[i]<-sqrt(Sigma[i,i])
  }


  a <- mvtnorm::rmvnorm(R, mean=mu, sigma=Sigma)


  for(i in 1:R){
    for(s in 1:S){
      lambda[i,s,1]<-exp(a[i,s])
      N[i,s,1]<-ifelse(zeroInflation[i,s,1]==0,
                       0,
                       extraDistr::rtpois(n=1, lambda=lambda[i,s,1],a=0))

      for(k in 2:K){
        lambda[i,s,k]<-exp(a[i,s]+phi[s]*log(N[i,s,k-1]+1))
        N[i,s,k]<-ifelse(zeroInflation[i,s,k]==0,
                         0,
                         extraDistr::rtpois(n=1, lambda=lambda[i,s,k],a=0))
      }
    }
  }


  for(i in 1:R){
    for(t in 1:T){
      for(s in 1:S){
        for(k in 1:K){
          Y[i,t,s,k]<-stats::rbinom(1, size=N[i,s,k], prob=p[i,s])
        }
      }
    }
  }

  ylist<-list("Y"=Y,
              "N"=N,
              "R"=R,
              "T"=T,
              "S"=S,
              "K"=K,
              "p"=p,
              "Sigma"=Sigma,
              "sigma"=sigma,
              "theta"=theta,
              "zeros"=zeroInflation,
              "lambda"=lambda,
              "a"=a,
              "mu"=mu,
              "correlation"=correlation,
              "phi"=phi,
              "muPhi"=muPhi,
              "varPhi"=sigmaPhi)
  return(ylist)
}



