pkgname <- "MultiNMix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "MultiNMix-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('MultiNMix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MNM")
### * MNM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MNM
### Title: Fit a Multi-Species N-Mixture Model (MNM) using Nimble
### Aliases: MNM

### ** Examples

# Example 1:
# Covariates must be of dimension (R, S, P1/P2). If covariates of an alternative dimension
# are used, they must first be coerced into the right format.
# If we have two abundance-covariates, one site-level covariate and one species-level
# covariate, they may be combined as follows:
R <- 10  # Number of sites
S <- 5   # Number of species
T<-5
Y <- array(sample(0:10, 100, replace = TRUE), dim = c(R, T, S))
covariate_1 <- runif(R)  # Site-level covariate
covariate_2 <- runif(S)  # Species-level covariate
# Expand covariate_1 to have S columns
expanded_covariate_1 <- matrix(rep(covariate_1, S), nrow = R, ncol = S)
# Expand covariate_2 to have R rows
expanded_covariate_2 <- t(matrix(rep(covariate_2, R), nrow = S, ncol = R))
# Combine into an array of dimensions (R, S, 2)
Xn <- array(c(expanded_covariate_1, expanded_covariate_2), dim = c(R, S, 2))
dim(Xn) # this is now in the correct format and can be used.
## Not run: result <- MNM(Y, Xn = Xn)
## Not run: print(result@summary)
#' data(birds_raw)
# Example 2: North American Breeding Bird Data
# Data must first be reformatted to an array of dimension (R,T,S,K)
R <- 24
T <- 10
S <- 20
K <- 6
# Ensure data is ordered consistently
birds_raw <- birds_raw[order(birds_raw$Route, birds_raw$Year, birds_raw$English_Common_Name), ]
# Create a 4D array with proper dimension
Y <- array(NA, dim = c(R, T, S, K))
# Map route, species, and year to indices
route_idx <- as.numeric(factor(birds_raw$Route))
species_idx <- as.numeric(factor(birds_raw$English_Common_Name))
year_idx <- as.numeric(factor(birds_raw$Year))
# Populate the array
stop_data <- as.matrix(birds_raw[, grep("^Stop", colnames(birds_raw))])
for (i in seq_len(nrow(birds_raw))) {
  Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
  }
  # Assign dimnames
  dimnames(Y) <- list(
    Route = sort(unique(birds_raw$Route)),
      Stop = paste0("Stop", 1:T),
        Species = sort(unique(birds_raw$English_Common_Name)),
          Year = sort(unique(birds_raw$Year))
          )
# Selecting only 5 bird species and 1 year for analysis:
Y<-Y[,,1:5,1]
## Not run: model<-MNM_fit(Y=Y, AR=FALSE, Hurdle=FALSE, iterations=5000, burnin=1000)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MNM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MNM_AR")
### * MNM_AR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MNM_AR
### Title: Fit a Multi-Species N-Mixture Model with AR-1 Using Nimble
### Aliases: MNM_AR

### ** Examples

# Example 1:
Y <- array(rpois(1000, lambda = 10), dim = c(10, 10, 5, 2))
Xp <- array(runif(500), dim = c(10, 5, 2, 3))
Xn <- array(runif(1000), dim = c(10, 5, 2, 4))

# Fit the AR-1 model
## Not run: result <- MNM_AR(Y = Y, Xp = Xp, Xn = Xn)

# Check fitted vs observed abundance
## Not run: plot(result@data, result@fitted_Y)

#' data(birds_raw)

# Example 2: North American Breeding Bird Data
# Data must first be reformatted to an array of dimension (R,T,S,K)
R <- 24
T <- 10
S <- 20
K <- 6
# Ensure data is ordered consistently
birds_raw <- birds_raw[order(birds_raw$Route, birds_raw$Year, birds_raw$English_Common_Name), ]

# Create a 4D array with proper dimension
Y <- array(NA, dim = c(R, T, S, K))

# Map route, species, and year to indices
route_idx <- as.numeric(factor(birds_raw$Route))
species_idx <- as.numeric(factor(birds_raw$English_Common_Name))
year_idx <- as.numeric(factor(birds_raw$Year))

# Populate the array
stop_data <- as.matrix(birds_raw[, grep("^Stop", colnames(birds_raw))])

for (i in seq_len(nrow(birds_raw))) {
  Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
  }

  # Assign dimnames
  dimnames(Y) <- list(
    Route = sort(unique(birds_raw$Route)),
      Stop = paste0("Stop", 1:T),
        Species = sort(unique(birds_raw$English_Common_Name)),
          Year = sort(unique(birds_raw$Year))
          )

# Selecting only 5 bird species  for analysis:
Y<-Y[,,1:5,]

## Not run: model<-MNM_fit(Y=Y, AR=TRUE, Hurdle=FALSE, iterations=10000, burnin=2000)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MNM_AR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MNM_Hurdle")
### * MNM_Hurdle

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MNM_Hurdle
### Title: Fit a Multi-Species N-Mixture Model with Hurdle Component using
###   Nimble
### Aliases: MNM_Hurdle

### ** Examples

# Example 1:
Y <- array(rpois(100, lambda = 5), dim = c(10, 5, 2))
Xp <- array(runif(100), dim = c(10, 2, 5))
Xn <- array(runif(100), dim = c(10, 2, 3))

## Not run: model <- MNM_Hurdle(Y = Y, Xp = Xp, Xn = Xn)

# Accessing results
## Not run: print(model@summary)

#' data(birds_raw)

# Example 2: North American Breeding Bird Data
# Data must first be reformatted to an array of dimension (R,T,S,K)
R <- 24
T <- 10
S <- 20
K <- 6
# Ensure data is ordered consistently
birds_raw <- birds_raw[order(birds_raw$Route, birds_raw$Year, birds_raw$English_Common_Name), ]

# Create a 4D array with proper dimension
Y <- array(NA, dim = c(R, T, S, K))

# Map route, species, and year to indices
route_idx <- as.numeric(factor(birds_raw$Route))
species_idx <- as.numeric(factor(birds_raw$English_Common_Name))
year_idx <- as.numeric(factor(birds_raw$Year))

# Populate the array
stop_data <- as.matrix(birds_raw[, grep("^Stop", colnames(birds_raw))])

for (i in seq_len(nrow(birds_raw))) {
  Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
  }

  # Assign dimnames
  dimnames(Y) <- list(
    Route = sort(unique(birds_raw$Route)),
      Stop = paste0("Stop", 1:T),
        Species = sort(unique(birds_raw$English_Common_Name)),
          Year = sort(unique(birds_raw$Year))
          )

# Selecting only 5 bird species and 1 year for analysis:
Y<-Y[,,1:5,1]

## Not run: model<-MNM_fit(Y=Y, AR=FALSE, Hurdle=TRUE, iterations=5000, burnin=1000)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MNM_Hurdle", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MNM_Hurdle_AR")
### * MNM_Hurdle_AR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MNM_Hurdle_AR
### Title: Fit a multi-species N-mixture (MNM) Hurdle-AR model using Nimble
### Aliases: MNM_Hurdle_AR

### ** Examples

# Example 1: Simulate data and fit the model
# Simulating example data
set.seed(42)
R <- 5  # Number of sites
T <- 10  # Number of replicates
S <- 3  # Number of species
K <- 2  # Number of time periods
P1 <- 2  # Number of detection covariates
P2 <- 3  # Number of abundance covariates

x<-simulateData(model="HurdleAR", R=R, T=T, ,S=S, K=K)
Xp <- array(runif(R * S * K * P1), dim = c(R, S, K, P1))
Xn <- array(runif(R * S * K * P2), dim = c(R, S, K, P2))
# Fit the MNM_Hurdle_AR model
## Not run: result <- MNM_Hurdle_AR(Y = x[["Y"]], Xp = Xp, Xn = Xn)
# Access results
## Not run: print(result@summary)

#' data(birds_raw)

# Example 2: North American Breeding Bird Data
# Data must first be reformatted to an array of dimension (R,T,S,K)
R <- 24
T <- 10
S <- 20
K <- 6
# Ensure data is ordered consistently
birds_raw <- birds_raw[order(birds_raw$Route, birds_raw$Year, birds_raw$English_Common_Name), ]

# Create a 4D array with proper dimension
Y <- array(NA, dim = c(R, T, S, K))

# Map route, species, and year to indices
route_idx <- as.numeric(factor(birds_raw$Route))
species_idx <- as.numeric(factor(birds_raw$English_Common_Name))
year_idx <- as.numeric(factor(birds_raw$Year))

# Populate the array
stop_data <- as.matrix(birds_raw[, grep("^Stop", colnames(birds_raw))])

for (i in seq_len(nrow(birds_raw))) {
  Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
  }

  # Assign dimnames
  dimnames(Y) <- list(
    Route = sort(unique(birds_raw$Route)),
      Stop = paste0("Stop", 1:T),
        Species = sort(unique(birds_raw$English_Common_Name)),
          Year = sort(unique(birds_raw$Year))
          )

# Selecting only 5 bird species for analysis:
Y<-Y[,,1:5,]

## Not run: model<-MNM_fit(Y=Y, AR=TRUE, Hurdle=TRUE, iterations=5000, burnin=1000)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MNM_Hurdle_AR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MNM_control")
### * MNM_control

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MNM_control
### Title: Generate Nimble Code for Multi-Species N-Mixture (MNM) Models
### Aliases: MNM_control

### ** Examples

# Example
# In order to implement scenarios involving nonlinear covariate
# effects or complex interaction terms,
# the user is invited to extract and modify the MNM Nimble code:
model<-MNM_control(model="MNM")
cat(model)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MNM_control", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MNM_fit")
### * MNM_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MNM_fit
### Title: Fit a Multi-Species N-Mixture (MNM) Model in Nimble
### Aliases: MNM_fit

### ** Examples

# Example 1: Fit a standard MNM model
Y <- array(data = rpois(60, lambda = 5), dim = c(3, 5, 4))  # Simulated counts
Xp <- array(data = rnorm(60), dim = c(3, 4, 2))  # Detection covariates
Xn <- array(data = rnorm(60), dim = c(3, 4, 2))  # Abundance covariates

## Not run: model <- MNM_fit(Y = Y, AR = FALSE, Hurdle = FALSE, Xp = Xp, Xn = Xn)

#' # Example 2: Fit an MNM model with AR-1 component
Y <- array(data = rpois(180, lambda = 5), dim = c(3, 5, 4, 3))  # Simulated counts
Xp <- array(data = rnorm(180), dim = c(3, 4, 3, 2))  # Detection covariates
Xn <- array(data = rnorm(180), dim = c(3, 4, 3, 2))  # Abundance covariates

## Not run: model <- MNM_fit(Y = Y, AR = TRUE, Hurdle = FALSE, Xp = Xp, Xn = Xn)

# Example 3: Fit an MNM model with user-specified prior distributions
Y <- array(data = rpois(60, lambda = 5), dim = c(3, 5, 4))  # Simulated counts
Xp <- array(data = rnorm(60), dim = c(3, 4, 2))  # Detection covariates
Xn <- array(data = rnorm(60), dim = c(3, 4, 2))  # Abundance covariates

## Not run: 
##D model <- MNM_fit(Y = Y, AR = FALSE, Hurdle = TRUE, Xp = Xp, Xn = Xn,
##D                           prior_detection_probability="dnorm(0.01,0.01)")
## End(Not run)
# Access traceplots and density plots
## Not run: tracePlot(y, "N[10, 1]")
## Not run: density(y, "N[10, 1]")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MNM_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("birds_raw")
### * birds_raw

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: birds_raw
### Title: Raw Birds Dataset - Subset of the North American Breeding Bird
###   Survey Dataset
### Aliases: birds_raw
### Keywords: datasets

### ** Examples

data(birds_raw)
head(birds_raw)

# Example: MNM
# Data must first be reformatted to an array of dimension (R,T,S,K)
R <- 24
T <- 10
S <- 20
K <- 6
# Ensure data is ordered consistently
birds <- birds_raw[order(birds_raw$Route, birds_raw$Year, birds_raw$English_Common_Name), ]

# Create a 4D array with proper dimension
Y <- array(NA, dim = c(R, T, S, K))

# Map route, species, and year to indices
route_idx <- as.numeric(factor(birds_raw$Route))
species_idx <- as.numeric(factor(birds_raw$English_Common_Name))
year_idx <- as.numeric(factor(birds_raw$Year))

# Populate the array
stop_data <- as.matrix(birds_raw[, grep("^Stop", colnames(birds))])

for (i in seq_len(nrow(birds))) {
  Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
  }

  # Assign dimnames
  dimnames(Y) <- list(
    Route = sort(unique(birds_raw$Route)),
      Stop = paste0("Stop", 1:T),
        Species = sort(unique(birds_raw$English_Common_Name)),
          Year = sort(unique(birds_raw$Year))
          )

# Selecting only 5 bird species and 1 year for analysis:
Y<-Y[,,1:5,1]

## Not run: model<-MNM_fit(Y=Y, AR=FALSE, Hurdle=FALSE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("birds_raw", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("density-MNM-method")
### * density-MNM-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: density,MNM-method
### Title: Density Plot Method for "MNM" Class
### Aliases: density,MNM-method

### ** Examples

# Calling the density function
## Not run: density(y, "N[10, 1]")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("density-MNM-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulateData")
### * simulateData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulateData
### Title: Simulate Data for Multi-Species N-Mixture Models
### Aliases: simulateData

### ** Examples

# Simulate data for a standard MNM model
data <- simulateData(model = "MNM", R = 10, S = 3, T = 5, prob = "all",
abundance = "small")

# Simulate data for a hurdle model
data <- simulateData(model = "Hurdle", R = 10, S = 3, T = 5, prob = "large",
abundance = "large", theta = 0.3)

# Simulate data for an autoregressive model
data <- simulateData(model = "AR", R = 10, S = 2, T = 5, K = 4, prob = "small",
abundance = "small")

# Simulate data for a hurdle autoregressive model
data <- simulateData(model = "HurdleAR", R = 10, S = 3, T = 5, K = 4, prob = "all",
abundance = "large", theta = 0.5)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulateData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tracePlot")
### * tracePlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tracePlot
### Title: Trace Plot Generic and Method for "MNM" Class
### Aliases: tracePlot tracePlot,MNM-method

### ** Examples

## Not run: 
##D # Assuming `y` is an object of class "MNM" with plots stored
##D tracePlot(y, "N[8, 1]")  # Generates a trace plot for parameter N[8, 1]
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tracePlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("validate_prior")
### * validate_prior

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: validate_prior
### Title: Validate User-Specified Prior for a Single Prior
### Aliases: validate_prior

### ** Examples

validate_prior("prior_mean", "dnorm(0, 0.001)")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("validate_prior", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
