
# MultiNMix: A Package for Multispecies N-Mixture Models


`MultinMix` is an R package designed for fitting Multispecies N-Mixture (MNM) Models (Mimnagh, Niamh, et al. (2022)), a powerful tool for estimating abundance and occurrence of multiple species in a hierarchical Bayesian framework.

### Features
- **Bayesian Modeling**: Fit hierarchical Bayesian MNM models using Nimble.
- **Customisable Priors**: Define prior distributions easily for each parameter.
- **Comprehensive Outputs**: Includes posterior summaries, convergence diagnostics, and model fit statistics (log-likelihood, AIC, BIC).
- **User-Friendly API**: Simple interface to specify data, initial values, and model parameters.
- **Visualisation**: Built-in methods for producing density plots and traceplots, for model diagnostics.
  

### Installation
To install the development version of `MultiNMix`, use the following commands in R:
```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")

devtools::install_github("niamhmimnagh/MultiNMix")
```

### Getting Started
Here is a quick example to get you started with `MultiNMix`:

```R
library(MultinMix)

# Example data
x <- simulateData(model = "MNM")
 R<-x$R
 T<-x$T
 S<-x$S
 K<-x$K

 Xp <- array(rnorm(R * S * 2), dim = c(R,  S, 2)) # creating 2 detection probability covariates
 Xn <- array(rnorm(R * S *3), dim = c(R, S,  3)) # creating 3 abundance covariates
 
# Fit 
fit <- MNM_fit(
  Y = species_counts,
  Xp = Xp,
  Xn = Xn,
  Hurdle=FALSE,
  AR = FALSE,
  iterations = 5000,  # Number of iterations
  burnin = 1000,  # Burn-in period
  thin = 10,  # Thinning interval
  prior_detection_probability="dnorm(0,0.01)" # user-defined normal prior distribution
)

# Summarize results
fit@summary

# Plot diagnostic results by specifying the model and the parameter
tracePlot(fit, param="N[8,1]")
density(fit, param="N[8,1]")

# A list of all available diagnostic plots can be found:
View(y@plot)
```

### Functions
#### Main Function
- `MNM_fit()`: Fits a Multispecies N-Mixture Model using specified data and parameters.

#### Utility Functions
- `tracePlot()`: Generates traceplots of monitored parameters.
- `density()`: Generates density plots of monitored parameters.
- `logLik()`: Extracts the log-likelihood of the model.
- `AIC()`, `BIC()`: Computes AIC and BIC values for model comparison.
- `check_convergence()`: Assesses model convergence using Gelman-Rubin diagnostics.

### Documentation
Detailed documentation and vignettes are available in the package. After installation, access them using:

```R
??MultiNMix
```

### Datasets
There are two datasets available in the package `birds` and the zero-inflated `birds_ZI`. Both are a subset of the North American Breeding Bird Survey dataset (https://www.pwrc.usgs.gov/BBS/). `birds` is a dataframe with 2,880 observations and 13 columns (R=24, T=10, S=20, K=6) while `birds_ZI` is a dataframe with 600 observations and 13 columns (R=15, T=10, S=10, K=4). 

In this vignette, we will show the `birds` dataset, the processing steps required and a worked example of it.

### The birds Dataset

 ```R
data(birds)
head(birds)
```

|  Route |  Year |  English_Common_Name |Stop 1|Stop 2|...|Stop 10|
| - | - | - | - | - | - |-|
|  001 |  2016 |  Mourning Dove |0|1|...|0|
|  007 |  2016 |  Mourning Dove |6|4|...|5|
|   009|  2016 |  Mourning Dove |0|0|...|0|


The `birds` dataset is currently a data frame of dimension (600, 10). It needs to be reformatted into an array of dimension (R=15, T=10, S=10, K=4) before it can be used with the `MultiNMix` functions.

```R
 # Data must first be reformatted to an array of dimension (R,T,S,K)
   R <- 15
   T <- 10
   S <- 10
   K <- 4

 # Ensure data is ordered consistently
   birds <- birds[order(birds$Route, birds$Year, birds$English_Common_Name), ]
  
 # Create a 4D array with proper dimension
   Y <- array(NA, dim = c(R, T, S, K))
  
 # Map route, species, and year to indices
   route_idx <- as.numeric(factor(birds$Route))
   species_idx <- as.numeric(factor(birds$English_Common_Name))
   year_idx <- as.numeric(factor(birds$Year))
  
 # Populate the array
   stop_data <- as.matrix(birds[, grep("^Stop", colnames(birds))])
  
   for (i in seq_len(nrow(birds))) {
     Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
     }
  
 # Assign dimnames
     dimnames(Y) <- list(
       Route = sort(unique(birds$Route)),
         Stop = paste0("Stop", 1:T),
           Species = sort(unique(birds$English_Common_Name)),
             Year = sort(unique(birds$Year)))
```
The function `MNM_fit` in the `MultiNMix` package allows for easy implementation of a multi-species N-mixture model using data of this format. 

```R
model<-MNM_fit((Y=Y, AR=FALSE, Hurdle=FALSE))
```

We can then access elements of the model as follows:

``` R
model@summary # outputs the mean estimate, standard deviation, standard error, 95% credible interval, effective sample size and gelman rubin statistic for each monitored variable

model@estimates$N # outputs the estimated mean  N

logLik(model) # estimates the log likelihood of the model

AIC(model)/BIC(model) # outputs the AIC or BIC values

tracePlot(model, param="N[1,1]") # outputs the traceplot of the N[1,1] parameter

density(model, param="N[1,1]") #outputs the density plot for the N[1,1] parameter

```

### Contributions
Contributions are welcome! If you encounter any issues or have suggestions for improvement, please submit a report or a pull request.

### References
Mimnagh, Niamh, et al. "Bayesian multi-species N-mixture models for unmarked animal communities." Environmental and Ecological Statistics 29.4 (2022): 755-778.



### Acknowledgements
`MultiNMix` was developed as part of research into multispecies abundance modeling. Special thanks to the creators of Nimble (r-nimble.org) for their invaluable tools in Bayesian modeling.

