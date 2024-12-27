
# MultiNMix: A Package for Multispecies N-Mixture Models


`MultinMix` is an R package designed for fitting Multispecies N-Mixture (MNM) Models, a powerful tool for estimating abundance and occurrence of multiple species in a hierarchical Bayesian framework.

### Features
- **Bayesian Modeling**: Fit hierarchical Bayesian MNM models using Nimble.
- **Customisable Priors**: Define prior distributions easily for each parameter.
- **Comprehensive Outputs**: Includes posterior summaries, convergence diagnostics, and model fit statistics (log-likelihood, AIC, BIC).
- **User-Friendly API**: Simple interface to specify data, initial values, and model parameters.
- **Visualisation**: Built-in methods for producing density plots and traceplots, for model diagnostics.
  

### Installation
To install the development version of \code{MultinMix}, use the following commands in R:
```
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")

devtools::install_github("yourusername/multinmix")
```

### Getting Started
Here is a quick example to get you started with `MultiNMix`:

```
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
  prior_detection_probability='dnorm(0,0.01)' # user-defined normal prior distribution
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

```
??MultiNMix
```

### Datasets
There are two datasets available in the package `birds_raw` and the processed format `birds`. Both are a subset of the North American Breeding Bird Survey dataset (https://www.pwrc.usgs.gov/BBS/), containing data collected at 24 routes in Michigan, USA.Each route has 10 stops, and the dataset includes counts for 20 bird species. `birds_raw` is a dataframe with 2,880 rows and 13 columns, and so requires processing into an array of dimension (R=24, T=10, S=20, K=6) before it can be used. `birds` has already been processed into the correct format, and so is ready to use out of the box.

In this vignette, we will show the `birds_raw` dataset, the processing steps required and a worked example of it.

### The birds_raw Dataset

 ```
data(birds_raw)
```
The birds_raw dataset has a large number of zero-counts. Of the 28,800 observations, 24,291 are zero-counts. 


### Contributions
Contributions are welcome! If you encounter any issues or have suggestions for improvement, please open an issue or submit a pull request on the \href{https://github.com/yourusername/multinmix}{GitHub repository}.



### Acknowledgements
`MultiNMix` was developed as part of research into multispecies abundance modeling. Special thanks to the creators of Nimble (r-nimble.org) for their invaluable tools in Bayesian modeling.

