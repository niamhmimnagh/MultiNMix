#' Birds Dataset -  Subset of the North American Breeding Bird Survey Dataset
#'
#' #' This dataset is a subset of the North American Breeding Bird Survey, containing data collected at 15 routes in Michigan, USA.
#' Each route has 10 stops, and the dataset includes counts for 10 bird species.
#'
#' @format A data frame with 600 rows and 13 columns:
#' \describe{
#'   \item{Route}{An identifier for the 15 survey sites.}
#'   \item{Year}{The year of the survey (numeric).}
#'   \item{English_Common_Name}{The English common name of the bird species surveyed (character).}
#'   \item{Stop1}{Count for replicate 1 at the site (numeric).}
#'   \item{Stop2}{Count for replicate 2 at the site (numeric).}
#'   \item{Stop3}{Count for replicate 3 at the site (numeric).}
#'   \item{Stop4}{Count for replicate 4 at the site (numeric).}
#'   \item{Stop5}{Count for replicate 5 at the site (numeric).}
#'   \item{Stop6}{Count for replicate 6 at the site (numeric).}
#'   \item{Stop7}{Count for replicate 7 at the site (numeric).}
#'   \item{Stop8}{Count for replicate 8 at the site (numeric).}
#'   \item{Stop9}{Count for replicate 9 at the site (numeric).}
#'   \item{Stop10}{Count for replicate 10 at the site (numeric).}
#' }
#'
#' @source  North American Breeding Bird Survey (\url{https://www.pwrc.usgs.gov/BBS/})
#' @examples
#' # MNM Model
#' #' # Example: Hurdle Model
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
#' \donttest{model<-MNM_fit(Y=Y, AR=FALSE, Hurdle=FALSE)}
#'
"birds"

#' Zero-Inflated Birds Dataset - Subset of the North American Breeding Bird Survey Dataset
#'
#' This dataset is a subset of the North American Breeding Bird Survey, containing data collected at 24 routes in Michigan, USA.
#' Each route has 10 stops, and the dataset includes counts for 20 bird species.
#'
#' @format A data frame with 2,880 rows and 13 columns:
#' \describe{
#'   \item{Route}{An identifier for the 24 survey sites.}
#'   \item{Year}{The year of the survey (numeric).}
#'   \item{English_Common_Name}{The English common name of the bird species surveyed (character).}
#'   \item{Stop1}{Count for replicate 1 at the site (numeric).}
#'   \item{Stop2}{Count for replicate 2 at the site (numeric).}
#'   \item{Stop3}{Count for replicate 3 at the site (numeric).}
#'   \item{Stop4}{Count for replicate 4 at the site (numeric).}
#'   \item{Stop5}{Count for replicate 5 at the site (numeric).}
#'   \item{Stop6}{Count for replicate 6 at the site (numeric).}
#'   \item{Stop7}{Count for replicate 7 at the site (numeric).}
#'   \item{Stop8}{Count for replicate 8 at the site (numeric).}
#'   \item{Stop9}{Count for replicate 9 at the site (numeric).}
#'   \item{Stop10}{Count for replicate 10 at the site (numeric).}
#' }
#'
#' @details
#' This dataset represents a subset of the North American Breeding Bird Survey. Data was collected in Michigan over six years,
#' with observations for 20 bird species recorded at 24 routes, each surveyed 10 times. The dataset is used
#' to study avian biodiversity and population trends.
#'
#' @source  North American Breeding Bird Survey (\url{https://www.pwrc.usgs.gov/BBS/})
#' @examples
#' data(birds_ZI)
#' head(birds_ZI)
#'
#' # Example: Hurdle Model
#' # Data must first be reformatted to an array of dimension (R,T,S,K)
#' R <- 24
#' T <- 10
#' S <- 20
#' K <- 6
#' # Ensure data is ordered consistently
#' birds_ZI <- birds_ZI[order(birds_ZI$Route, birds_ZI$Year, birds_ZI$English_Common_Name), ]
#'
#' # Create a 4D array with proper dimension
#' Y <- array(NA, dim = c(R, T, S, K))
#'
#' # Map route, species, and year to indices
#' route_idx <- as.numeric(factor(birds_ZI$Route))
#' species_idx <- as.numeric(factor(birds_ZI$English_Common_Name))
#' year_idx <- as.numeric(factor(birds_ZI$Year))
#'
#' # Populate the array
#' stop_data <- as.matrix(birds_ZI[, grep("^Stop", colnames(birds))])
#'
#' for (i in seq_len(nrow(birds))) {
#'   Y[route_idx[i], , species_idx[i], year_idx[i]] <- stop_data[i, ]
#'   }
#'
#'   # Assign dimnames
#'   dimnames(Y) <- list(
#'     Route = sort(unique(birds_ZI$Route)),
#'       Stop = paste0("Stop", 1:T),
#'         Species = sort(unique(birds_ZI$English_Common_Name)),
#'           Year = sort(unique(birds_ZI$Year))
#'           )
#'
#' # Selecting only 5 bird species  for analysis:
#' Y<-Y[,,1:5,]
#'
#' \donttest{model<-MNM_fit(Y=Y, AR=TRUE, Hurdle=TRUE)}

"birds_ZI"

