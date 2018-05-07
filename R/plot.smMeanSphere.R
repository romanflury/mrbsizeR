#' Plotting of scale-dependent features on a sphere.
#'
#' Scale-dependent features are plotted using differences of smooths at
#' neighboring scales. The features are summarized by their posterior mean. 
#' Continental lines are added to the plots.
#' 
#' \code{x} corresponds to the \code{smmean}-part of the
#' output of \code{\link{mrbsizeRsphere}}.
#'
#' @param x List containing the posterior mean of all differences of smooths.
#' @param lon Vector containing the longitudes of the data points.
#' @param lat Vector containing the latitudes of the data points.
#' @param color.pallet The color pallet to be used for plotting scale-dependent
#'    features. 
#' @param turnOut Logical. Should the output images be turned 90 degrees 
#'     counter-clockwise?
#' @param title Vector containing one string per plot. The required 
#'     number of titles is equal to \code{length(mrbOut$smMean)}. If no \code{title} 
#'     is passed, defaults are used. 
#' @param ... Further graphical parameters can be passed.
#' @return Plots of the differences of smooths are created.
#' @examples
#' # Artificial spherical sample data
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(2000), nrow = 200)
#' sampleData[50:65, ] <- sampleData[50:65, ] + 5
#' lon <- seq(-180, 180, length.out = 20)
#' lat <- seq(-90, 90, length.out = 10)
#' 
#' # mrbsizeRsphere analysis
#' mrbOut <- mrbsizeRsphere(posteriorFile = sampleData, mm = 20, nn = 10, 
#'                          lambdaSmoother = c(0.1, 1), prob = 0.95)
#'                            
#' # Posterior mean of the differences of smooths
#' plot(x = mrbOut$smMean, lon = lon, lat = lat,
#'      color = fields::tim.colors(), turnOut = FALSE) 
#'
plot.smMeanSphere <- function(x, lon, lat, color.pallet = fields::tim.colors(), 
                              turnOut = TRUE, title, ...) {
  
  if (length(lon) * length(lat) != dim(x[[1]])[1] * dim(x[[1]])[2]) {
    stop("Longitude and latitude are non-comformable to the dimension of x")
  }
  
  if (methods::hasArg(title) == TRUE && length(title) != length(x)) {
    stop("Number of titles must be equal to the number of plots")
  }
  
  graphics::par(mfrow = c(ceiling(length(x) / 2), 2))
  graphics::par(mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
  
  if (methods::hasArg(title) == FALSE) {
    namesVec <- names(x)
    namesVec <- sapply(strsplit(namesVec, "_"), "[", 1)
    namesVec <- gsub("NA", "infinity", namesVec)
  }
  
  for(i in 1:length(x)){
    if (methods::hasArg(title) == FALSE) {
      if (i < length(x)) {
        mainTxt <- as.expression(
          substitute(
            paste(E, "(", z[y], "|y) for ", lambda, "-range = [", x, " - ", w, "]"), 
            list(x = as.name(namesVec[i]), y = i, w = as.name(namesVec[i + 1]))
          )
        )
      } else {
        mainTxt <- as.expression(
          substitute(
            paste(E, "(", z[y], "|y) for ", lambda, "-range = [", x, "]"), 
            list(x = as.name(namesVec[i]), y = i)
          )
        )
      }
    } else {
      mainTxt <- title[i]
    }
    
    if(turnOut) {
      tmp <- turnmat(Re(x[[i]]))
      lon.tmp <- lat
      lat.tmp <- lon
    } else {
      tmp <- Re(x[[i]])
      lon.tmp <- lon
      lat.tmp <- lat
    }

    fields::image.plot(lon.tmp, lat.tmp, tmp, col = color.pallet, main = mainTxt, 
                       xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.main = 1, 
                       bty = "n",...)
    maps::map("world", add = TRUE)
  }
}


