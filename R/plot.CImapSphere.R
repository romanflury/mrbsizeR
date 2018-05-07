#' Plotting of simultaneous credible intervals on a sphere.
#'
#' Maps with simultaneous credible intervals for all differences of smooths
#' at neighboring scales \eqn{z_{i}} are plotted. Continental lines are added. 
#'
#' The default colors of the maps have the following meaning:
#' \itemize{
#' \item \strong{Blue}: Credibly positive pixels.
#' \item \strong{Red}: Credibly negative pixels.
#' \item \strong{Grey}: Pixels that are not credibly different from zero.
#' }
#' \code{x} corresponds to the \code{ciout}-part of the
#'     output of \code{\link{mrbsizeRsphere}}.
#'
#' @param x List containing the simultaneous credible intervals of all
#'     differences of smooths.
#' @param lon Vector containing the longitudes of the data points.
#' @param lat Vector containing the latitudes of the data points.
#' @param color Vector of length 3 containing the colors to be used in the 
#'     credibility maps. The first color represents the credibly negative pixels, 
#'     the second color the pixels that are not credibly different from zero
#'     and the third color the credibly positive pixels. 
#' @param turnOut Logical. Should the output images be turned 90 degrees 
#'     counter-clockwise?
#' @param title Vector containing one string per plot. The required 
#'     number of titles is equal to \code{length(mrbOut$ciout)}. If no \code{title} 
#'     is passed, defaults are used. 
#' @param ... Further graphical parameters can be passed. 
#' @return Plots of simultaneous credible intervals for all differences of
#'     smooths are created.
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
#'      color = fields::tim.colors()) 
#' 
#' # Credibility analysis using simultaneous credible intervals
#' plot(x = mrbOut$ciout, lon = lon, lat = lat)
#'
plot.CImapSphere <- function(x, lon, lat, color = c("firebrick1", "gainsboro", "dodgerblue3"), 
                            turnOut = FALSE, title, ...) {
  
  if (length(color) != 3) {
    stop("'color' must be a vector containing exactly three colors!")
  }
  
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
            paste("CI-Map for ", z[y], " with ", lambda, "-range = [", x, " - ", w, "]"), 
            list(x = as.name(namesVec[i]), y = i, w = as.name(namesVec[i + 1]))
          )
        )
      } else {
        mainTxt <- as.expression(
          substitute(
            paste("CI-Map for ", z[y], " with ", lambda, "-range = [", x, "]"), 
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
    graphics::image(lon.tmp, lat.tmp, tmp, col = color, main = mainTxt, 
                    xaxt = "n", yaxt = "n", cex.main = 1, bty = "n", ...)
    maps::map("world", add = TRUE)
  }
}
