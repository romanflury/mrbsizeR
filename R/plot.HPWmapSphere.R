#' Plotting of pointwise and highest pointwise probabilities on a sphere.
#'
#' Maps with pointwise (PW) probabilities and/or highest pointwise (HPW) 
#' probabilities of all differences of smooths at neighboring scales are plotted. 
#' Continental lines are added. 
#'
#' The default colors of the maps have the following meaning:
#' \itemize{
#' \item \strong{Blue}: Credibly positive pixels.
#' \item \strong{Red}: Credibly negative pixels.
#' \item \strong{Grey}: Pixels that are not credibly different from zero.
#' }
#' \code{x} corresponds to the \code{hpout}-part of the
#' output of \code{\link{mrbsizeRsphere}}.
#'
#' @param x List containing the pointwise (PW) and highest pointwise (HPW)
#'     probabilities of all differences of smooths.
#' @param lon Vector containing the longitudes of the data points.
#' @param lat Vector containing the latitudes of the data points.
#' @param plotWhich Which probabilities shall be plotted? \code{HPW}, \code{PW}
#'     or \code{Both}?
#' @param color Vector of length 3 containing the colors to be used in the 
#'     credibility maps. The first color represents the credibly negative pixels, 
#'     the second color the pixels that are not credibly different from zero
#'     and the third color the credibly positive pixels. 
#' @param turnOut Logical. Should the output images be turned 90 degrees 
#'     counter-clockwise?
#' @param title Vector containing one string per plot. The required 
#'     number of titles is equal to \code{length(mrbOut$hpout)}. If no \code{title} 
#'     is passed, defaults are used. 
#' @param ... Further graphical parameters can be passed.
#' @return Plots of pointwise and/or highest pointwise probabilities for all
#'     differences of smooths are created.
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
#' # Credibility analysis using pointwise (PW) maps
#' plot(x = mrbOut$hpout, lon = lon, lat = lat, plotWhich = "PW")
#' 
#' # Credibility analysis using highest pointwise probability (HPW) maps
#' plot(x = mrbOut$hpout, lon = lon, lat = lat, plotWhich = "HPW")
#'
plot.HPWmapSphere <- function(x, lon, lat, plotWhich = "Both", 
                             color = c("firebrick1", "gainsboro", "dodgerblue3"), 
                             turnOut = FALSE, title, ...) {
  
  if (length(color) != 3) {
    stop("'color' must be a vector containing exactly three colors!")
  }
  
  if (length(lon) * length(lat) != dim(x[[1]]$pw)[1] * dim(x[[1]]$pw)[2]) {
    stop("Longitude and latitude are non-comformable to the dimension of x")
  }

  if (methods::hasArg(title) == TRUE && length(title) != length(x)) {
    stop("Number of titles must be equal to the number of plots")
  }
  
  if (methods::hasArg(title) == FALSE) {
    namesVec <- names(x)
    namesVec <- sapply(strsplit(namesVec, "_"), "[", 1)
    namesVec <- gsub("NA", "infinity", namesVec)
  }
  
  if(plotWhich == "HPW") {
    graphics::par(mfrow = c(ceiling(length(x) / 2), 2))
    graphics::par(mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
    for(i in 1:length(x)){
      if (methods::hasArg(title) == FALSE) {
        if (i < length(x)) {
          mainTxt <- as.expression(
            substitute(
              paste("HPW-Map for ", z[y], " with ", lambda, "-range = [", x, " - ", w, "]"), 
              list(x = as.name(namesVec[i]), y = i, w = as.name(namesVec[i + 1]))
            )
          )
        } else {
          mainTxt <- as.expression(
            substitute(
              paste("HPW-Map for ", z[y], " with ", lambda, "-range = [", x, "]"), 
              list(x = as.name(namesVec[i]), y = i)
            )
          )
        }
      } else {
        mainTxt <- title[i]
      }
      
      if(turnOut) {
        tmp <- turnmat(x[[i]]$hpw)
        lon.tmp <- lat
        lat.tmp <- lon
      } else {
        tmp <- x[[i]]$hpw
        lon.tmp <- lon
        lat.tmp <- lat
      }
      graphics::image(lon.tmp, lat.tmp, tmp, col = color, main = mainTxt, 
                      xaxt = "n", yaxt = "n", cex.main = 1, bty = "n", ...)
      maps::map("world", add = TRUE)
    }
  } 
  
  if(plotWhich == "PW") {
    graphics::par(mfrow = c(ceiling(length(x) / 2), 2))
    graphics::par(mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
    for(i in 1:length(x)){
      if (methods::hasArg(title) == FALSE) {
        if (i < length(x)) {
          mainTxt <- as.expression(
            substitute(
              paste("PW-Map for ", z[y], " with ", lambda, "-range = [", x, " - ", w, "]"), 
              list(x = as.name(namesVec[i]), y = i, w = as.name(namesVec[i + 1]))
            )
          )
        } else {
          mainTxt <- as.expression(
            substitute(
              paste("PW-Map for ", z[y], " with ", lambda, "-range = [", x, "]"), 
              list(x = as.name(namesVec[i]), y = i)
            )
          )
        }
      } else {
        mainTxt <- title[i]
      }
      
      if(turnOut) {
        tmp <- turnmat(x[[i]]$pw)
        lon.tmp <- lat
        lat.tmp <- lon
      } else {
        tmp <- x[[i]]$pw
        lon.tmp <- lon
        lat.tmp <- lat
      }
      graphics::image(lon.tmp, lat.tmp, tmp, col = color, main = mainTxt, 
                      xaxt = "n", yaxt = "n", cex.main = 1, bty = "n", ...)
      maps::map("world", add = TRUE)
    }
  }
  
  if(plotWhich == "Both") {
    graphics::par(mfrow = c(ceiling(length(x) / 2), 2))
    graphics::par(mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
    for(i in 1:length(x)){
      if (methods::hasArg(title) == FALSE) {
        if (i < length(x)) {
          mainTxt <- as.expression(
            substitute(
              paste("HPW-Map for ", z[y], " with ", lambda, "-range = [", x, " - ", w, "]"), 
              list(x = as.name(namesVec[i]), y = i, w = as.name(namesVec[i + 1]))
            )
          )
        } else {
          mainTxt <- as.expression(
            substitute(
              paste("HPW-Map for ", z[y], " with ", lambda, "-range = [", x, "]"), 
              list(x = as.name(namesVec[i]), y = i)
            )
          )
        }
      } else {
        mainTxt <- title[i]
      }
      
      if(turnOut) {
        tmp <- turnmat(x[[i]]$hpw)
        lon.tmp <- lat
        lat.tmp <- lon
      } else {
        tmp <- x[[i]]$hpw
        lon.tmp <- lon
        lat.tmp <- lat
      }
      graphics::image(lon.tmp, lat.tmp, tmp, col = color, main = mainTxt, 
                      xaxt = "n", yaxt = "n", cex.main = 1, bty = "n", ...)
      maps::map("world", add = TRUE)
    }
    readline(prompt = "Press [enter] for the PW maps!")
    
    graphics::par(mfrow = c(ceiling(length(x) / 2), 2))
    graphics::par(mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
    for(i in 1:length(x)){
      if (methods::hasArg(title) == FALSE) {
        if (i < length(x)) {
          mainTxt <- as.expression(
            substitute(
              paste("PW-Map for ", z[y], " with ", lambda, "-range = [", x, " - ", w, "]"), 
              list(x = as.name(namesVec[i]), y = i, w = as.name(namesVec[i + 1]))
            )
          )
        } else {
          mainTxt <- as.expression(
            substitute(
              paste("PW-Map for ", z[y], " with ", lambda, "-range = [", x, "]"), 
              list(x = as.name(namesVec[i]), y = i)
            )
          )
        }
      } else {
        mainTxt <- title[i]
      }
      
      if(turnOut) {
        tmp <- turnmat(x[[i]]$pw)
        lon.tmp <- lat
        lat.tmp <- lon
      } else {
        tmp <- x[[i]]$pw
        lon.tmp <- lon
        lat.tmp <- lat
      }
      graphics::image(lon.tmp, lat.tmp, tmp, col = color, main = mainTxt, 
                      xaxt = "n", yaxt = "n", cex.main = 1, bty = "n", ...) 
      maps::map("world", add = TRUE)
    }
  }
}
