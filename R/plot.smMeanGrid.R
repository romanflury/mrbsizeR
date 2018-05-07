#' Plotting of scale-dependent features.
#'
#' Scale-dependent features are plotted using differences of smooths at
#' neighboring scales. The features are summarized by their posterior mean.
#'
#' \code{x} corresponds to the \code{smmean}-part of the
#' output of \code{\link{mrbsizeRgrid}}.
#'
#' @param x List containing the posterior mean of all differences of smooths.
#' @param color.pallet The color pallet to be used for plotting scale-dependent
#'    features.
#' @param turnOut Logical. Should the output images be turned 90 degrees 
#'     counter-clockwise?
#' @param title Vector containing one string per plot. The required 
#'     number of titles is equal to \code{length(mrbOut$smMean)}. If no \code{title} 
#'     is passed, defaults are used. 
#' @param aspRatio Adjust the aspect ratio of the plots. The default \code{aspRatio = 1}
#'     produces square plots.
#' @param ... Further graphical parameters can be passed.
#' @return Plots of the differences of smooths are created.
#' @examples
#' # Artificial sample data
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(100), nrow = 10)
#' sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5
#' 
#' # Generate samples from multivariate t-distribution
#' tSamp <- rmvtDCT(object = sampleData, lambda = 0.2, sigma = 6, nu0 = 15,
#'                    ns = 1000)  
#'  
#' # mrbsizeRgrid analysis
#' mrbOut <- mrbsizeRgrid(posteriorFile = tSamp$sample, mm = 10, nn = 10, 
#'                        lambdaSmoother = c(1, 1000), prob = 0.95)
#' 
#' # Posterior mean of the differences of smooths
#' plot(x = mrbOut$smMean, turnOut = TRUE) 
#'
plot.smMeanGrid <- function(x, color.pallet = fields::tim.colors(), turnOut = TRUE, 
                            title, aspRatio = 1, ...) {
  
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

  for(i in 1:length(x)) {
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
      graphics::image(turnmat(x[[i]]), col = color.pallet,
                      main = mainTxt, xaxt = "n", yaxt = "n", cex.main = 1, 
                      asp = aspRatio, bty = "n", ...)
    } else {
      graphics::image(x[[i]], col = color.pallet,
                      main = mainTxt, xaxt = "n", yaxt = "n", cex.main = 1, 
                      asp = aspRatio, bty = "n", ...)
    }
  }
}


