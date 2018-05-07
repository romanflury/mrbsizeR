#' Plot of simultaneous credible intervals.
#'
#' Maps with simultaneous credible intervals for all differences of smooths
#' at neighboring scales \eqn{z_{i}} are plotted.
#'
#' The default colors of the maps have the following meaning:
#' \itemize{
#' \item \strong{Blue}: Credibly positive pixels.
#' \item \strong{Red}: Credibly negative pixels.
#' \item \strong{Grey}: Pixels that are not credibly different from zero.
#' }
#' \code{x} corresponds to the \code{ciout}-part of the output 
#' of \code{\link{mrbsizeRgrid}}.
#'
#' @param x List containing the simultaneous credible intervals for all
#'     differences of smooths.
#' @param color Vector of length 3 containing the colors to be used in the 
#'     credibility maps. The first color represents the credibly negative pixels, 
#'     the second color the pixels that are not credibly different from zero
#'     and the third color the credibly positive pixels. 
#' @param turnOut Logical. Should the output images be turned 90 degrees 
#'     counter-clockwise?
#' @param title Vector containing one string per plot. The required 
#'     number of titles is equal to \code{length(mrbOut$ciout)}. If no \code{title} 
#'     is passed, defaults are used. 
#' @param aspRatio Adjust the aspect ratio of the plots. The default \code{aspRatio = 1}
#'     produces square plots.
#' @param ... Further graphical parameters can be passed. 
#' @return Plots of simultaneous credible intervals for all differences of
#'     smooths are created.
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
#' # Credibility analysis using simultaneous credible intervals
#' plot(x = mrbOut$ciout, turnOut = TRUE)
#' 
plot.CImapGrid <- function(x, color = c("firebrick1", "gainsboro", "dodgerblue3"), turnOut = TRUE, 
                           title, aspRatio = 1, ...){ 

  if (length(color) != 3) {
    stop("'color' must be a vector containing exactly three colors!")
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
      
    if(turnOut){
      graphics::image(turnmat(x[[i]]), col = color, breaks = c(-0.1, 0.4, 0.6, 1.1),
                      main = mainTxt, xaxt = "n", yaxt = "n", cex.main = 1, 
                      asp = aspRatio, bty = "n", ...)
    } else {
      graphics::image(x[[i]], col = color, breaks = c(-0.1, 0.4, 0.6, 1.1),
                      main = mainTxt, xaxt = "n", yaxt = "n", cex.main = 1, 
                      asp = aspRatio, bty = "n", ...)
    }
  }
}
