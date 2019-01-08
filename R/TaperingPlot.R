#' Plot of tapering functions.
#'
#' Tapering functions corresponding to the smoothing levels in \code{lambdaSmoother}
#' are drawn. This plot helps to assess if the chosen smoothing levels are
#' appropriate.
#'
#' The tapering functions of the smoothing levels chosen should be generally
#' approximately disjoint. This will produce features which are somewhat orthogonal.
#' With orthogonal features, it is likely that each difference of smooths
#' corresponds to a different pattern in the input image.
#'
#' Sometimes, not all patterns of the input image can be extracted using smoothing
#' levels whose tapering functions are disjoint. It might then be necessary to
#' include additional smoothing levels, and the disjointedness might not be satisfied
#' anymore. The selection of appropriate smoothing levels with this method
#' therefore requires some user interaction. Still, choosing disjoint tapering
#' functions for finding appropriate smoothing levels is a good starting point.
#'
#' Better results could be obtained if the structure of the posterior mean of
#' the input is also taken into account. If the posterior mean is available,
#' it can be added with the argument \code{Xmu} and moving averages of the absolute
#' values of the signal-dependent tapering functions are drawn. \code{\link{MinLambda}} 
#' offers a more formal approach of optimizing the disjointedness of the tapering 
#' functions and can help finding appropriate smoothing levels when the input 
#' signal is taken into account.
#'
#' @param lambdaSmoother Vector consisting of the smoothing levels to be used.
#' @param mm Number of rows of the original input object.
#' @param nn Number of columns of the original input object.
#' @param Xmu If availabe, posterior mean of the input object.
#' @param returnseq instead of plotting return the tapering sequences.
#' @param ... Further graphical parameters can be passed.
#' @return Plots of the tapering functions for all differences of smooths
#'     at neighboring scales are created.
#' @examples
#' # Signal-independent tapering function plot for a 30-by-10 object with 
#' # the smoothing parameter sequence [0, 1, 10, 1000, inf]: 
#' 
#' TaperingPlot(lambdaSmoother = c(1, 10, 1000), mm = 30, nn = 10)
#' 
#' 
#' # Signal-dependent tapering function plot for a 30-by-10 object with 
#' # the smoothing parameter sequence [0, 1, 10, 1000, inf]: 
#' 
#' set.seed(987)
#' xmuExample <- c(stats::rnorm(300))
#' TaperingPlot(lambdaSmoother = c(1, 10, 1000), mm = 30, nn = 10, 
#'              Xmu = xmuExample)
#'
#'
TaperingPlot <- function(lambdaSmoother, mm, nn, Xmu, returnseq = FALSE, ...){

  #---------- Initiating variables --------------------------#

  if (min(lambdaSmoother) != 0) {
    lambdaSmoother <- c(0, sort(lambdaSmoother, decreasing = FALSE))
  }

  N <- mm * nn
  MU <- t(eigenLaplace(mm, nn))^2
  D <- sort(MU)
  indx <- order(MU)

  if (missing(Xmu)){
    fac <- rep(1, N)
  } else {
    if (mm * nn != length(Xmu)) {
      stop("Dimensions (mm * nn) do not fit the number of locations (length(Xmu))")
    }
    fac <- matrix(dctMatrix(mm) %*% matrix(Xmu, nrow = mm) %*% t(dctMatrix(nn)), nrow = N)
    fac <- fac[indx]
  }

  lenL <- length(lambdaSmoother) - 1
  val1 <- rep(0, N)
  val2 <- rep(0, N)
  x <- seq(1, N, 1)
  diff <- matrix(NA, nrow = lenL + 1, ncol = N)


  #--------- Defining the tapering functions and drawing the plot ---------#
  for (l in 1:lenL){
    lfirst <- lambdaSmoother[l]
    lnext <- lambdaSmoother[l+1]

    for (i in 2:N){
      val2[i] <- 1 / (1 + lfirst * D[i]) * fac[i]
      val1[i] <- 1 / (1 + lnext * D[i]) * fac[i]
    }

    diff[l, ] <- val2 - val1
  }
  diff[lenL + 1, ] <- val1

  # If signal-dependent, use running mean
  if (missing(Xmu) == FALSE){
    for (i in 1:dim(diff)[1]){
      for (j in 1:dim(diff)[2]){
        diff[i, j] <- abs(diff[i, j]) * j
      }
      diff[i, ] <- stats::filter(diff[i, ], rep(1, 5 + 20 * (lenL - i + 1)), sides = 2)
      diff[which(is.na(diff))] <- 0
    }
  }

  if (returnseq) {
    return(diff)
  }

  # Plot all tapering functions except of the last one
  for (i in 1:lenL){
    if (i == 1){
      graphics::par(mfrow = c(1, 1), mar = c(4.1, 4.1, 2.1, 9.1), xpd = TRUE)
      graphics::plot(x, diff[1, ], type = "l", lty = 1, col = 1, lwd = 2, 
                     ylim = c(min(diff, na.rm = TRUE), max(diff, na.rm = TRUE)),
                     xlim = c(1, N), xlab = "Eigenvalue index", 
                     ylab = "Tapering function value", bty = "L", 
                     cex.axis = 0.8, cex.lab = 0.8, log = "x", ...)
    } else {
      graphics::lines(x, diff[i, ], type = "l", lty = 1, col = i, lwd = 2)
    }
  }


  # Add missing Tapering Functions ((x - inf) and inf)
  graphics::lines(x, diff[lenL + 1, ], type = "l", lty = 1, col = lenL + 1, lwd = 2)
  graphics::segments(x0 = 1, y0 = 0, x1 = 1, y1 = max(diff), col = lenL + 2, lwd = 2)
  graphics::segments(x0 = 1, y0 = 0, x1 = N, y1 = 0, col = lenL + 2, lwd = 2)
  graphics::box()


  # Add Legend
  legendlabel <- vector(mode = "expression", length = lenL + 2)
  for (i in 1:lenL){
    legendlabel[i] <- paste(lambdaSmoother[i], " - ", lambdaSmoother[i + 1])
  }
  legendlabel[lenL + 1] <- substitute(expression(paste(x, " - ", infinity)), list(x = lambdaSmoother[lenL + 1]))[2]
  legendlabel[lenL + 2] <- expression(infinity)

  graphics::legend("right", legend = c(legendlabel), lty = 1, col = seq(1, lenL + 2, 1), lwd = 2, bty = "n",
                   inset = c(-0.4, 0), cex = 0.8, title = "Smoothing range")
}


