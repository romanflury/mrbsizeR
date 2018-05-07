#' Plot of objective function for finding appropriate smoothing parameters. 
#'
#' The objective function \eqn{G} is plotted on a grid. The minimum is 
#' indicated with a white point. 
#' 
#' When minimizing over 2 \eqn{\lambda}'s, one plot is generated: (\eqn{\lambda_2} vs \eqn{\lambda_3}).
#' With 3 \eqn{\lambda}'s, 3 plots are generated: \eqn{\lambda_2} vs. \eqn{\lambda_3}, 
#' \eqn{\lambda_2} vs. \eqn{\lambda_4} and \eqn{\lambda_3} vs. \eqn{\lambda_4}.
#'
#' @param x Output of function \code{\link{MinLambda}}.
#' @param ... Further graphical parameters can be passed.
#' @return Plot of \eqn{G} on a grid.
#' @examples
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(100), nrow = 10)
#' sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5
#' 
#' # Minimization of two lambdas on a 20-by-20-grid
#' minLamOut <- MinLambda(Xmu = c(sampleData), mm = 10, nn = 10, 
#'                         nGrid = 20, nLambda = 3)
#' 
#' # Plot of the objective function
#' plot(x = minLamOut)
#'
plot.minLambda <- function(x, ...) {
  
  n <- length(x$lambda)
  axisLabels <- floor(log10(range(x$lambda)))
  pow <- 10^seq(axisLabels[1], axisLabels[2])
  pow <- format(pow, scientific = TRUE)
  
  
  if (length(x$minind) == 2) {
    graphics::par(mfcol = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
    graphics::image(1:(n - 1), 2:n, x$G[-n, -1], xlab = expression(lambda[2]), 
          ylab = expression(lambda[3]), xaxt = "n", yaxt = "n", col = grDevices::rainbow(64),
          cex.lab = 1.2, ...)
  
    graphics::axis(1, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::axis(2, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::box()
    graphics::points(x$minind[1], x$minind[2], col = "white", pch = 19, cex = 1.2)
    invisible()
    
  } else {
    graphics::par(mfcol = c(1, 3), mar = c(5.1, 4.1, 4.1, 2.1))
    graphics::image(1:n, 1:n, x$G[ , , x$minind[3]], xlab = expression(lambda[2]),
                    ylab = expression(lambda[3]), xaxt = "n", yaxt = "n", col = grDevices::rainbow(64),
                    zlim = range(x$G, na.rm = TRUE), cex.lab = 1.2, ...)
    
    graphics::axis(1, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::axis(2, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::box()
    graphics::points(x$minind[1], x$minind[2], col = "white", pch = 19, cex = 1.2)
    
    graphics::image(1:n, 1:n, x$G[ , x$minind[2], ], xlab = expression(lambda[2]),
                    ylab = expression(lambda[4]), xaxt = "n", yaxt = "n", col = grDevices::rainbow(64),
                    zlim = range(x$G, na.rm = TRUE), cex.lab = 1.2, ...)
    graphics::axis(1, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::axis(2, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::box()
    graphics::points(x$minind[1], x$minind[3], col = "white", pch = 19, cex = 1.2)
    
    graphics::image(1:n, 1:n, x$G[x$minind[1], , ], xlab = expression(lambda[3]),
                    ylab = expression(lambda[4]), xaxt = "n", yaxt = "n", col = grDevices::rainbow(64),
                    zlim = range(x$G, na.rm = TRUE), cex.lab = 1.2, ...)
    graphics::axis(1, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::axis(2, at = seq(1, n, length.out = length(pow)), pow, cex.axis = 0.9, ...)
    graphics::box()
    graphics::points(x$minind[2], x$minind[3], col = "white", pch = 19, cex = 1.2)
    invisible()
  }
}



