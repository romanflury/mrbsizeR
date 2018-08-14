#' Numerical optimization for finding appropriate smoothing levels.
#'
#' Numerical optimization of an objective function \eqn{G} is carried out to find 
#' appropriate signal-dependent smoothing levels (\eqn{\lambda}'s). This is easier 
#' than visual inspection via the signal-dependent tapering function in \code{\link{TaperingPlot}}.
#'
#' As signal-dependent tapering functions are quiet irregular, it is hard to 
#' find appropriate smoothing values only by visual inspection of the tapering
#' function plot. A more formal approach is the numerical optimization of an 
#' objective function. 
#' 
#' Optimization can be carried out with 2 or 3 smoothing parameters. As the 
#' smoothing parameters 0 and \eqn{\infty} are always added, this results
#' in a mrbsizeR analysis with 4 or 5 smoothing parameters. 
#'
#' Sometimes, not all features of the input object can be extracted using the 
#' smoothing levels proposed by \code{MinLambda}. It might then be necessary to
#' include additional smoothing levels. 
#' 
#' \code{\link{plot.minLambda}} creates a plot of the objective function \eqn{G} 
#' on a grid. The minimum is indicated with a white point. The minimum values of 
#' the \eqn{\lambda}'s can be extracted from the output of \code{MinLambda}, 
#' see examples. 
#'
#' @param Xmu Posterior mean of the input object as a vector.
#' @param mm Number of rows of the original input object.
#' @param nn Number of columns of the original input object.
#' @param nGrid Size of grid where objective function is evaluated (nGrid-by-nGrid).
#' This argument is ignorded if a sequence \code{lambda} is specified.
#' @param nLambda Number of lambdas to minimize over. Possible arguments: 2 (default) or 3. 
#' @param lambda \eqn{\lambda}-sequence which is used for optimization. If nothing is provided, \cr
#'     \code{lambda <- 10^seq(-3, 10, len = nGrid)} is used for data on a grid and \cr
#'     \code{lambda <- 10^seq(-6, 1, len = nGrid)} is used for spherical data.
#' @param sphere \code{TRUE} or \code{FALSE}: Is the input object defined on a sphere?
#' @return A list with 3 objects: 
#' @return \code{G} Value of objective function \eqn{G}.
#' @return \code{lambda} Evaluated smoothing parameters \eqn{\lambda}.
#' @return \code{minind} Index of minimal \eqn{\lambda}'s. \code{lambda}[\code{minind}] 
#'     gives the minimal values.
#' @examples
#' # Artificial sample data
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(100), nrow = 10)
#' sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5
#' 
#' # Minimization of two lambdas on a 20-by-20-grid
#' minlamOut <- MinLambda(Xmu = c(sampleData), mm = 10, nn = 10, 
#'                        nGrid = 20, nLambda = 2)
#' 
#' # Minimal lambda values
#' minlamOut$lambda[minlamOut$minind]
#' 
MinLambda <- function(Xmu, mm, nn, nGrid, nLambda = 2, lambda, sphere = FALSE){

  if (mm * nn != length(Xmu)) {
    stop("Dimensions (mm * nn) do not fit the number of locations (length(Xmu))")
  }
  
  #----------- Calculation of Eigenvalues ----------------------------------#  

  if (sphere == FALSE) {
    eigMu <- eigenLaplace(mm, nn)
    D <- eigMu^2
    fac <- c(dctMatrix(mm) %*% matrix(Xmu, nrow = mm) %*% t(dctMatrix(nn)))
    if (methods::hasArg(lambda) == FALSE) {
      lambda <- 10^seq(-3, 10, len = nGrid)
    } else {
      lambda <- unique(lambda)
      nGrid <- length(lambda)
    }
  } else {
    eigOut <- eigenQsphere((180 / mm), (180 - 180 / mm), mm, nn)
    eigMu <- eigOut$eigval
    eigVec <- eigOut$eigvec
    D <- t(sort(eigMu))
    indx <- sort(eigMu, index.return = TRUE)$ix
    D[1] <- 0
    eigVec <- eigVec[ , indx]
    fac <- c(t(Xmu) %*% eigVec)
    if (methods::hasArg(lambda) == FALSE) {
      lambda <- 10^seq(-6, 1, len = nGrid)
    } else {
      lambda <- unique(lambda)
      nGrid <- length(lambda)
    }
  }

  #---------- Initiating variables ---------------------------------------#

  N <- mm * nn

  minimum <- 10^11
  if (identical(nLambda, 2) | identical(nLambda, 3)) {
    G <- array(NA, c(rep(length(lambda), nLambda)))
  } else{
    stop("Optimization can only be done over 2 or 3 lambdas")
  }
  lambda1 <- fac
  lambda1[1] <- 0
  LambdaMat <- initLambdaMat(nGrid, N, lambda, fac, D)

  #--------- Computing the value of G on the grid for nLambda = 2 --------------#
  if (identical(nLambda, 2)) {
    min2LambdaListoutput <- min2Lambda(nGrid, LambdaMat, lambda1, minimum)
    G <- apply(min2LambdaListoutput$G, c(1,2), function(x) { if(identical(x, 0)) { x <- NA } else { x <- x}})

    minil <- lambda[min2LambdaListoutput$mini]
    minjl <- lambda[min2LambdaListoutput$minj]
    #--------- Computing the value of G on the grid for nLambda = 3 --------------#  
  } else {

    tmp <- sapply(0:(nGrid-2), function(i) {
      lambda2 <- LambdaMat[ ,i+1]
      min3LambaListoutput <- min3Lambda(i, nGrid, LambdaMat, lambda1, lambda2, minimum)
      if (min3LambaListoutput$minimum < minimum) {
        minimum <- min3LambaListoutput$minimum
        mini <- min3LambaListoutput$mini
        minj <- min3LambaListoutput$minj
        mink <- min3LambaListoutput$mink
      }

      return(list(Gtmp = min3LambaListoutput$G, mini = min3LambaListoutput$mini, minj = min3LambaListoutput$minj, mink = min3LambaListoutput$mink))
    })

    for(i in 1:(nGrid-2)) {
      G[i, , ] <- apply(tmp[,i]$Gtmp, c(1,2), function(x) { if(identical(x, 0)) { x <- NA } else { x <- x}})
    }

    minil <- lambda[tmp$mini]
    minjl <- lambda[tmp$minj]
    minkl <- lambda[tmp$mink]
  }

  minlamout <- list(G = G, lambda = lambda, 
                    minind = which(G == min(G, na.rm = TRUE), arr.ind = TRUE))
  class(minlamout) <- append(class(minlamout), "minLambda")
  return(minlamout)
}

