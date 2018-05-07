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
    if (methods::hasArg(lambda) == FALSE) lambda <- 10^seq(-3, 10, len = nGrid)
  } else {
    eigOut <- eigenQsphere((180 / mm), (180 - 180 / mm), mm, nn)
    eigMu <- eigOut$eigval
    eigVec <- eigOut$eigvec
    D <- t(sort(eigMu))
    indx <- sort(eigMu, index.return = TRUE)$ix
    D[1] <- 0
    eigVec <- eigVec[ , indx]
    fac <- c(t(Xmu) %*% eigVec)
    if (methods::hasArg(lambda) == FALSE) lambda <- 10^seq(-6, 1, len = nGrid)
  }

  #---------- Initiating variables ---------------------------------------#

  N <- mm * nn

  minimum <- 10^11
  if (identical(nLambda, 2) | identical(nLambda, 3)) {
    G <- array(NA, c(rep(nGrid, nLambda)))
  } else{
    stop("Optimization can only be done over 2 or 3 lambdas")
  }
  lambda1 <- fac
  lambda1[1] <- 0
  LambdaMat <- matrix(0, nrow = N, ncol = nGrid)

  for (i in 1:nGrid) {
    l1 <- lambda[i]
    lambda2 <- fac / (1 + l1 * D)
    lambda2[1] <- 0
    LambdaMat[ ,i] <- lambda2
  } 


  #--------- Computing the value of G on the grid for nLambda = 2 --------------#
  if (identical(nLambda, 2)) {
   
    for (i in 1:(nGrid-1)) {
      lambda2 <- LambdaMat[ ,i]
      
      for (j in (i+1):nGrid) {
        lambda3 <- LambdaMat[ ,j]
        diff12 <- lambda1 - lambda2
        diff12 <- diff12 / sqrt(sum(diff12 * diff12))
        
        diff23 <- lambda2 - lambda3
        diff23 <- diff23 / sqrt(sum(diff23 * diff23))
        
        diff34 <- lambda3
        diff34 <- diff34 / sqrt(sum(diff34 * diff34))
                           
        val <- abs(sum(diff12 * diff23)) + abs(sum(diff23 * diff34)) + 
          abs(sum(diff12 * diff34))
        
        G[i, j] <- val
        
        if (val < minimum) {
          minimum <- val
          mini <- i
          minj <- j
        }
      }
    } 
    minil <- lambda[mini]
    minjl <- lambda[minj]
    #--------- Computing the value of G on the grid for nLambda = 3 --------------#  
  } else {
    
    for (i in 1:(nGrid - 2)) {
      lambda2 <- LambdaMat[ ,i]
      
      for (j in (i + 1):(nGrid - 1)) {
        lambda3 <- LambdaMat[ ,j]
        
        for (k in (j + 1):nGrid) {
          lambda4 <- LambdaMat[ ,k]
          
          diff12 <- lambda1 - lambda2
          diff12 <- diff12 / sqrt(sum(diff12 * diff12))
          
          diff23 <- lambda2 - lambda3
          diff23 <- diff23 / sqrt(sum(diff23 * diff23))
          
          diff34 <- lambda3 - lambda4
          diff34 <- diff34 / sqrt(sum(diff34 * diff34))
          
          diff45 <- lambda4
          diff45 <- diff45 / sqrt(sum(diff45 * diff45))
          
          val <- abs(sum(diff12 * diff23)) + abs(sum(diff12 * diff34)) + 
            abs(sum(diff12 * diff45)) + abs(sum(diff23 * diff34)) + 
            abs(sum(diff23 * diff45)) + abs(sum(diff34 * diff45))
          
          G[i, j, k] <- val
          
          if (val < minimum) {
            minimum <- val
            mini <- i
            minj <- j
            mink <- k
          }
        }
      }
    }
  } 
  minlamout <- list(G = G, lambda = lambda, 
                    minind = which(G == min(G, na.rm = TRUE), arr.ind = TRUE))
  class(minlamout) <- append(class(minlamout), "minLambda")
  return(minlamout)
}

