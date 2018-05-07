#' Computation of simultaneous credible intervals.
#'
#' Simultaneous credible intervals for all differences of smooths at neighboring
#' scales \eqn{z_{i}} are computed.
#'
#' \code{CImap} is an internal function of \code{\link{mrbsizeRgrid}} and is usually
#' not used independently. The output can be analyzed with the plotting function 
#' \code{\link{plot.CImapGrid}}.
#'
#' @param smoothVec Differences of smooths at neighboring scales.
#' @param mm Number of rows of the original input object.
#' @param nn Number of columns of the original input object.
#' @param prob Credibility level for the posterior credibility analysis. 
#'     By default \code{prob = 0.95}.
#' @return An array with simultaneous credible intervals \code{VmapCI} and
#'     the dimensions of the original input object, \code{mm} and \code{nn}.
#' @examples
#' # Artificial sample data: 10 observations (5-by-2 object), 10 samples
#' set.seed(987)
#' sampleData  <- matrix(stats::rnorm(100), nrow = 10)
#' sampleData [4:6, ] <- sampleData [4:6, ] + 5
#'
#' # Calculation of the simultaneous credible intervals
#' CImap(smoothVec = sampleData , mm = 5, nn = 2, prob = 0.95)
#'
CImap <- function(smoothVec, mm, nn, prob = 0.95){

  #-------------Initiating variables---------#
  N <- dim(smoothVec)[1]
  ns <- dim(smoothVec)[2]
  smmean <- rep(0, N)
  smsd <- rep(0, N)
  Z <- rep(0, min(ns, N))
  if (floor(N / 5000) == 0) {
    Nloop <- 1
  } else {
    Nloop <- floor(N / 5000)
  }

  VmapCI <- rep(0.5, N)


  #--------------Defining mean and std--------#
  for (j in 1:Nloop) {
    lb <- (j - 1) * 5000 + 1
    if (N < 5000) {
      ub <- N
    } else {
      ub <- j * 5000
    }
    smooth <- Re(smoothVec[lb:ub, ])
    smmean[lb:ub] <- rowSums(smooth) / ns
    smsd[lb:ub] <- apply(smooth, 1, stats::sd)
  }


  if (N > 5000 && (N / 5000) > floor(N / 5000)) {
    lb <- Nloop * 5000 + 1
    ub <- N
    smooth <- Re(smoothVec[lb:ub, ])
    smmean[lb:ub] <- rowSums(smooth) / ns
    smsd[lb:ub] <- apply(smooth, 1, stats::sd)
  }


  #-------------- Defining Delta-----------------#
  for (i in 1:min(ns,N)) {
    resid <- abs((Re(smoothVec[ ,i]) - smmean) / smsd)
    maks <- max(resid)
    I <- which.max(resid)
    Z[i] <- maks
  }


  if (sum(is.na(Z)) == length(Z)) {
    delta <- 0
  } else {
    delta <- stats::quantile(Z, prob)
  }


  #---------------Coloring CI-maps---------------#
  for (i in 1:N) {
    if (Re(smmean[i] - delta * smsd[i]) > 0){
      VmapCI[i] <- 1
    } else if (Re(smmean[i] + delta * smsd[i]) < 0) {
      VmapCI[i] <- 0
    }
  }
  ciout <- array(VmapCI, c(mm, nn))
  return(ciout)
}

