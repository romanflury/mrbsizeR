#' Multiresolution analysis of random signals.
#'
#' \code{mrbsizeRgrid} is the interface of the scale space multiresolution method
#' for data on a regular grid. Here, the differences of smooths as well as the posterior 
#' credibility analysis are computed. The output can be analyzed with the plotting
#' functions \code{\link{plot.smMeanGrid}}, \code{\link{plot.CImapGrid}} and 
#' \code{\link{plot.HPWmapGrid}}.
#'
#' \code{mrbsizeRgrid} conducts two steps of the scale space multiresolution analysis:
#' \enumerate{
#'     \item Extraction of scale-dependent features from the reconstructed signal.
#'         This is done by smoothing at different smoothing levels and taking the
#'         difference of smooths at neighboring scales.
#'     \item Posterior credibility analysis of the differences of smooths created.
#'         Three different methods are applied: Pointwise probabilities (see \code{\link{HPWmap}}),
#'         highest pointwise probabilities (see \code{\link{HPWmap}}) and simultaneous 
#'         credible intervals (see \code{\link{CImap}}).
#' }
#' The signal can be reconstructed using the build-in multivariate t-distribution
#' sampling \code{\link{rmvtDCT}}. It is also possible to provide samples 
#' generated with other methods, see the parameter \code{posteriorFile} and the
#' examples.
#' 
#' For further information and examples, see the vignette.
#'
#' @param posteriorFile Matrix with posterior samples as column vectors. 
#' @param mm Number of rows of the original object.
#' @param nn Number of columns of the original object.
#' @param lambdaSmoother Vector consisting of the smoothing levels to be used.
#' @param prob Credibility level for the posterior credibility analysis.
#' @param smoothOut Should the differences of smooths at neighboring scales be 
#'     returned as output (FALSE by default)?
#' @return A list containing the following sublists:
#' @return \code{smMean} Posterior mean of all differences of smooths created.
#' @return \code{hpout} Pointwise (PW) and highest pointwise (HPW) probabilities 
#'     of all differences of smooths created.
#' @return \code{ciout} Simultaneous credible intervals (CI) of all differences of
#'     smooths created.
#' @return \code{smoothSamples} Samples of differences of smooths at neighboring scales, 
#'     as column vectors.
#' @examples
#' # Artificial sample data
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(100), nrow = 10)
#' sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5
#' 
#' # Generate samples from multivariate t-distribution
#' tSamp <- rmvtDCT(object = sampleData, lambda = 0.2, sigma = 6, nu0 = 15,
#'                   ns = 1000)  
#'  
#' # mrbsizeRgrid analysis
#' mrbOut <- mrbsizeRgrid(posteriorFile = tSamp$sample, mm = 10, nn = 10, 
#'                        lambdaSmoother = c(1, 1000), prob = 0.95)
#'
mrbsizeRgrid <- function(posteriorFile, mm, nn, lambdaSmoother, prob = 0.95, smoothOut = FALSE) {

  if (mm * nn != nrow(posteriorFile)) {
    stop("Dimensions (mm * nn) do not match the number of locations (nrow(posteriorFile))")
  }
  
  #------------------ Sampling ---------------------------------------------#

  sampleVek <- posteriorFile
  muVek <- rowMeans(posteriorFile)

  #------------------Initiating variables-----------------------------------#

  N <- mm * nn
  ns <- ncol(posteriorFile)
  
  dct_mtx_mm <- dctMatrix(mm)
  dct_mtx_nn <- dctMatrix(nn)
  t_dct_mtx_mm <- t(dctMatrix(mm))
  t_dct_mtx_nn <- t(dctMatrix(nn))

  if (lambdaSmoother[1] != 0) {
    lambdaSmoother <- c(0, lambdaSmoother)
  }

  nl <- length(lambdaSmoother)

  smMean <- list()
  hpout <- list()
  ciout <- list()
  smoothSamples <- list()

  #------------------Drawing posterior mean ----------------------------------#
  mu <- matrix(muVek, nrow = mm)

  smoothNew <- sampleVek
  smoothVec <- sampleVek

  DCTsample <- matrix(NA, nrow = dim(sampleVek)[1], ncol = dim(sampleVek)[2])
  eigMu <- matrix(eigenLaplace(mm, nn)^2, nrow = mm)

  for (s in 1:ns) {
    DCTsample[ ,s] <- dct_mtx_mm %*% matrix(sampleVek[ , s], nrow = mm) %*% t_dct_mtx_nn
  }


  # loop over smoothing values
  for (b in 1:(nl+2)){


  #--------------- Smoothing the sampled objects #-------------------#
    if (b > nl) {
      smoothVec <- matrix(rep(apply(sampleVek, 2, mean), N), ncol = ns, byrow = TRUE)  # Infinty
    } else {
      beta <- lambdaSmoother[b]
      if (beta > 0) {
        for (s in 1:ns){
          c <- 1 / (1 + beta * eigMu) * DCTsample[ , s]
          sm <- t_dct_mtx_mm %*% c %*% dct_mtx_nn
          smoothVec[ , s] <- sm
        }
      } else {
        smoothVec <- sampleVek
      }
    }


    #----------------- Saving the smooths #-------------------------------#
    smoothOld <-  smoothNew
    smoothNew <-  smoothVec
    smoothVec <-  smoothOld - smoothNew

    if (b == (nl + 2)) {
      smoothVec <-  smoothNew
    }
    smoothMean <- apply(smoothVec, 1, mean)

    #-------------------Saving the mrbsizeR-maps -----------------------#
    # First difference of smooths (empty matrix) is not needed
    if (b > 1){
      smMean[[b-1]] <- matrix(smoothMean, nrow = mm)
      names(smMean)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      hpout[[b-1]] <- HPWmap(smoothVec, mm, nn, prob = prob)
      names(hpout)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      ciout[[b-1]] <- CImap(smoothVec, mm, nn, prob = prob)
      names(ciout)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      
      if(smoothOut == TRUE) {
        smoothSamples[[b-1]] <- smoothVec
        names(smoothSamples)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      }
    }

  }
  class(smMean) <- append(class(smMean), "smMeanGrid")
  class(hpout) <- append(class(hpout), "HPWmapGrid") 
  class(ciout) <- append(class(ciout), "CImapGrid") 

  if(smoothOut == TRUE){
    return(list(smMean = smMean, hpout = hpout, ciout = ciout, smoothSamples = smoothSamples))
  } else {
    return(list(smMean = smMean, hpout = hpout, ciout = ciout))
  }
}




