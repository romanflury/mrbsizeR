#' Multiresolution analysis of random signals for spherical data.
#' 
#' \code{mrbsizeRSphere} is the interface of the scale space multiresolution method
#' for spherical data. Here, the differences of smooths as well as the posterior 
#' credibility analysis are computed. The output can be analyzed with the plotting
#' functions \code{\link{plot.smMeanSphere}}, \code{\link{plot.CImapSphere}} 
#' and \code{\link{plot.HPWmapSphere}}.
#'
#' In contrast to \code{mrbsizeRgrid}, \code{mrbsizeRsphere} does not conduct 
#' Bayesian signal reconstruction via sampling from a posterior distribution. 
#' Samples of the posterior distribution have to be provided instead. 
#'
#' For further information and examples, see \code{\link{mrbsizeRgrid}} and 
#' the vignette.
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
#'    of all differences of smooths created.
#' @return \code{ciout} Simultaneous credible intervals (CI) of all differences of
#'     smooths created.
#' @return \code{smoothSamples} Samples of differences of smooths at neighboring scales, 
#'     as column vectors.
#' @examples
#' # Artificial spherical sample data
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(2000), nrow = 200)
#' sampleData[50:65, ] <- sampleData[50:65, ] + 5
#' 
#' # mrbsizeRsphere analysis
#' mrbOut <- mrbsizeRsphere(posteriorFile = sampleData, mm = 10, nn = 20,  
#'                           lambdaSmoother = c(1, 1000), prob = 0.95)
mrbsizeRsphere <- function (posteriorFile, mm, nn, lambdaSmoother, prob = 0.95,
                            smoothOut = FALSE) {

  if (mm * nn != nrow(posteriorFile)) {
    stop("Dimensions (mm * nn) do not match the number of locations (nrow(posteriorFile))")
  }
  
  #------------------ Initializing variables -------------------------------#
  
  N <- mm * nn
  ns <- ncol(posteriorFile)
  
  sampleVek <- posteriorFile
  sampleMat <- array(posteriorFile, c(mm, nn, ns))
  sm <- array(rep(0, (N * ns)), c(mm, nn, ns))
  

  #------------------Drawing posterior mean  -------------------------------#
    
  mu <- apply(sampleMat, c(1, 2), mean)

    
  #---------- Computing variables needed for smoothing ---------------------#
    
  phi <- seq(180 / mm, 180 - (180 / mm), len = mm)
  deltaPhi <- pi / 180 * (phi[2] - phi[1])
  Ufft <- array(NA, c(mm, nn, ns))
  for(i in 1:dim(sampleMat)[3]){
    Ufft[ , , i] <- t(apply(sampleMat[ , , i], MARGIN  = 1, FUN = stats::fft))
  }
  
  for(i in 1:dim(Ufft)[3]) {
    Ufft[ , , i] <- ifftshift(Ufft[ , , i], dimension = 2)
  }
  
  sinPhi2 <- -1 / sin(phi * pi / 180)^2
  
  CNeven <- tridiag(-2 * rep(1, mm), rep(1, mm - 1), rep(1, mm - 1)) / deltaPhi^2
  c2 <- diag(sinPhi2)
  CNeven[1, 1] <- CNeven[1, 1] + 1 / deltaPhi^2
  CNeven[mm, min(mm, nn)] <- CNeven[mm, min(mm, nn)] + 1 / deltaPhi^2
  
  
  #------------  Smoothing ------------------------------------------------#
        
  smoothVec <- sampleVek
  smoothNew <- smoothVec
      
  if (lambdaSmoother[1] != 0) {
    lambdaSmoother <- c(0, lambdaSmoother)
  }
        
  nl <- length(lambdaSmoother)
  
  smMean <- list() 
  hpout <- list() 
  ciout <- list()
  smoothSamples <- list()
  
 
  for (b in 2:(nl + 2)) {
    if (b > nl) {
      smoothVec <- matrix(rep(apply(sampleVek, 2, mean), N), ncol = ns, byrow = TRUE) 
    } else { 
      beta <- lambdaSmoother[b]
      for (i in 1:nn) {
        CNevenMM <- CNeven + c2 * (i - (1 + nn / 2))^2
        invSeven <- diag(1, mm, mm) + beta * t(CNevenMM) %*% CNevenMM
        Uffi <- matrix(Ufft[ , i, ], nrow = mm)
        sm[ , i, ] <- solve(invSeven, Uffi)
      }
      
      smooth <- array(NA, c(mm, nn, ns))
      for(i in 1:dim(sm)[3]) {
        sm[ , , i] <- fftshift(sm[ , , i], dimension = 2)
        smooth[ , , i] <- t(apply(sm[ , , i], MARGIN = 1, FUN = stats::fft, inverse = TRUE) / max(nrow(sm), ncol(sm))) 
      }
      
      for (i in 1:ns){
        smoothVec[ , i] <- c(smooth[ , , i])
      }
    }
        
    smoothOld <- smoothNew
    smoothNew <- smoothVec
    smoothVec <- smoothOld - smoothNew
        
    if (b == (nl+2)) {
      smoothVec <- smoothNew
    }
        
    smoothMean <- apply(smoothVec, 1, mean)  

    #-------------------Saving the MRBSiZer-maps -----------------------#

    if (b > 1){
      smMean[[b-1]] <- matrix(smoothMean, nrow = mm)
      names(smMean)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      hpout[[b-1]] <- HPWmap(smoothVec, mm, nn, prob = 0.95)
      names(hpout)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      ciout[[b-1]] <- CImap(smoothVec, mm, nn, prob = 0.95)
      names(ciout)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      
      if (smoothOut == TRUE){
        smoothSamples[[b-1]] <- smoothVec
        names(smoothSamples)[[b-1]] <- paste(lambdaSmoother[b-1], "_", lambdaSmoother[b], sep = "")
      }
    }
  }
  class(smMean) <- append(class(smMean), "smMeanSphere")
  class(hpout) <- append(class(hpout), "HPWmapSphere") 
  class(ciout) <- append(class(ciout), "CImapSphere") 
  
  if (smoothOut == TRUE){
    return(list(smMean = smMean, hpout = hpout, ciout = ciout, smoothSamples = smoothSamples))
  } else {
    return(list(smMean = smMean, hpout = hpout, ciout = ciout)) 
  }
}
        
        
