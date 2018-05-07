#' Computation of pointwise and highest pointwise probabilities.
#'
#' Pointwise (PW) probabilities and highest pointwise (HPW) probabilities of
#' all differences of smooths at neighboring scales are computed.
#'
#' \code{HPWmap} is an internal function of \code{\link{mrbsizeRgrid}} and is usually
#' not used independently. The output can be analyzed with the plotting function 
#' \code{\link{plot.HPWmapGrid}}.
#'
#' @param smoothVec Differences of smooths at neighboring scales.
#' @param mm Number of rows of the original input image.
#' @param nn Number of columns of the original input image.
#' @param prob Credibility level for the posterior credibility analysis
#' @return List with two arrays:
#' \itemize{
#' \item \code{pw}: Pointwise probabilities (\code{VmapPW}) including the dimensions
#'     of the original input image, \code{mm} and \code{nn}.
#' \item \code{hpw}: Highest pointwise probabilities (\code{VmapHPW}) including
#'     the dimensions of the original input image, \code{mm} and \code{nn}.
#' }
#' @examples
#' # Artificial sample data: 10 observations (5-by-2 object), 10 samples
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(100), nrow = 10)
#' sampleData[4:6, ] <- sampleData[4:6, ] + 5
#'
#' # Calculation of the simultaneous credible intervals
#' HPWmap(smoothVec = sampleData, mm = 5, nn = 2, prob = 0.95)
#'
HPWmap <- function(smoothVec, mm, nn, prob = 0.95) {

  #-----------------------Initiating variables-----------#
  N <- dim(smoothVec)[1]
  ns <- dim(smoothVec)[2]
  reds <- rep.int(1L, N)
  blues <- rep.int(0L, N)
  color <- rep.int(0L, N)
  VmapHPW <- rep(0.5, N)
  VmapPW <- rep(0.5, N)
  VMaximum <- rep.int(0L, N)


  #----------------------- PW-maps ----------------#
  for (i in 1:N) {
    SmoothIJ <- smoothVec[i, ]

    redi <- sum(Re(SmoothIJ) < 0)
    bluesi <- sum(Re(SmoothIJ) > 0)

    reds[i] <- redi
    blues[i] <- bluesi

    if (redi > bluesi) {
      VMaximum[i] <- redi / ns
      color[i] <- 0
    } else if (redi < bluesi) {
      VMaximum[i] <- bluesi / ns
      color[i] <- 1
    }

    if (redi >=  prob * ns) {
      VmapPW[i]  <-  0
    } else if (bluesi >= prob * ns) {
      VmapPW[i]  <-  1
    }

  }


  #----------------------- HPW-maps ----------------#
  Start <- N
  Ord <- sort(-VMaximum)
  ind <- order(-VMaximum)

  #----pixels with VMaximum = 1 are handled separately to accelerate computation----#

  for (q in 1:N) {
    if (Ord[q] > -1) {
      Start <- q
      break
    } else {
      if (color[ind[q]] == 1) {
        VmapHPW[ind[q]] <- 1
      } else {
        VmapHPW[ind[q]] <- 0
      }
    }
  }


  #--------------------rest of the pixels-----------------------------------#
  Correct <- rep.int(1L, ns)

  for (r in Start:N) {
    if (Ord[r] > -prob) {
      break
    }

    index <- ind[r]
    indVec <- color[index]
    NewCorrect <- ((Re(smoothVec[index, ]) > 0) == indVec) # indVec is 0 (FALSE) or 1 (TRUE)

    Correct <-  NewCorrect & Correct
    evalProb <- sum(Correct) / ns

    if ((evalProb) < prob) { # Stop when joint probability is < prob
      break
    } else {                  #-Coloring the HPW-maps
      if (color[index] == 1) {
        VmapHPW[index] <- 1
      } else {
        VmapHPW[index] <- 0
      }
    }
  }
  
  hpout <- list(pw = array(VmapPW, c(mm, nn)),
                hpw = array(VmapHPW, c(mm, nn)))
  return(hpout)

}
