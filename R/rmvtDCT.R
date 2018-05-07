#' Sampling from marginal posterior multivariate t-distribution.
#'
#' Samples from a marginal posterior multivariate t-distribution with
#' normal-inverse-chi-squared-prior are generated.
#'
#' An eigenvalue decomposition is used for sampling. To speed up computations,
#' a 2D discrete cosine transform (DCT) has been implemented, see \code{\link{dctMatrix}}.
#' The output is a list containing
#' \enumerate{
#'     \item Samples of the marginal posterior of the input as column vectors.
#'     \item The mean of the marginal posterior of the input as a vector.
#'}
#'
#' @param object Observed object, as \code{matrix}.
#' @param lambda Scaling parameter (\eqn{\lambda}) of the normal-inverse-chi-squared-prior.
#' @param sigma Square root of the \eqn{\sigma_{0}^{2}} parameter of the 
#'     normal-inverse-chi-squared-prior.
#' @param nu0 Degrees of freedom (\eqn{\nu_{0}}) of the normal-inverse-chi-square-prior.
#' @param ns Number of samples that should be generated.
#' @return A list containing the following elements:
#' @return \code{sample} Samples of the marginal posterior of the input.
#' @return \code{mu} Mean of the marginal posterior of the input.
#' @examples
#' # Artificial sample data
#' set.seed(987)
#' sampleData <- matrix(stats::rnorm(100), nrow = 10)
#' sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5
#' 
#' # Sampling from a multivariate t-distribution
#' t_dist_samp <- rmvtDCT(object = sampleData, lambda = 1, sigma = 10,
#'                        nu0 = 50, ns = 1000)
rmvtDCT <- function(object, lambda, sigma, nu0, ns){

  #---------------------------Initiating variables--------------------------#
  Y <- object
  mm <- dim(Y)[1]
  nn <- dim(Y)[2]
  N <- mm * nn

  dct_mat_mm <- dctMatrix(mm)
  dct_mat_nn <- dctMatrix(nn)
  t_dct_mat_mm <- t(dctMatrix(mm))
  DCTy <- dct_mat_mm %*% Y %*% t(dct_mat_nn)
  DCTy <- as.vector(DCTy)

  y <- as.vector(Y)
  sigma20 <- sigma^2
  deg_freed <- nu0 + N - 1
  eigMu <- eigenLaplace(mm, nn)
  eigMu <- eigMu^2
  sample <- array(0, dim = c(N, ns))


  #--------------------------Drawing samples--------------------------------#
  c <- 1 / t((1 + lambda %*% eigMu)) * DCTy
  d <- t_dct_mat_mm %*% matrix(c, nrow = mm) %*% dct_mat_nn

  mu <- as.vector(d)
  prop <- sqrt((t(y) %*% y - t(y) %*% mu + nu0 %*% sigma20) %/% deg_freed)
  ev <- sqrt((1 / (1 + lambda %*% eigMu)))

  for (i in 1:ns) {
    randi <- ev * stats::rnorm(N)
    b <- t_dct_mat_mm %*% matrix(randi, nrow = mm) %*% dct_mat_nn
    u <- stats::rchisq(1, deg_freed)
    factor <- sqrt(deg_freed / u)
    x <- mu + prop %*% factor %*% as.vector(b)
    sample[, i] <- x
  }

  return(list(sample = sample, mu = mu))
}

