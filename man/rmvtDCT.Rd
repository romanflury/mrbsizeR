% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmvtDCT.R
\name{rmvtDCT}
\alias{rmvtDCT}
\title{Sampling from marginal posterior multivariate t-distribution.}
\usage{
rmvtDCT(object, lambda, sigma, nu0, ns)
}
\arguments{
\item{object}{Observed object, as \code{matrix}.}

\item{lambda}{Scaling parameter (\eqn{\lambda}) of the normal-inverse-chi-squared-prior.}

\item{sigma}{Square root of the \eqn{\sigma_{0}^{2}} parameter of the 
normal-inverse-chi-squared-prior.}

\item{nu0}{Degrees of freedom (\eqn{\nu_{0}}) of the normal-inverse-chi-square-prior.}

\item{ns}{Number of samples that should be generated.}
}
\value{
A list containing the following elements:

\code{sample} Samples of the marginal posterior of the input.

\code{mu} Mean of the marginal posterior of the input.
}
\description{
Samples from a marginal posterior multivariate t-distribution with
normal-inverse-chi-squared-prior are generated.
}
\details{
An eigenvalue decomposition is used for sampling. To speed up computations,
a 2D discrete cosine transform (DCT) has been implemented, see \code{\link{dctMatrix}}.
The output is a list containing
\enumerate{
    \item Samples of the marginal posterior of the input as column vectors.
    \item The mean of the marginal posterior of the input as a vector.
}
}
\examples{
# Artificial sample data
set.seed(987)
sampleData <- matrix(stats::rnorm(100), nrow = 10)
sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5

# Sampling from a multivariate t-distribution
t_dist_samp <- rmvtDCT(object = sampleData, lambda = 1, sigma = 10,
                       nu0 = 50, ns = 1000)
}
