context("rcpp")

test_that("compare previous dctMatrix and rcpp implementation", {
  n <- 10

  r_dctMatrix <- function(n){
    dctMat <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i == 1) {
          dctMat[i, j] <- 1 / sqrt(n)
        }
        if (i > 1) {
          dctMat[i, j] <- sqrt(2 / n) * cos((2 * (j - 1) + 1) * (i - 1) * pi / (2 * n))
        }
      }
    }
    return(dctMat)
  }

  expect_equal(r_dctMatrix(n), mrbsizeR::dctMatrix(n))
})


test_that("compare previous eigenLaplace and rcpp implementation", {
  r_eigenLaplace <- function(mm, nn) {
    N <- mm * nn
    Mu <- rep(0, N)
    In <- seq(0, mm - 1, 1)
    Im <- seq(0, nn - 1, 1)
    In <- t(In)
    Im <- t(Im)
    mu <- 2 - 2 * cos(pi * Im / nn)
    lambda <- 2 - 2 * cos(pi * In / mm)
    for (k in 1:nn) {
      for (s in 1:mm) {
        Mu[(k - 1) * mm + s] <- lambda[s] + mu[k]
      }
    }
    return(Mu)
  }

  expect_equal(r_eigenLaplace(20,20), mrbsizeR::eigenLaplace(20,20))
})


test_that("compare previous initLambdaMat loops and rcpp implementation", {
  nGrid <- 100
  N <- 100
  fac <- rnorm(N)
  D <- rnorm(N)
  LambdaMat <- matrix(0, nrow = N, ncol = nGrid)
  lambda <- 1:nGrid

  for (i in 1:nGrid) {
    l1 <- lambda[i]
    lambda2 <- fac / (1 + l1 * D)
    lambda2[1] <- 0
    LambdaMat[ ,i] <- lambda2
  }

  expect_equal(LambdaMat, mrbsizeR:::initLambdaMat(nGrid, N, lambda, fac, D))
})


test_that("compare previous min2Lambda loops and rcpp implementation", {
  nGrid <- 100
  N <- 100
  fac <- rnorm(N)
  D <- rnorm(N)
  LambdaMat <- matrix(0, nrow = N, ncol = nGrid)
  lambda <- 1:nGrid

  for (i in 1:nGrid) {
    l1 <- lambda[i]
    lambda2 <- fac / (1 + l1 * D)
    lambda2[1] <- 0
    LambdaMat[ ,i] <- lambda2
  }

  minimum <- 10^11
  nLambda <- 2
  G <- array(NA, c(rep(length(lambda), nLambda)))
  lambda1 <- fac
  lambda1[1] <- 0
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

  minimum <- 10^11
  min2LambdaListoutput <- mrbsizeR:::min2Lambda(nGrid, LambdaMat, lambda1, minimum)
  min2LambdaListoutput$G <- apply(min2LambdaListoutput$G, c(1,2), function(x) { if(identical(x, 0)) { x <- NA } else { x <- x}})
  expect_equal(G, min2LambdaListoutput$G)
  expect_equal(mini, min2LambdaListoutput$mini)
  expect_equal(minj, min2LambdaListoutput$minj)
})





test_that("compare previous min3Lambda loops and rcpp implementation", {
  nGrid <- 100
  N <- 100
  fac <- rnorm(N)
  D <- rnorm(N)
  LambdaMat <- matrix(0, nrow = N, ncol = nGrid)
  lambda <- 1:nGrid
  
  for (i in 1:nGrid) {
    l1 <- lambda[i]
    lambda2 <- fac / (1 + l1 * D)
    lambda2[1] <- 0
    LambdaMat[ ,i] <- lambda2
  }
  
  minimum <- 10^11
  nLambda <- 3
  G <- array(NA, c(rep(length(lambda), nLambda)))
  lambda1 <- fac
  lambda1[1] <- 0

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

  mini <- 0
  minj <- 0
  mink <- 0
  minimum <- 10^11
  test <- sapply(0:(nGrid-2), function(i) {
    lambda2 <- LambdaMat[ ,i+1]
    min3LambaListoutput <- mrbsizeR:::min3Lambda(i, nGrid, LambdaMat, lambda1, lambda2, minimum)
    if (min3LambaListoutput$minimum < minimum) {
      minimum <- min3LambaListoutput$minimum
      mini <- min3LambaListoutput$mini
      minj <- min3LambaListoutput$minj
      mink <- min3LambaListoutput$mink
    }

    min3LambaListoutput$G <- apply(min3LambaListoutput$G, c(1,2), function(x) { if(identical(x, 0)) { x <- NA } else { x <- x}})
    return(list(Gtmp = min3LambaListoutput$G, mini = min3LambaListoutput$mini, minj = min3LambaListoutput$minj, mink = min3LambaListoutput$mink))
  })

  Gtest <- array(NA, c(rep(length(lambda), nLambda)))
  for(i in 1:(nGrid-2)) {
    Gtest[i, , ] <- apply(test[,i]$Gtmp, c(1,2), function(x) { if(identical(x, 0)) { x <- NA } else { x <- x}})
  }

  all.equal(G, Gtest)
})
