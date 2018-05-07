#' Generate eigenvalues of precision matrix Q on the surface of a sphere.
#'
#' The eigenvalues of the precision matrix Q with dimension (\code{mm}, \code{nn}) 
#' and polar angle limits \code{phimin}, \code{phimax} are calculated.
#'
#' The corresponding function for data on a grid is \code{ \link{eigenLaplace}}.
#'
#' @param phimin Polar angle minimum. 
#' @param phimax Polar angle maximum.
#' @param mm Number of rows of precision matrix Q.
#' @param nn Number of columns of precision matrix Q.
#' @return A list containing 2 elements: 
#' \itemize{
#'     \item eigval Row vector containing the eigenvalues of Q. 
#'     \item eigvec Matrix containing the eigenvectors of Q as columns. 
#' }
#' @examples
#' eig_out <- eigenQsphere(phimin = 180/10, phimax = 180 - 180/10, mm = 10, nn = 20)
#' 
eigenQsphere <- function(phimin, phimax, mm, nn) {
  
  #------------ Computing variables needed later --------------------------#
    
  phi <- seq(phimin, phimax, len = mm)
  deltaPhi <- pi / 180 * (phi[2] - phi[1])
  
  sinPhi2 <- -1 / sin(phi* pi / 180)^2
  N <- mm * nn
  
  c1 <- tridiag(-2 * rep(1, mm), rep(1, mm-1), rep(1, mm-1)) / deltaPhi^2

  c3 <- diag(sinPhi2)
  CNeven <- c1
  CNeven[1, 1] <- CNeven[1, 1] + 1 / deltaPhi^2
  CNeven[mm, min(mm, nn)] <- CNeven[mm, min(mm, nn)] + 1 / deltaPhi^2
  
  Ctotal <- matrix(NA, nrow = mm * nn, ncol = mm * nn)
  eigen <- rep(0, mm * nn)
  eigv <- matrix(0, nrow = mm * nn, ncol = mm * nn)
  
  #------Computing the eigenvalues and eigenvectors -----------#
    
  for (i in (1:nn)) {
    Ceven <- CNeven + c3 * (i - (1 + nn / 2))^2
    Ctotal[((i - 1) * mm + 1):(i * mm),((i - 1) * mm + 1):(i * mm)] <- Ceven
    v <- apply(eigen(t(Ceven) %*% Ceven)$vectors, 2, rev)[ , ncol(eigen(t(Ceven) %*% Ceven)$vectors):1]
    e <- rev(eigen(t(Ceven) %*% Ceven)$values)
    eigen[((i - 1) * mm + 1):(mm * i)] <- e
    eigv[((i - 1) * mm + 1):(i * mm), ((i - 1) * mm + 1):(i * mm)] <- v  
  }
  

  for (i in (1:N)) {
    eigvmtrx <- fftshift(matrix(eigv[ ,i], nrow = mm), dimension = 2)
    eigv[ ,i] <- eigvmtrx[ , ]
  }
  
  eigV <- kronecker(Conj(dftMatrix(nn)) / sqrt(nn), diag(1, mm, mm)) %*% eigv
  indx <- sort(eigen, index.return = TRUE, decreasing = FALSE)$ix
  eigval <- sort(eigen, index.return = TRUE, decreasing = FALSE)$x
  eigV <- eigV[ , indx]
  
  
  #------------- Computing eigenvectors with real values ------#
    
  eigVr <- Re(eigV)
  eigVc <- Im(eigV)
  
  for (i in 1:(nn * mm)) {
    eigVr[ , i] <- eigVr[ , i] / sqrt(t(eigVr[ , i]) %*% eigVr[ , i])
  }
  
  eigV <- eigVr
  ei2 <- matrix(NA, nrow = 2, ncol = (N+1))
  
  ind <- 1
  i <- 1   

  while (i < (N + 1)) {
    valueTest <- 0
    j <- i
    while (valueTest < (300 * .Machine$double.eps)) {
      j <- j + 1
      if (j < (N + 1)) {
        valueTest <- abs(eigval[i] - eigval[j])
      } else {
        valueTest <- 1
      }
    }
    
    if ((j - i) == 2) {
      if (is.complex(eigV[ , j]) == FALSE) {
        eigV[ , i] <- eigVc[ , i] / sqrt(t(eigVc[ , i]) %*% eigVc[ , i])
      } else {
        eigV[ , j] <- eigVc[ , j] / sqrt(t(eigVc[ , j]) %*% eigVc[ , j])
      }
    } else if ((j - i) > 2) {   
        ei2[ , ind] <- c(i, j)
        ind <- ind + 1
    }

    i <- j
  }

  ei2 <- ei2[,!colSums(!is.finite(ei2))]
  
  if(ncol(ei2) >= 1){
    for (i in 1:ncol(ei2)) {
      sqq <- t(ei2[ , i])
      s <- sqq[1]
      j <- sqq[2] - 1
      difff <- j - s + 1
      matrxx <- matrix(0, N, (2 * difff))
      
      for (k in 1:difff) {
        matrxx[ , (1 + (k - 1) * 2):(k * 2)] <- c(eigVr[ , (s + (k - 1))], eigVc[ , (s + (k - 1))])
      }
      
      tmp <- svd(matrxx) 
      U <- tmp$u 
      V <- tmp$v 
      S <- diag(tmp$d)
      
      eigV[ , s:j] <- U[ , 1:difff] 
    }
  } 

  for (i in 1:(nn * mm)) {
    eigV[ , i] <- eigV[ , i] / sqrt(t(eigV[ , i]) %*% eigV[ , i])
  }

  return(list(eigval = eigval, eigvec = eigV))
}



