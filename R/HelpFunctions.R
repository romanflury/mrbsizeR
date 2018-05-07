#' Turn matrix 90 degrees counter-clockwise.
#'
#' Help function to turn matrix 90 degrees counter-clockwise.
#' \code{turnmat} is used as an internal function in some of the plotting 
#' functions. Depending on how the data is stored, it can be necessary to turn 
#' the matrices 90 degrees counter-clockwise for being able to plot them correctly.
#'
#' @param x Matrix to be turned.
#' @return Matrix \code{x}, turned by 90 degrees counter-clockwise.
#' @examples
#' set.seed(987)
#' sampleMat <- matrix(stats::rnorm(100), nrow = 10)
#' sampleMatTurn <- turnmat(x = sampleMat)
#'
turnmat <- function(x) {
  t(x)[, nrow(x):1]
}



#' Swap the quadrants or halves of a 2d matrix.
#' 
#' \code{fftshift} is an R equivalent to the Matlab function \code{fftshift} 
#' applied on matrices. For more information about \code{fftshift} see 
#' the Matlab documentation.
#' 
#' It is possible to swap the halves or the quadrants of the input matrix. 
#' Halves can be swapped along the rows (\code{dimension = 1}) or along the 
#' columns (\code{dimension = 2}). When swapping the quadrants, \code{fftshift} 
#' swaps the first quadrant with the third and the second quadrant with the fourth
#' (\code{dimension = -1}).
#' 
#' @param inputMatrix Matrix to be swapped.
#' @param dimension Which swap should be performed? 
#' \itemize{
#'     \item \code{1}: swap halves along the rows.
#'     \item \code{2}: swap halves along the columns. 
#'     \item \code{-1}: swap first quadrant with third and second quadrant with fourth. 
#' }
#' @return Swapped matrix.
#' @examples
#' set.seed(987) 
#' sampleMat <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5)
#' 
#' # Swap halves along the rows:
#' fftshift(sampleMat, dimension = 1)
#' 
#' # Swap halves along the columns:
#' fftshift(sampleMat, dimension = 2)
#' 
#' # Swap first quadrant with third and second quadrant with fourth:
#' fftshift(sampleMat, dimension = -1)
#' 

# Function inspired by http://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r
fftshift <- function(inputMatrix, dimension = -1) {
  rows <- dim(inputMatrix)[1]
  cols <- dim(inputMatrix)[2]

  swapUpDown <- function(inputMatrix) {
    rowsHalf <- ceiling(rows / 2)
    return(rbind(inputMatrix[((rowsHalf + 1):rows), (1:cols)], inputMatrix[(1:rowsHalf), (1:cols)]))
  }

  swapLeftRight <- function(inputMatrix) {
    colsHalf <- ceiling(cols / 2)
    return(cbind(inputMatrix[1:rows, ((colsHalf + 1):cols)], inputMatrix[1:rows, 1:colsHalf]))
  }

  if (dimension == -1) {
    inputMatrix <- swapUpDown(inputMatrix)
    return(swapLeftRight(inputMatrix))
  }
  else if (dimension == 1) {
    return(swapUpDown(inputMatrix))
  }
  else if (dimension == 2) {
    return(swapLeftRight(inputMatrix))
  }
  else {
    stop("Invalid dimension parameter")
  }
}



#' Inverse FFT shift of a 2d matrix.
#' 
#' \code{ifftshift} is an R equivalent to the Matlab function \code{ifftshift} 
#' applied on matrices. For more information about \code{ifftshift} see 
#' the Matlab documentation.
#' 
#' \code{ifftshift} is the inverse function to \code{\link{fftshift}}. For more
#' information see the details of \code{\link{fftshift}}
#' 
#' @param inputMatrix Matrix to be swapped.
#' @param dimension Which swap should be performed? 
#' \itemize{
#'     \item \code{1}: swap halves along the rows.
#'     \item \code{2}: swap halves along the columns. 
#'     \item \code{-1}: swap first quadrant with third and second quadrant with fourth. 
#' }
#' @return Swapped matrix.
#' @examples
#' set.seed(987)
#' sampleMat <- matrix(sample(1:10, size = 25, replace = TRUE), nrow = 5)
#' 
#' # Swap halves along the rows:
#' ifftshift(sampleMat, dimension = 1)
#' 
#' # Swap halves along the columns:
#' ifftshift(sampleMat, dimension = 2)
#' 
#' # Swap first quadrant with third and second quadrant with fourth:
#' ifftshift(sampleMat, dimension = -1)
#' 

# Function inspired by http://stackoverflow.com/questions/38230794/how-to-write-fftshift-and-ifftshift-in-r
ifftshift <- function(inputMatrix, dimension = -1) {

  rows <- dim(inputMatrix)[1]
  cols <- dim(inputMatrix)[2]

  swapUpDown <- function(inputMatrix) {
    rowsHalf <- floor(rows/2)
    return(rbind(inputMatrix[((rowsHalf+1):rows), (1:cols)], inputMatrix[(1:rowsHalf), (1:cols)]))
  }

  swapLeftRight <- function(inputMatrix) {
    colsHalf <- floor(cols/2)
    return(cbind(inputMatrix[1:rows, ((colsHalf+1):cols)], inputMatrix[1:rows, 1:colsHalf]))
  }
  
  if (dimension == -1) {
    inputMatrix <- swapLeftRight(inputMatrix)
    return(swapUpDown(inputMatrix))
  }
  else if (dimension == 1) {
    return(swapUpDown(inputMatrix))
  }
  else if (dimension == 2) {
    return(swapLeftRight(inputMatrix))
  }
  else {
    stop("Invalid dimension parameter")
  }
}



#' Generate a tridiagonal matrix.
#' 
#' Generate a tridiagonal matrix with \code{upperDiag} as superdiagonal, 
#' \code{lowerDiag} as subdiagonal and \code{mainDiag} as diagonal. 
#' 
#' @param mainDiag Diagonal of tridiagonal matrix.
#' @param upperDiag Superdiagonal of tridiagonal matrix. Must have length \code{length(mainDiag) - 1}.
#' @param lowerDiag Subdiagonal of tridiagonal matrix. Must have length \code{length(mainDiag) - 1}.
#' 
#' @return Tridiagonal matrix. 
#' @examples 
#' set.seed(987)
#' mainDiag <- sample(100:110, size = 6, replace = TRUE)
#' upperDiag <- sample(10:20, size = 5, replace = TRUE)
#' lowerDiag <- sample(1:10, size = 5, replace = TRUE)
#' 
#' tridiag(mainDiag, upperDiag, lowerDiag)  
#'  

# Function inspired by http://stackoverflow.com/questions/28974507/efficient-creation-of-tridiagonal-matrices/28974577
tridiag <- function (mainDiag, upperDiag, lowerDiag) {

  n <- max(length(mainDiag), length(upperDiag) + 1, length(lowerDiag) + 1)
  matOut <- matrix(0, n, n)
  diag(matOut) <- mainDiag
  indx <- seq.int(n-1)
  matOut[cbind(indx + 1, indx)] <- lowerDiag
  matOut[cbind(indx, indx + 1)] <- upperDiag
  return(matOut)
}
