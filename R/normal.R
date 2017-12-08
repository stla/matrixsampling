#' Matrix normal sampler
#'
#' Samples the matrix normal distribution.
#'
#' @param n sample size, a positive integer
#' @param M mean matrix, without constraints
#' @param U columns covariance matrix, a positive semidefinite matrix of order equal
#' to \code{nrow(M)}
#' @param V rows covariance matrix, a positive semidefinite matrix of order equal
#' to \code{ncol(M)}
#' @param checkSymmetry logical, whether to check the symmetry of \code{U} and \code{V}
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' m <- 3
#' p <- 2
#' M <- matrix(1, m, p)
#' U <- toeplitz(m:1)
#' V <- toeplitz(p:1)
#' MNsims <- rmatrixnormal(10000, M, U, V)
#' dim(MNsims) # 3 2 10000
#' apply(MNsims, 1:2, mean) # approximates M
#' vecMNsims <- t(apply(MNsims, 3, function(X) c(t(X))))
#' round(cov(vecMNsims)) # approximates kronecker(U,V)
rmatrixnormal <- function(n, M, U, V, checkSymmetry=TRUE){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  Uroot <- matrixroot(U, matrixname = "U", symmetric=!checkSymmetry)
  Vroot <- matrixroot(V, matrixname = "V", symmetric=!checkSymmetry)
  m <- nrow(Uroot)
  p <- nrow(Vroot)
  M <- as.matrix(M)
  if(m != nrow(M) || p != ncol(M)){
    stop("Incorrect dimensions")
  }
  out <- array(NA_real_, dim=c(m,p,n))
  Z <- array(rnorm(m*p*n), dim=c(m,p,n))
  for(i in 1:n){
    out[,,i] <- M + Uroot %*% Z[,,i] %*% Vroot
  }
  out
}

