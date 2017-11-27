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
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm
#'
#' @examples
#' m <- 3
#' p <- 4
#' M <- matrix(1, m, p)
#' U <- toeplitz(m:1)
#' V <- toeplitz(p:1)
#' MNsims <- rmatrixnormal(10000, M, U, V)
#' dim(MNsims) # 3 4 10000
#' apply(MNsims, 1:2, mean) # approximates M
#' vecMNsims <- t(apply(MNsims, 3, function(X) c(t(X))))
#' round(cov(vecMNsims)) # approximates kronecker(U,V)
rmatrixnormal <- function(n, M, U, V){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  Uroot <- matrixroot(U, matrixname = "U")
  Vroot <- matrixroot(V, matrixname = "V")
  m <- nrow(Uroot)
  p <- nrow(Vroot)
  if(isScalar(M)){
    M <- as.matrix(M)
  }
  if(m != nrow(M) || p != ncol(M)){
    stop("Incorrect dimensions")
  }
  out <- array(NA_real_, dim=c(m,p,n))
  for(i in 1:n){
    out[,,i] <- M + Uroot %*% matrix(rnorm(m*p), m, p) %*% Vroot
  }
  out
}
