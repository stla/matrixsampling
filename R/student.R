#' Matrix t sampler
#'
#' Samples the matrix t-distribution.
#'
#' @param n sample size, a positive integer
#' @param nu degrees of freedom, a positive number
#' @param M mean matrix, without constraints
#' @param U columns covariance matrix, a positive semidefinite matrix of order equal
#' to \code{nrow(M)}
#' @param V rows covariance matrix, a positive semidefinite matrix of order equal
#' to \code{ncol(M)}
#' @param checkSymmetry logical, whether to check the symmetry of \code{U} and \code{V}
#' @param keep logical, whether to return an array with class \pkg{\link[keep]{keep}}
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm
#' @importFrom keep as.karray
#'
#' @note When \code{p=1} and \code{V=nu} or when \code{m=1} and \code{U=nu}, the
#' distribution is the multivariate t-distribution.
#'
#' @examples
#' nu <- 4
#' m <- 2
#' p <- 3
#' M <- matrix(1, m, p)
#' U <- toeplitz(m:1)
#' V <- toeplitz(p:1)
#' Tsims <- rmatrixt(10000, nu, M, U, V)
#' dim(Tsims) # 2 3 10000
#' apply(Tsims, 1:2, mean) # approximates M
#' vecTsims <- t(apply(Tsims, 3, function(X) c(t(X))))
#' round(cov(vecTsims), 1) # approximates 1/(nu-2) * kronecker(U,V)
#' ## the `keep` class is nice when m=1 or p=1:
#' Tsims <- rmatrixt(2, nu, M=1:3, U=diag(3), V=1)
#' Tsims[,,1] # dimensions 3 1
#' # without `keep`, dimensions are lost:
#' rmatrixt(2, nu, M=1:3, U=diag(3), V=1, keep=FALSE)[,,1]
rmatrixt <- function(n, nu, M, U, V, checkSymmetry=TRUE, keep=TRUE){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isRealScalar(nu) || nu <= 0){
    stop("`nu` must be a positive number")
  }
  Vroot <- matrixroot(V, matrixname="V", symmetric=!checkSymmetry)
  U <- as.matrix(U)
  m <- nrow(U)
  p <- nrow(Vroot)
  M <- as.matrix(M)
  if(m != nrow(M) || p != ncol(M)){
    stop("Incorrect dimensions")
  }
  IW <- rinvwishart(n, nu+m-1, U, checkSymmetry=checkSymmetry)
  out <- array(NA_real_, dim=c(m,p,n))
  Z <- array(rnorm(m*p*n), dim=c(m,p,n))
  for(i in 1:n){
    out[,,i] <- M + crossprod(chol(IW[,,i]), Z[,,i] %*% Vroot)
  }
  if(keep){
    as.karray(out)
  }else{
    out
  }
}

#' Matrix inverted-t sampler
#'
#' Samples the matrix inverted-t distribution.
#'
#' @param n sample size, a positive integer
#' @param nu degrees of freedom, any positive number or an
#' integer strictly greater than \code{1-nrow(M)}
#' @param M mean matrix, without constraints
#' @param U columns covariance matrix, a positive semidefinite matrix of order equal
#' to \code{nrow(M)}
#' @param V rows covariance matrix, a positive semidefinite matrix of order equal
#' to \code{ncol(M)}
#' @param checkSymmetry logical, whether to check the symmetry of \code{U} and \code{V}
#' @param keep logical, whether to return an array with class \pkg{\link[keep]{keep}}
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm
#' @importFrom keep as.karray
#'
#' @examples
#' nu <- 0
#' m <- 2
#' p <- 3
#' M <- matrix(1, m, p)
#' U <- toeplitz(m:1)
#' V <- toeplitz(p:1)
#' ITsims <- rmatrixit(10000, nu, M, U, V)
#' dim(ITsims) # 2 3 10000
#' apply(ITsims, 1:2, mean) # approximates M
#' vecITsims <- t(apply(ITsims, 3, function(X) c(t(X))))
#' round(cov(vecITsims),2) # approximates 1/(nu+m+p-1) * kronecker(U,V)
rmatrixit <- function(n, nu, M, U, V, checkSymmetry=TRUE, keep=TRUE){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  Uroot <- matrixroot(U, matrixname="U", symmetric=!checkSymmetry)
  Vroot <- matrixroot(V, matrixname="V", symmetric=!checkSymmetry)
  m <- nrow(Uroot)
  p <- nrow(Vroot)
  M <- as.matrix(M)
  if(m != nrow(M) || p != ncol(M)){
    stop("Incorrect dimensions")
  }
  if(!isRealScalar(nu) || nu <= 1-m){
    stop("`nu` must be a number strictly greater than `1-nrow(M)`")
  }
  W <- rwishart_I(n, nu+m-1, m)
  out <- array(NA_real_, dim=c(m,p,n))
  Z <- array(rnorm(m*p*n), dim=c(m,p,n))
  if(m>1){
    for(i in 1:n){
      out[,,i] <- M + Uroot %*%
        forwardsolve(t(chol(W[,,i] + tcrossprod(Z[,,i]))), diag(m)) %*% Z[,,i] %*% Vroot
    }
  }else{
    for(i in 1:n){
      out[,,i] <- M + Uroot %*%
        forwardsolve(t(chol(W[,,i] + crossprod(Z[,,i]))), diag(m)) %*% Z[,,i] %*% Vroot
    }
  }
  if(keep){
    as.karray(out)
  }else{
    out
  }
}
