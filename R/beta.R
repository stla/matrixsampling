#' Matrix Beta sampler
#'
#' Samples a matrix Beta distribution.
#'
#' @param n sample size, a positive integer
#' @param p dimension, a positive integer
#' @param a,b parameters of the distribution, positive numbers with constraints given
#' in Details
#' @param Theta1 numerator noncentrality parameter, a positive semidefinite real
#' matrix of order \code{p}; setting it to \code{NULL} (default) is
#' equivalent to setting it to the zero matrix
#' @param Theta2 denominator noncentrality parameter, a positive semidefinite real
#' matrix of order \code{p}; setting it to \code{NULL} (default) is
#' equivalent to setting it to the zero matrix
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#'
#' @details Parameters \code{a} and \code{b} are positive numbers that satisfy the
#' following constraints:
#' \itemize{
#' \item if both \code{Theta1} and \code{Theta2} are the null matrix,
#' \code{a+b >= p/2}; if \code{a < p/2}, it must be half an integer;
#' if \code{b < p/2}, it must be half an integer
#' \item if \code{Theta1} is not the null matrix, \code{a >= p/2};
#' if \code{b < p/2}, it must be half an integer
#' \item if \code{Theta2} is the null matrix, \code{b >= p/2};
#' if \code{a < p/2}, it must be half an integer}
#'
#' @examples
#' Bsims <- rmatrixbeta(10000, 3, 1, 1)
#' dim(Bsims) # 3 3 10000
rmatrixbeta <- function(n, p, a, b, Theta1=NULL, Theta2=NULL){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isPositiveInteger(p)){
    stop("`p` must be a positive integer")
  }
  if(!isRealScalar(a) || a <= 0){
    stop("`a` must be a positive number")
  }
  if(!isRealScalar(b) || b <= 0){
    stop("`b` must be a positive number")
  }
  if(isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    if(2*a+2*b < p){
      stop("`a` and `b` must satisfy `a+b >= p/2`")
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(!isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(isNullOrZeroMatrix(Theta1) && !isNullOrZeroMatrix(Theta2)){
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p, Theta2)
  }else{
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p, Theta2)
  }
  out <- array(NA_real_, dim=c(p,p,n))
  for(i in 1:n){
    A <- invsqrtm(W1[,,i] + W2[,,i])
    out[,,i] <- A %*% W1[,,i] %*% A
  }
  out
}

#' Matrix Beta II sampler
#'
#' Samples a matrix Beta type II distribution.
#'
#' @param n sample size, a positive integer
#' @param p dimension, a positive integer
#' @param a,b parameters of the distribution, positive numbers with constraints given
#' in Details
#' @param Theta1 numerator noncentrality parameter, a positive semidefinite real
#' matrix of order \code{p}; setting it to \code{NULL} (default) is
#' equivalent to setting it to the zero matrix
#' @param Theta2 denominator noncentrality parameter, a positive semidefinite real
#' matrix of order \code{p}; setting it to \code{NULL} (default) is
#' equivalent to setting it to the zero matrix
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#'
#' @details Parameters \code{a} and \code{b} are positive numbers that satisfy the
#' following constraints:
#' \itemize{
#' \item in any case, \code{b >= p/2}
#' \item if \code{Theta1} is not the null matrix, \code{a >= p/2}
#' \item if \code{Theta1} is the null matrix and \code{a < p/2}, then \code{a}
#' must be half an integer}
#'
#' @examples
#' Bsims <- rmatrixbetaII(10000, 3, 1, 1.5)
#' dim(Bsims) # 3 3 10000
rmatrixbetaII <- function(n, p, a, b, Theta1=NULL, Theta2=NULL){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isPositiveInteger(p)){
    stop("`p` must be a positive integer")
  }
  if(!isRealScalar(a) || a <= 0){
    stop("`a` must be a positive number")
  }
  if(!isRealScalar(b) || b <= 0){
    stop("`b` must be a positive number")
  }
  if(isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    if(2*b < p){
      stop("`b` must satisfy `b >= p/2`")
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(!isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(isNullOrZeroMatrix(Theta1) && !isNullOrZeroMatrix(Theta2)){
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p, Theta2)
  }else{
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p, Theta2)
  }
  out <- array(NA_real_, dim=c(p,p,n))
  for(i in 1:n){
    A <- invsqrtm(W2[,,i])
    out[,,i] <- A %*% W1[,,i] %*% A
  }
  out
}
