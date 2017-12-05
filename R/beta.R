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
#' \code{a+b > (p-1)/2}; if \code{a <= (p-1)/2}, it must be half an integer;
#' if \code{b <= (p-1)/2}, it must be half an integer
#' \item if \code{Theta1} is the null matrix, \code{a > (p-1)/2} and \code{a}
#' must be half an integer if \code{a < p-1/2};
#' if \code{b <= (p-1)/2}, it must be half an integer
#' \item if \code{Theta2} is the null matrix, \code{b > (p-1)/2} and \code{b}
#' must be half an integer if \code{b < p-1/2};
#' if \code{a <= (p-1)/2}, it must be half an integer}
#'
#' @note The matrix variate Beta distribution is usually defined only for
#' \eqn{a > (p-1)/2} and \eqn{b > (p-1)/2}. In this case, a random matrix \eqn{U}
#' generated from this distribution satisfies \eqn{0 < U < I}.
#' For an half integer \eqn{a \le (p-1)/2}, it satisfies \eqn{0 \le U < I}
#' and \eqn{rank(U)=2a}.
#' For an half integer \eqn{b \le (p-1)/2}, it satisfies \eqn{0 < U \le I}
#' and \eqn{rank(I-U)=2b}.
#'
#' @examples
#' Bsims <- rmatrixbeta(10000, 3, 1, 1)
#' dim(Bsims) # 3 3 10000
rmatrixbeta <- function(n, p, a, b, Theta1=NULL, Theta2=NULL, def=2){
  def <- match.arg(as.character(def), choices=1:2)
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
    if(2*a+2*b <= p-1){
      stop("`a` and `b` must satisfy `a+b > (p-1)/2`")
    }
    W1root <- rwishart_chol_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p)
    out <- array(NA_real_, dim=c(p,p,n))
    for(i in 1:n){
      out[,,i] <- W1root[,,i] %*%
        chol2inv(chol(tcrossprod(W1root[,,i]) + W2[,,i])) %*% t(W1root[,,i])
    }
    return(out)
  }else if(!isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(isNullOrZeroMatrix(Theta1) && !isNullOrZeroMatrix(Theta2)){
    W2 <- rwishart_I(n, 2*b, p, Theta2)
    if(def==2){
      W1root <- rwishart_chol_I(n, 2*a, p)
      out <- array(NA_real_, dim=c(p,p,n))
      for(i in 1:n){
        out[,,i] <- W1root[,,i] %*%
          chol2inv(chol(tcrossprod(W1root[,,i]) + W2[,,i])) %*% t(W1root[,,i])
      }
    }
    W1 <- rwishart_I(n, 2*a, p)
  }else{
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p, Theta2)
  }
  out <- array(NA_real_, dim=c(p,p,n))
  if(def==1){
    for(i in 1:n){
      A <- invsqrtm(W1[,,i] + W2[,,i])
      out[,,i] <- A %*% W1[,,i] %*% A
    }
  }else{
    for(i in 1:n){
      A <- chol(W1[,,i])
      out[,,i] <- A %*% chol2inv(chol(W1[,,i] + W2[,,i])) %*% t(A)
    }
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
#' \item in any case, \code{b > (p-1)/2}
#' \item if \code{Theta1} is the null matrix and \code{a <= (p-1)/2}, then
#' \code{a} must be half an integer
#' \item if \code{Theta1} is not the null matrix, \code{a > (p-1)/2} and
#' \code{a} must be half an integer if \code{a < p-1/2}
#' \item if \code{Theta2} is not the null matrix,
#' \code{b} must be half an integer if \code{b < p-1/2}}
#'
#' @note The matrix variate Beta distribution of type II is usually defined only for
#' \eqn{a > (p-1)/2} and \eqn{b > (p-1)/2}. In this case, a random matrix \eqn{V}
#' generated from this distribution satisfies \eqn{V > 0}.
#' For an half integer \eqn{a \le (p-1)/2}, it satisfies \eqn{V \ge 0} and
#' \eqn{rank(V)=2a}.
#'
#' @examples
#' Bsims <- rmatrixbetaII(10000, 3, 1, 1.5)
#' dim(Bsims) # 3 3 10000
rmatrixbetaII <- function(n, p, a, b, Theta1=NULL, Theta2=NULL, def=1){
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
  if(2*b <= p-1){
    stop("`b` must satisfy `b > (p-1)/2`")
  }
  if(isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
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
  if(def==1){
    for(i in 1:n){
      A <- invsqrtm(W2[,,i])
      out[,,i] <- A %*% W1[,,i] %*% A
    }
  }else{
    for(i in 1:n){
      A <- matrixroot(W1[,,i]) # con dans le cas Theta1=0
      out[,,i] <- A %*% chol2inv(chol(W2[,,i])) %*% A # con dans le cas Theta2=0
    }
  }
  out
}
