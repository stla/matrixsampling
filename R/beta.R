#' Matrix Beta sampler
#'
#' Samples a matrix Beta (type I) distribution.
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
#' @param def \code{1} or \code{2}, the definition used; see Details
#' @param checkSymmetry logical, whether to check the symmetry of \code{Theta1}
#' and \code{Theta2}
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#'
#' @details A Beta random matrix \eqn{U} is defined as follows.
#' Take two independent Wishart random matrices
#' \ifelse{html}{\out{S<sub>1</sub> ~ W<sub>p</sub>(2a,I<sub>p</sub>,&Theta;<sub>1</sub>)}}{\eqn{S_1 \sim \mathcal{W}_p(2a,I_p,\Theta_1)}}
#' and
#' \ifelse{html}{\out{S<sub>2</sub> ~ W<sub>p</sub>(2b,I<sub>p</sub>,&Theta;<sub>2</sub>)}}{\eqn{S_2 \sim \mathcal{W}_p(2b,I_p,\Theta_2)}}.
#' \itemize{
#' \item \strong{definition 1}:
#' \ifelse{html}{\out{U = (S<sub>1</sub>+S<sub>2</sub>)<sup>-&frac12;</sup>S<sub>1</sub>(S<sub>1</sub>+S<sub>2</sub>)<sup>-&frac12;</sup>}}{\eqn{U = {(S_1+S_2)}^{-\frac12}S_1{(S_1+S_2)}^{-\frac12}}}
#' \item \strong{definition 2}:
#' \ifelse{html}{\out{U = S<sub>1</sub><sup>&frac12;</sup>(S<sub>1</sub>+S<sub>2</sub>)<sup>-1</sup>S<sub>1</sub><sup>&frac12;</sup>}}{\eqn{U=S_1^\frac12{(S_1+S_2)}^{-1}S_1^\frac12}}
#' }
#' In the central case, the two definitions yield the same distribution.
#' Under definition 2, the Beta distribution is related to the Beta type II
#' distribution by
#' \ifelse{html}{\out{U ~ V(I+V)<sup>-1</sup>}}{\eqn{U \sim V{(I+V)}^{-1}}}.
#'
#' Parameters \code{a} and \code{b} are positive numbers that satisfy the
#' following constraints:
#' \itemize{
#' \item if both \code{Theta1} and \code{Theta2} are the null matrix,
#' \code{a+b > (p-1)/2}; if \code{a < (p-1)/2}, it must be half an integer;
#' if \code{b < (p-1)/2}, it must be half an integer
#' \item if \code{Theta1} is not the null matrix, \code{a >= (p-1)/2};
#' if \code{b < (p-1)/2}, it must be half an integer
#' \item if \code{Theta2} is not the null matrix, \code{b >= (p-1)/2};
#' if \code{a < (p-1)/2}, it must be half an integer}
#'
#' @note The matrix variate Beta distribution is usually defined only for
#' \eqn{a > (p-1)/2} and \eqn{b > (p-1)/2}. In this case, a random matrix \eqn{U}
#' generated from this distribution satisfies \eqn{0 < U < I}.
#' For an half integer \eqn{a \le (p-1)/2}, it satisfies \eqn{0 \le U < I}
#' and \eqn{rank(U)=2a}.
#' For an half integer \eqn{b \le (p-1)/2}, it satisfies \eqn{0 < U \le I}
#' and \eqn{rank(I-U)=2b}.
#'
#' @section Warning:
#' Definition 2 requires the calculation of the square root of
#' \ifelse{html}{\out{S<sub>1</sub> ~ W<sub>p</sub>(2a,I<sub>p</sub>,&Theta;<sub>1</sub>)}}{\eqn{S_1 \sim \mathcal{W}_p(2a,I_p,\Theta_1)}}
#' (see Details). While \ifelse{html}{\out{S<sub>1</sub>}}{\eqn{S_1}} is always
#' positive semidefinite in theory, it could happen that the simulation of
#' \ifelse{html}{\out{S<sub>1</sub>}}{\eqn{S_1}} is not positive semidefinite,
#' especially when \code{a} is small. In this case the calculation of the square root
#' will return \code{NaN}.
#'
#' @examples
#' Bsims <- rmatrixbeta(10000, 3, 1, 1)
#' dim(Bsims) # 3 3 10000
rmatrixbeta <- function(n, p, a, b, Theta1=NULL, Theta2=NULL, def=1, checkSymmetry=TRUE){
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
    if(2*a < p-1 && floor(2*a) != 2*a){
      stop("`a < (p-1)/2`, it must be half an integer")
    }
    if(2*b < p-1 && floor(2*b) != 2*b){
      stop("`b < (p-1)/2`, it must be half an integer")
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(!isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    if(2*a < p-1){
      stop("`a` must be greater than `(p-1)/2`")
    }
    if(2*b < p-1 && floor(2*b) != 2*b){
      stop("`b < (p-1)/2`, it must be half an integer")
    }
    W1 <- rwishart_I(n, 2*a, p, Theta1, checkSymmetry=checkSymmetry)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(isNullOrZeroMatrix(Theta1) && !isNullOrZeroMatrix(Theta2)){
    if(2*a < p-1 && floor(2*a) != 2*a){
      stop("`a < (p-1)/2`, it must be half an integer")
    }
    if(2*b < p-1){
      stop("`b` must be greater than `(p-1)/2`")
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p, Theta2, checkSymmetry=checkSymmetry)
  }else{
    if(2*a < p-1){
      stop("`a` must be greater than `(p-1)/2`")
    }
    if(2*b < p-1){
      stop("`b` must be greater than `(p-1)/2`")
    }
    W1 <- rwishart_I(n, 2*a, p, Theta1, checkSymmetry=checkSymmetry)
    W2 <- rwishart_I(n, 2*b, p, Theta2, checkSymmetry=checkSymmetry)
  }
  out <- array(NA_real_, dim=c(p,p,n))
  if(def==1){
    for(i in 1:n){
      A <- invsqrtm(W1[,,i] + W2[,,i])
      out[,,i] <- A %*% W1[,,i] %*% A
    }
  }else{
    for(i in 1:n){
      A <- sqrtm(W1[,,i]) # not cholesky
      out[,,i] <- A %*% chol2inv(chol(W1[,,i] + W2[,,i])) %*% A
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
#' @param def \code{1} or \code{2}, the definition used; see Details
#' @param checkSymmetry logical, whether to check the symmetry of \code{Theta1}
#' and \code{Theta2}
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#'
#' @details A Beta type II random matrix \eqn{V} is defined as follows.
#' Take two independent Wishart random matrices
#' \ifelse{html}{\out{S<sub>1</sub> ~ W<sub>p</sub>(2a,I<sub>p</sub>,&Theta;<sub>1</sub>)}}{\eqn{S_1 \sim \mathcal{W}_p(2a,I_p,\Theta_1)}}
#' and
#' \ifelse{html}{\out{S<sub>2</sub> ~ W<sub>p</sub>(2b,I<sub>p</sub>,&Theta;<sub>2</sub>)}}{\eqn{S_2 \sim \mathcal{W}_p(2b,I_p,\Theta_2)}}.
#' \itemize{
#' \item \strong{definition 1}:
#' \ifelse{html}{\out{V = S<sub>2</sub><sup>-&frac12;</sup>S<sub>1</sub>S<sub>2</sub><sup>-&frac12;</sup>}}{\eqn{V = S_2^{-\frac12}S_1S_2^{-\frac12}}}
#' \item \strong{definition 2}:
#' \ifelse{html}{\out{V = S<sub>1</sub><sup>&frac12;</sup>S<sub>2</sub><sup>-1</sup>S<sub>1</sub><sup>&frac12;</sup>}}{\eqn{V = S_1^\frac12 S_2^{-1}S_1^\frac12}}
#' }
#' In the central case, the two definitions yield the same distribution.
#' Under definition 2, the Beta type II distribution is related to the Beta
#' distribution by
#' \ifelse{html}{\out{V ~ U(I-U)<sup>-1</sup>}}{\eqn{V \sim U{(I-U)}^{-1}}}.
#'
#' Parameters \code{a} and \code{b} are positive numbers that satisfy the
#' following constraints:
#' \itemize{
#' \item in any case, \code{b > (p-1)/2}
#' \item if \code{Theta1} is the null matrix and \code{a < (p-1)/2}, then
#' \code{a} must be half an integer
#' \item if \code{Theta1} is not the null matrix, \code{a >= (p-1)/2}}
#'
#' @note The matrix variate Beta distribution of type II is usually defined only for
#' \eqn{a > (p-1)/2} and \eqn{b > (p-1)/2}. In this case, a random matrix \eqn{V}
#' generated from this distribution satisfies \eqn{V > 0}.
#' For an half integer \eqn{a \le (p-1)/2}, it satisfies \eqn{V \ge 0} and
#' \eqn{rank(V)=2a}.
#'
#' @section Warning:
#' The issue described in the \strong{Warning} section of \code{\link{rmatrixbeta}}
#' also concerns \code{rmatrixbetaII}.
#'
#' @examples
#' Bsims <- rmatrixbetaII(10000, 3, 1, 1.5)
#' dim(Bsims) # 3 3 10000
rmatrixbetaII <- function(n, p, a, b, Theta1=NULL, Theta2=NULL, def=1, checkSymmetry=TRUE){
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
  if(2*b <= p-1){
    stop("`b` must satisfy `b > (p-1)/2`")
  }
  if(isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    if(2*a < p-1 && floor(2*a) != 2*a){
      stop("`a < (p-1)/2`, it must be half an integer")
    }
    if(def==2){
      out <- array(NA_real_, dim=c(p,p,n))
      W2root <- rwishart_chol_I(n, 2*b, p, upper=TRUE)
      if(2*a > p-1){
        W1root <- rwishart_chol_I(n, 2*a, p)
        for(i in 1:n){
          out[,,i] <- W1root[,,i] %*% chol2inv(W2root[,,i]) %*% t(W1root[,,i])
        }
      }else{
        W1 <- rwishart_I(n, 2*a, p)
        for(i in 1:n){
          W1root <- sqrtm(W1[,,i])
          out[,,i] <- W1root %*% chol2inv(W2root[,,i]) %*% W1root
        }
      }
      return(out)
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(!isNullOrZeroMatrix(Theta1) && isNullOrZeroMatrix(Theta2)){
    if(2*a < p-1){
      stop("`a` must be greater than `(p-1)/2`")
    }
    W1 <- rwishart_I(n, 2*a, p, Theta1, checkSymmetry=checkSymmetry)
    if(def==2){
      W2root <- rwishart_chol_I(n, 2*b, p, upper=TRUE)
      out <- array(NA_real_, dim=c(p,p,n))
      for(i in 1:n){
        W1root <- sqrtm(W1[,,i]) # could use chol, but sqrtm better if not >0
        out[,,i] <- W1root %*% chol2inv(W2root[,,i]) %*% W1root
      }
      return(out)
    }
    W2 <- rwishart_I(n, 2*b, p)
  }else if(isNullOrZeroMatrix(Theta1) && !isNullOrZeroMatrix(Theta2)){
    if(2*a <= p-1 && floor(2*a) != 2*a){
      stop("`a <= (p-1)/2`, it must be half an integer")
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p, Theta2, checkSymmetry=checkSymmetry)
  }else{
    if(2*a < p-1){
      stop("`a` must be greater than `(p-1)/2`")
    }
    W1 <- rwishart_I(n, 2*a, p, Theta1, checkSymmetry=checkSymmetry)
    W2 <- rwishart_I(n, 2*b, p, Theta2, checkSymmetry=checkSymmetry)
  }
  out <- array(NA_real_, dim=c(p,p,n))
  if(def==1){
    for(i in 1:n){
      A <- invsqrtm(W2[,,i])
      out[,,i] <- A %*% W1[,,i] %*% A
    }
  }else{
    for(i in 1:n){
      A <- sqrtm(W1[,,i]) # not chol if Theta2 is not scalar
      out[,,i] <- A %*% chol2inv(chol(W2[,,i])) %*% A
    }
  }
  out
}
