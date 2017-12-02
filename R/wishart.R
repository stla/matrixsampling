#' Sampling Cholesky factor of a Wishart matrix
#'
#' Samples the lower triangular Cholesky factor of a Wishart random matrix.
#'
#' @param n sample size, a positive integer
#' @param nu degrees of freedom, a positive number at least equal to the dimension
#' (the order of \code{Sigma})
#' @param Sigma scale matrix, a positive definite real matrix
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm rchisq
#'
#' @examples
#' nu <- 4
#' p <- 3
#' Sigma <- diag(p)
#' Wsims <- rwishart_chol(10000, nu, Sigma)
#' dim(Wsims) # 3 3 10000
#' Wsims[,,1]
rwishart_chol <- function(n, nu, Sigma){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isRealScalar(nu)){
    stop("`nu` must be a number")
  }
  Sigma_chol <- chol(Sigma)
  if(!isSymmetricMatrix(Sigma)){
    stop("`Sigma` is not symmetric")
  }
  p <- nrow(Sigma_chol)
  if(nu <= p-1){
    stop(sprintf("`nu+1` (%s) must be greater than the dimension (%s)", nu+1, p))
  }
  out <- array(NA_real_, dim=c(p, p, n))
  for(i in 1:n){
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
    Z[lower.tri(Z)] <- rnorm(p*(p-1)/2)
    out[,,i] <- crossprod(Sigma_chol, Z)
  }
  out
}

rwishart_root <- function(n, nu, Sigma, Sigma_root=NULL, check=TRUE){
  # samples R such that RR'~ W(nu, Sigma)
  # nu > p-1 only
  # Sigma_root Sigma_root' = Sigma
  # Sigma is ignored is Sigma_root is provided
  if(is.null(Sigma_root)){
    p <- ifelse(isScalar(Sigma), 1L, nrow(Sigma))
  }else{
    p <- ifelse(isScalar(Sigma_root), 1L, nrow(Sigma_root))
  }
  if(check){
    if(nu <= p-1){
      stop(sprintf("`nu+1` (%s) must be greater the dimension (%s)", nu+1, p))
    }
    if(is.null(Sigma_root)){
      Sigma_root <- matrixroot(Sigma)
    }else{
      if(!isSquareRealMatrix(Sigma_root)){
        stop("`Sigma_root` must be a square real matrix")
      }
    }
  }else{
    if(is.null(Sigma_root)){
      Sigma_eig <- eigen(Sigma, symmetric = TRUE)
      Sigma_root <- Sigma_eig$vectors %*% (sqrt(Sigma_eig$values) * t(Sigma_eig$vectors))
    }
  }
  out <- array(NA_real_, dim=c(p, p, n))
  for(i in 1:n){
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
    Z[lower.tri(Z)] <- rnorm(p*(p-1)/2)
    out[,,i] <- crossprod(Sigma_root, Z)
  }
  out
}

rwishart_root_I <- function(n, nu, p){
  # samples R such that RR' ~ W(nu, I_p)
  # nu > p-1 only
  out <- array(NA_real_, dim=c(p, p, n))
  for(i in 1:n){
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
    Z[lower.tri(Z)] <- rnorm(p*(p-1)/2)
    out[,,i] <- Z
  }
  out
}

#' Wishart sampler
#'
#' Samples a Wishart distribution.
#'
#' @param n sample size, a positive integer
#' @param nu degrees of freedom, a positive number;
#' if \code{nu <= p-1} where \code{p} is the dimension (the order of \code{Sigma}),
#' must be an integer;
#' in the noncentral case (i.e. when \code{Theta} is not the null matrix), \code{nu}
#' must satisfy these constraints:
#' \itemize{
#' \item \code{nu >= p}
#' \item if \code{nu <= 2*p-1}, it must be an integer}
#' @param Sigma scale matrix, a positive semidefinite real matrix
#' @param Theta noncentrality parameter, a positive semidefinite real matrix of
#' same order as \code{Sigma}; setting it to \code{NULL} (default) is
#' equivalent to setting it to the zero matrix
#'
#' @details A sampled Wishart matrix is always positive semidefinite.
#' If \code{nu > p-1}, it is positive definite.
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm rchisq
#'
#' @examples
#' nu <- 4
#' p <- 3
#' Sigma <- toeplitz(p:1)
#' Theta <- diag(p)
#' Wsims <- rwishart(10000, nu, Sigma, Theta)
#' dim(Wsims) # 3 3 10000
#' apply(Wsims, 1:2, mean) # approximately nu*Sigma+Theta
rwishart <- function(n, nu, Sigma, Theta=NULL){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isRealScalar(nu)){
    stop("`nu` must be a positive number")
  }
  Sigma_root <- matrixroot(Sigma)
  p <- nrow(Sigma_root)
  if(isNullOrZeroMatrix(Theta)){
    if(nu > p-1){
      R <- rwishart_root(n, nu, Sigma_root=Sigma_root, check=FALSE)
      out <- array(NA_real_, dim=dim(R))
      for(i in 1:n){
        out[,,i] <- tcrossprod(R[,,i])
      }
    }else{
      if(!isPositiveInteger(nu)){
        stop(sprintf("`nu` (%s) must be an integer when it is lower than `p-1` (%s)", nu, p-1))
      }
      out <- array(NA_real_, dim=c(p,p,n))
      for(i in 1:n){
        out[,,i] <- tcrossprod(Sigma_root %*% matrix(rnorm(nu*p), p, nu))
      }
    }
    return(out)
  }else{
    if(!(isScalar(Sigma) && isScalar(Theta)) && !identical(dim(Sigma), dim(Theta))){
      stop("`Sigma` and `Theta` must have the same dimension")
    }
    if(nu <= p-1){
      stop(
        sprintf(
          "`nu` (%s) must be greater than `p-1` (%s) in the noncentral case", nu, p-1))
    }
    Theta_root <- matrixroot(Theta, matrixname = "Theta")
    W <- array(NA_real_, dim=c(p,p,n))
    if(nu > 2*p-1){
      WrootI <- rwishart_root_I(n, nu-p, p)
      for(i in 1:n){
        Z <- matrix(rnorm(p*p), p, p)
        W[,,i] <- (Theta_root + Sigma_root %*% Z) %*%
          (Theta_root + t(Z) %*% Sigma_root) +
          tcrossprod(Sigma_root %*% WrootI[,,i])
      }
    }else if(nu>p){
      if(floor(nu) != nu){
        stop("In the noncentral case, `nu` must be an integer if `nu<2*p-1`")
      }
      for(i in 1:n){
        Z <- matrix(rnorm(p*p), p, p)
        W[,,i] <- (Theta_root + Sigma_root %*% Z) %*%
          (Theta_root + t(Z) %*% Sigma_root) +
          Sigma_root %*%
          tcrossprod(matrix(rnorm((nu-p)*p), p, nu-p)) %*% t(Sigma_root)
      }
    }else{ # nu=p
      for(i in 1:n){
        W[,,i] <- tcrossprod(Theta_root + Sigma_root %*% matrix(rnorm(p*p), p, p))
      }
    }
    return(W)
  }
}

rwishart_I <- function(n, nu, p, Theta=NULL){
  # samples W(nu, Ip, Theta)
  if(is.null(Theta)){
    if(nu > p-1){
      R <- rwishart_root_I(n, nu, p)
      out <- array(NA_real_, dim=c(p,p,n))
      for(i in 1:n){
        Z <- matrix(0, p, p)
        diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
        Z[lower.tri(Z)] <- rnorm(p*(p-1)/2)
        out[,,i] <- tcrossprod(Z)
      }
    }else{
      if(!isPositiveInteger(nu)){
        stop("`nu` must be an integer")
      }
      out <- array(NA_real_, dim=c(p,p,n))
      for(i in 1:n){
        out[,,i] <- tcrossprod(matrix(rnorm(nu*p), p, nu))
      }
    }
    return(out)
  }else{
    if(!(p==1 && isScalar(Theta)) && !all(c(p,p)==dim(Theta))){
      stop("`Theta` must have dimension p x p")
    }
    if(nu <= p-1){
      stop(
        sprintf(
          "`nu` (%s) must be greater than `p-1` (%s) in the noncentral case", nu, p-1))
    }
    Theta_root <- matrixroot(Theta, matrixname = "Theta")
    W <- array(NA_real_, dim=c(p,p,n))
    if(nu > 2*p-1){
      WrootI <- rwishart_root_I(n, nu-p, p)
      for(i in 1:n){
        Z <- matrix(rnorm(p*p), p, p)
        W[,,i] <- (Theta_root + Z) %*% (Theta_root + t(Z)) +
          tcrossprod(WrootI[,,i])
      }
    }else if(nu>p){
      if(!isPositiveInteger(nu)){
        stop("`nu` must be an integer")
      }
      for(i in 1:n){
        Z <- matrix(rnorm(p*p), p, p)
        W[,,i] <- (Theta_root + Z) %*% (Theta_root + t(Z)) +
          tcrossprod(matrix(rnorm((nu-p)*p), p, nu-p))
      }
    }else{
      for(i in 1:n){
        W[,,i] <- tcrossprod(Theta_root + matrix(rnorm(p*p), p, p))
      }
    }
    return(W)
  }
}

#' Inverse-Wishart sampler
#'
#' Samples the inverse-Wishart distribution.
#'
#' @param n sample size, a positive integer
#' @param nu degrees of freedom, must satisfy \code{nu > p-1}, where \code{p} is
#' the dimension (the order of \code{Omega})
#' @param Omega scale matrix, a positive definite real matrix
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm rchisq
#'
#' @details The inverse-Wishart distribution with scale matrix
#' \ifelse{html}{\out{&Omega;}}{\eqn{\Omega}} is
#' defined as the inverse of the Wishart distribution with scale matrix
#' \ifelse{html}{\out{&Omega;<sup>-1</sup>}}{\eqn{\Omega^{-1}}}
#' and same number of degrees of freedom.
#'
#' @examples
#' nu <- 6
#' p <- 3
#' Omega <- toeplitz(p:1)
#' IWsims <- rinvwishart(10000, nu, Omega)
#' dim(IWsims) # 3 3 10000
#' apply(IWsims, 1:2, mean) # approximately Omega/(nu-p-1)
rinvwishart <- function(n, nu, Omega){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isRealScalar(nu)){
    stop("`nu` must be a positive number")
  }
  if(!isSymmetricMatrix(Omega)){
    stop("`Sigma` is not symmetric")
  }
  Omegainv_chol <- chol(chol2inv(chol(Omega)))
  p <- nrow(Omegainv_chol)
  if(nu <= p-1){
    stop(
      sprintf(
        "`nu+1` (%s) must be greater than the dimension (%s)", nu+1, p))
  }
  out <- array(NA_real_, dim=c(p, p, n))
  for(i in 1:n){
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
    Z[lower.tri(Z)] <- rnorm(p*(p-1)/2)
    out[,,i] <- chol2inv(crossprod(Z, Omegainv_chol))
  }
  out
}
