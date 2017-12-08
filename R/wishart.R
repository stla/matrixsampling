#' Sampling Cholesky factor of a Wishart matrix
#'
#' Samples the lower triangular Cholesky factor of a Wishart random matrix.
#'
#' @param n sample size, a positive integer
#' @param nu degrees of freedom, a number strictly greater than \code{p-1},
#' where \code{p} is the dimension (the order of \code{Sigma})
#' @param Sigma scale matrix, a positive definite real matrix
#' @param epsilon a number involved in the algorithm only if it positive; its role
#' is to guarantee the invertibility of the sampled matrices; see Details
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm rchisq
#'
#' @details The argument \code{epsilon} is a threshold whose role is to guarantee
#' that the algorithm samples invertible matrices.
#' The matrices sampled by the algorithm are theoretically invertible.
#' However, because of numerical precision, they are not always invertible when
#' \code{nu} is close to \code{p-1}, i.e. when \code{nu-p+1} is small. In this case,
#' the simulations of chi-squared distributions involved in the algorithm can
#' generate zero values or values close to zero, yielding the non-invertibility
#' of the sampled matrices. These values are replaced with \code{epsilon} if they are
#' smaller than \code{epsilon}.
#'
#' @examples
#' nu <- 4
#' p <- 3
#' Sigma <- diag(p)
#' Wsims <- rwishart_chol(10000, nu, Sigma)
#' dim(Wsims) # 3 3 10000
#' Wsims[,,1]
rwishart_chol <- function(n, nu, Sigma, epsilon=0){
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
    stop(sprintf("`nu` (%s) must be greater than `p-1` (%s)", nu, p-1))
  }
  out <- array(NA_real_, dim=c(p, p, n))
  chisims <- matrix(NA_real_, p, n)
  for(k in 1:p){
    chisims[k,] <- sqrt(rchisq(n, nu-k+1))
  }
  if(epsilon>0){
    chisims <- pmax(chisims, epsilon)
  }
  Y <- matrix(rnorm(p*(p-1)/2*n), p*(p-1)/2, n)
  Z <- matrix(0, p, p)
  lowertri <- lower.tri(Z)
  for(i in 1:n){
    diag(Z) <- chisims[,i]
    Z[lowertri] <- Y[,i]
    out[,,i] <- crossprod(Sigma_chol, Z)
  }
  out
}

rwishart_root <- function(n, nu, Sigma, Sigma_root=NULL, check=TRUE,
                          epsilon=0){
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
  chisims <- matrix(NA_real_, p, n)
  for(k in 1:p){
    chisims[k,] <- sqrt(rchisq(n, nu-k+1))
  }
  if(epsilon>0){
    chisims <- pmax(chisims, epsilon)
  }
  out <- array(NA_real_, dim=c(p, p, n))
  Y <- matrix(rnorm(p*(p-1)/2*n), p*(p-1)/2, n)
  Z <- matrix(0, p, p)
  lowertri <- lower.tri(Z)
  for(i in 1:n){
    diag(Z) <- chisims[,i]
    Z[lowertri] <- Y[,i]
    out[,,i] <- crossprod(Sigma_root, Z)
  }
  out
}

rwishart_chol_I <- function(n, nu, p, epsilon=0){
  # samples R such that RR' ~ W(nu, I_p)
  # nu > p-1 only
  chisims <- matrix(NA_real_, p, n)
  for(k in 1:p){
    chisims[k,] <- sqrt(rchisq(n, nu-k+1))
  }
  if(epsilon>0){
    chisims <- pmax(chisims, epsilon)
  }
  out <- array(NA_real_, dim=c(p, p, n))
  Y <- matrix(rnorm(p*(p-1)/2*n), p*(p-1)/2, n)
  Z <- matrix(0, p, p)
  lowertri <- lower.tri(Z)
  for(i in 1:n){
    diag(Z) <- chisims[,i]
    Z[lowertri] <- Y[,i]
    out[,,i] <- Z
  }
  out
}

rwishart_AA_Im <- function(n, nu, m, Theta, epsilon=0){
  d <- nrow(Theta)
  rChi2 <- if(epsilon>0 && m==d && nu>d-1){
    function(n, df, ncp){
      pmax(rchisq(n, df, ncp), epsilon)
    }
  }else{
    rchisq
  }
  ec <- extendedCholesky(Theta[-1L,-1L])
  L <- ec$L; Ctilde <- ec$Ctilde; P <- ec$P
  r <- nrow(L)
  Pi <- cbind(c(1,rep(0,d-1L)), rbind(0,P))
  xtilde <- Pi %*% Theta %*% t(Pi)
  u1 <- c(forwardsolve(L, diag(r)) %*% xtilde[1L, 2L:(r+1L)])
  U11 <- rChi2(n, df=nu-r, ncp=max(0, xtilde[1L,1L]-sum(u1^2)))
  U <- sweep(matrix(rnorm(r*n),r,n), 1L, u1, "+")
  B <- t(Pi) %*% cbind(c(1,rep(0,d-1L)), rbind(0, Ctilde))
  Wsims <- array(NA_real_, dim=c(d, d, n))
  Y <- matrix(0, d, d)
  for(i in 1:n){
    Y[1L:(r+1L), 1L:(r+1L)] <- cbind(c(U11[i] + sum(U[,i]^2), U[,i]),
                                     rbind(U[,i], diag(r)))
    Wsims[,,i] <- B %*% Y %*% t(B)
  }
  for(k in (1L+seq_len(m-1L))){
    p <- 1L:d; p[k] <- 1L; p[1L] <- k
    Wsims <- Wsims[p,p,,drop=FALSE] # don't drop if n=1
    for(i in 1L:n){
      ec <- extendedCholesky(Wsims[-1L,-1L,i])
      L <- ec$L; Ctilde <- ec$Ctilde; P <- ec$P
      r <- nrow(L)
      Pi <- cbind(c(1,rep(0,d-1L)), rbind(0,P))
      xtilde <- Pi %*% Wsims[,,i] %*% t(Pi)
      u1 <- c(forwardsolve(L, diag(r)) %*% xtilde[1L, 2L:(r+1L)])
      U11 <- rChi2(1L, df=nu-r, ncp=max(0, xtilde[1L,1L]-sum(u1^2)))
      U <- rnorm(r) + c(forwardsolve(L, diag(r)) %*% xtilde[1L, 2L:(r+1L)])
      B <- t(Pi) %*% cbind(c(1,rep(0,d-1L)), rbind(0, Ctilde))
      Y <- matrix(0, d, d)
      Y[1L:(r+1L), 1L:(r+1L)] <- cbind(c(U11 + sum(U^2), U), rbind(U, diag(r)))
      Wsims[,,i] <- (B %*% Y %*% t(B))[p,p]
    }
  }
  Wsims
}

rwishart_AA <- function(n, nu, Sigma, Theta, epsilon=0){
  ec <- extendedCholesky(Sigma)
  L <- ec$L; Ctilde <- ec$Ctilde; P <- ec$P
  theta <- t(P) %*% Ctilde
  thetainv <- forwardsolve(Ctilde, diag(nrow(P))) %*% P
  Y <- rwishart_AA_Im(n, nu, nrow(L), thetainv%*%Theta%*%t(thetainv), epsilon)
  array(apply(Y, 3L, function(x) theta%*%x%*%t(theta)), dim=dim(Y))
}

#' Wishart sampler
#'
#' Samples a Wishart distribution.
#'
#' @param n sample size, a positive integer
#' @param nu degrees of freedom, a positive number;
#' if \code{nu < p-1} where \code{p} is the dimension (the order of \code{Sigma}),
#' must be an integer;
#' in the noncentral case (i.e. when \code{Theta} is not the null matrix), \code{nu}
#' must satisfy \code{nu >= p-1}
#' @param Sigma scale matrix, a positive semidefinite real matrix
#' @param Theta noncentrality parameter, a positive semidefinite real matrix of
#' same order as \code{Sigma}; setting it to \code{NULL} (default) is
#' equivalent to setting it to the zero matrix
#' @param epsilon a number involved in the algorithm only if it positive; its role
#' is to guarantee the invertibility of the sampled matrices; see Details
#' @param checkSymmetry logical, whether to check the symmetry of \code{Sigma}
#' and \code{Theta}
#'
#' @note A sampled Wishart matrix is always positive semidefinite.
#' It is positive definite if \code{nu > p-1} and \code{Sigma} is positive
#' definite, in theory (see Details).
#'
#' In the noncentral case, i.e. when \code{Theta} is not null, the Ahdida & Alfonsi
#' algorithm is used if \code{nu} is not an integer and \code{p-1 < nu < 2p-1}, or
#' if \code{nu = p-1}. The simulations are slower in this case.
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension (see example).
#' @export
#' @importFrom stats rnorm rchisq
#'
#' @details The argument \code{epsilon} is a threshold whose role is to guarantee
#' that the algorithm samples invertible matrices when \code{nu > p-1} and
#' \code{Sigma} is positive definite.
#' They are theoretically invertible under this condition, but due to
#' numerical precision, they are not always invertible when
#' \code{nu} is close to \code{p-1}, i.e. when \code{nu-p+1} is small. In this case,
#' the chi-squared distributions involved in the algorithm can
#' generate zero values or values close to zero, yielding the non-invertibility
#' of the sampled matrices. These values are replaced with \code{epsilon} if they are
#' smaller than \code{epsilon}.
#'
#' @references A. Ahdida & A. Alfonsi. Exact and high-order discretization schemes
#' for Wishart processes and their affine extensions.
#' \emph{The Annals of Applied Probability} \strong{23}, 2013, 1025-1073.
#'
#' @examples
#' nu <- 4
#' p <- 3
#' Sigma <- toeplitz(p:1)
#' Theta <- diag(p)
#' Wsims <- rwishart(10000, nu, Sigma, Theta)
#' dim(Wsims) # 3 3 10000
#' apply(Wsims, 1:2, mean) # approximately nu*Sigma+Theta
#' # the epsilon argument:
#' Wsims_det <- apply(rwishart(10000, nu=p-1+0.001, Sigma), 3, det)
#' length(which(Wsims_det < .Machine$double.eps))
#' Wsims_det <- apply(rwishart(10000, nu=p-1+0.001, Sigma, epsilon=1e-8), 3, det)
#' length(which(Wsims_det < .Machine$double.eps))
rwishart <- function(n, nu, Sigma, Theta=NULL, epsilon=0, checkSymmetry=TRUE){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isRealScalar(nu)){
    stop("`nu` must be a positive number")
  }
  p <- ifelse(isScalar(Sigma), 1L, nrow(Sigma))
  if(isNullOrZeroMatrix(Theta)){
    Sigma_root <- matrixroot(Sigma, symmetric=!checkSymmetry)
    if(nu > p-1){
      R <- rwishart_root(n, nu, Sigma_root=Sigma_root, check=FALSE, epsilon=epsilon)
      out <- array(NA_real_, dim=dim(R))
      for(i in 1:n){
        out[,,i] <- tcrossprod(R[,,i])
      }
    }else{
      if(!isPositiveInteger(nu)){
        stop(sprintf("`nu` (%s) must be an integer when it is lower than `p-1` (%s)", nu, p-1))
      }
      out <- array(NA_real_, dim=c(p,p,n))
      Z <- array(rnorm(p*nu*n), dim=c(p,nu,n))
      for(i in 1:n){
        out[,,i] <- tcrossprod(Sigma_root %*% Z[,,i])
      }
    }
    return(out)
  }else{
    if(!(isScalar(Sigma) && isScalar(Theta)) && !identical(dim(Sigma), dim(Theta))){
      stop("`Sigma` and `Theta` must have the same dimension")
    }
    if(p == 1){
      W <- array(NA_real_, dim=c(p,p,n))
      W[1,1,] <- c(Sigma)*rchisq(n, df=nu, ncp=c(Theta))
      return(W)
    }
    if(nu < p-1){
      stop(
        sprintf(
          "`nu` (%s) must be greater than `p-1` (%s) in the noncentral case", nu, p-1))
    }
    W <- array(NA_real_, dim=c(p,p,n))
    if(nu > 2*p-1){
      Sigma_root <- matrixroot(Sigma, symmetric=!checkSymmetry)
      Theta_root <- matrixroot(Theta, matrixname="Theta", symmetric=!checkSymmetry)
      WrootI <- rwishart_chol_I(n, nu-p, p, epsilon)
      Z <- array(rnorm(p*p*n), dim=c(p,p,n))
      for(i in 1:n){
        W[,,i] <- (Theta_root + Sigma_root %*% Z[,,i]) %*%
          (Theta_root + t(Z[,,i]) %*% Sigma_root) +
          tcrossprod(Sigma_root %*% WrootI[,,i])
      }
    }else if(floor(nu) == nu && nu != p-1){ # nu is an integer >= p
      Sigma_root <- matrixroot(Sigma, symmetric=!checkSymmetry)
      Theta_root <- matrixroot(Theta, matrixname="Theta", symmetric=!checkSymmetry)
      if(nu != p){
        Z <- array(rnorm(p*p*n), dim=c(p,p,n))
        Y <- array(rnorm(p*(nu-p)*n), dim=c(p,nu-p,n))
        for(i in 1:n){
          W[,,i] <- (Theta_root + Sigma_root %*% Z[,,i]) %*%
            (Theta_root + t(Z[,,i]) %*% Sigma_root) +
            Sigma_root %*%
            tcrossprod(Y[,,i]) %*% t(Sigma_root)
        }
      }else{ # nu=p
        for(i in 1:n){
          W[,,i] <- tcrossprod(Theta_root + Sigma_root %*% matrix(rnorm(p*p), p, p))
        }
      }
    }else{ # nu is not an integer or nu = p-1
      invisible(isSymmetricPositive(Sigma, symmetric=!checkSymmetry) &&
                  isSymmetricPositive(Theta, matrixname="Theta", symmetric=!checkSymmetry))
      W <- rwishart_AA(n, nu, Sigma, Theta, epsilon)
    }
    return(W)
  }
}

rwishart_I <- function(n, nu, p, Theta=NULL, epsilon=0, checkSymmetry=TRUE){
  # samples W(nu, Ip, Theta)
  if(is.null(Theta)){
    if(nu > p-1){
      chisims <- matrix(NA_real_, p, n)
      for(k in 1:p){
        chisims[k,] <- sqrt(rchisq(n, nu-k+1))
      }
      if(epsilon>0){
        chisims <- pmax(chisims, epsilon)
      }
      out <- array(NA_real_, dim=c(p,p,n))
      Z <- matrix(0, p, p)
      lowertri <- lower.tri(Z)
      Y <- matrix(rnorm(p*(p-1)/2*n), p*(p-1)/2, n)
      for(i in 1:n){
        diag(Z) <- chisims[,i]
        Z[lowertri] <- Y[,i]
        out[,,i] <- tcrossprod(Z)
      }
    }else{
      if(!isPositiveInteger(nu)){
        stop("`nu` must be an integer")
      }
      out <- array(NA_real_, dim=c(p,p,n))
      Y <- array(rnorm(p*nu*n), dim=c(p,nu,n))
      for(i in 1:n){
        out[,,i] <- tcrossprod(Y[,,i])
      }
    }
    return(out)
  }else{
    if(!(p==1 && isScalar(Theta)) && !all(c(p,p)==dim(Theta))){
      stop("`Theta` must have dimension `p x p`")
    }
    if(p == 1){
      W <- array(NA_real_, dim=c(1,1,n))
      W[1,1,] <- rchisq(n, df=nu, ncp=c(Theta))
      return(W)
    }
    if(nu < p-1){
      stop(
        sprintf(
          "`nu` (%s) must be greater than `p-1` (%s) in the noncentral case", nu, p-1))
    }
    W <- array(NA_real_, dim=c(p,p,n))
    if(nu > 2*p-1){
      Theta_root <- matrixroot(Theta, matrixname="Theta", symmetric=!checkSymmetry)
      WrootI <- rwishart_chol_I(n, nu-p, p, epsilon)
      Z <- array(rnorm(p*p*n), dim=c(p,p,n))
      for(i in 1:n){
        W[,,i] <- (Theta_root + Z[,,i]) %*% (Theta_root + t(Z[,,i])) +
          tcrossprod(WrootI[,,i])
      }
    }else if(floor(nu) == nu && nu != p-1){ # nu is an integer >= p
      Theta_root <- matrixroot(Theta, matrixname="Theta", symmetric=!checkSymmetry)
      if(nu != p){
        Z <- array(rnorm(p*p*n), dim=c(p,p,n))
        Y <- array(rnorm(p*(nu-p)*n), dim=c(p,nu-p,n))
        for(i in 1:n){
          W[,,i] <- (Theta_root + Z[,,i]) %*% (Theta_root + t(Z[,,i])) +
            tcrossprod(Y[,,i])
        }
      }else{ # nu=p
        Z <- array(rnorm(p*p*n), dim=c(p,p,n))
        for(i in 1:n){
          W[,,i] <- tcrossprod(Theta_root + Z[,,i])
        }
      }
    }else{ # nu is not an integer or nu = p-1
      invisible(isSymmetricPositive(Theta, matrixname="Theta", symmetric=!checkSymmetry))
      W <- rwishart_AA_Im(n, nu, p, Theta, epsilon)
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
#' @param epsilon threshold to force invertibility in the algorithm; see Details
#' @param checkSymmetry logical, whether to check the symmetry of \code{Omega};
#' if \code{FALSE}, only the upper triangular part of \code{Omega} is used
#'
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
#'     The argument \code{epsilon} is a threshold whose role is to guarantee
#' the invertibility of the sampled Wishart distributions.
#' See Details in \code{\link{rwishart}}.
#'
#' @examples
#' nu <- 6
#' p <- 3
#' Omega <- toeplitz(p:1)
#' IWsims <- rinvwishart(10000, nu, Omega)
#' dim(IWsims) # 3 3 10000
#' apply(IWsims, 1:2, mean) # approximately Omega/(nu-p-1)
#' # the epsilon argument:
#' IWsims <- tryCatch(rinvwishart(10000, nu=p+0.001, Omega),
#'                    error=function(e) e)
#' IWsims <- tryCatch(rinvwishart(10000, nu=p+0.001, Omega, epsilon=1e-8),
#'                    error=function(e) e)
rinvwishart <- function(n, nu, Omega, epsilon=0, checkSymmetry=TRUE){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  if(!isRealScalar(nu)){
    stop("`nu` must be a positive number")
  }
  if(checkSymmetry && !isSymmetricMatrix(Omega)){
    stop("`Omega` is not symmetric")
  }
  Omegainv_chol <- chol(chol2inv(chol(Omega)))
  p <- nrow(Omegainv_chol)
  if(nu <= p-1){
    stop(
      sprintf(
        "`nu+1` (%s) must be greater than the dimension (%s)", nu+1, p))
  }
  chisims <- matrix(NA_real_, p, n)
  for(k in 1:p){
    chisims[k,] <- sqrt(rchisq(n, nu-k+1))
  }
  if(epsilon>0){
    chisims <- pmax(chisims, epsilon)
  }
  out <- array(NA_real_, dim=c(p, p, n))
  Y <- matrix(rnorm(p*(p-1)/2*n), p*(p-1)/2, n)
  Z <- matrix(0, p, p)
  lowertri <- lower.tri(Z)
  for(i in 1:n){
    diag(Z) <- chisims[,i]
    Z[lowertri] <- Y[,i]
    out[,,i] <- chol2inv(crossprod(Z, Omegainv_chol))
  }
  out
}
