isScalar <- function(x){
  (is.complex(x) || is.numeric(x)) && (length(x) == 1L)
}

isRealScalar <- function(x){
  is.numeric(x) && (length(x) == 1L)
}

isPositiveInteger <- function(x){
  isRealScalar(x) && (floor(x) == x)
}

isSymmetricMatrix <- function(Sigma){
  isTRUE(all.equal.numeric(Sigma, t(Sigma), tolerance=100*.Machine$double.eps))
}

isSquareRealMatrix <- function(M){
  (isScalar(M) || (is.matrix(M) && (nrow(M) == ncol(M)))) && is.numeric(M)
}

isSymmetricPositive <- function(Sigma, matrixname="Sigma", symmetric=FALSE){
  if(!symmetric && !isSymmetricMatrix(Sigma)){
    stop(sprintf("`%s` must be a symmetric positive matrix - it is not symmetric",
                 matrixname), call. = FALSE)
  }
  Sigma_eig <- eigen(Sigma, symmetric = TRUE)
  if(any(Sigma_eig$values < 0)){
    stop(sprintf("`%s` is symmetric but not positive", matrixname), call. = FALSE)
  }
  TRUE
}

matrixroot <- function(Sigma, symmetric=FALSE, matrixname="Sigma"){
  if(isScalar(Sigma)){
    if(Sigma >= 0){
      return(as.matrix(sqrt(Sigma)))
    }else{
      stop(sprintf("`%s` is not positive", matrixname), call. = FALSE)
    }
  }
  p <- nrow(Sigma)
  if(p != ncol(Sigma)){
    stop(sprintf("`%s` must be a symmetric positive matrix - it is not square",
                 matrixname), call. = FALSE)
  }
  if(!symmetric && !isSymmetricMatrix(Sigma)){
    stop(sprintf("`%s` must be a symmetric positive matrix - it is not symmetric",
                 matrixname), call. = FALSE)
  }
  if(is.complex(Sigma)){
    stop(sprintf("`%s` has complex entries", matrixname), call. = FALSE)
  }
  Sigma_eig <- eigen(Sigma, symmetric = TRUE)
  if(any(Sigma_eig$values < 0)){
    stop(sprintf("`%s` is symmetric but not positive", matrixname), call. = FALSE)
  }
  Sigma_eig$vectors %*% (sqrt(Sigma_eig$values) * t(Sigma_eig$vectors))
}

sqrtm <- function(Sigma){ # square root for positive symmetric Sigma
  Sigma_eig <- eigen(Sigma, symmetric = TRUE)
  Sigma_eig$vectors %*% (sqrt(Sigma_eig$values) * t(Sigma_eig$vectors))
}

invsqrtm <- function(Sigma){ # inverse square root for positive symmetric Sigma
  Sigma_eig <- eigen(Sigma, symmetric = TRUE)
  Sigma_eig$vectors %*% (1/sqrt(Sigma_eig$values) * t(Sigma_eig$vectors))
}

isZeroMatrix <- function(M){
  if(isScalar(M)){
    M <- as.matrix(M)
  }
  is.matrix(M) && isTRUE(all.equal.numeric(M, matrix(0, nrow(M), ncol(M))))
}

isNullOrZeroMatrix <- function(M){
  is.null(M) || isZeroMatrix(M)
}

extendedCholesky <- function(S){ # does not check S >= 0
  C <- suppressWarnings(chol(S, pivot=TRUE))
  d <- nrow(C)
  P <- matrix(0, d, d)
  P[cbind(1L:d, attr(C,"pivot"))] <- 1
  r <- attr(C, "rank")
  return(list(L = t(C[seq_len(r), seq_len(r), drop=FALSE]),
              Ctilde = cbind(t(C[seq_len(r), , drop=FALSE]),
                             rbind(matrix(0, r, d-r), diag(d-r))),
              P = P))
}
