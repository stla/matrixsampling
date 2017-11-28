isScalar <- function(x){
  (is.complex(x) || is.numeric(x)) && (length(x) == 1)
}

isRealScalar <- function(x){
  is.numeric(x) && (length(x) == 1)
}

isPositiveInteger <- function(x){
  isRealScalar(x) && (floor(x) == x)
}

isSymmetricMatrix <- function(Sigma){
  isTRUE(all.equal.numeric(Sigma, t(Sigma), tolerance=100*.Machine$double.eps))
}

matrixroot <- function(Sigma, matrixname="Sigma"){
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
  if(!isSymmetricMatrix(Sigma)){
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
