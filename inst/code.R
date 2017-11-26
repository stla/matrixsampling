isScalar <- function(x){
  (is.complex(x) || is.numeric(x)) && (length(x) == 1)
}
isRealScalar <- function(x){
  is.numeric(x) && (length(x) == 1)
}
isPositiveInteger <- function(x){
  isRealScalar(x) && (floor(x) == x)
}
# not used yet:
sqrtm <- function(M){ # for symmetric M, returns R such that RR=M
  M.eig <- eigen(M, symmetric = TRUE)
  M.eig$vectors %*% (sqrt(M.eig$values) * t(M.eig$vectors))
}

is_squareMatrix <- function(M){
  isScalar(M) || (nrow(M) == ncol(M))
}
is_squareRealMatrix <- function(M){
  is_squareMatrix(M) && is.numeric(M)
}

checkSigma <- function(Sigma){
  if(isScalar(Sigma)){
    if(Sigma >= 0){
      return(TRUE)
    }else{
      stop("`Sigma` is not positive")
    }
  }
  p <- nrow(Sigma)
  if(p != ncol(Sigma)){
    stop("`Sigma` must be a symmetric positive matrix - it is not square")
  }
  if(!isTRUE(all.equal.numeric(Sigma, t(Sigma), tol=100*.Machine$double.eps))){
    stop("`Sigma` must be a symmetric positive matrix - it is not symmetric")
  }
  if(is.complex(Sigma)){
    stop("`Sigma` has complex entries")
  }
  Sigma_eig <- eigen(Sigma, symmetric = TRUE)
  if(any(is.complex(Sigma_eig$values)) || any(Sigma_eig$values < 0)){
    stop("`Sigma` is symmetric but not positive")
  }
  TRUE
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
  if(!isTRUE(all.equal.numeric(Sigma, t(Sigma), tol=100*.Machine$double.eps))){
    stop(sprintf("`%s` must be a symmetric positive matrix - it is not symmetric",
                 matrixname), call. = FALSE)
  }
  if(is.complex(Sigma)){
    stop(sprintf("`%s` has complex entries", matrixname), call. = FALSE)
  }
  Sigma_eig <- eigen(Sigma, symmetric = TRUE)
  if(any(is.complex(Sigma_eig$values)) || any(Sigma_eig$values < 0)){
    stop(sprintf("`%s` is symmetric but not positive", matrixname), call. = FALSE)
  }
  Sigma_eig$vectors %*% (sqrt(Sigma_eig$values) * t(Sigma_eig$vectors))
}

isZeroMatrix <- function(M){
  if(isScalar(M)){
    M <- as.matrix(M)
  }
  isTRUE(all.equal.numeric(M, matrix(0, nrow(M), ncol(M))))
}

#### Wishart ####
rwishart_chol <- function(n, nu, Sigma){
  # returns lower triangular Cholesky
  # nonsingular only
  Sigma_chol <- chol(Sigma)
  if(any(Sigma != t(Sigma))){
    stop("`Sigma` is not symmetric")
  }
  p <- nrow(Sigma_chol)
  if(nu < p){
    stop(sprintf("`nu` (%s) must be at least equal to the dimension (%s)", nu, p))
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
  # ! returns C such that CC'=Wishart
  # nu >= p only
  # Sigma_root Sigma_root' = Sigma
  if(is.null(Sigma_root)){
    p <- ifelse(isScalar(Sigma), 1L, nrow(Sigma))
  }else{
    p <- ifelse(isScalar(Sigma_root), 1L, nrow(Sigma_root))
  }
  if(check){
    if(nu < p){
      stop(sprintf("`nu` (%s) must be at least equal to the dimension (%s)", nu, p))
    }
    if(is.null(Sigma_root)){
      Sigma_root <- matrixroot(Sigma)
    }else{
      if(!is_squareRealMatrix(Sigma_root)){
        stop("`Sigma_root` must be a square real matrix")
      }
      checkSigma(Sigma_root %*% t(Sigma_root)) # inutile de check symétrique et positive - ça l'est !
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
  # ! returns C such that CC'=Wishart(nu, I_p)
  # nu >= p only
  out <- array(NA_real_, dim=c(p, p, n))
  for(i in 1:n){
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
    Z[lower.tri(Z)] <- rnorm(p*(p-1)/2)
    out[,,i] <- Z
  }
  out
}

rwishart <- function(n, nu, Sigma, Theta=NULL){
  Sigma_root <- matrixroot(Sigma)
  p <- nrow(Sigma_root)
  # check n ? check nu ?
  if(is.null(Theta) || isZeroMatrix(Theta)){
    if(nu >= p){
      R <- rwishart_root(n, nu, Sigma_root=Sigma_root, check=FALSE)
      out <- array(NA_real_, dim=dim(R))
      # out <- structure(apply(R, 3, dim=dim(R)))
      for(i in 1:n){
        out[,,i] <- tcrossprod(R[,,i])
      }
    }else{
      if(!isPositiveInteger(nu)){
        stop("`nu` must be an integer")
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
    if(nu < p){
      stop(
        sprintf("`nu` must be greater than the dimension (%s) in the noncentral case",
                p))
    }
    Theta_root <- matrixroot(Theta, matrixname = "Theta")
    W <- array(NA_real_, dim=c(p,p,n))
    if(nu >= 2*p){
      WrootI <- rwishart_root_I(n, nu-p, p)
      for(i in 1:n){
        Z <- matrix(rnorm(p*p), p, p)
        W[,,i] <- (Theta_root + Sigma_root %*% Z) %*%
          (Theta_root + t(Z) %*% Sigma_root) +
          tcrossprod(Sigma_root %*% WrootI[,,i])
      }
    }else if(nu>p){
      if(!isPositiveInteger(nu)){
        stop("`nu` must be an integer")
      }
      for(i in 1:n){
        Z <- matrix(rnorm(p*p), p, p)
        W[,,i] <- (Theta_root + Sigma_root %*% Z) %*%
          (Theta_root + t(Z) %*% Sigma_root) +
          Sigma_root %*%
          tcrossprod(matrix(rnorm((nu-p)*p), p, nu-p)) %*% t(Sigma_root)
      }
    }else{
      for(i in 1:n){
        W[,,i] <- tcrossprod(Theta_root + Sigma_root %*% matrix(rnorm(p*p), p, p))
      }
    }
    return(W)
  }
}
# rwishart2 <- function(n, nu, Sigma){
#   R <- rwishart_chol(n, nu, Sigma)
#   array(apply(R, 3, tcrossprod), dim=dim(R))
# }
# rwishart3 <- function(n, nu, Sigma){
#   R <- rwishart_chol(n, nu, Sigma)
#   plyr::aaply(R, 3, tcrossprod, .parallel = TRUE)
# }
# library(microbenchmark)
# nu <- 16
# p <- 13
# Sigma <- rwishart(1, p, diag(p))[,,1]
# nsims <- 1000
# microbenchmark( # f1 f2 kif-kif; f3 lent
#   f1 = rwishart(nsims, nu, Sigma),
#   f2 = rwishart2(nsims, nu, Sigma),
#   f3 = rwishart3(nsims, nu, Sigma)
# )

rwishart_I <- function(n, nu, p, Theta=NULL){
  # check n ? check nu ?
  if(is.null(Theta)){
    if(nu >= p){
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
    if(nu < p){
      stop(
        sprintf("`nu` must be greater than the dimension (%s) in the noncentral case",
                p))
    }
    Theta_root <- matrixroot(Theta, matrixname = "Theta")
    W <- array(NA_real_, dim=c(p,p,n))
    if(nu >= 2*p){
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
        W[,,i] <- tcrossprod(Theta_root + Sigma_root %*% matrix(rnorm(p*p), p, p))
      }
    }
    return(W)
  }
}

# check moments - ok singular
nu <- 6
p <- 3
S <- rwishart(1, p+20, diag(p))[,,1]/(p+20)
#S <- cbind(c(3,2,0),c(2,2,0), c(0,0,0))
nsims <- 150000
sims <- rwishart(nsims, nu, S)
# mean
apply(sims, 1:2, mean)
nu*S
# variances
apply(sims, c(1,2), var)
nu * (S^2 + tcrossprod(diag(S))) # nu * (S[i,j]^2 + S[i,i]*S[j,j])
# covariances
vecsims <- apply(sims, 3, c)
round(cov(t(vecsims))-2*nu*kronecker(S,S),1)

# check determinant
detsims <- apply(sims, 3, det)
chi2sims <- rchisq(nsims, nu)
for(i in seq_len(p-1)){
  chi2sims <- chi2sims*rchisq(nsims, nu-i)
}
curve(ecdf(detsims/det(S))(x), from=0, to=200)
curve(ecdf(chi2sims)(x), add=TRUE, col="red")

# check trace - ok for nu < p
tr <- function(x) sum(diag(x))
mean(apply(sims, 3, tr))
nu*tr(S) # ok singular
mean(apply(sims, 3, function(x) tr(solve(S)%*%x)*tr(x)))
nu*2*(nu*p/2+1)*tr(S)

# characteristic function - ok singular
# of the trace
z <- .04
mean(apply(sims, 3, function(x) exp(1i*z*tr(x))))
complexplus::Det(diag(p) - 2*1i*z*S)^(-nu/2)
#
Z <- matrix(0.04, p, p)
mean(apply(sims, 3, function(x) exp(tr(1i*Z%*%x))))
complexplus::Det(diag(p) - 2*1i*Z%*%S)^(-nu/2)

# check singular case
Sigma <- cbind(c(3,2,0),c(2,2,0), c(0,0,0))
nu <- 2
nsims <- 5000
sims <- rwishart(nsims, nu, Sigma)
r <- 2 # rank Sigma
sims2 <- numeric(nsims)
SigmaPlus <- MASS::ginv(Sigma)
for(i in 1:nsims){
  y <- rnorm(3)
  sims2[i] <- (t(y) %*% SigmaPlus %*% y) / (t(y) %*% MASS::ginv(sims[,,i]) %*% y)
}
curve(ecdf(sims2)(x), to=10)
curve(pchisq(x, nu-r+1), add=TRUE, col="red")


# check U-distribution
p <- 3
Sigma <- rwishart(1, p, diag(p))[,,1] # arbitray U matrix
nu1 <- 6 # degrees of freedom
nu2 <- 4
nsims <- 50000
U <- rwishart(nsims, nu1, Sigma)
V <- rwishart(nsims, nu2, Sigma)
sims <- apply(U, 3, det)/apply(U + V, 3, det)
betasims <- rbeta(nsims, nu1/2, nu2/2)
for(i in seq_len(p-1)){
  betasims <- betasims*rbeta(nsims, (nu1-i)/2, nu2/2)
}
curve(ecdf(sims)(x))
curve(ecdf(betasims)(x), add=TRUE, col="red")

#### inverse Wishart ####
rinvwishart <- function(n, nu, Sigma){ # only for Sigma>0 and nu>=p
  p <- ifelse(isScalar(Sigma), 1L, nrow(Sigma))
  if(nu < p){
    stop("xxx")
  }
  # if(p != ncol(Sigma)){
  #   stop("`Sigma` must be a symmetric positive matrix - it is not square")
  # }
  # if(any(Sigma != t(Sigma))){
  #   stop("`Sigma` must be a symmetric positive matrix - it is not symmetric")
  # }
  # Sigma_eig <- eigen(Sigma, symmetric = TRUE)
  # if(any(is.complex(Sigma_eig$values)) || any(Sigma_eig$values < 0)){
  #   stop("`Sigma` is not positive")
  # }
  Sigmainv_chol <- chol(chol2inv(chol(Sigma)))
  out <- array(NA_real_, dim=c(p, p, n))
  for(i in 1:n){
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, nu:(nu-p+1)))
    Z[lower.tri(Z)] <- rnorm(p*(p-1)/2)
    out[,,i] <- chol2inv(t(crossprod(Sigmainv_chol, Z)))
  }
  out
  # structure(
  #   apply(structure(apply(W, 3, chol), dim=dim(W)), 3, chol2inv),
  #   dim=dim(W))
}

nu <- 5
p <- 3
S <- rwishart(1, p, diag(p))[,,1]
nsims <- 100000
sims <- rinvwishart(nsims, nu, S)
apply(sims, c(1,2), mean)
S / (nu-p-1)

#### noncentral ####
# rwishartnc <- function(n, nu, Sigma, Theta){ # nu >= p+1
#   p <- ifelse(isScalar(Sigma), 1L, nrow(Sigma))
#   if(!(isScalar(Sigma) && isScalar(Theta)) && !identical(dim(Sigma), dim(Theta))){
#     stop("`Sigma` and `Theta` must have the same dimension")
#   }
#   Theta_root <- matrixroot(Theta) # pb: matrix name in error message
#   SigmaRoot <- matrixroot(Sigma) # t(chol(Sigma))
#   W1 <- W2 <- array(NA_real_, dim=c(p,p,n))
#   if(nu >= 2*p){
#     W <- rwishart_root(n, nu-p, Sigma_root=diag(p), check=FALSE)
#     # un peu con ça fait produit avec Sigma_root=diag(p)
#     for(i in 1:n){
#       Z <- matrix(rnorm(p*p), p, p)
#       W1[,,i] <- (Theta_root + SigmaRoot %*% Z) %*% (Theta_root + t(Z) %*% t(SigmaRoot))
#       W2[,,i] <- tcrossprod(SigmaRoot %*% W[,,i])
#     }
#   }else{
#     for(i in 1:n){
#       Z <- matrix(rnorm(p*p), p, p)
#       W1[,,i] <- (Theta_root + SigmaRoot %*% Z) %*% (Theta_root + t(Z) %*% t(SigmaRoot))
#       W <- tcrossprod(matrix(rnorm((nu-p)*p), p, nu-p))
#       W2[,,i] <- SigmaRoot %*% W %*% t(SigmaRoot)
#     }
#   }
#   W1 + W2
# }

p <- 3
nu <- 6.5
Sigma <- rwishart(1, p, diag(p))[,,1]
Theta <- rwishart(1, p, diag(p))[,,1]
nsims <- 50000
sims <- rwishart(nsims, nu, Sigma, Theta)
apply(sims, 1:2, mean)
Theta + nu*Sigma
# characteristic function - faux dans Gupta Nagar
# formule trouvée dans http://www.tandfonline.com/doi/abs/10.1080/00949658208810570
Z <- matrix(0.04, p, p)
mean(apply(sims, 3, function(x) exp(tr(1i*Z%*%x))))
complexplus::Det(diag(p) - 2*1i*Z%*%Sigma)^(-nu/2) *
  exp(1i*tr(solve(diag(p) - 2*1i*Z%*%Sigma) %*% Z %*% Theta))
#  exp(tr(-Theta/2 + 1/2*solve(diag(p) - 2*1i*Z%*%Sigma) %*% Theta))

#### Beta matrix ####
invsqrtm <- function(M){ # for symmetric M
  M.eig <- eigen(M, symmetric = TRUE)
  M.eig$vectors %*% (1/sqrt(M.eig$values) * t(M.eig$vectors))
}
rmatrixbeta <- function(n, p, a, b, Theta1=NULL, Theta2=NULL){
  # check n, p, a, b
  # remplacer is.null par is.null or isZeroMatrix
  if(is.null(Theta1) && is.null(Theta2)){
    if(2*a+2*b < p){
      stop("`a` and `b` must satisfy `a+b >= p/2`")
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(!is.null(Theta1) && is.null(Theta2)){
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(is.null(Theta1) && !is.null(Theta2)){
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
rmatrixbetaII <- function(n, p, a, b, Theta1=NULL, Theta2=NULL){
  # check n, p, a, b
  if(is.null(Theta1) && is.null(Theta2)){
    if(2*b < p){
      stop("`b` must satisfy `b >= p/2`") # et si j'utilise la transfo de BetaI à BetaII ?
    }
    W1 <- rwishart_I(n, 2*a, p)
    W2 <- rwishart_I(n, 2*b, p) # ou rwishart_root_I et chol2inv(chol(W2)) ?
  }else if(!is.null(Theta1) && is.null(Theta2)){
    W1 <- rwishart_I(n, 2*a, p, Theta1)
    W2 <- rwishart_I(n, 2*b, p)
  }else if(is.null(Theta1) && !is.null(Theta2)){
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

n1 <- 16
n2 <- 3
p <- 3
Theta <- rwishart(1, p, diag(p))[,,1]
nsims <- 20000
sims <- rmatrixbeta(nsims, p, n1, n2, Theta2=Theta)
detsims <- apply(sims, 3, det)

mgamma <- function(p, a){ # multivariate Gamma
  exp(log(pi) * (p * (p - 1)/4) + sum(lgamma(a + (1-(1:p))/2)))
}
h1F1 <- function(alpha, beta, R){
  p <- nrow(R)
  sims <- apply(rmatrixbeta(20000, p, alpha, beta-alpha), 3,
                function(x) exp(tr(R%*%x)))
  mean(sims) # est-ce bien "2*" dans Gupta-Nagar ?
  # non, donc il y a une erreur dans E[det(Beta non central)]
  # non c'est bon, G-N considère Beta(n1/2, n2/2, ncp)
}
h <- 1
mgamma(p,n1+h)*mgamma(p,n1+n2)/mgamma(p,n1+n2+h)/mgamma(p,n1) *
  exp(tr(-Theta/2)) *
  h1F1(n1+n2, n1+n2+h, Theta/2)
mean(detsims^h)

# check p=1
n1 <- 16
n2 <- 3
p <- 1
Theta <- 1
nsims <- 20000
sims <- rmatrixbeta(nsims, p, n1, n2, Theta1=Theta)[1,1,]
curve(ecdf(sims)(x))
curve(pbeta(x, n1, n2, ncp=Theta), add=TRUE, col="red")
#
sims <- rmatrixbeta(nsims, p, n1, n2, Theta2=Theta)[1,1,]
chi1sims <- rchisq(nsims, 2*n1)
chi2sims <- rchisq(nsims, 2*n2, ncp=Theta)
curve(ecdf(sims)(x))
curve(ecdf(chi1sims/(chi1sims+chi2sims))(x), add=TRUE, col="green")
curve(1-pbeta(1-x, n2, n1, ncp=Theta), add=TRUE, col="blue")


#### matrix normal ####
rmatrixnormal <- function(n, M, U, V){
  if(!isPositiveInteger(n)){
    stop("`n` must be a positive integer")
  }
  tcholU <- t(chol(U)) # prendre matrixroot pour autoriser semidéfinie positive
  cholV <- chol(V)
  m <- nrow(tcholU)
  p <- nrow(cholV)
  if(isScalar(M)){
    M <- as.matrix(M)
  }
  if(m != nrow(M) || p != ncol(M)){
    stop("Incorrect dimensions")
  }
  if(any(U != t(U))){
    stop("`U` is not symmetric")
  }
  if(any(V != t(V))){
    stop("`V` is not symmetric")
  }
  out <- array(NA_real_, dim=c(m,p,n))
  for(i in 1:n){
    out[,,i] <- M + tcholU %*% matrix(rnorm(m*p), m, p) %*% cholV
  }
  out
}

#### matrix Student ####
rmatrixt <- function(n, nu, M, U, V){
  if(!isRealScalar(nu)){
    stop("`nu` must be a positive number")
  }
  cholV <- chol(V) # matrixroot ?
  if(any(V != t(V))){
    stop("`V` is not symmetric")
  }
  m <- ifelse(isScalar(U), 1L, nrow(U))
  p <- nrow(cholV)
  if(isScalar(M)){
    M <- as.matrix(M)
  }
  if(m != nrow(M) || p != ncol(M)){
    stop("Incorrect dimensions")
  }
  IW <- rinvwishart(n, nu+m-1, U)
  out <- array(NA_real_, dim=c(m,p,n))
  for(i in 1:n){
    out[,,i] <- M + crossprod(chol(IW[,,i]), matrix(rnorm(m*p), m, p) %*% chol(V))
  }
  out
}

#### checks matrix t ####
m <- 3; p <- 5
M <- matrix(1, m, p)
U <- rwishart(1, m, diag(m))[,,1] # arbitray U matrix
V <- rwishart(1, p, diag(p))[,,1] # arbitrary V matrix
Uinv <- solve(U)
Vinv <- solve(V) # inverse V
nu <- 6 # degrees of freedom

# ratio of determinants #
nsims <- 20000
sims <- rmatrixt(nsims, nu, M, U, V)
X <- sweep(sims, 1:2, M, "-")
detsims <- apply(X, 3,
                 function(x) 1/det(diag(m) + Uinv %*% x %*% Vinv %*% t(x)))
# simulations of product of m Beta distributions
k <- (nu + (m-1))/2
sims2 <- rbeta(nsims, k, p/2)
for(i in seq_len(m-1)){
  sims2 <- sims2*rbeta(nsims, k-i/2, p/2)
}
k <- (nu + p - 1)
sims3 <- rbeta(nsims, k/2, m/2)
for(i in seq_len(p-1)){
  sims3 <- sims3 * rbeta(nsims, (k-i)/2, m/2)
}
# compare cumulative distribution functions
curve(ecdf(detsims)(x), from=0, to=1)
curve(ecdf(sims2)(x), col="red", add=TRUE, lty=2, lwd=3)
curve(ecdf(sims3)(x), col="green", add=TRUE, lty=2, lwd=3)

## covariance matrix
Cov <- cov(t(apply(sims, 3, function(x) c(t(x)))))
k <- nu-2
round((1/k * kronecker(U,V))[1:3,1:3],2)
round(Cov[1:3,1:3], 2)

## compare with matrix beta
Uroot <- sqrtm(U)
bsims <- array(apply(X, 3,
                     function(x) Uroot %*% solve(U + x %*% Vinv %*% t(x)) %*% Uroot),
               dim = c(m,m,nsims))
bsims2 <- rmatrixbeta(nsims, m, (nu + m - 1)/2, p/2)
Z <- matrix(0.04, m, m)
mean(apply(bsims, 3, function(x) exp(tr(1i*Z%*%x))))
mean(apply(bsims2, 3, function(x) exp(tr(1i*Z%*%x))))
## compare with matrix beta II
Urootinv <- invsqrtm(U)
bsims <- array(apply(X, 3,
                     function(x) Urootinv %*% x %*% Vinv %*% t(x) %*% Urootinv),
               dim = c(m,m,nsims))
bsims2 <- rmatrixbetaII(nsims, m, p/2, (nu+m-1)/2)
mean(apply(bsims, 3, function(x) exp(tr(1i*Z%*%x))))
mean(apply(bsims2, 3, function(x) exp(tr(1i*Z%*%x))))


#### inverted matrix Student ####
rmatrixit <- function(n, nu, M, U, V){
  if(!isRealScalar(nu)){
    stop("`nu` must be a positive number")
  }
  Uroot <- matrixroot(U)
  Vroot <- matrixroot(V)
  m <- nrow(Uroot)
  p <- nrow(Vroot)
  if(isScalar(M)){
    M <- as.matrix(M)
  }
  if(m != nrow(M) || p != ncol(M)){
    stop("Incorrect dimensions")
  }
  W <- rwishart_I(n, nu+m-1, m)
  out <- array(NA_real_, dim=c(m,p,n))
  for(i in 1:n){
    Z <- matrix(rnorm(m*p), m, p)
    out[,,i] <- M + Uroot %*%
      forwardsolve(t(chol(W[,,i] + Z%*%t(Z))), diag(m)) %*% Z %*% Vroot
  }
  out
}

#### checks inverted t ####
# thm 5.2.3
nsims <- 50000
p <- 3
n1 <- 5
n2 <- 4
Sigma <- rwishart(1, p, diag(p))[,,1]
X <- rmatrixnormal(nsims, matrix(0, p, n1), Sigma, diag(n1))
S2 <- rwishart(nsims, n2, Sigma)
Z <- array(NA_real_, dim=c(p,n1,nsims))
for(i in 1:nsims){
  Z[,,i] <- invsqrtm(S2[,,i]+tcrossprod(X[,,i]))%*%X[,,i]
}
IT <- rmatrixit(nsims, n2-p+1, matrix(0, p, n1), diag(p), diag(n1))
z <- 0.04*diag(p) #matrix(0.04, p, p)
mean(apply(Z, 3, function(x) exp(tr(1i*z%*%x))))
mean(apply(IT, 3, function(x) exp(tr(1i*z%*%x))))
apply(Z, 1:2, var)
apply(IT, 1:2, var)
curve(ecdf(Z[1,1,]+Z[1,2,])(x), from=-1.5, to=1.5)
curve(ecdf(IT[1,1,]+IT[1,2,])(x), add=TRUE, col="red")

