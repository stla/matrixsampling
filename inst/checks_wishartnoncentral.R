
# characteristic function
tr <- function(x) sum(diag(x))
phi <- function(Z, nu, Sigma, Theta){
  p <- nrow(Sigma)
  complexplus::Det(diag(p) - 2*1i*Z%*%Sigma)^(-nu/2) *
    exp(1i*tr(solve(diag(p) - 2*1i*Z%*%Sigma) %*% Z %*% Theta))
}

# THM 3.5.1 - error: Theta = MM'
  # => 3.5.2 also wrong
p <- 2
n <- 3
Sigma <- diag(c(2,1))  #rwishart(1,n,diag(p))[,,1] #toeplitz(p:1)
M <- matrix(rpois(p*n,2),p,n)
nsims <- 50000
Gsims <- rmatrixnormal(nsims, M, Sigma, diag(n))
Wsims <- array(apply(Gsims, 3, tcrossprod), dim=c(p,p,nsims))
# Theta <- solve(Sigma)%*%M%*%t(M)
# Theta <- t(chol(solve(Sigma))) %*% M %*% t(M) %*% chol(solve(Sigma))
# Theta <- matrixsampling:::invsqrtm(Sigma) %*% M %*% t(M) %*% matrixsampling:::invsqrtm(Sigma)
# Theta <- apply(Wsims, 1:2, mean) - n*Sigma
# Theta <- matrixsampling:::matrixroot(Sigma) %*% M %*% t(M) %*% matrixsampling:::matrixroot(Sigma)
Theta <- M%*%t(M)
S <- 2*diag(p) + matrix(0.1, p, p)
Z <- 1/2i*(diag(p)-S) %*% solve(Sigma)
phi(Z, n, Sigma, Theta)
mean(apply(Wsims, 3, function(x) exp(tr(1i*Z%*%x))))
round(apply(Wsims, 1:2, mean), 3)
n*Sigma + Theta

