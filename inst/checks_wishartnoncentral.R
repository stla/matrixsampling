
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
Theta <- M%*%t(M)
S <- 2*diag(p) + matrix(0.1, p, p)
Z <- 1/2i*(diag(p)-S) %*% solve(Sigma)
phi(Z, n, Sigma, Theta)
mean(apply(Wsims, 3, function(x) exp(tr(1i*Z%*%x))))
round(apply(Wsims, 1:2, mean), 3)
n*Sigma + Theta

## invertible when Theta is not ? yes
p <- 3; nu <- p+3
detsims <- apply(rwishart(10000, nu, Sigma=diag(p), Theta=tcrossprod(1:p)),
                 3, det)
any(detsims < .Machine$double.eps)

# p-1 <= nu < 2p-1
nsims <- 50000
p <- 4
nu <- 4 # 3 # 3.5
Sigma <- diag(p) # toeplitz(p:1)
Theta <- tcrossprod(c(1, rep(0,p-1))) #  toeplitz(p:1)
Wsims <- rwishart(nsims, nu, Sigma, Theta)
round(apply(Wsims, 1:2, mean), 3)
nu*Sigma + Theta
Z <- 0.04*diag(p) + matrix(0.01, p, p)
phi(Z, nu, Sigma, Theta)
mean(apply(Wsims, 3, function(x) exp(tr(1i*Z%*%x))))
# for Theta = (theta, 0, ...) and Sigma=Id
detsims <- apply(Wsims, 3, function(x) det(x, symmetric=TRUE))
chi2sims <- rchisq(10000, nu, 1)
for(i in 1:(p-1)){
  chi2sims <- chi2sims*rchisq(10000, nu-i)
}
curve(ecdf(detsims)(x), to=quantile(detsims, 0.95))
curve(ecdf(chi2sims)(x), add=TRUE, col="red")

