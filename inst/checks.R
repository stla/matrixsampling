#### checks Wishart ####

## check moments - ok singular
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

## check determinant
detsims <- apply(sims, 3, det)
chi2sims <- rchisq(nsims, nu)
for(i in seq_len(p-1)){
  chi2sims <- chi2sims*rchisq(nsims, nu-i)
}
curve(ecdf(detsims/det(S))(x), from=0, to=200)
curve(ecdf(chi2sims)(x), add=TRUE, col="red")

## check trace - ok for nu < p
tr <- function(x) sum(diag(x))
mean(apply(sims, 3, tr))
nu*tr(S) # ok singular
mean(apply(sims, 3, function(x) tr(solve(S)%*%x)*tr(x)))
nu*2*(nu*p/2+1)*tr(S)

## check characteristic function - ok singular
# of the trace
z <- .04
mean(apply(sims, 3, function(x) exp(1i*z*tr(x))))
complexplus::Det(diag(p) - 2*1i*z*S)^(-nu/2)
# of the matrix
Z <- matrix(0.04, p, p)
mean(apply(sims, 3, function(x) exp(tr(1i*Z%*%x))))
complexplus::Det(diag(p) - 2*1i*Z%*%S)^(-nu/2)

## check singular case
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


## check U-distribution
p <- 3
Sigma <- rwishart(1, p, diag(p))[,,1]
nu1 <- 6
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

#### checks inverse Wishart ####
## check moment
nu <- 5
p <- 3
S <- rwishart(1, p, diag(p))[,,1]
nsims <- 100000
sims <- rinvwishart(nsims, nu, S)
apply(sims, c(1,2), mean)
S / (nu-p-1)

#### checks noncentral Wishart ####
p <- 3
nu <- 6.5
Sigma <- rwishart(1, p, diag(p))[,,1]
Theta <- rwishart(1, p, diag(p))[,,1]
nsims <- 50000
sims <- rwishart(nsims, nu, Sigma, Theta)

## check mean
apply(sims, 1:2, mean)
Theta + nu*Sigma
## check characteristic function - faux dans Gupta Nagar
# formule trouvée dans http://www.tandfonline.com/doi/abs/10.1080/00949658208810570
Z <- matrix(0.04, p, p)
mean(apply(sims, 3, function(x) exp(tr(1i*Z%*%x))))
complexplus::Det(diag(p) - 2*1i*Z%*%Sigma)^(-nu/2) *
  exp(1i*tr(solve(diag(p) - 2*1i*Z%*%Sigma) %*% Z %*% Theta))

#### checks Beta matrix ####
n1 <- 16
n2 <- 3
p <- 3
Theta <- rwishart(1, p, diag(p))[,,1]
nsims <- 20000
sims <- rmatrixbeta(nsims, p, n1, n2, Theta2=Theta)

## check moments of determinant
detsims <- apply(sims, 3, det)

mgamma <- function(p, a){ # multivariate Gamma
  exp(log(pi) * (p * (p - 1)/4) + sum(lgamma(a + (1-(1:p))/2)))
}
h1F1 <- function(alpha, beta, R){ # matricial 1F1 function
  p <- nrow(R)
  sims <- apply(rmatrixbeta(20000, p, alpha, beta-alpha), 3,
                function(x) exp(tr(R%*%x)))
  mean(sims)
}
h <- 1
mgamma(p,n1+h)*mgamma(p,n1+n2)/mgamma(p,n1+n2+h)/mgamma(p,n1) *
  exp(tr(-Theta/2)) *
  h1F1(n1+n2, n1+n2+h, Theta/2)
mean(detsims^h)

## check p=1
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


#### checks matrix normal ####

#### checks matrix Student ####
m <- 3; p <- 5
M <- matrix(1, m, p)
U <- rwishart(1, m, diag(m))[,,1] # arbitray U matrix
V <- rwishart(1, p, diag(p))[,,1] # arbitrary V matrix
Uinv <- solve(U)
Vinv <- solve(V) # inverse V
nu <- 6 # degrees of freedom

## check U-distribution
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

## check covariance matrix
Cov <- cov(t(apply(sims, 3, function(x) c(t(x)))))
k <- nu-2
round((1/k * kronecker(U,V))[1:3,1:3],2)
round(Cov[1:3,1:3], 2)

## check transformation to matrix beta
Uroot <- matrixsampling:::matrixroot(U)
bsims <- array(apply(X, 3,
                     function(x) Uroot %*% solve(U + x %*% Vinv %*% t(x)) %*% Uroot),
               dim = c(m,m,nsims))
bsims2 <- rmatrixbeta(nsims, m, (nu + m - 1)/2, p/2)
Z <- matrix(0.04, m, m)
mean(apply(bsims, 3, function(x) exp(tr(1i*Z%*%x))))
mean(apply(bsims2, 3, function(x) exp(tr(1i*Z%*%x))))

## check transformation to matrix beta II
Urootinv <- matrixsampling:::invsqrtm(U)
bsims <- array(apply(X, 3,
                     function(x) Urootinv %*% x %*% Vinv %*% t(x) %*% Urootinv),
               dim = c(m,m,nsims))
bsims2 <- rmatrixbetaII(nsims, m, p/2, (nu+m-1)/2)
mean(apply(bsims, 3, function(x) exp(tr(1i*Z%*%x))))
mean(apply(bsims2, 3, function(x) exp(tr(1i*Z%*%x))))


#### checks inverted matrix Student ####
## thm 5.2.3 G & N
nsims <- 50000
p <- 3
n1 <- 5
n2 <- 4
Sigma <- rwishart(1, p, diag(p))[,,1]
X <- rmatrixnormal(nsims, matrix(0, p, n1), Sigma, diag(n1))
S2 <- rwishart(nsims, n2, Sigma)
Z <- array(NA_real_, dim=c(p,n1,nsims))
for(i in 1:nsims){
  Z[,,i] <- matrixsampling:::invsqrtm(S2[,,i]+tcrossprod(X[,,i]))%*%X[,,i]
}
IT <- rmatrixit(nsims, n2-p+1, matrix(0, p, n1), diag(p), diag(n1))
z <- 0.04*diag(p) #matrix(0.04, p, p)
mean(apply(Z, 3, function(x) exp(tr(1i*z%*%x))))
mean(apply(IT, 3, function(x) exp(tr(1i*z%*%x))))
apply(Z, 1:2, var)
apply(IT, 1:2, var)
curve(ecdf(Z[1,1,]+Z[1,2,])(x), from=-1.5, to=1.5)
curve(ecdf(IT[1,1,]+IT[1,2,])(x), add=TRUE, col="red")

## checks moments
nu <- -1 ## pourquoi je peux inverser dans rmatrixit ?? car >0 + >=0
m <- 3
p <- 2
M <- matrix(1, m, p)
U <- toeplitz(m:1)
V <- toeplitz(p:1)
ITsims <- rmatrixit(10000, nu, M, U, V)
apply(ITsims, 1:2, mean) # approximates M
vecITsims <- t(apply(ITsims, 3, function(X) c(t(X))))
round(cov(vecITsims) - 1/(nu+m+p-1) * kronecker(U,V), 1)

## check p = m = 1
nu <- 3
m <- 1
p <- 1
M <- matrix(0, m, p)
U <- toeplitz(m:1)
V <- toeplitz(p:1)
ITsims <- rmatrixit(10000, nu, M, U, V)
vecITsims <- t(apply(ITsims, 3, function(X) c(t(X))))
curve(ecdf(vecITsims[1,])(x), from=-1, to =1)
chi2sims <- rchisq(10000, nu)
normsims <- rnorm(10000)
curve(ecdf(normsims/sqrt(chi2sims + normsims^2))(x), add=TRUE, col="red")
# le carré est une Beta(1/2, nu/2) - THM 5.2.4 G&N
# distribution of correlation ?
library(mvtnorm)
corsims <- numeric(10000)
for(i in 1:10000){
  corsims[i] <- cor(rmvnorm(nu+2, sigma = diag(2)))[1,2]
}
curve(ecdf(corsims)(x), add=TRUE, col="blue")
