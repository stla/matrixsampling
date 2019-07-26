tr <- function(x) sum(diag(x))
mgamma <- function(p, a){ # multivariate Gamma
  exp(log(pi) * (p * (p - 1)/4) + sum(lgamma(a + (1-(1:p))/2)))
}
h1F1 <- function(alpha, beta, R){
  p <- nrow(R)
  if(alpha < beta){
    sims <- apply(rmatrixbeta(20000, p, alpha, beta-alpha, def=1), 3,
                  function(x) exp(tr(R%*%x)))
  }else{ # Constantine
    sims <- apply(rwishart(20000, beta/2, R), 3, det)
    t <- alpha - beta
    mean(sims^t)/mgamma(p, t+beta/2)*mgamma(beta/2)/exp(-tr(R))
  }
  mean(sims)
}

h1F1(2,4,0.5*diag(1))
gsl::hyperg_1F1(2,4,0.5)


n <- 2
A <- diag(n) +1
h1F1(1/2, n/2, A)
#
nsims <- 10000
sims <- numeric(nsims)
for(i in 1:nsims){
  z <- rnorm(n)
  sims[i] <- f(z/crossprod(z))
}
mean(sims)

p <- 3
a <- 4; b <- 5; Theta <- toeplitz(p:1)
nsims <- 20000
sims <- rmatrixbeta(nsims, p, a, b, Theta2=Theta)
detsims <- apply(sims, 3, det)
h <- 1
mgamma(p,a+h)*mgamma(p,a+b)/mgamma(p,a+b+h)/mgamma(p,a) *
  exp(tr(-Theta/2)) *
  h1F1(a+b, a+b+h, Theta/2)
mean(detsims^h)


