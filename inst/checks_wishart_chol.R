# characteristic function
tr <- function(x) sum(diag(x))
phi <- function(Z, nu, Sigma){
  p <- nrow(Sigma)
  complexplus::Det(diag(p) - 2*1i*Z%*%Sigma)^(-nu/2)
}

## nu = p-1 - marche pas bien
p <- 3
nu <- 2
nsims <- 150000
csims <- matrixsampling:::rwishart_chol_I(nsims, nu, p)
wsims <- array(apply(csims, 3, tcrossprod), dim=c(p,p,nsims))
apply(wsims, 1:2, mean)
Z <- 0.5*diag(p) + matrix(0.1, p, p)
mean(apply(wsims, 3, function(x) exp(1i*tr(Z*x))))
phi(Z, nu, diag(p))

Sigma <- toeplitz(p:1)
csims <- rwishart_chol(nsims, nu, Sigma)
wsims <- array(apply(csims, 3, tcrossprod), dim=c(p,p,nsims))
apply(wsims, 1:2, mean); nu*Sigma
Z <- 0.25*diag(p) + matrix(0.1, p, p)
mean(apply(wsims, 3, function(x) exp(1i*tr(Z*x))))
phi(Z, nu, Sigma)


