## BetaI
nsims <- 50000
p <- 3
a <- 2; b <- 2.5
Usims <- rmatrixbeta(nsims, p, a, b)

## check U-distribution for BetaI
detsims <- apply(Usims, 3, det)
curve(ecdf(detsims)(x), to=quantile(detsims, 0.9))
e <- 2*a; h <- 2*b
betasims <- rbeta(nsims, e/2, h/2)
for(i in seq_len(p-1)){
  betasims <- betasims * rbeta(nsims, (e-i)/2, h/2)
}
curve(ecdf(betasims)(x), add=TRUE, col="red")

# BetaI noncentral rank one
theta <- 3
Theta <- diag(c(theta, rep(0,p-1)))
Usims <- rmatrixbeta(nsims, p, a, b, Theta2=Theta, def=2)
n <- 2*a; m <- 2*b
betasims <- 1/(1 + rchisq(nsims, p, theta)/rchisq(nsims, n-p+m))
for(i in seq_len(m-1)){
  betasims <- betasims * rbeta(nsims, (n-p+m-i)/2, p/2)
}
detsims <- apply(Usims, 3, det)
curve(ecdf(detsims)(x), to=quantile(detsims,0.8))
curve(ecdf(betasims)(x), add=TRUE, col="red")
Usims <- rmatrixbeta(nsims, p, a, b, Theta1=Theta)
betasims <- rbeta(nsims, n/2, m/2, theta)
for(i in seq_len(p-1)){
  betasims <- betasims * rbeta(nsims, (n-i)/2, m/2)
}
detsims <- apply(Usims, 3, det)
curve(ecdf(detsims)(x), to=quantile(detsims,0.8))
curve(ecdf(betasims)(x), add=TRUE, col="red")

## BetaII
nsims <- 50000
p <- 3
a <- 2; b <- 2.5
Vsims <- rmatrixbetaII(nsims, p, a, b)

## check U'-distribution for BetaI
detsims <- apply(Vsims, 3, det)
curve(ecdf(detsims)(x), to=quantile(detsims, 0.9))
e <- 2*a; h <- 2*b
betasims <- brr::rbeta2(nsims, e/2, h/2, scale=1)
for(i in seq_len(p-1)){
  betasims <- betasims * brr::rbeta2(nsims, (e-i)/2, (h-i)/2, scale=1)
}
curve(ecdf(betasims)(x), add=TRUE, col="red")

## check transfo BetaI <-> BetaII
nsims <- 50000
p <- 4
a <- 2; b <- 2
Theta1 <- toeplitz(p:1) #diag(p:1) # diag(p) # toeplitz(p:1)
Usims <- rmatrixbeta(nsims, p, a, b, Theta2=Theta1) # !!! pas terrible quand c'est Theta1 et a,b petits!
Vsims <- rmatrixbetaII(nsims, p, a, b, Theta2=Theta1)
Usims2 <- array(apply(Vsims, 3,
                      function(V){
                        # A <- matrixsampling:::invsqrtm(diag(p)+V)
                        # A%*%V%*%A
                        chol2inv(chol(diag(p)+V))%*%V
                      }), dim=dim(Vsims))
Vsims2 <- array(apply(Usims, 3,
                      function(U){
                        # A <- matrixsampling:::invsqrtm(diag(p)-U)
                        # A%*%U%*%A
                        U%*%chol2inv(chol(diag(p)-U))
                      }), dim=dim(Vsims))
curve(ecdf(Usims[1,1,]+Usims[1,2,])(x))
curve(ecdf(Usims2[1,1,]+Usims2[1,2,])(x), add=TRUE, col="red")
curve(ecdf(Vsims[1,1,]+Vsims[1,2,])(x), from=-1, to=40)
curve(ecdf(Vsims2[1,1,]+Vsims2[1,2,])(x), add=TRUE, col="red")
detU <- apply(Usims, 3, det)
detU2 <- apply(Usims2, 3, det)
curve(ecdf(detU)(x), to=quantile(detU, 0.9))
curve(ecdf(detU2)(x), add=TRUE, col=2)
detV <- apply(Vsims, 3, function(V) det(diag(p)+V))
detV2 <- apply(Vsims2, 3, function(V) det(diag(p)+V))
curve(ecdf(detV)(x), to=quantile(detV, 0.9))
curve(ecdf(detV2)(x), add=TRUE, col=2)
trU <- apply(Usims, 3, tr)
trU2 <- apply(Usims2, 3, tr)
curve(ecdf(trU)(x), to=quantile(trU, 0.9))
curve(ecdf(trU2)(x), add=TRUE, col=2)

Z <- 0.4*diag(p) + matrix(0.2,p,p)
mean(apply(Usims, 3, function(x) exp(1i*tr(Z%*%x))))
mean(apply(Usims2, 3, function(x) exp(1i*tr(Z%*%x))))
z <- seq(0.01, 4, length.out = 20)
phi1 <- sapply(z, function(z) mean(apply(Usims, 3, function(x) exp(1i*tr((z*diag(p)+matrix(z,p,p))%*%x)))))
phi2 <- sapply(z, function(z) mean(apply(Usims2, 3, function(x) exp(1i*tr((z*diag(p)+matrix(z,p,p))%*%x)))))
plot(Re(phi1), type="o", pch=19)
lines(Re(phi2), col="red", type="o", pch=19)
plot(Im(phi1), type="l")
lines(Im(phi2), col="red")


#### noncentral ####
n1 <- 2; n2 <- 2
p <- 3
Theta <- tcrossprod(1:p) # rwishart(1, p, diag(p))[,,1]
nsims <- 100000
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
n1 <- 16; n2 <- 3
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


#### simulates Beta II by conditional approach ####
p <- 3
a <- (p-1)/2+0.501; b <- (p-1)/2+0.501
p2 <- 1; p1 <- p-1
b-p2/2; p1/2 # pas poss
a-p1/2; p2/2 # pas poss
nsims <- 10000
V11 <- rmatrixbetaII(nsims, p1, a, b-p2/2)
V221 <- rmatrixbetaII(nsims, p2, a-p1/2, b)
V <- array(NA_real_, dim=c(p,p,nsims))
for(i in 1:nsims){
  V21 <-
    matrix(rmatrixt(1, 2*a+2*b-p+1, matrix(0, p2, p1),
                    diag(p2)+V221[,,i],
                    V11[,,i]*(diag(p1)+V11[,,i]), checkSymmetry = FALSE)[,,1], p2, p1)
  V22 <- V221[,,i] + V21%*%chol2inv(chol(V11[,,i]))%*%t(V21)
  V[,,i] <- cbind(rbind(V11[,,i], V21), rbind(t(V21), V22))
}

detsims <- apply(V, 3, det)
curve(ecdf(detsims)(x), to=quantile(detsims, 0.9))
e <- 2*a
h <- 2*b
betasims <- brr::rbeta2(nsims, e/2, h/2, scale=1)
for(i in seq_len(p-1)){
  betasims <- betasims * brr::rbeta2(nsims, (e-i)/2, (h-i)/2, scale=1)
}
curve(ecdf(betasims)(x), add=TRUE, col="red")


#### the two definitions ####
Phi <- function(sims, z = seq(0.01, 4, length.out = 20)){
  sapply(z,
         function(z) mean(apply(sims, 3, function(x) exp(1i*tr((z*diag(p)+matrix(z,p,p))%*%x)))))
}
VtoU <- function(Vsims){
  array(apply(Vsims, 3,
              function(V){
                chol2inv(chol(diag(p)+V))%*%V
              }), dim=dim(Vsims))
}
#
nsims <- 50000
p <- 4
a <- 2; b <- 2
Theta <- diag(c(5,5,5,5)) #tcrossprod(1:p)
U1 <- rmatrixbeta(nsims, p, a, b, def=1, Theta2=Theta)
V1 <- rmatrixbetaII(nsims, p, a, b, def=1, Theta2=Theta)
U2 <- rmatrixbeta(nsims, p, a, b, def=2, Theta2=Theta)
V2 <- rmatrixbetaII(nsims, p, a, b, def=2, Theta2=Theta)
curve(ecdf(U1[3,2,]+U1[1,2,]*exp(U1[1,2,]))(x))
curve(ecdf(U2[3,2,]+U2[1,2,]*exp(U2[1,2,]))(x), add=TRUE, col="red")
phiU1 <- Phi(U1)
phiU2 <- Phi(U2)
phiV1 <- Phi(VtoU(V1))
phiV2 <- Phi(VtoU(V2))
plot(Re(phiU1), type="o", pch=19)
lines(Re(phiU2), col="red", type="o", pch=19)
lines(Re(phiV1), col="green", type="o", pch=19)
lines(Re(phiV2), col="blue", type="o", pch=19)
plot(Im(phiU1), type="o", pch=19)
lines(Im(phiU2), col="red", type="o", pch=19)
lines(Im(phiV1), col="green", type="o", pch=19)
lines(Im(phiV2), col="blue", type="o", pch=19)
curve(ecdf(V1[4,4,])(x), to=quantile(V1[4,4,], 0.8))
curve(ecdf(V2[4,4,])(x), add=TRUE, col="red")

# Clsion:
# noncentral: tout ok
# Theta2 : différence def1 def2 sauf Theta2 = a*Id ;
#   transfo ok pour chaque def - PAS SÛR pour def1, légère différence !
# Theta1 : def1 def2 pareil pour U (PAS SÛR, légère différence), différence pour V
  # transfo ok pour def2, pas pour def1
# => def 3 pour type I: Asoo's definition U = V(I+V)^{-1} où V ~ BetaII def I


V1transfo <- array(apply(V1, 3,
                         function(V){
                           chol2inv(chol(diag(p)+V))%*%V
                         }), dim=dim(Vsims))
U1transfo <- array(apply(U1, 3,
                         function(U){
                           U%*%chol2inv(chol(diag(p)-U))
                         }), dim=dim(Vsims))
V2transfo <- array(apply(V2, 3,
                         function(V){
                           chol2inv(chol(diag(p)+V))%*%V
                         }), dim=dim(Vsims))
curve(ecdf(V1transfo[4,4,])(x))
curve(ecdf(U1[4,4,])(x), add=TRUE, col="red")
curve(ecdf(V2transfo[4,4,])(x), add=TRUE, col="blue")
