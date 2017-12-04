d <- 6
nsims <- 100000
nu <- d-1 + 0.001
Sigma <- matrixsampling::rwishart(1, d, diag(d))[,,1]
Theta <- matrixsampling::rwishart(1, 1, diag(d))[,,1]

detsims <-
  apply(matrixsampling:::rwishart_I(nsims, nu, d, epsilon=.Machine$double.eps^0.45), 3,
        function(x) det(x, symmetric=TRUE))
length(which(detsims<.Machine$double.eps))
length(which(detsims<sqrt(.Machine$double.eps)))

detsims <- apply(matrixsampling:::rwishart_root_I(nsims, nu, d), 3,
                 function(x) prod(diag(x)))
length(which(detsims<.Machine$double.eps))
length(which(detsims<sqrt(.Machine$double.eps)))
detsims <- apply(matrixsampling:::rwishart_root_I(nsims, nu, d), 3,
                 function(x) det(x, symmetric=TRUE))
length(which(detsims<.Machine$double.eps))
length(which(detsims<sqrt(.Machine$double.eps)))

detsims <- apply(matrixsampling:::rwishart_chol(nsims, nu, diag(d)), 3,
                 function(x) prod(diag(x)))
length(which(detsims<.Machine$double.eps))
length(which(detsims<sqrt(.Machine$double.eps)))

detsims <-
  apply(matrixsampling:::rwishart(nsims, nu, Sigma, epsilon=.Machine$double.eps^0.25), 3,
        function(x) det(x, symmetric=TRUE))
length(which(detsims<.Machine$double.eps))
length(which(detsims<sqrt(.Machine$double.eps)))

winvsims <- matrixsampling::rinvwishart(nsims, nu, Sigma)
