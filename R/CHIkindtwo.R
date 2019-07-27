#' Sampler of the matrix variate type I confluent hypergeometric kind two
#' distribution
#'
#' Samples the matrix variate type I confluent hypergeometric kind two
#' distribution.
#'
#' @param n sample size, a positive integer
#' @param nu shape parameter, a positive number; if \code{nu < (p-1)/2},
#' then \code{nu} must be a half integer
#' @param alpha,beta shape parameters; \code{alpha > nu + (p-1)/2},
#' \code{beta < nu + 1}
#' @param theta scale parameter, a positive number
#' @param p dimension (order of the sampled matrices), an integer greater than
#' or equal to one
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension.
#' @export
#'
#' @references A. K. Gupta & D. K. Nagar. \emph{Matrix Variate Distributions}.
#' Chapman & Hall/CRC, Boca Raton (2000).
#'
#' @examples
#' nu <- 5; alpha <- 13; beta <- 4; theta <- 2; p <- 2
#' sims <- rmatrixCHIkind2(50000, nu, alpha, beta, theta, p = 2)
#' # simulations of the trace
#' trsims <- apply(sims, 3, function(X) sum(diag(X)))
#' mean(trsims)
#' p * theta * nu * (nu+(p+1)/2-beta) / (alpha-nu-(p+1)/2)
rmatrixCHIkind2 <- function(n, nu, alpha, beta, theta = 1, p){
  stopifnot(nu > 0, alpha > nu+(p-1)/2, beta < nu+1)
  Beta <- rmatrixbetaII(n, p, nu+(p+1)/2-beta, alpha-nu, def = 1)
  W <- rwishart_I(n, 2*nu, p)
  CH <- array(NA_real_, dim = c(p,p,n))
  for(i in 1:n){
    S <- sqrtm(Beta[,,i])
    CH[,,i] <- S %*% W[,,i] %*% S
  }
  theta/2 * CH
}
