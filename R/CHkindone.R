#' Sampler of the matrix variate confluent hypergometric kind one distribution
#'
#' Samples the matrix variate confluent hypergometric kind one distribution.
#'
#' @param n sample size, a positive integer
#' @param nu shape parameter, a positive number; if \code{nu < (p-1)/2},
#' where \code{p} is the dimension (the order of \code{Sigma}), then \code{nu}
#' must be a half integer
#' @param alpha,beta shape parameters with the following constraints:
#' \code{b = a} or \code{b > a > nu + (p-1)/2}
#' @param theta scale parameter, a positive number
#' @param Sigma scale matrix, a symmetric positive definite matrix, or
#' \code{NULL} for the identity matrix of order \code{p}
#' @param p if \code{Sigma} is \code{NULL}, this sets \code{Sigma} to the
#' identity matrix of order \code{p}; ignored if \code{Sigma} is not \code{NULL}
#' @param checkSymmetry logical, whether to check that \code{Sigma} is a
#' symmetric positive definite matrix
#'
#' @return A numeric three-dimensional array;
#' simulations are stacked along the third dimension.
#' @export
#'
#' @note For \code{alpha = beta}, this is the matrix variate Gamma distribution
#' with parameters \code{nu}, \code{theta}, \code{Sigma}.
#'
#' @references Gupta & al.
#' Properties of Matrix Variate Confluent Hypergeometric Function Distribution.
#' \emph{Journal of Probability and Statistics}
#' vol. 2016, Article ID 2374907, 12 pages, 2016.
#'
#' @examples nu <- 5; alpha <- 10; beta <- 12; theta <- 2; p <- 3; Sigma <- toeplitz(3:1)
#' CHsims <- rmatrixCHkind1(10000, nu, alpha, beta, theta, Sigma)
#' # simulations of the trace
#' sims <- apply(CHsims, 3, function(X) sum(diag(X)))
#' mean(sims)
#' theta * nu * (nu-beta+(p+1)/2) / (nu-alpha+(p+1)/2) * sum(diag(Sigma))
rmatrixCHkind1 <-
  function(n, nu, alpha, beta, theta = 1, Sigma = NULL, p, checkSymmetry = TRUE){
    if(alpha == beta){
      return(rmatrixgamma(n, nu, theta, Sigma, p, checkSymmetry))
    }
    if(is.null(Sigma) && missing(p)){
      stop("Enter `Sigma` or `p`")
    }
    if(!is.null(Sigma) && !(isRealScalar(Sigma) || nrow(Sigma) == ncol(Sigma))){
      stop("`Sigma` is not square")
    }
    if(checkSymmetry && !is.null(Sigma)){
      . <- isSymmetricPositive(Sigma)
    }
    SigmaIsId <- is.null(Sigma) || isTRUE(all.equal(Sigma, diag(nrow(Sigma))))
    if(!is.null(Sigma)) p <- nrow(Sigma)
    stopifnot(nu > 0, alpha-nu > (p-1)/2, beta > alpha, theta > 0)
    Beta <- rmatrixbeta(n, p, alpha-nu, beta-alpha, def = 1)
    if(is.null(Sigma)){
      W <- rwishart_I(n, 2*nu, p)
    }else{
      W <- rwishart(n, 2*nu, Sigma, checkSymmetry = FALSE)
    }
    CH <- array(NA_real_, dim = c(p,p,n))
    for(i in 1:n){
      S <- invsqrtm(Beta[,,i])
      CH[,,i] <- S %*% W[,,i] %*% S
    }
    theta/2 * CH
  }
