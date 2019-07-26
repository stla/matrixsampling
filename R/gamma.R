#' Matrix Gamma sampler
#'
#' Samples a matrix Gamma distribution.
#'
#' @param n sample size, a positive integer
#' @param nu shape parameter, a positive number
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
#' @details This is the distribution of
#' \ifelse{html}{\out{&theta;/2&times;S}}{\eqn{\frac{\theta}{2}S}} where
#' \ifelse{html}{\out{S ~ W<sub>p</sub>(2&nu;,&Sigma;)}}{\eqn{S \sim \mathcal{W}_p(2\nu,\Sigma)}}.
#'
#' @references Gupta & al.
#' Properties of Matrix Variate Confluent Hypergeometric Function Distribution.
#' \emph{Journal of Probability and Statistics}
#' vol. 2016, Article ID 2374907, 12 pages, 2016.
#'
#' @examples nu <- 3; theta <- 4; Sigma <- toeplitz(2:1)
#' Gsims <- rmatrixgamma(10000, nu, theta, Sigma)
#' apply(Gsims, c(1,2), mean) # should be nu * theta * Sigma
#' nu * theta * Sigma
rmatrixgamma <- function(n, nu, theta, Sigma = NULL, p, checkSymmetry = TRUE){
  stopifnot(nu > 0, theta > 0)
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
  if(SigmaIsId){
    if(!is.null(Sigma)) p <- nrow(Sigma)
    theta/2 * rwishart_I(n, 2*nu, p)
  }else{
    theta/2 * rwishart(n, 2*nu, Sigma, checkSymmetry = FALSE)
  }
}












