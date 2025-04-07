#' Simulate Factor Model Data with Structural Breaks (Duan 2023)
#'
#' Generates high-dimensional time series data from a dynamic factor model with optional structural breaks
#' in the factor loadings. This simulation design follows Duan et al. (2023) and allows for several types
#' of change-point behaviour in the factor structure.
#'
#' @param Time An integer specifying the length of the time series to be returned (post-burn-in).
#' @param p An integer giving the number of observed variables.
#' @param type An integer vector specifying the type of change at each break:
#' \describe{
#'   \item{0}{No change.}
#'   \item{1}{Rank reduction (drop in one factor loading).}
#'   \item{2}{Rotational change in the loadings.}
#'   \item{3}{Linear transformation (e.g., upper triangular transformation).}
#'   \item{4}{Completely new random loadings (increase in rank).}
#' }
#' @param dep Logical. If \code{TRUE}, induces temporal dependence in the factors and idiosyncratic errors.
#' @param seed Optional integer to set a reproducible random seed.
#' @param theta_coef A numeric vector (default based on \code{Time}) specifying relative change-point locations, as proportions of \code{Time}.
#'
#' @return A list with components:
#' \describe{
#'   \item{x}{A \code{p × Time} matrix of simulated observations.}
#'   \item{r}{An integer giving the final number of latent factors after all breaks.}
#'   \item{m}{A vector of change-point locations.}
#'   \item{type}{A vector indicating the type of change applied at each change-point.}
#' }
#'
#' @details
#' The function generates factors \code{f} and idiosyncratic errors \code{e}, optionally with AR(1) structure
#' (if \code{dep = TRUE}). Factor loadings are generated initially and then transformed at each break according
#' to the specified \code{type}. The signal \code{chi} is constructed by projecting the factors through the
#' piecewise-constant loading matrices, and noise is added to obtain the final observed matrix \code{x}.
#'
#' The number of change-points is determined by the length of \code{type}, and their locations are based on
#' \code{theta_coef}, unless \code{type = 0}, which gives a no-change scenario.
#'
#' @importFrom mvtnorm rmvnorm
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- duan_dgp(Time = 300, p = 50, type = c(1, 2, 4), dep = TRUE, seed = 123)
#' image(sim$x)
#' }


duan_dgp <- function(Time, p, type, dep = FALSE, seed = NULL, theta_coef = c(max(0.5-10*log(Time)/Time,0.42),0.5,0.5+6*sqrt(Time)/Time)){
  if (!is.null(seed)) set.seed(seed)

  burnin <- Time
  r <- r0 <- 3
  if(type[1] == 0){
    k0 <- Time
    brks <- c(0, Time)
  } else{
    #k0 <- floor(Time * 1:length(type) / (length(type) + 1))
    k0 <- c(floor(Time*theta_coef[1]), floor(Time*theta_coef[2]), floor(Time*theta_coef[3]))
    brks <- c(0, k0, Time)
  }

  beta <- c(0, .3)[2]
  Omega <- toeplitz(beta^(1:p - 1))

  if(dep){
    rho <- .7
    alpha <- .3
  } else rho <- alpha <- 0

  f <- matrix(rnorm(r0 * (Time + burnin)), nrow = r0)
  e <- t(mvtnorm::rmvnorm(Time + burnin, sigma = Omega))
  for(tt in 2:dim(f)[2]){
    f[, tt] <- rho * f[, tt - 1] + f[, tt]
    e[, tt] <- alpha * e[, tt - 1] + e[, tt]
  }
  f <- f[, -(1:burnin), drop = FALSE]
  e <- e[, -(1:burnin)]

  lam <- array(0, dim = c(p, r0, length(k0) + 1))
  lam[,, 1] <- matrix(rnorm(r0 * p), nrow = p) / sqrt(r0)

  for(kk in 1:length(k0)){
    if(type[kk] == 1){
      C <- diag(rep(1, r0)); C[r0, r0] <- 0
      lam[,, kk + 1] <- lam[,, 1] %*% C
      r <- r + 0
    } else if(type[kk] == 2){
      C <- matrix(0, r0, r0)
      diag(C) <- c(.5, 1, 1.5)[1:r0]
      C[lower.tri(C, diag = FALSE)] <- rnorm(r0 * (r0 - 1)/2)
      lam[,, kk + 1] <- lam[,, 1] %*% C
      r <- r + 0
    } else if(type[kk] == 3){
      m <- 1
      C <- matrix(c(1, 0, 0, 2, 1, 0, 3, 2, m), byrow = TRUE, nrow = r0)
      lam[,, kk + 1] <- lam[,, 1] %*% C
      r <- r + 0
    } else if(type[kk] == 4){
      lam[,, kk + 1] <- matrix(rnorm(r0 * p), nrow = p) / sqrt(r0)
      r <- r + r0
    }
  }

  chi <- matrix(0, nrow = p, ncol = Time)
  for(jj in 1:(length(brks) - 1)){
    int <- (brks[jj] + 1):brks[jj + 1]
    chi[, int] <- lam[,, jj] %*% f[, int, drop = FALSE]
  }

  x <- chi + e

  return(list(x = x, r = r, m = k0, type = type))

}
