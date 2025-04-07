#' Estimate Long-Run Covariance Matrix
#'
#' Computes localised estimates of the long-run covariance matrix of the vectorised outer products \eqn{F_t F_t^\top - I_r},
#' using a rolling window HAC-type approach. This is used in standardising test statistics in change point detection algorithms.
#'
#' @param F_hat A matrix of dimension \eqn{r \times T}, where each column represents observations of \eqn{r} latent factors at a given time.
#' @param Time An integer specifying the number of time points (defaults to \code{ncol(F_hat)}).
#' @param window_size Optional integer giving the size of the local rolling window used to estimate time-varying covariances.
#'                    Defaults to \code{floor(Time^{1/4})}.
#'
#' @return A list with components:
#' \describe{
#'   \item{V_full}{The full estimated long-run covariance matrix of dimension \eqn{p \times p}, where \eqn{p = r(r+1)/2}.}
#'   \item{V_diag}{The diagonal matrix formed from the diagonal entries of \code{V_full}.}
#'   \item{V_non}{An identity matrix of dimension \eqn{p \times p}, used as a naive alternative.}
#' }
#'
#' @details
#' The function computes, for each time point \eqn{t}, the vectorised difference \eqn{F_t F_t^\top - I_r}, using the \code{vech}
#' operator (half-vectorisation with diagonal). Then, for each \eqn{t}, a sample covariance matrix is computed over a local window
#' of neighbouring time points. The full long-run covariance matrix is obtained by averaging over all local covariance matrices.
#'
#' This method avoids global assumptions and accommodates time variation in second-order structure.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' F_hat <- matrix(rnorm(3 * 300), nrow = 3)
#' V_list <- long_run_V(F_hat)
#' str(V_list$V_full)
#' }


long_run_V <- function(F_hat, Time = ncol(F_hat), window_size = NULL) {
  # HAC
  r <- nrow(F_hat)
  p <- r * (r + 1) / 2

  if (is.null(window_size)) {
    window_size <- floor(Time^(1/4))
  }

  vech <- function(A) as.vector(A[upper.tri(A, diag = TRUE)])
  diff_series <- matrix(NA, nrow = p, ncol = Time)
  for (t in 1:Time) {
    F_t <- F_hat[, t, drop = FALSE]
    diff_matrix <- F_t %*% t(F_t) - diag(r)
    diff_series[, t] <- vech(diff_matrix)
  }

  local_cov_matrices <- array(0, dim = c(p, p, Time))
  for (t in 1:Time) {
    win_idx <- seq(max(1, t - floor(window_size / 2)),
                   min(Time, t + floor(window_size / 2)))
    if (length(win_idx) > 1) {
      local_cov_matrices[,,t] <- cov(t(diff_series[, win_idx, drop = FALSE]))
    } else {
      local_cov_matrices[,,t] <- matrix(0, nrow = p, ncol = p)
    }
  }

  V_full <- apply(local_cov_matrices, c(1, 2), mean, na.rm = TRUE)
  V_diag <- diag(diag(V_full))
  V_non <- diag(p)

  return(list(V_full = V_full, V_diag = V_diag, V_non = V_non))
}

