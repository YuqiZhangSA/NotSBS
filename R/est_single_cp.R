#' Estimate a Single change point within a Given Interval
#'
#' This function computes a CUSUM-based statistic for detecting a single change point
#' within the interval \code{[st, ed]} in a multivariate time series matrix.
#' It compares empirical covariance matrices before and after each candidate point \eqn{k}, adjusted for long-run variance or local estimates.
#'
#' @param F_hat A matrix of dimension \eqn{r \times T}, representing the observed time series data with \eqn{r} variables and \eqn{T} time points.
#' @param st An integer specifying the start index of the interval.
#' @param ed An integer specifying the end index of the interval.
#' @param trim A positive integer specifying the trimming parameter to avoid estimating near the boundaries. Requires \code{ed - st > 2 * trim}.
#' @param V An optional positive-definite covariance matrix for standardising the test statistic. If \code{NULL}, a local estimate is used.
#' @param window_size An optional integer indicating the local window size used when estimating the covariance matrix for standardisation (only relevant if \code{V = NULL}).
#'
#' @return A data frame with one row and columns:
#' \describe{
#'   \item{est.cp}{Estimated location of the change point within \code{[st, ed]}.}
#'   \item{val}{The associated CUSUM statistic value.}
#'   \item{st}{The interval start index.}
#'   \item{ed}{The interval end index.}
#'   \item{trim}{The trimming parameter used.}
#' }
#'
#' @details
#' For each candidate point \eqn{k \in (st + trim, ed - trim)}, the function computes the difference between
#' the empirical covariance matrices before and after \eqn{k}. The Frobenius norm (or Mahalanobis distance)
#' of the vectorised difference is used as a CUSUM score.
#'
#' If \code{V} is supplied, it is used for standardisation. Otherwise, the function estimates a local covariance
#' matrix around each \eqn{k} using a moving window and applies the corresponding inverse for scaling.
#'
#' Boundary points are forced to have score zero to prevent spurious maxima near edges.
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' F_hat <- matrix(rnorm(200), nrow = 5)
#' find_single_cp_std(F_hat, st = 1, ed = 40, trim = 5)
#' }
#'
#' @export


find_single_cp_std <- function(F_hat, st, ed, trim, V = NULL, window_size = NULL) {
  stopifnot(ed - st > 2 * trim)
  T_time <- ncol(F_hat)

  if (is.null(window_size)) {
    window_size <- floor(T_time^(1/4))
  }

  k_range <- (st + trim + 1):(ed - trim)
  scores <- numeric(length(k_range))

  vech <- function(A) as.matrix(A[upper.tri(A, diag = TRUE)])

  M_ske_list <- vector("list", length(k_range))
  for (i in seq_along(k_range)) {
    k <- k_range[i]

    F_left <- F_hat[, (st + 1):k, drop = FALSE]
    F_right <- F_hat[, (k + 1):ed, drop = FALSE]

    Gamma_F_left <- cov(t(F_left))
    Gamma_F_right <- cov(t(F_right))

    diff_cov <- Gamma_F_right - Gamma_F_left

    M_ske <- sqrt(((k - st) * (ed - k)) / (ed - st)) * vech(diff_cov)
    M_ske_list[[i]] <- M_ske
  }

  if (!is.null(V)) {
    V_inv <- solve(V)
    for (i in seq_along(k_range)) {
      scores[i] <- sqrt(abs(t(M_ske_list[[i]]) %*% V_inv %*% M_ske_list[[i]]))
    }
  } else {
    scores <- sapply(seq_along(k_range), function(i) {
      k <- k_range[i]
      win_idx <- seq(max(1, k - floor(window_size / 2)),
                     min(T_time, k + floor(window_size / 2)))

      r <- nrow(F_hat)
      p <- r * (r + 1) / 2
      diff_series <- matrix(NA, nrow = p, ncol = length(win_idx))
      for (j in seq_along(win_idx)) {
        t_idx <- win_idx[j]
        F_t <- F_hat[, t_idx, drop = FALSE]
        diff_matrix <- F_t %*% t(F_t) - diag(r)
        diff_series[, j] <- vech(diff_matrix)
      }

      local_cov <- if (length(win_idx) > 1) cov(t(diff_series)) else diag(1, p)
      V_local_inv <- tryCatch(solve(local_cov), error = function(e) diag(1, nrow(local_cov)))
      sqrt(abs(t(M_ske_list[[i]]) %*% V_local_inv %*% M_ske_list[[i]]))
    })
  }

  scores[1] <- 0
  scores[length(k_range)] <- 0

  best_index <- which.max(scores)
  best_k <- k_range[best_index]
  best_score <- scores[best_index]

  return(data.frame(est.cp = best_k, val = best_score, st = st, ed = ed, trim = trim))
}
