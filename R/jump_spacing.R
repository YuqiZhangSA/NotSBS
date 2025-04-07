#' Jump Magnitudes and Segment Spacings
#'
#' Computes the operator norm of changes in covariance structure between adjacent time segments,
#' as well as the spacing (length) between estimated change points.
#'
#' @param X_t A \code{p × T} data matrix, where \code{p} is the number of variables and \code{T} is the number of time points.
#' @param theta A numeric vector of estimated change points (values in \eqn{1:(T-1)}). No need to include \code{0} or \code{T}.
#' @param Time The total number of time points \code{T}. Used to close the final segment.
#'
#' @return A list containing:
#' \describe{
#'   \item{jump_magnitudes}{A numeric vector of spectral norms of covariance differences between adjacent segments.}
#'   \item{min_jump}{The minimum jump magnitude.}
#'   \item{max_jump}{The maximum jump magnitude.}
#'   \item{spacings}{A numeric vector of segment lengths between adjacent change points.}
#'   \item{min_spacing}{The minimum segment length.}
#' }
#'
#' @details
#' The function divides the series into segments using the supplied change points \code{theta}, then computes
#' the empirical covariance matrix within each segment. The spectral norm of the difference between successive
#' covariances is recorded as a measure of "jump" across change points. Segment spacings are defined as the
#' number of time points between consecutive change points.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' X_t <- matrix(rnorm(50 * 200), nrow = 50)
#' theta <- c(60, 130)
#' result <- jump_spacing(X_t, theta, Time = 200)
#' result$jump_magnitudes
#' result$min_spacing
#' }


jump_spacing <- function(X_t, theta, Time) {
  theta <- c(0, theta, Time)

  num_segments <- length(theta) - 1
  segment_covariances <- vector("list", num_segments)
  jump_magnitudes <- numeric(num_segments - 1)
  spacings <- numeric(num_segments)

  for (j in 1:num_segments) {
    t_start <- theta[j] + 1
    t_end <- theta[j + 1]

    X_segment <- X_t[, t_start:t_end]

    cov_X_j <- cov(t(X_segment))
    segment_covariances[[j]] <- cov_X_j

    spacings[j] <- t_end - t_start + 1

    if (j > 1) {
      cov_diff <- segment_covariances[[j]] - segment_covariances[[j - 1]]
      jump_magnitudes[j - 1] <- norm(cov_diff, type = "2")
    }
  }

  spacings <- spacings[-length(spacings)]
  min_spacing <- min(spacings)
  min_jump <- if (length(jump_magnitudes) > 0) min(jump_magnitudes) else NA
  max_jump <- if (length(jump_magnitudes) > 0) max(jump_magnitudes) else NA

  return(list(
    jump_magnitudes = jump_magnitudes,
    min_jump = min_jump,
    max_jump = max_jump,
    spacings = spacings,
    min_spacing = min_spacing
  ))
}
