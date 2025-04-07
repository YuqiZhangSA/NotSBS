#' Align Estimated change points with Ground Truth
#'
#' Matches estimated change points to true ones within specified tolerance-based alignment zones
#' and computes signed deviations. Each estimated point can be matched to at most one true point.
#'
#' @param estimated_points A numeric vector of estimated change point locations.
#' @param true_points A numeric vector of ground truth change point locations.
#' @param Time An integer specifying the total length of the time series.
#' @param tolerance A non-negative integer specifying the margin before the first and after the last true change point. Default is \code{50}.
#'
#' @return A numeric vector of length equal to the number of true change points. Each entry gives the signed deviation
#' between the matched estimated and true change point. If no estimated point is matched to a true one, \code{NA} is returned for that position.
#'
#' @details
#' The interval \code{[1, Time]} is partitioned into alignment zones based on the midpoints between adjacent true change points.
#' For each true change point, the algorithm searches for unmatched estimated points within the corresponding zone,
#' and selects the closest one (in absolute distance). This matching is one-to-one and greedy.
#'
#' @examples
#' \dontrun{
#' true_cp <- c(100, 200, 300)
#' est_cp <- c(95, 205, 290)
#' align_change_points(est_cp, true_cp, Time = 400)
#' }
#'
#' @export

align_change_points <- function(estimated_points, true_points, Time, tolerance = 50) {
  estimated_points <- sort(estimated_points)
  true_points <- sort(true_points)

  n_true <- length(true_points)
  zones <- c(tolerance, (true_points[-n_true] + true_points[-1]) / 2, Time - tolerance)

  result <- rep(NA, n_true)
  used <- rep(FALSE, length(estimated_points))

  for (i in seq_len(n_true)) {
    zone_start <- zones[i]
    zone_end <- zones[i + 1]

    in_zone <- which(estimated_points >= zone_start & estimated_points <= zone_end & !used)

    if (length(in_zone) > 0) {
      deviations <- abs(estimated_points[in_zone] - true_points[i])
      best_index <- in_zone[which.min(deviations)]

      result[i] <- estimated_points[best_index] - true_points[i]
      used[best_index] <- TRUE
    } else {
      result[i] <- NA
    }
  }

  return(result)
}
