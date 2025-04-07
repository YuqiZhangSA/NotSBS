#' Compute Simulation Statistics for change point Detection
#'
#' Evaluates accuracy and deviation metrics of estimated change points across multiple simulations,
#' comparing them against ground truth. Reports how often the number of detected change points matches
#' or deviates from the true number, as well as detection rate and localisation accuracy.
#'
#' @param Time An integer specifying the time series length.
#' @param theta A numeric vector of true change point locations.
#' @param detected_cp_list A list of numeric vectors, where each vector contains estimated change points for a simulation run.
#' @param trim A non-negative numeric value for the zone tolerance in \code{\link{align_change_points}}. Usually \code{trim = log(Time)}.
#' @param num_simulations An integer giving the number of simulation runs (should match \code{length(detected_cp_list)}).
#' @param method_name A character string for labelling the output (e.g., the name of the detection method).
#'
#' @return A list of numeric vectors (length = number of true change points), where each element contains
#' the aligned deviations between detected and true change points across simulations.
#'
#' @details
#' For each true change point, the function uses \code{\link{align_change_points}} to match the closest
#' estimated point (within a tolerance zone) from each simulation. It then:
#' \itemize{
#'   \item Computes the deviation from truth (signed difference).
#'   \item Calculates detection rate (proportion of simulations with a match).
#'   \item Computes accuracy: whether the deviation falls within \eqn{2 \log(T)}.
#'   \item Summarises frequency of detecting too few or too many change points.
#' }
#'
#' Summary statistics are printed to the console, including detection and accuracy rates per change point
#' and frequencies of \code{m_est - m} values (\code{0}, \code{±1}, \code{≥2}, \code{≤-2}).
#'
#' @seealso \code{\link{align_change_points}}
#'
#' @examples
#' \dontrun{
#' theta <- c(80, 160, 240)
#' cp_list <- replicate(100, sort(sample(1:300, 3)), simplify = FALSE)
#' compute_stats(Time = 300, theta = theta, detected_cp_list = cp_list,
#'               trim = log(300), num_simulations = 100, method_name = "TestMethod")
#' }
#'
#' @export

compute_stats <- function(Time, theta, detected_cp_list, trim = NULL, num_simulations, method_name) {
  m_true <- length(theta)
  deviations <- vector("list", m_true)
  accuracy_counts_per_cp <- numeric(m_true)
  m_est_diff <- integer(num_simulations)

  for (sim in 1:num_simulations) {
    detected_cp_sorted <- detected_cp_list[[sim]]
    aligned_deviations <- align_change_points(detected_cp_sorted, theta, Time, trim)
    for (i in 1:m_true) {
      deviations[[i]] <- c(deviations[[i]], aligned_deviations[i])
    }

    for (j in 1:m_true) {
      # Accuracy use 2*log(Time)
      found_closest <- !is.na(aligned_deviations[j]) && abs(aligned_deviations[j]) <= 2 * log(Time)
      if (found_closest) {
        accuracy_counts_per_cp[j] <- accuracy_counts_per_cp[j] + 1
      }
    }
  }

  m_est_diff <- sapply(detected_cp_list, length) - m_true
  freq_m_diff_0 <- sum(m_est_diff == 0)
  freq_m_diff_p1 <- sum(m_est_diff == 1)
  freq_m_diff_m1 <- sum(m_est_diff == -1)
  freq_m_diff_ge2 <- sum(m_est_diff >= 2)
  freq_m_diff_le2 <- sum(m_est_diff <= -2)

  mean_deviation <- numeric(m_true)
  median_deviation <- numeric(m_true)
  sd_deviation <- numeric(m_true)
  detection_rate <- numeric(m_true)

  for (i in 1:m_true) {
    deviations_cp <- deviations[[i]]
    valid_deviations <- deviations_cp[!is.na(deviations_cp)]
    detection_rate[i] <- sum(!is.na(deviations_cp)) / num_simulations
  }

  accuracy_per_cp <- accuracy_counts_per_cp / num_simulations

  cat("Method:", method_name, "\n")
  cat("Frequency of m_est - m:\n")
  cat("  Value = 0:", freq_m_diff_0 / num_simulations, "\n")
  cat("  Value = +1:", freq_m_diff_p1 / num_simulations, "\n")
  cat("  Value = -1:", freq_m_diff_m1 / num_simulations, "\n")
  cat("  Value >= 2:", freq_m_diff_ge2 / num_simulations, "\n")
  cat("  Value <= -2:", freq_m_diff_le2 / num_simulations, "\n\n")

  for (i in 1:m_true) {
    cat("Change Point", i, " (True CP =", theta[i], "):\n")
    cat("  Detection Rate:", round(detection_rate[i] * 100, 2), "%\n")
    cat("  Accuracy for CP", i, ":", round(accuracy_per_cp[i] * 100, 2), "%\n\n")
  }

  return(deviations)
}
