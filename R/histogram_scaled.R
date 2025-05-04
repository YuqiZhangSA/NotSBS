#' Visualise Distribution of Estimated change points across Simulations
#'
#' Produces a histogram of scaled estimated change points across multiple simulations.
#'
#' @param detected_cp_list A list of numeric vectors, where each vector contains the estimated change points from a single simulation.
#' @param theta A numeric vector of true change point locations.
#' @param Time An integer specifying the time series length.
#' @param p An integer specifying the data dimension (not used internally but retained for consistency with other simulation functions).
#' @param num_simulations An integer giving the number of simulations (typically equal to \code{length(detected_cp_list)}).
#' @param method_name Optional character string to be included in the plot title, identifying the method used.
#'
#' @return No return value. A histogram is generated as a side effect.
#'
#' @details
#' The function visualises the empirical distribution of estimated change points after scaling to the \eqn{[0, 1]} interval.
#' The histogram uses bin width \eqn{2 \log(T)/T} to reflect a theoretical detection tolerance.
#'
#' True change point locations (scaled) are marked with red dotted vertical lines, and orange dashed lines denote the tolerance bands
#' of width \eqn{2 \log(T)/T} around each true change point.
#'
#' @examples
#' \dontrun{
#' # Simulated output from 100 repetitions
#' detected_list <- replicate(100, sort(sample(1:300, 3)), simplify = FALSE)
#' true_cp <- c(75, 150, 225)
#' visualise_stats(detected_list, theta = true_cp, Time = 300, p = 50,
#'                 num_simulations = 100, method_name = "Example Method")
#' }
#'
#' @export

visualise_stats <- function(detected_cp_list, theta, Time, p, num_simulations, method_name = "") {
  all_cps <- unlist(detected_cp_list)

  scaled_cps   <- all_cps / Time
  scaled_theta <- theta / Time

  hist(
    scaled_cps,
    breaks = seq(0, 1, by = 2*(log(Time)/Time)),
    main   = paste(
      "Distribution of Scaled Estimated Change Points\n",
      method_name
    ),
    xlab   = "Scaled Estimated Change Points",
    col    = "skyblue",
    border = "black",
    xlim   = c(0, 1),
    ylim   = c(0, num_simulations),
    las    = 1
  )

  abline(v = scaled_theta, col = "red", lty = "dotted", lwd = 2)
  for (cp in scaled_theta) {
    abline(v = cp - 2* log(Time)/Time, col = "orange", lty = "dashed")
    abline(v = cp + 2* log(Time)/Time, col = "orange", lty = "dashed")
  }


  grid(nx = NULL, ny = NULL, col = "grey70", lty = "dotted")

  #box()
}
