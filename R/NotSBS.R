#' Change Point Detection by NotSBS
#'
#' This function detects change points in high-dimensional time series data using a Narrowest-over-threshold principle under seeded binary segmentation
#' (NotSBS) approach based on a standardised CUSUM statistic.
#'
#' @param x A numeric matrix of dimension \eqn{p \times T}, representing the observed data, where \eqn{p} is the number of series
#'          and \eqn{T} the number of time points.
#' @param type An optional numeric vector used to determine the number of change points if \code{m} is not explicitly provided.
#'             If \code{type = 0}, the method assumes no change point.
#' @param m An optional integer specifying the expected number of change points (used when \code{method = "oracle"}).
#'          If not supplied, it is inferred from \code{type}.
#' @param trim A positive integer specifying the trimming parameter. Only intervals of length greater than \code{2 * trim}
#'             are considered for change point detection.
#' @param threshold A numeric value specifying the fixed threshold. Must be provided when \code{method = "fixed"}.
#' @param method A character string indicating the selection method. One of \code{"fixed"} or \code{"oracle"}.
#' @param V.diag Logical. If \code{TRUE}, only the diagonal of the long-run covariance matrix is used for standardisation.
#' @param lrv Logical. If \code{TRUE}, long-run variance is estimated using HAC weights.
#' @param lbd Optional numeric value specifying the minimum interval length for seeding.
#'
#' @return A data frame containing the estimated change points and corresponding statistics. The columns include:
#' \describe{
#'   \item{est.cp}{Estimated location of the change point.}
#'   \item{val}{CUSUM statistic value at the detected point.}
#'   \item{st}{Start index of the interval.}
#'   \item{ed}{End index of the interval.}
#'   \item{trim}{Trimming parameter used.}
#' }
#' If \code{method = "oracle"}, additional columns are returned:
#' \describe{
#'   \item{selected_threshold}{Threshold value that produced exactly \code{m} change points.}
#'   \item{next_highest_threshold}{Next largest threshold.}
#'   \item{no_change_threshold}{Maximum threshold that leads to no change point detection.}
#' }
#'
#' @details
#' This function operates in two modes:
#' \itemize{
#'   \item \strong{Fixed}: Retains change points where the CUSUM statistic exceeds a user-specified threshold.
#'   \item \strong{Oracle}: Iteratively searches for the highest threshold yielding exactly \code{m} non-overlapping change points.
#' }
#' The underlying CUSUM statistics are computed via \code{\link{cusum.fts}}, and overlap between candidate points is removed
#' using interval exclusion.
#'
#' @seealso \code{\link{cusum.fts}}, \code{\link{find_single_cp}}
#'
#' @examples
#' \dontrun{
#' # Simulated matrix with time series signals
#' p <- 100; T <- 400
#' F_hat <- matrix(rnorm(p * T), nrow = p)
#'
#' # Fixed threshold method
#' results_fixed <- NotSBS(x = F_hat, m = 3, type = c(4,1,2), trim = 12,
#'                         threshold = 5, method = "fixed",
#'                         V.diag = TRUE, lrv = TRUE, lbd = 20)
#'
#' # Oracle method
#' results_oracle <- NotSBS(x = F_hat, m = 3, type = c(4,1,2), trim = 12,
#'                          method = "oracle",
#'                          V.diag = TRUE, lrv = TRUE, lbd = 20)
#' }
#'
#' @export

NotSBS <- function(x, r = NULL, type = NULL, m = 3, trim, threshold = NULL,
                   method = c("fixed", "oracle"),
                   V.diag = TRUE, lrv = TRUE, lbd = NULL) {

  method <- match.arg(method)
  Time <- ncol(x)

  if (is.null(trim)) {
    stop("Must supply 'trim'.")
  }

  results <- cusum.fts(x, V.diag = V.diag, lrv = lrv, trim = trim, lbd = lbd)
  results <- results[order(results$ed - results$st, decreasing = FALSE), ]

  # Method!
  if (method == "fixed") {

    if (is.null(threshold)) {
      stop("Must supply 'threshold'.")
    }

    results <- results[results$val >= threshold + 5e-3, ]

    cid <- 0
    selected_indices <- integer()
    while (cid < nrow(results)) {
      cid <- cid + 1
      selected_indices <- c(selected_indices, cid)

      if (cid < nrow(results)) {
        tmp <- results[(cid + 1):nrow(results), ]
        tmp <- tmp[!(tmp$st < results$est.cp[cid] & tmp$ed >= results$est.cp[cid]), ]
        results <- rbind(results[1:cid, ], tmp)
      }
    }

    results <- results[selected_indices, ]
    results <- results[order(results$est.cp), ]


  } else if (method == "oracle") {
    if (is.null(m)){
      if (identical(type, 0)) { m <- 0 } else { m <- length(type) }}

    thds <- sort(results$val, decreasing = TRUE)
    selected_threshold <- NA
    next_highest_threshold <- NULL
    no_change_threshold <- max(thds, na.rm = TRUE)

    selected_solution <- NULL

    for (i in seq_along(thds)) {
      thd <- thds[i]
      aux <- results[results$val >= thd, ]
      cid <- 0

      while (cid < m && cid < nrow(aux)) {
        cid <- cid + 1
        if (cid < nrow(aux)) {
          tmp <- aux[(cid + 1):nrow(aux), ]
          tmp <- tmp[!(tmp$st < aux$est.cp[cid] & tmp$ed >= aux$est.cp[cid]), ]
          aux <- rbind(aux[1:cid, ], tmp)
        }
      }

      if (cid == m) {
        selected_threshold <- thd
        selected_solution <- aux
        S <- selected_solution$est.cp

        # === Updated logic for next_highest_threshold ===
        for (j in (i + 1):length(thds)) {
          candidate_thd <- thds[j]
          candidate_intervals <- results[results$val == candidate_thd, ]

          for (k in seq_len(nrow(candidate_intervals))) {
            intvl <- candidate_intervals[k, ]
            st <- intvl$st
            ed <- intvl$ed
            if (all(S - st < log(Time) | ed - S < log(Time))) {
              next_highest_threshold <- candidate_thd
              break
            }
          }

          if (!is.null(next_highest_threshold)) break
        }


        results <- selected_solution
        break
      }
    }

    results$selected_threshold <- selected_threshold
    results$next_highest_threshold <- next_highest_threshold
    results$no_change_threshold <- no_change_threshold
    results <- results[order(results$est.cp), ]
  }

  rownames(results) <- NULL
  return(results)
}

