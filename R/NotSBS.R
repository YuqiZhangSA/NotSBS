#' Change Point Detection via NotSBS
#'
#' This function performs change point detection on a time series by utilising a matrix of estimated signals,
#' \code{F_hat}, and applying either a fixed or an oracle method for threshold selection.
#'
#' @param F_hat A numeric matrix of signals where each column corresponds to a time point.
#' @param m For the oracle method, an integer specifying the expected number of change points.
#'          Ignored if \code{method} is not \code{"oracle"}. Note: if \code{m} is not provided,
#'          internal logic may adjust it based on the supplied \code{type} (if available).
#' @param beta A numeric value intended to inform the seeding of intervals; currently not actively used
#'             in the interval selection.
#' @param trim A positive numeric value specifying the trimming length; only intervals with length greater
#'             than \code{2 * trim} will be considered.
#' @param threshold A numeric threshold for the fixed method. Must be supplied if \code{method = "fixed"}.
#' @param method A character string specifying the threshold selection method. Choices are \code{"fixed"},
#'               \code{"oracle"} and \code{"none"}.
#' @param V_shap (Optional) A character string indicating the form of the variance estimator to use.
#'               Choices are \code{"diag"}, \code{"full"}, or \code{"non"}. If not supplied, no variance
#'               correction is applied.
#' @param lbd A numeric value specifying the minimum length for the seeded intervals (often set as a multiple
#'              of \code{trim}).
#'
#' @return A data frame containing the detected change points along with associated statistics. For the oracle
#'         method, additional columns \code{selected_threshold}, \code{next_highest_threshold}, and
#'         \code{no_change_threshold} are appended.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming F_hat is a matrix with time series data as columns:
#' results_fixed <- NotSBS_not(F_hat = F_hat, m = 3, type = c(4,1,2), beta = 0.5, trim = 10,
#'                             threshold = 5, method = "fixed",
#'                             V_shap = "diag", lbd = 20)
#'
#' results_oracle <- NotSBS_not(F_hat = F_hat, m = 3, type = c(4,1,2), beta = 0.5, trim = 10,
#'                              method = "oracle",
#'                              V_shap = "diag", lbd = 20)
#' }

NotSBS_not <- function(F_hat, m = NULL, type = NULL, beta = NULL, trim = NULL, threshold = NULL,
                       method = c("fixed", "oracle", "none"),
                       V_shap = NULL, lbd = NULL) {
  method <- match.arg(method)

  Time <- ncol(F_hat)

  if (is.null(trim)) {
    stop("Must supply 'trim'.")
  }

  if (!is.null(V_shap)) {
    V_shap <- match.arg(V_shap, choices = c("diag", "full", "non"))
    long_run_V_est <- long_run_V(F_hat, Time)

    if (V_shap == "full") {
      V <- long_run_V_est$V_full
    } else if (V_shap == "diag") {
      V <- long_run_V_est$V_diag
    } else if (V_shap == "non") {
      V <- long_run_V_est$V_non
    }
  } else {
    V <- NULL
  }


  # seeded intervals with min length = 2*trim + 2
  #intervals <- seeded_intervals(Time, a = 0.5^(beta), minl = 2 * trim + 2)
  #intervals <- seeded_intervals(Time, minl = 2 * trim + 2)
  # perform best at (6 * trim + 2)
  intervals <- seeded_intervals(Time, minl = lbd)

  results <- data.frame()
  for (i in seq_len(nrow(intervals))) {
    st <- intervals$st[i]
    ed <- intervals$ed[i]

    if ((ed - st) > 2 * trim) {
      res <- find_single_cp_std(F_hat, st, ed, trim, V)
      results <- rbind(results, res)
    }
  }

  results <- results[results$est.cp > trim & results$est.cp < (Time - trim), ]
  results <- results[order(results$ed - results$st, decreasing = FALSE), ]

  #if (nrow(results) == 0) {
  # warning("No intervals are long enough; no candidate change points found.")
  #return(results)
  #}

  # Method!
  if (method == "fixed") {

    if (is.null(threshold)) {
      stop("Must supply 'threshold'.")
    }

    #results <- results[results$val >= threshold, ]
    # note that we change it to lager but not equal for "exceeding" the NH
    results <- results[results$val >= threshold + 5e-3, ]
    # Change ordering to be consistent with the oracle method:
    #results <- results[order(results$val, decreasing = TRUE), ]

    #if (nrow(results) == 0) {
    # warning("No intervals exceed the given threshold.")
    #return(results)
    #}

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
    if (is.null(m)) {
      if (identical(type, 0)) { m <- 0 } else { m <- length(type) }
    }

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
