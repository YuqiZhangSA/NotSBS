#' Multiple Change Point Estimation within Seeded Intervals
#'
#' This function detects multiple change points in the second moment structure of a high-dimensional time series matrix
#' using a standardised CUSUM statistic.
#'
#' @param x A matrix of dimension \eqn{p \times T}, representing the observed data with \eqn{p} variables and \eqn{T} dimensionality.
#' @param r An optional integer specifying the number of factors. If \code{NULL}, a default data-driven choice is used.
#' @param V.diag Logical. If \code{TRUE}, only the diagonal elements of the long-run covariance matrix \eqn{V} are used for standardisation.
#' @param lrv Logical. If \code{TRUE}, long-run variance estimation is used with HAC weights up to lag \code{n}.
#' @param n Optional bandwidth for HAC estimation. Default is \code{floor(T^{0.25})}.
#' @param trim A positive integer specifying the trimming parameter to avoid change point estimation near boundaries.
#' @param lbd An optional minimum segment length. If \code{NULL}, a default data-adaptive value is used.
#'
#' @return A data frame with one row per seeded interval and columns:
#' \describe{
#'   \item{est.cp}{Estimated location of the change point within the interval.}
#'   \item{val}{The associated CUSUM statistic value.}
#'   \item{st}{The interval start index.}
#'   \item{ed}{The interval end index.}
#'   \item{trim}{The trimming parameter used.}
#' }
#'
#' @seealso \code{\link{find_single_cp}} for the CUSUM computation in a single given interval.
#'
#' @examples
#' \dontrun{
#' # Simulate a data matrix with a change in factor covariance
#' p <- 100; T <- 400
#' x <- matrix(rnorm(p * T), nrow = p)
#' cusum.fts(x, trim = 12)
#' }
#'
#' @export



cusum.fts <- function(x, r = NULL,
                      V.diag = TRUE, lrv = TRUE,
                      n = NULL,
                      trim, lbd = NULL) {
  p <- nrow(x)
  Time <- ncol(x)
  if (is.null(r)) {r <- median(abc.factor.number(x)$r[4:6])}
  if (is.null(lbd)) {lbd <- round(dim(x)[2]^(max(2/5, 1 - min(1, log(dim(x)[1])/log(dim(x)[2])))) * log(dim(x)[2])^1.1)}
  if (is.null(n)) {n <- floor(dim(x)[2]^0.25)}
  #if (is.null(trim)) trim <- 2 * round(log(Time))
  d <- r * (r + 1) / 2

  #sv <- svd(x, nu = 0, nv = r); f_hat <- t(sv$v) * sqrt(Time)
  sv <- svd(x, nu =r, nv = 0); f_hat <- (t(sv$u) %*% x / sqrt(p))
  # sv <- svd(x, nu = r, nv = 0); svals <- sv$d[1:r]; f_hat <- sweep(t(sv$u) %*% x, 1, svals, "/") * sqrt(Time)

  ind <- lower.tri(diag(r), diag = TRUE)
  covF <- f_hat %*% t(f_hat) / Time
  mean_vec <- covF[ind]
  ff <- matrix(0, nrow = d, ncol = Time)
  II <- diag(1, r)[ind]

  for (t in seq_len(Time)) {
    ftft <- f_hat[,t] %*% t(f_hat[,t])
    #ff[,t] <- (ftft[ind]-II)
    ff[,t] <- ftft[ind] - mean_vec
  }

  V <- ff %*% t(ff) / Time # Gamma(0)
  if (lrv && n >= 1) {
    for (ell in 1:n) {
      tmp <- ff[, 1:(Time-ell), drop=FALSE] %*% t(ff[, (1:(Time-ell))+ell, drop=FALSE]) / Time
      w   <- 1 - ell/(n+1)
      V   <- V + w*(tmp + t(tmp))
    }
  }

  flag <- FALSE
  if (V.diag) {
    dd <- diag(V)
    if (min(dd) <= 0) {
      # replace negative or zero by smallest positive
      posmin <- min(dd[dd>0])
      dd     <- pmax(dd, posmin)
      flag   <- TRUE
    }
    Vhalf_inv <- diag(1/sqrt(dd))
  } else {
    eig <- eigen(V, symmetric = TRUE)
    dd  <- eig$values
    if (min(dd) <= 0) {
      posmin <- min(dd[dd>0])
      dd     <- pmax(dd, posmin)
      flag   <- TRUE
    }
    Vhalf_inv <- diag(1/sqrt(dd), d) %*% t(eig$vectors)
  }

  D <- Vhalf_inv %*% ff    # d x Time

  intervals <- seeded_intervals(Time, minl = lbd)
  results <- data.frame()
  for (i in seq_len(nrow(intervals))) {
    st <- intervals$st[i]
    ed <- intervals$ed[i]
    if ((ed - st) > 2 * trim) {
      res <- find_single_cp(D, st, ed, trim)
      results <- rbind(results, res)
    }
  }

  return(results)
}





find_single_cp <- function(D, st, ed, trim) {
  Ttot <- ncol(D)

  st_i <- max(1, as.integer(st))
  ed_i <- min(Ttot, as.integer(ed))
  len  <- ed_i - st_i + 1
  stopifnot(len > 2 * trim + 1)

  stat <- numeric(len)
  total <- rowSums(D[, st_i:ed_i, drop = FALSE])
  lsum  <- numeric(nrow(D))


  for (k in seq_len(len - 1)) {
    lsum <- lsum + D[, st_i + k - 1]
    rsum <- total - lsum
    nL   <- k
    nR   <- len - k
    diff <- rsum/nR - lsum/nL
    stat[k] <- sqrt((nL * nR) / len) * sqrt(sum(diff^2))
  }

  stat[seq_len(trim)]         <- 0
  stat[(len - trim + 1):len] <- 0

  best_rel <- which.max(stat)
  best_abs_i <- st_i + best_rel - 1

  data.frame(
    est.cp = best_abs_i,
    val    = stat[best_rel],
    st     = st,
    ed     = ed,
    trim   = trim
  )
}

