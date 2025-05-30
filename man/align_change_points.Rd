% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align_cp.R
\name{align_change_points}
\alias{align_change_points}
\title{Align Estimated Change Points with True Change Points}
\usage{
align_change_points(estimated_points, true_points, Time, tolerance = 50)
}
\arguments{
\item{estimated_points}{A numeric vector of estimated change point locations.}

\item{true_points}{A numeric vector of ground truth change point locations.}

\item{Time}{An integer specifying the total length of the time series.}

\item{tolerance}{A non-negative integer specifying the margin before the first and after the last true change point. Default is \code{50}.}
}
\value{
A numeric vector of length equal to the number of true change points. Each entry gives the signed deviation
between the matched estimated and true change point. If no estimated point is matched to a true one, \code{NA} is returned for that position.
}
\description{
Matches estimated change points to true ones within specified tolerance-based alignment zones
and computes signed deviations. Each estimated point can be matched to at most one true point.
}
\details{
The interval \code{[1, Time]} is partitioned into alignment zones based on the midpoints between adjacent true change points.
For each true change point, the algorithm searches for unmatched estimated points within the corresponding zone,
and selects the closest one (in absolute distance). This matching is one-to-one and greedy.
}
\examples{
\dontrun{
true_cp <- c(100, 200, 300)
est_cp <- c(95, 205, 290)
align_change_points(est_cp, true_cp, Time = 400)
}

}
