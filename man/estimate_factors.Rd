% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_factor.R
\name{estimate_factors}
\alias{estimate_factors}
\title{Estimate Latent Factors via Principal Component Analysis}
\usage{
estimate_factors(X, r)
}
\arguments{
\item{X}{A matrix of dimension \eqn{p \times T}, where \eqn{p} is the number of observed variables
and \eqn{T} is the number of time points. Columns are assumed to represent time.}

\item{r}{An integer specifying the number of latent factors to extract.}
}
\value{
A list containing:
\describe{
\item{F_t}{A matrix of estimated factors of dimension \eqn{r \times T}, obtained by projecting the data onto the leading eigenvectors.}
\item{R}{A matrix of estimated factor loadings of dimension \eqn{p \times r}, scaled by \eqn{\sqrt{p}}.}
\item{V_X}{The matrix of leading \code{r} eigenvectors (dimension \eqn{p \times r}) of the sample covariance matrix.}
\item{M_X}{A diagonal matrix containing the top \code{r} eigenvalues of the sample covariance matrix.}
}
}
\description{
Performs a principal component decomposition of the covariance matrix of observed high-dimensional data
and returns estimated latent factors and loadings based on the leading \code{r} eigencomponents.
}
\details{
This function implements a standard principal component estimator for approximate factor models.
The sample covariance matrix of the transposed data \eqn{X} is decomposed via eigendecomposition,
and the first \code{r} components are retained to form factor estimates \code{F_t} and loadings \code{R}.
}
\examples{
\dontrun{
X <- matrix(rnorm(1000), nrow = 20)  # 20 variables, 50 time points
result <- estimate_factors(X, r = 3)
str(result)
}

}
