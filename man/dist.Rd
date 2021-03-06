% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist.R
\name{dist}
\alias{dist}
\alias{dist.default}
\alias{dist.sgd}
\alias{dist.gmd}
\title{Distance Matrix for Data Sets}
\usage{
dist(x, ...)

\method{dist}{default}(x, method = "euclidean", diag = FALSE,
  upper = FALSE, p = 2)

\method{dist}{sgd}(x, ref_mean = 0, ref_precision = 1)

\method{dist}{gmd}(x, rule = 2, squared = TRUE)
}
\arguments{
\item{x}{Either a numeric matrix, a data frame, a \code{\link[stats]{dist}}
object, an \code{\link{sgd}} object or a \code{\link{gmd}} object.}

\item{method}{the distance measure to be used.  This must be one of
    \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
    \code{"canberra"}, \code{"binary"} or \code{"minkowski"}.
    Any unambiguous substring can be given.}

\item{diag}{logical value indicating whether the diagonal of the
    distance matrix should be printed by \code{print.dist}.}

\item{upper}{logical value indicating whether the upper triangle of the
    distance matrix should be printed by \code{print.dist}.}

\item{p}{The power of the Minkowski distance.}
}
\value{
\code{dist} returns an object of class \code{"dist"}.

  The lower triangle of the distance matrix stored by columns in a
  vector, say \code{do}. If \code{n} is the number of
  observations, i.e., \code{n <- attr(do, "Size")}, then
  for \eqn{i < j \le n}, the dissimilarity between (row) i and j is
  \code{do[n*(i-1) - i*(i-1)/2 + j-i]}.
  The length of the vector is \eqn{n*(n-1)/2}, i.e., of order \eqn{n^2}.

  The object has the following attributes (besides \code{"class"} equal
  to \code{"dist"}):
  \item{Size}{integer, the number of observations in the dataset.}
  \item{Labels}{optionally, contains the labels, if any, of the
    observations of the dataset.}
  \item{Diag, Upper}{logicals corresponding to the arguments \code{diag}
    and \code{upper} above, specifying how the object should be printed.}
  \item{call}{optionally, the \code{\link{call}} used to create the
    object.}
  \item{method}{optionally, the distance method used; resulting from
    \code{\link{dist}()}, the (\code{\link{match.arg}()}ed) \code{method}
    argument.}
}
\description{
This function computes and returns the distance matrix computed by using the
specified distance measure to compute the pairwise distances between data
points in a data set.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: This is the \code{\link[stats]{dist}} function of the
\pkg{stats} package. We refer the user to the corresponding documentation
for more details on the available distances and examples.

\item \code{sgd}: Implementation for Single Gaussian Data (stored in objects
of class \code{\link{sgd}}).

\item \code{gmd}: Implementation for Gaussian Mixture Data (stored in objects
of class \code{\link{gmd}}).
}}

\examples{
x <- sgd(
  c(mean =  0, precision = 1  ),
  c(mean =  3, precision = 0.5),
  c(mean = -1, precision = 2  )
)
dist(x)

N <- 100
M <- 4
w <- matrix(runif(N * M), N, M)
w <- w / rowSums(w)
s <- tidyr::crossing(
  observation = paste0("O", 1:N),
  component = paste0("C", 1:M)
) \%>\%
dplyr::mutate(mixing = as.numeric(t(w)))
d <- tibble::tibble(
  component = paste0("C", 1:M),
  mean = numeric(M),
  precision = 1:M
)
y <- gmd(s, d)
dist(y)
}
