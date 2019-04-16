#' Distance Matrix for Data Sets
#'
#' This function computes and returns the distance matrix computed by using the
#' specified distance measure to compute the pairwise distances between data
#' points in a data set.
#'
#' @param x Either a numeric matrix, a data frame, a \code{\link[stats]{dist}}
#'   object, an \code{\link{sgd}} object or a \code{\link{gmd}} object.
#' @inheritParams stats::dist
#'
#' @inherit stats::dist return
#'
#' @examples
#' x <- sgd(
#'   c(mean =  0, precision = 1  ),
#'   c(mean =  3, precision = 0.5),
#'   c(mean = -1, precision = 2  )
#' )
#' dist(x)
#'
#' N <- 100
#' M <- 4
#' w <- matrix(runif(N * M), N, M)
#' w <- w / rowSums(w)
#' s <- tidyr::crossing(
#'   observation = paste0("O", 1:N),
#'   component = paste0("C", 1:M)
#' ) %>%
#' dplyr::mutate(mixing = as.numeric(t(w)))
#' d <- tibble::tibble(
#'   component = paste0("C", 1:M),
#'   mean = numeric(M),
#'   precision = 1:M
#' )
#' y <- gmd(s, d)
#' dist(y)
#' @export
dist <- function(x, ...) {
  UseMethod("dist", x)
}

#' @describeIn dist This is the \code{\link[stats]{dist}} function of the
#'   \pkg{stats} package. We refer the user to the corresponding documentation
#'   for more details on the available distances and examples.
#' @export
dist.default <- stats::dist

#' @describeIn dist Implementation for Single Gaussian Data (stored in objects
#'   of class \code{\link{sgd}}).
#' @export
dist.sgd <- function(x, ref_mean = 0, ref_precision = 1) {
  nu0 <- as_sgm(c(mean = ref_mean, precision = ref_precision))
  n <- nrow(x)
  d <- numeric(n * (n - 1) / 2)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      d[n*(i-1) - i*(i-1)/2 + j-i] <- .dist_sgm(x$data[[i]], x$data[[j]], ref_mean, ref_precision)
    }
  }
  attributes(d) <- NULL
  attr(d, "Labels") <- x$point
  attr(d, "Size") <- n
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}

#' @describeIn dist Implementation for Gaussian Mixture Data (stored in objects
#'   of class \code{\link{gmd}}).
#' @export
dist.gmd <- function(x, rule = 2, squared = TRUE) {
  x <- unfold_gmd(x)
  rule = ghRules[[rule]]
  d <- GetSquaredDistanceMatrix(x$data, rule$x, rule$w);
  if (!squared) d <- sqrt(d)
  attributes(d) <- NULL
  nm <- x$observation
  attr(d, "Labels") <- nm
  attr(d, "Size") <- nrow(x)
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}
