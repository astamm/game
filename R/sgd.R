#' Constructor for Single Gaussian Data
#'
#' @param ... A number of named length-2 numeric vectors providing the
#'   \code{mean} and \code{precision} of each Gaussian distribution in the data
#'   set.
#'
#' @return An object of class \code{\link{sgd}} which is effectively a list of
#'   named length-2 numeric vectors providing the \code{mean} and
#'   \code{precision} of each Gaussian distribution in the data set.
#' @export
#'
#' @examples
#' x <- sgd(
#'   c(mean =  0, precision = 1  ),
#'   c(mean =  3, precision = 0.5),
#'   c(mean = -1, precision = 2  )
#' )
sgd <- function(...) {
  dots <- rlang::list2(...)
  as_sgd(dots)
}

as_sgd <- function(x) {
  UseMethod("as_sgd", x)
}

as_sgd.list <- function(x) {
  stopifnot(!any(purrr::map_lgl(x, ~ is.null(names(.x)))))
  if (all(purrr::map_lgl(x, ~ "weight" %in% names(.x)))) {
    w <- purrr::map_dbl(x, "weight")
    stopifnot(all(w >= 0) && sum(w) <= 1)
  } else {
    n <- length(x)
    w <- rep(1 / n, n)
  }
  x <- x %>%
    purrr::map(as_sgm) %>%
    .list2data() %>%
    dplyr::mutate(weight = w)
  class(x) <- c("sgd", class(x))
  x
}

#' @export
as_sgd.tbl_df <- function(x) {
  p <- length(x)
  n <- nrow(x)
  nm <- names(x)
  if (p == 2) { # Assuming it is point and data
    stopifnot(all(c("point", "data") %in% nm))
    stopifnot(all(x$data$precision > 0))
    x$weight <- rep(1 / n, n)
  } else if (p == 3) { # Assuming it is point, data and weight
    stopifnot(all(c("point", "data", "weight") %in% nm))
    stopifnot(all(x$data$precision > 0))
    print(paste("POS:", all(x$weight >= 0)))
    print(paste("SUM:", sum(x$weight) <= 1, sum(x$weight)))
    stopifnot(all(x$weight >= 0) && sum(x$weight) <= 1)
  } else if (p == 4) { # Assuming it is point, mean, precision and weight
    stopifnot(all(c("point", "mean", "precision", "weight") %in% nm))
    x <- tidyr::nest(x, mean, precision)
    stopifnot(all(x$data$precision > 0))
    print(paste("POS:", all(x$weight >= 0)))
    print(paste("SUM:", sum(x$weight) <= 1, sum(x$weight)))
    stopifnot(all(x$weight >= 0) && sum(x$weight) <= 1)
  }
  class(x) <- c("sgd", class(x))
  x
}

is_sgd <- function(x) {
  "sgd" %in% class(x)
}

#' Mean for Single Gaussian Data
#'
#' @param x An object of class \code{\link{sgd}}.
#' @param ref_mean A scalar giving the mean of the Gaussian distribution used as
#'   reference for the Bayes space (default: \code{0}).
#' @param ref_precision A scalar giving the precision (inverse of the variance)
#'   of the Gaussian distribution used as reference for the Bayes space
#'   (default: \code{1}). This should be a strictly positive value.
#'
#' @return A named length-2 numeric vector providing the \code{mean} and
#'   \code{precision} of the mean Gaussian distribution of the data set.
#' @export
#'
#' @examples
#' x <- sgd(
#'   c(mean =  0, precision = 1  ),
#'   c(mean =  3, precision = 0.5),
#'   c(mean = -1, precision = 2)
#' )
#' mean(x)
mean.sgd <- function(x, ref_mean = 0, ref_precision = 1, ...) {
  stopifnot(ref_precision > 0)
  n <- nrow(x)
  stopifnot(n > 0)
  means <- purrr::map_dbl(x$data, "mean")
  precisions <- purrr::map_dbl(x$data, "precision")
  weights <- x$weight
  sum_weights <- sum(weights)
  out_precision <- sum(weights * precisions) + (1 - sum_weights) * ref_precision
  out_mean <- (sum(weights * precisions * means) + (1 - sum_weights) * ref_precision * ref_mean) / out_precision
  print(x)
  print(out_mean)
  print(out_precision)
  as_sgm(c(mean = out_mean, precision = out_precision))
}

#' Plot for Single Gaussian Data
#'
#' @param x An object of class \code{\link{sgd}}.
#' @param scale An integer that defines the range of \code{x} values for
#'   plotting as a multiple of the maximum observed standard deviation (default:
#'   \code{3L}).
#' @param centered A boolean specifying whether individual densities should be
#'   centered prior to plotting (default: \code{FALSE}).
#'
#' @return Invisibly returns the input object.
#' @export
#'
#' @examples
#' x <- sgd(
#'   c(mean =  0, precision = 1  ),
#'   c(mean =  3, precision = 0.5),
#'   c(mean = -1, precision = 2  )
#' )
#' plot(x, centered = FALSE)
#' plot(x, centered = TRUE, xlab = "values")
plot.sgd <- function(x, scale = 3L, centered = FALSE, add_mean = FALSE) {
  max_sd <- 1 / sqrt(min(x$precision))
  if (centered) x$mean <- 0
  lb <- min(x$mean) - scale * max_sd
  ub <- max(x$mean) + scale * max_sd
  p <- ggplot2::ggplot(tibble::tibble(x = c(lb, ub)), ggplot2::aes(x = x)) +
    ggplot2::theme_bw()
  n <- nrow(x)
  cols <- rainbow(n)
  u <- seq(1e-4, 1 - 1e-4, length.out = 2^9 + 1)
  for (i in 1:n) {
    m <- x$mean[i]
    s <- 1 / sqrt(x$precision[i])
    c <- cols[i]
    t <- tibble::tibble(
      x = qnorm(p = u, mean = m, sd = s),
      y = dnorm(x = x, mean = m, sd = s)
    )
    p <- p + ggplot2::geom_line(data = t, ggplot2::aes(y = y), col = c)
  }
  if (add_mean) {
    x <- mean(x)
    p <- p +
      ggplot2::stat_function(
        fun = dnorm,
        args = list(mean = x[["mean"]], sd = 1 / sqrt(x[["precision"]])),
        col = "black",
        size = 1.2
      )
  }
  p
}
