#' Sample Mean for Gaussian Mixture Data
#'
#' @param x An object of class \code{\link{gmd}}.
#' @param native Boolean specyfing whether the resulting mean should be mapped
#'   back as a density (default: \code{TRUE}) or if its crl-transform should be
#'   return instead.
#'
#' @return The sample mean of the input Gaussian mixture data in density space
#'   (if \code{native==TRUE}) or in clr-space otherwise.
#' @export
#'
#' @examples
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
#' x <- gmd(s, d)
#' mean(x)
mean.gmd <- function(x, trim = 0, na.rm = TRUE, K = NULL, rule = 2) {
  x <- unfold_gmd(x)
  if (is.null(K)) K <- get_mean_raw_moment(x$data, 0, rule)
  function(t) GetMean(t, x$data, trim) / K
}

#' @export
get_mean_raw_moment <- function(x, order = 0, rule = 2) {
  x <- unfold_gmd(x)
  rule <- ghRules[[rule]]
  GetMeanRawMoment(x$data, order, rule$x, rule$w)
}
