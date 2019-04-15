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
mean.gmd <- function(x, trim = 0, na.rm = TRUE, log = TRUE, normalize = FALSE) {
#   tmp <- unfold_gmd(x)
#   out <- function(t) {
#     tmp$data %>%
#       purrr::map(~ GetLogDensityWRTLebesgue(
#         inputValues = t,
#         meanValues = .x$mean,
#         precisionValues = .x$precision,
#         mixingValues = .x$mixing
#       )) %>%
#       purrr::transpose() %>%
#       purrr::simplify_all() %>%
#       purrr::map_dbl(mean, trim = trim, na.rm = na.rm)
#   }
#   if (log) return(out)
#   function(t) {exp(out(t))}
# }
#
# mean2 <- function(x, trim = 0, na.rm = TRUE, log = TRUE, normalize = FALSE) {
  tmp <- unfold_gmd(x)
  out <- function(t) GetMean(t, tmp$data, trim)
  if (log) return(out)
  f <- function(t) {exp(out(t))}
  K <- integrate(f, lower = -350, upper = 350)$value
  function(t) {f(t) / K}
}

#' @export
mean_expval <- function(x, rule = 2) {
  x <- unfold_gmd(x)
  rule <- ghRules[[rule]]
  GetMeanExpectedValue(x$data, rule$x, rule$w)
}

#' @export
mean_k <- function(x, rule = 2) {
  x <- unfold_gmd(x)
  rule <- ghRules[[rule]]
  GetMeanNormalizationFactor(x$data, rule$x, rule$w)
}
