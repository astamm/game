#' Constructor for Gaussian Mixture Data
#'
#' @param sample A data frame with 3 variables: \code{observation} which labels
#'   each observation, \code{component} which labels each mixture component and
#'   \code{mixing} which gives the weight of each component within each
#'   observation.
#' @param dictionary A data frame with 3 variables: \code{component} which
#'   labels each mixture component and \code{mean} and \code{precision} which
#'   respectively give the mean and the precision (inverse variance) of the
#'   Gaussian distribution of the corresponding component.
#'
#' @return An object of class \code{\link{gmd}}, which is effectively a list
#'   with 2 components: \code{sample} and \code{dictionary}, as described by the
#'   input arguments.
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
gmd <- function(sample, dictionary, compute_support = TRUE, compute_centring = TRUE) {
  stopifnot(is.data.frame(sample) && is.data.frame(dictionary))
  stopifnot(length(sample) == 3L && length(dictionary) == 3L)
  stopifnot(all(c("observation", "component", "mixing") %in% names(sample)))
  stopifnot(all(c("component", "mean", "precision") %in% names(dictionary)))
  stopifnot(all(dictionary$precision > 0))
  stopifnot(all(sample$mixing >= 0))
  sum_weights <- sample %>%
    dplyr::group_by(observation) %>%
    dplyr::summarise(sum_weights = sum(mixing)) %>%
    pull(sum_weights)
  stopifnot(all(abs(sum_weights - 1) < sqrt(.Machine$double.eps)))
  res <- list(
    sample = dplyr::filter(sample, mixing > 0),
    dictionary = dictionary
  )
  as_gmd(res, compute_support, compute_centring)
}

as_gmd <- function(x, compute_support = TRUE, compute_centring = TRUE) {
  if (compute_support)
    x$support <- support_gmd(x)
  if (compute_centring) {
    x$centring <- x %>%
      unfold_gmd() %>%
      dplyr::transmute(
        centring = purrr::map_dbl(data, ~ centring(
          .x$mean,
          .x$precision,
          .x$mixing,
          x$support$ref_mean,
          x$support$ref_precision
        ))
      ) %>%
      dplyr::pull(centring)
  }
  class(x) <- c("gmd", class(x))
  x
}

is_gmd <- function(x) {
  "gmd" %in% class(x)
}

unfold_gmd <- function(x) {
  res <- x$sample %>%
    dplyr::left_join(x$dictionary, by = "component") %>%
    tidyr::nest(-observation)
  if ("centring" %in% names(x))
    res <- dplyr::mutate(res, centring = x$centring)
  res
}

align_gmd <- function(x, warping = "shift", target_mean = 0, target_precision = 1) {
  x <- x %>%
    unfold_gmd() %>%
    dplyr::mutate(data = purrr::map(
      data,
      align_gmm,
      warping = warping,
      target_mean = target_mean,
      target_precision = target_precision
    )) %>%
    tidyr::unnest() %>%
    tidyr::unite(component, observation, component, remove = FALSE)
  samp <- dplyr::distinct(x, observation, component, mixing)
  dict <- dplyr::distinct(x, component, mean, precision)
  as_gmd(list(sample = samp, dictionary = dict))
}

support_gmd <- function(x, width = 5.85) {
  m <- x %>%
    unfold_gmd() %>%
    dplyr::mutate(m = purrr::map_dbl(data, ~ sum(.x$mixing * .x$mean))) %>%
    dplyr::pull(m) %>%
    mean()
  x <- x$dictionary %>%
    dplyr::mutate(
      aj = mean - width / sqrt(precision),
      bj = mean + width / sqrt(precision)
    )
  a <- min(x$aj)
  b <- max(x$bj)
  sl <- max(b - m, m - a)
  lb <- m - sl
  ub <- m + sl
  v <- (ub - lb)^2 / 12
  list(
    ref_mean = m, ref_precision = 1 / v,
    ref_lb = lb, ref_ub = ub
  )
}
