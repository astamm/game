as_sgm <- function(x) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) %in% 2:3)
  stopifnot(!is.null(names(x)) && all(c("mean", "precision") %in% names(x)))
  stopifnot(x[["precision"]] > 0)
  x <- x %>%
    tibble::enframe() %>%
    tidyr::spread(name, value) %>%
    dplyr::select(mean, precision)
  class(x) <- c("sgm", class(x))
  x
}

.dist_sgm <- function(nu1, nu2, m0 = 0, t0 = 1) {
  m1 <- nu1[["mean"]]
  m2 <- nu2[["mean"]]
  d1 <- m1 - m0
  d2 <- m2 - m0
  t1 <- nu1[["precision"]]
  t2 <- nu2[["precision"]]
  sqrt(((t1 - t2) / t0)^2 / 2 + (t1 * d1 - t2 * d2)^2 / t0)
}
