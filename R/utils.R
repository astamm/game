.aggregate <- function(data, labels, fun) {
  k <- dplyr::n_distinct(labels)
  1:k %>%
    purrr::map(~ data[labels == .x]) %>%
    purrr::map_dbl(fun)
}

.find_mixings <- function(x, y, dictionary) {
  # Some input checks
  N <- length(x)
  stopifnot(length(y) == N)
  # Compute X dictionary matrix
  X <- tidyr::crossing(x, dictionary) %>%
    dplyr::mutate(
      value = purrr::pmap_dbl(
        list(x, mean, precision),
        ~ dnorm(x = ..1, mean = ..2, sd = 1 / sqrt(..3))
      )
    ) %>%
    dplyr::select(x, component, value) %>%
    tidyr::spread(component, value) %>%
    dplyr::select(-x) %>%
    as.matrix()
  # Solve non negative LS with sum-to-one constraint using the quadprog library
  D <- t(X) %*% X
  Rinv <- solve(chol(D))
  M <- nrow(dictionary)
  C <- cbind(rep(1, M), diag(M))
  b <- c(1, rep(0, M))
  d <- t(y) %*% X
  sol <- quadprog::solve.QP(
    Dmat = Rinv,
    dvec = d,
    Amat = C,
    bvec = b,
    meq = 1,
    factorized = TRUE
  )
  stopifnot(1 %in% sol$iact)
  mixings <- sol$solution
  mixings[sort(sol$iact)[-1] - 1] <- 0
  mixings <- mixings / sum(mixings)
  mixings
}

plot_gmd <- function(x, resolution = 1) {
  dict <- x$dictionary
  df <- x$sample %>%
    dplyr::left_join(x$dictionary, by = "component") %>%
    tidyr::nest(-observation)
  grid <- .dictionary_support(dict, resolution = resolution)
  y <- tidyr::crossing(grid, df) %>%
    dplyr::mutate(
      value = purrr::map2_dbl(grid, data, ~ {
        dgmm(.x, .y$mean, .y$precision, .y$mixing)
      })
    ) %>%
    dplyr::select(grid, observation, value) %>%
    tidyr::spread(observation, value) %>%
    dplyr::select(-grid) %>%
    as.matrix()
  matplot(grid, y, type = "l")
}

integrate_wrt_gaussian <- function(f, ref_mean, ref_precision, rule = ghRules[[20]]) {
  weights <- rule$w
  nodes <- ref_mean + sqrt(2) * rule$x / sqrt(ref_precision)
  sum(weights * f(nodes), na.rm = TRUE) / piPowerOneOverTwo
}

integrate_wrt_mixture <- function(f, ref_means, ref_precisions, ref_mixings, rule = ghRules[[20]]) {
  values <- purrr::map2_dbl(ref_means, ref_precisions, integrate_wrt_gaussian, f = f, rule = rule)
  sum(ref_mixings * values, na.rm = TRUE)
}

evaluate_at_nodes <- function(f, ref_means, ref_precisions, rule = ghRules[[20]]) {
  purrr::map2(ref_means, ref_precisions, ~ f(.x + sqrt(2) * rule$x / sqrt(.y)))
}

#' @export
aligned_dist <- function(f, g, method = "L-BFGS-B", rule = 2) {
  rule <- ghRules[[rule]]
  rprecisions <- c(f$precision, g$precision)
  rmixings <- c(f$mixing, g$mixing) / 2
  k <- length(rprecisions)
  workvec <- double(0)


  mf <- sum(f$mixing * f$mean)
  mg <- sum(g$mixing * g$mean)
  x0 <- mf - mg

  cost <- function(x) {
    gmeans <- g$mean + x[1]
    rmeans <- c(f$mean, gmeans)
    GetSquaredDistance(f$mean, f$precision, f$mixing, gmeans, g$precision, g$mixing, rmeans, rprecisions, rmixings, k, rule$x, rule$w, workvec)
  }
  # cost(x0)
  optim(par = x0, fn = cost, method = method)
}
