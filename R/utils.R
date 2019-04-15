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

frm_old <- function(x, max_iter = 10) {
  # Compute mean in B^2 (the mean is independent from the reference measure)
  xmean <- mean(x)
  # Find a Gaussian approx by mean of means and minimum of precisions
  tmp <- x %>%
    unfold_gmd() %>%
    dplyr::mutate(
      m = purrr::map_dbl(data, ~ sum(.x$mixing * .x$mean)),
      p = purrr::map2_dbl(data, m, ~ (sum(.x$mixing * (1 / .x$precision + .x$mean^2) - .y^2))^(-1))
    )
  ref_means <- mean(tmp$m)
  ref_precisions <- min(tmp$p)
  ref_mixings <- 1
  idx <- 0
  old_indices <- NULL
  for (i in 1:max_iter) {
    writeLines(paste("- Index of reference measure:", idx))
    # Get clr of xmean wrt ref
    xmean_ctr <- centring_f(
      f = xmean,
      ref_means = ref_means,
      ref_precisions = ref_precisions,
      ref_mixings = ref_mixings,
      already_log = TRUE
    )
    writeLines("  * Mean centring done.")
    xmean_clr <- clr_f(
      f = xmean,
      center = xmean_ctr,
      ref_means = ref_means,
      ref_precisions = ref_precisions,
      ref_mixings = ref_mixings,
      already_log = TRUE
    )
    writeLines("  * Mean clr done.")
    # Evaluate mean clr at nodes
    xmean_nodes <- evaluate_at_nodes(
      f = xmean_clr,
      ref_means = ref_means,
      ref_precisions = ref_precisions
    )
    writeLines("  * Mean nodes done.")
    # Compute distances to xmean in B^2(ref)
    centrings <- tmp$data %>%
      purrr::map(~ function(t) GetLogDensityWRTLebesgue(t, .x$mean, .x$precision, .x$mixing)) %>%
      purrr::map_dbl(
      centring_f,
      ref_means = ref_means,
      ref_precisions = ref_precisions,
      ref_mixings = ref_mixings,
      already_log = TRUE
    )
    writeLines("  * Centrings done.")
    clrs <- tmp$data %>%
      purrr::map(~ function(t) GetLogDensityWRTLebesgue(t, .x$mean, .x$precision, .x$mixing)) %>%
      purrr::map2(
        centrings,
        clr_f,
        ref_means = ref_means,
        ref_precisions = ref_precisions,
        ref_mixings = ref_mixings,
        already_log = TRUE
      )
    writeLines("  * Clrs done.")
    node_evals <- clrs %>%
      purrr::map(
        evaluate_at_nodes,
        ref_means = ref_means,
        ref_precisions = ref_precisions
      )
    writeLines("  * Node evaluations done.")
    distances <- purrr::map_dbl(
      node_evals,
      dist_new2,
      y = xmean_nodes,
      ref_means = ref_means,
      ref_precisions = ref_precisions,
      ref_mixings = ref_mixings
    )
    writeLines("  * Distances done.")
    if (sort(which.min(distances))[1] %in% old_indices) break()
    # Pick observed mixture with smallest distance as reference
    idx <- sort(which.min(distances))[1]
    old_indices <- c(old_indices, idx)
    ref_means <- tmp$data[[idx]]$mean
    ref_precisions <- tmp$data[[idx]]$precision
    ref_mixings <- tmp$data[[idx]]$mixing

  }
  list(
    index = idx,
    distances = distances,
    ref_means = ref_means,
    ref_precisions = ref_precisions,
    ref_mixings = ref_mixings,
    centrings = centrings,
    clrs = clrs,
    node_evals = node_evals, # Individual clrs evaluated at nodes for each Gaussian component of the reference Gaussian mixture
    xmean_nodes = xmean_nodes, # Mean clr evaluated at nodes for each Gaussian component of the reference Gaussian mixture
    xmean = xmean
  )
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
