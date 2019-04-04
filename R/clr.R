centring <- function(means, precisions, mixings, ref_mean = 0, ref_precision = 1, rule = ghRules[[20]]) {
  # Restrict mixture to actual non-null components
  means <- means[mixings > 0]
  precisions <- precisions[mixings > 0]
  mixings <- mixings[mixings > 0]

  weights <- rule$w
  nodes <- ref_mean + sqrt(2) * rule$x / sqrt(ref_precision)
  sum(weights * GetLogDensityWRTGaussian(
    nodes,
    means,
    precisions,
    mixings,
    ref_mean,
    ref_precision
  ), na.rm = TRUE) / sqrt(pi)
}

centring_f <- function(f, ref_means = 0, ref_precisions = 1, ref_mixings = 1, already_log = FALSE, rule = ghRules[[20]]) {
  integrand <- function(t) {
    log_f <- if (already_log) f(t) else log(f(t))
    log_ref <- GetLogDensityWRTLebesgue(
      inputValues = t,
      meanValues = ref_means,
      precisionValues = ref_precisions,
      mixingValues = ref_mixings
    )
    log_f - log_ref
  }
  integrate_wrt_mixture(integrand, ref_means, ref_precisions, ref_mixings)
}

clr <- function(x, means, precisions, mixings, center, ref_mean = 0, ref_precision = 1) {
  # Restrict mixture to actual non-null components
  means <- means[mixings > 0]
  precisions <- precisions[mixings > 0]
  mixings <- mixings[mixings > 0]
  GetCenteredLogRatio(x, means, precisions, mixings, ref_mean, ref_precision, center)
}

clr_f <- function(f, center, ref_means = 0, ref_precisions = 1, ref_mixings = 1, already_log = FALSE) {
  function(t) {
    log_f <- if (already_log) f(t) else log(f(t))
    log_ref <- GetLogDensityWRTLebesgue(
      inputValues = t,
      meanValues = ref_means,
      precisionValues = ref_precisions,
      mixingValues = ref_mixings
    )
    log_f - log_ref - center
  }
}

clr_inv <- function(f, ref_mean = 0, ref_precision = 1, rule = ghRules[[20]]) {
  weights <- rule$w
  nodes <- ref_mean + sqrt(2) * rule$x / sqrt(ref_precision)
  x <- f(nodes)
  xmax <- max(x)
  K <- sum(weights * exp(x - xmax), na.rm = TRUE) / sqrt(pi)
  function(t) exp(f(t) - xmax + dnorm(x = t, mean = ref_mean, sd = 1 / sqrt(ref_precision), log = TRUE)) / K
}

dist_gmm <- function(x, y, cx, cy, ref_mean = 0, ref_precision = 1, rule = ghRules[[20]]) {
  isXgmm <- FALSE
  if (tibble::is_tibble(x)) {
    if (missing(cx))
      stop("The first input is a GMM model. The associated centring must be specified.")
    isXgmm <- TRUE
  }
  isYgmm <- FALSE
  if (tibble::is_tibble(y)) {
    if (missing(cy))
      stop("The second input is a GMM model. The associated centring must be specified.")
    isYgmm <- TRUE
  }
  integrand <- function(t) {
    clr1 <- if (isXgmm)
      clr(t, x$mean, x$precision, x$mixing, cx, ref_mean, ref_precision)
    else
      x(t)

    clr2 <- if (isYgmm)
      clr(t, y$mean, y$precision, y$mixing, cy, ref_mean, ref_precision)
    else
      y(t)

    (clr1 - clr2)^2
  }
  weights <- rule$w
  nodes <- ref_mean + sqrt(2) * rule$x / sqrt(ref_precision)
  sqrt(sum(weights * integrand(nodes), na.rm = TRUE)) / piPowerOneOverFour
}

dist_new <- function(x, y, ref_means = 0, ref_precisions = 1, ref_mixings = 1, rule = ghRules[[20]]) {
  x <- evaluate_at_nodes(x, ref_means, ref_precisions)
  values <- purrr::map2_dbl(x, y, ~ sum(rule$w * (.x - .y)^2, na.rm = TRUE))
  sqrt(sum(ref_mixings * values, na.rm = TRUE)) / piPowerOneOverFour
}

dist_new2 <- function(x, y, ref_means = 0, ref_precisions = 1, ref_mixings = 1, rule = ghRules[[20]]) {
  # x and y are clr evaluations at nodes for each Gaussian component within the reference Gaussian mixture

  # Compute the squared L2 distance between clr's w.r.t. each Gaussian component
  values <- purrr::map2_dbl(x, y, ~ sum(rule$w * (.x - .y)^2, na.rm = TRUE))

  # Aggregate with mixing weights of reference Gaussian mixture to get final distance
  sqrt(sum(ref_mixings * values, na.rm = TRUE)) / piPowerOneOverFour
}

norm_b2 <- function(x, ref_means = 0, ref_precisions = 1, ref_mixings = 1, rule = ghRules[[20]]) {
  # x is clr evaluations at nodes for each Gaussian component within the reference Gaussian mixture

  # Compute the squared L2 distance between clr's w.r.t. each Gaussian component
  values <- purrr::map_dbl(x, ~ sum(rule$w * .x^2, na.rm = TRUE))

  # Aggregate with mixing weights of reference Gaussian mixture to get final distance
  sqrt(sum(ref_mixings * values, na.rm = TRUE)) / piPowerOneOverFour
}
