#' Affine Alignment for Gaussian Mixture Models
#'
#' @param x A Gaussian mixture stored in an object of class \code{\link{gmm}}.
#' @param warping A character specifying the family of warping functions to
#'   consider, among \code{shift}, \code{dilation} and \code{affine} (default).
#'   Any other choice leads to no alignment
#' @param ref_mean The reference mean to align to when shift or affine warping functions are used (default: \code{0}).
#' @param ref_precision The reference precisions to align to when dilation or affine warping functions are used (default: \code{1}).
#' @param dictionary A data frame with the 3 variables \code{component}, \code{mean} and \code{precision} specifying a dictionary of Gaussian measures onto which projecting the resulting aligned Gaussian mixture. If not provided (default), this step is skipped.
#'
#' @return A Gaussian mixture in an object of class \code{\link{gmm}} aligned on the provided reference.
#' @export
#'
#' @examples
#' x <- gmm(rep(0, 4), 1:4, rep(0.25, 4))
#' dict <- data.frame(component = paste0("C", 1:4), mean = rep(0, 4), precision = 1:4)
#' align_gmm(x, dictionary = dict, warping = "shift")
align_gmm <- function(x, warping = "affine", target_mean = 0, target_precision = 1, dictionary = NULL) {
  if (!(warping %in% c("shift", "dilation", "affine"))) return(x)
  m <- sum(x$mixing * x$mean)
  a <- 1
  if (warping != "shift")
    a <- sqrt(target_precision * (sum(x$mixing * (1 / x$precision + x$mean^2)) - m^2))
  b <- 0
  if (warping != "dilation")
    b <- m - a * target_mean
  means <- (x$mean - b) / a
  precisions <- x$precision * a^2
  mixings <- x$mixing

  if (is.null(dictionary)) {
    res <- tibble::tibble(
      component = x$component,
      mean = means,
      precision = precisions,
      mixing = mixings
    )
    return(res)
  }

  M <- nrow(dictionary)
  grid <- .dictionary_support(dictionary, 1)
  # Compute X dictionary matrix
  X <- tidyr::crossing(grid, dictionary) %>%
    dplyr::mutate(
      value = purrr::pmap_dbl(
        list(grid, mean, precision),
        ~ dnorm(x = ..1, mean = ..2, sd = 1 / sqrt(..3))
      )
    ) %>%
    dplyr::select(grid, component, value) %>%
    tidyr::spread(component, value) %>%
    dplyr::select(-grid) %>%
    as.matrix()
  # Solve non negative LS with sum-to-one constraint using the quadprog library
  y <- dgmm(grid, means, precisions, mixings)
  means <- dictionary$mean
  precisions <- dictionary$precision
  D <- t(X) %*% X
  d <- t(y) %*% X
  Rinv <- solve(chol(D))
  if (warping == "shift") {
    C <- cbind(rep(1, M), means, diag(M))
    b <- c(1, ref_mean, rep(0, M))
  } else if (warping == "dilation") {
    C <- cbind(rep(1, M), 1 / precisions + means^2, diag(M))
    b <- c(1, 1 / ref_precision + ref_mean^2, rep(0, M))
  } else {
    C <- cbind(rep(1, M), means, 1 / precisions + means^2, diag(M))
    b <- c(1, ref_mean, 1 / ref_precision + ref_mean^2, rep(0, M))
  }
  eqcon <- 2 + as.numeric(warping == "affine")
  sol <- quadprog::solve.QP(
    Dmat = Rinv,
    dvec = d,
    Amat = C,
    bvec = b,
    meq = eqcon,
    factorized = TRUE
  )
  stopifnot(all(1:eqcon %in% sol$iact))
  mixings <- sol$solution
  mixings[sort(sol$iact)[-(1:eqcon)] - eqcon] <- 0
  mixings <- mixings / sum(mixings)
  gmm(mean = means, precision = precisions, mixing = mixings)
}
