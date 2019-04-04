#' K-Means Clustering
#'
#' This function performs k-means clustering of the data points in a data set.
#'
#' @param x A numeric matrix where each row is a data point or an object that
#'   can be coerced to such a matrix (such as a numeric vector or a data frame
#'   with all numeric columns), an \code{\link{sgd}} object or a
#'   \code{\link{gmd}} object.
#' @inheritParams stats::kmeans
#' @param k The number of clusters to look for (default: \code{2L}).
#'
#' @return An object of class \code{"kmeans"} which as a \code{print} and a
#'   \code{fitted} methods. It is a list with at least the following components:
#' \describe{
#'   \item{\code{cluster}}{A vector of integers (among \code{1:k}) indicating
#'   the cluster to which each point is allocated.}
#'   \item{\code{centers}}{A matrix of cluster centres.}
#'   \item{\code{totss}}{The total sum of squares.}
#'   \item{\code{withinss}}{Vector of within-cluster sum of squares, one
#'   component per cluster.}
#'   \item{\code{tot.withinss}}{Total within-cluster sum of squares.}
#'   \item{\code{betweenss}}{The between-cluster sum of squares.}
#'   \item{\code{size}}{The number of points in each cluster.}
#'   \item{\code{iter}}{The number of (outer) iterations.}
#'   \item{\code{ifault}}{integer: indicator of a possible algorithm problem –
#'   for experts.}
#' }
#'
#' @examples
#' x <- sgd(
#'   c(mean =  0, precision = 1  ),
#'   c(mean =  3, precision = 0.5),
#'   c(mean = -1, precision = 2  )
#' )
#' kmeans(x)
#'
#' N <- 100
#' M <- 4
#' w <- matrix(runif(N * M), N, M)
#' w <- w / rowSums(w)
#' samp <- tidyr::crossing(
#'   observation = paste0("O", 1:N),
#'   component = paste0("C", 1:M)
#' ) %>%
#' dplyr::mutate(mixing = as.numeric(t(w)))
#' dict <- tibble::tibble(
#'   component = paste0("C", 1:M),
#'   mean = numeric(M),
#'   precision = 1:M
#' )
#' x <- gmd(samp, dict)
#' kx <- kmeans(x)
#' @export
kmeans <- function(x, ...) {
  UseMethod("kmeans", x)
}

#' @describeIn kmeans This is the \code{\link[stats]{kmeans}} function of the
#'   \pkg{stats} package. We refer the user to the corresponding documentation
#'   for more details on the available algorithms and examples.
#' @export
kmeans.default <- stats::kmeans

#' @describeIn kmeans Implementation for Single Gaussian Data (stored in objects
#'   of class \code{\link{sgd}}).
#' @export
kmeans.sgd <- function(x, k = 2, iter.max = 50L) {
  n <- nrow(x)
  stopifnot(k <= n)
  totss <- sum(purrr::map_dbl(x$data, .dist_sgm, mean(x))^2)
  # Initialization
  d <- dist(x)
  h <- stats::hclust(d, method = "ward.D2")
  c <- stats::cutree(h, k)
  cweights <- .aggregate(x$weight, c, sum)
  x$weight <- x$weight / cweights[c]
  centroids <- 1:k %>%
    purrr::map(~ mean(x[c == .x, ]))
  stop_loop <- FALSE
  iter <- 0
  while (!stop_loop) {
    iter <- iter + 1
    # Membership assignment
    tmp <- x$data %>%
      purrr::map(function(data_point) {
        centroids %>%
          purrr::map_dbl(~ .dist_sgm(
            nu1 = data_point,
            nu2 = .x,
            m0 = .x$mean,
            t0 = .x$precision
          ))
      })
    c_up <- purrr::map_int(tmp, which.min)
    withinss <- purrr::map_dbl(tmp, min)^2 %>%
      .aggregate(c_up, sum)
    stop_loop <- (sum(c_up - c) == 0) || (iter > iter.max)

    # x$weight <- x$weight * cweights[c]
    # cweights <- .aggregate(x$weight, c_up, sum)
    # x$weight <- x$weight / cweights[c_up]

    c <- c_up
    # Centroid update
    centroids <- 1:k %>%
      purrr::map(~ mean(
        x[c == .x, ],
        ref_mean = 0, # centroids[[.x]][["mean"]],
        ref_precision = 1 # centroids[[.x]][["precision"]]
      ))
  }
  tot.withinss <- sum(withinss)
  res <- list(
    cluster = c,
    centers = tibble::tibble(
      center = 1:k,
      data = centroids,
      weight = cweights
    ),
    totss = totss,
    withinss = withinss,
    tot.withinss = tot.withinss,
    betweenss = totss - tot.withinss,
    size = table(c),
    iter = iter,
    ifault = as.integer(iter >= iter.max)
  )
  class(res) <- c("kmeans", class(res))
  res
}

#' @describeIn kmeans Implementation for Gaussian Mixture Data (stored in objects
#'   of class \code{\link{gmd}}).
#' @export
kmeans.gmd <- function(x, k = 2, iter.max = 50L, d = NULL, info = NULL, rule = ghRules[[20]], alpha = 0.5) {
  xf <- unfold_gmd(x)
  n <- nrow(xf)
  stopifnot(k <= n)

  if (is.null(info)) info <- frm(x)
  if (is.null(d)) d <- dist(x, info$index)

  totss <- 0 # sum(info$distances^2)

  # Initialization
  writeLines("***** INITIALIZATION *****")
  cl <- cluster::pam(x = d, k = k)
  c <- cl$clustering - 1
  # refs <- cl$id.med - 1
  refs <- rep(info$index - 1, k)
  print(table(c + 1))

  # Now iterating
  stop_loop <- FALSE
  iter <- 0

  while (!stop_loop) {
    iter <- iter + 1

    writeLines(paste("***** ITERATION", iter, "*****"))

    writeLines("- Assigning observations to clusters...")
    tmp <- PerformReassignment(c, xf$data, refs, rule$x, rule$w, alpha)
    c_up <- tmp$memberships
    distances <- tmp$distances

    withinss <- .aggregate(distances, c_up + 1, sum)
    stop_loop <- (sum(c_up - c) == 0) || (iter > iter.max)
    c <- c_up

    writeLines("- Updating reference measures...")
    # refs <- ComputeNewReferences(refs, xf$data, c, rule$x, rule$w);

    print(table(c + 1))
    print(sum(withinss))
  }

  # Wrap-up
  c <- c + 1
  tot.withinss <- sum(withinss)
  centroids <- sort(unique(c)) %>%
    purrr::map(~ mean(x[c == .x], log = TRUE, trim = alpha))

  # Build output
  res <- list(
    cluster = c,
    centers = tibble::tibble(
      center = sort(unique(c)),
      data = centroids
    ),
    dictionary = x$dictionary,
    totss = totss,
    withinss = withinss,
    tot.withinss = tot.withinss,
    betweenss = totss - tot.withinss,
    size = table(c),
    iter = iter,
    ifault = as.integer(iter >= iter.max)
  )
  class(res) <- c("kmeans", class(res))
  res
}
# kmeans.gmd <- function(x, k = 2, iter.max = 50L, d = NULL, info = NULL, rule = ghRules[[20]]) {
#   xf <- unfold_gmd(x)
#   n <- nrow(xf)
#   stopifnot(k <= n)
#
#   if (is.null(info)) info <- frm(x)
#   if (is.null(d)) d <- dist(x, info$index)
#
#   totss <- 0 # sum(info$distances^2)
#
#   # Initialization
#   writeLines("***** INITIALIZATION *****")
#   cl <- cluster::pam(x = d, k = k)
#   c <- cl$clustering
#   print(table(c))
#
#   # Now iterating
#   stop_loop <- FALSE
#   iter <- 0
#
#   while (!stop_loop) {
#     iter <- iter + 1
#
#     writeLines(paste("***** ITERATION", iter, "*****"))
#
#     tmp <- sort(unique(c)) %>%
#       purrr::map(~ GetSquaredDistancesToMean(xf$data, c == .x, rule$x, rule$w)) %>%
#       purrr::transpose() %>%
#       purrr::simplify_all()
#
#     c_up <- tmp %>%
#       purrr::map_int(which.min) %>%
#       forcats::as_factor() %>%
#       forcats::fct_relabel(~ as.character(seq_along(as.character(.)))) %>%
#       as.numeric()
#
#     withinss <- purrr::map_dbl(tmp, min)^2 %>%
#       .aggregate(c_up, sum)
#     stop_loop <- (sum(c_up - c) == 0) || (iter > iter.max)
#     c <- c_up
#
#     print(table(c))
#     print(sum(withinss))
#   }
#
#   # Wrap-up
#   tot.withinss <- sum(withinss)
#   centroids <- sort(unique(c)) %>%
#     purrr::map(~ function(t) {
#       GetMean(t, xf$data[c == .x])
#     })
#
#   # Build output
#   res <- list(
#     cluster = c,
#     centers = tibble::tibble(
#       center = sort(unique(c)),
#       data = centroids
#     ),
#     dictionary = x$dictionary,
#     totss = totss,
#     withinss = withinss,
#     tot.withinss = tot.withinss,
#     betweenss = totss - tot.withinss,
#     size = table(c),
#     iter = iter,
#     ifault = as.integer(iter >= iter.max)
#   )
#   class(res) <- c("kmeans", class(res))
#   res
# }
