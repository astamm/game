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
kmeans.gmd <- function(x, k = 2, iter.max = 50L, d2 = NULL, method = "ward.D", rule = 2, shift = FALSE, avoid_mean_computation = FALSE) {
  if (shift) x <- align_gmd(x)
  xf <- unfold_gmd(x)
  n <- nrow(xf)
  stopifnot(k <= n)

  xc <- x
  if (shift) {
    mu <- get_mean_raw_moment(x = x, order = 1, rule = rule)
    xc$dictionary$mean <- x$dictionary$mean + mu
  }

  # if (is.null(d2) || shift) d2 <- dist(xc, rule = rule, squared = TRUE)

  if (avoid_mean_computation) {
    totss <- 1:n %>%
      purrr::map_dbl(dist_to_centroid, j = 1, d2 = d2, m = rep(1, n)) %>%
      sum()
  } else {
    totss <- xc %>%
      get_squared_distances_to_mean(cluster = 1, memberships = rep(1, n), rule = rule) %>%
      sum()
  }

  # Initialization
  writeLines("***** INITIALIZATION *****")
  cl <- initialize_kmeans(d2, k = k, method = method)
  memberships <- cl$memberships

  print(table(memberships))
  print(totss)

  # Now iterating
  stop_loop <- FALSE
  iter <- 0

  while (!stop_loop) {
    iter <- iter + 1

    writeLines(paste("***** ITERATION", iter, "*****"))

    if (avoid_mean_computation) {
      if (shift) {
        d2matrices <- 1:k %>%
          purrr::map(~ {
            mu <- get_mean_raw_moment(x = x[memberships == .x], order = 1, rule = rule)
            xc <- x
            xc$dictionary$mean <- x$dictionary$mean + mu
            dist(xc, rule = rule)
          })
        dsts <- 1:n %>%
          purrr::map(~ purrr::imap_dbl(d2matrices, dist_to_centroid, i = .x, m = memberships))
      } else {
        dsts <- 1:n %>%
          purrr::map(~ purrr::map_dbl(1:k, dist_to_centroid, i = .x, d2 = d2, m = memberships))
      }
    } else {
      dsts <- 1:k %>%
        purrr::map(~ {
          if (shift)
            mu <- get_mean_raw_moment(x = x[memberships == .x], order = 1, rule = rule)
          else
            mu <- 0
          xc <- x
          xc$dictionary$mean <- x$dictionary$mean + mu
          get_squared_distances_to_mean(
            x = xc,
            cluster = .x,
            memberships = memberships,
            rule = rule
          )}) %>%
        purrr::transpose() %>%
        purrr::simplify_all()
    }

    umemberships <- purrr::map_int(dsts, which.min)
    withinss <- dsts %>% purrr::map_dbl(min) %>% .aggregate(umemberships, sum)
    stop_loop <- (sum(umemberships - memberships) == 0) || (iter > iter.max)

    memberships <- umemberships

    print(table(memberships))
    print(sum(withinss))
  }

  # Wrap-up
  tot.withinss <- sum(withinss)
  centroids <- sort(unique(memberships)) %>%
    purrr::map(~ mean(x[memberships == .x], rule = rule))

  # Build output
  res <- list(
    cluster = memberships,
    centers = tibble::tibble(
      center = sort(unique(memberships)),
      data = centroids
    ),
    dictionary = x$dictionary,
    totss = totss,
    withinss = withinss,
    tot.withinss = tot.withinss,
    betweenss = totss - tot.withinss,
    size = table(memberships),
    iter = iter,
    ifault = as.integer(iter >= iter.max)
  )
  class(res) <- c("kmeans", class(res))
  res
}

# initialize_kmeans <- function(d2, k = 1, R = 1) {
#   n <- attr(d2, "Size")
#   init_best <- 0
#   wss_best <- 0
#   for (i in 1:R) {
#     init <- sample.int(n, size = 1)
#     while (length(init) < k) {
#       dsts <- apply(d2[, init], 1, min)
#       wss <- sum(dsts)
#       pr <- dsts / wss
#       init <- c(init, sample.int(n, size = 1, prob = pr))
#     }
#     if (wss < wss_best || i == 1) {
#       init_best <- init
#       wss_best <- wss
#     }
#   }
#   list(id = init_best, wss = wss_best)
# }

initialize_kmeans <- function(d2, k = 1, method = "ward.D") {
  if (method == "pam") {
    cl <- cluster::pam(x = sqrt(d2), k = k)
    return(list(
      id = cl$id.med,
      wss = sum(apply(d2[, cl$id.med], 1, min)),
      memberships = cl$clustering
    ))
  }
  if (method == "ward.D")
    cl <- hclust(d2, method = method)
  else
    cl <- hclust(sqrt(d2), method = method)
  memberships <- cutree(cl, k = k)
  med <- integer(k)
  for (i in 1:k) {
    idx <- which(memberships == i)
    dsts <- rowSums(d2[idx, idx])
    med[i] <- idx[which.min(dsts)]
  }
  med
  list(
    id = med,
    wss = sum(apply(d2[, med], 1, min)),
    memberships = memberships
  )
}
