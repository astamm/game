#' Subset Operator for Gaussian Mixture Data
#'
#' @param x An object of class \code{\link{gmd}}.
#' @param i A valid expression giving instructions for subsetting observations.
#' @param j Unused.
#' @param ... Unused.
#' @param drop Unused.
#'
#' @return An object of the same class as the input object with only the subset
#'   of observations identified by the expression in argument \code{i}.
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
#' x[2]
#' x[2:4]
#' x[c(1, 3)]
#' x[-1]
#' x["O2"]
#' x[c("O1", "O3")]
`[.gmd` <- function(x, i, j, ..., drop = TRUE) {
  if (missing(i)) return(x)

  x$sample <- tidyr::nest(x$sample, -observation)
  obs <- x$sample$observation

  if (is.character(i))
    x$sample <- dplyr::filter(x$sample, observation %in% i)
  else if (is.logical(i))
    x$sample <- dplyr::filter(x$sample, i)
  else
    x$sample <- dplyr::slice(x$sample, i)

  if ("centring" %in% names(x)) {
    if (is.character(i))
      x$centring <- x$centring[obs %in% i]
    else
      x$centring <- x$centring[i]
  }

  x$sample <- tidyr::unnest(x$sample)

  x
}
