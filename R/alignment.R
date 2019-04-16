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
align_gmm <- function(x, warping = "affine", target_mean = 0, target_precision = 1) {
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

  tibble::tibble(
    component = x$component,
    mean = means,
    precision = precisions,
    mixing = mixings
  )
}
