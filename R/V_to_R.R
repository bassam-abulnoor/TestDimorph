#' @title Correlation matrix extraction
#' @description Returns a correlation matrix from a variance covariance matrix
#' @param V Variance-covariance matrix
#' @return Correlation matrix
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname V_to_R
<<<<<<< HEAD
#' @keywords internal
=======
>>>>>>> 6dea5f88f050c15e544b721a245e69cf92b48861

V_to_R <- function(V) {
  # Returns a correlation matrix from a variance-covariance matrix
  V <- as.matrix(V)
  d <- diag(1 / sqrt(diag(V)))
  R <- d %*% V %*% d
  return(R)
}
