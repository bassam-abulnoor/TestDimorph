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
#' @keywords internal

V_to_R <- function(V) {
  # Returns a correlation matrix from a variance-covariance matrix
  V <- as.matrix(V)
  d <- diag(1 / sqrt(diag(V)))
  R <- d %*% V %*% d
  return(R)
}
