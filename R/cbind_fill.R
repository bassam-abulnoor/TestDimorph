#' @title cbind_fill
#' @description cbind columns with unequal length
#' @param ... columns with unequal length
#' @return a data frame
#' @rdname cbind_fill
#' @keywords internal
cbind_fill <- function(...) {
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x) {
    as.data.frame(rbind(x, matrix(data = NA, n - nrow(x), ncol(x))))
  }))
}
