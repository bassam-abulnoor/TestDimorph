#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname to_params
<<<<<<< HEAD
#' @keywords internal
<<<<<<< HEAD
=======
#' @export
>>>>>>> 6dea5f88f050c15e544b721a245e69cf92b48861
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> 49bab1ae1b9ab0ccdc09b835509661bcaf96caeb
=======
#' @export
>>>>>>> 2a7c9b7ceb62f20bb9d8a0ffb70230424317dd2e

>>>>>>> 6dea5f88f050c15e544b721a245e69cf92b48861
>>>>>>> 3ffea2d74b01d4945f77061b801bd9783c34951e
<<<<<<< HEAD
>>>>>>> 49bab1ae1b9ab0ccdc09b835509661bcaf96caeb
=======
>>>>>>> 49bab1ae1b9ab0ccdc09b835509661bcaf96caeb
to_params <- function(x) {
  V <- x$V
  n.r <- NROW(V)
  Pops <- x$means[1:(n.r / 2), 1]
  Traits <- colnames(x$means)[-(1:3)]

  # Get means

  M.mu <- as.matrix(x$means[1:(n.r / 2), -(1:3)])
  row.names(M.mu) <- Pops
  colnames(M.mu) <- Traits
  F.mu <- as.matrix(x$means[(n.r / 2 + 1):n.r, -(1:3)])
  row.names(F.mu) <- Pops
  colnames(F.mu) <- Traits

  # Get sample sizes
  N <- x$means[, 3]
  m <- N[1:(n.r / 2)]
  names(m) <- Pops
  f <- N[(n.r / 2 + 1):n.r]
  names(f) <- Pops
  V.mat <- to_symm_diag(V[1, ])
  V.pooled <- V.mat * (N[1] - 1)
  n.t <- NCOL(V.mat)
  sdev <- matrix(NA, ncol = n.t, nrow = n.r)
  sdev[1, ] <- sqrt(diag(V.mat))
  for (i in 2:n.r) {
    V.mat <- to_symm_diag(V[i, ])
    sdev[i, ] <- sqrt(diag(V.mat))
    V.pooled <- V.pooled + (N[i] - 1) * V.mat
  }
  V.pooled <- V.pooled / (sum(N) - n.r)
  d <- diag(1 / sqrt(diag(V.pooled)))
  R.res <- round(V_to_R(V.pooled), 4)
  M.sdev <- round(sdev[1:(n.r / 2), ], 2)
  colnames(M.sdev) <- Traits
  row.names(M.sdev) <- Pops
  F.sdev <- round(sdev[(n.r / 2 + 1):n.r, ], 2)
  colnames(F.sdev) <- Traits
  row.names(F.sdev) <- Pops
  list(
    R.res = R.res,
    M.mu = M.mu,
    F.mu = F.mu,
    m = m,
    f = f,
    M.sdev = M.sdev,
    F.sdev = F.sdev
  )
}
