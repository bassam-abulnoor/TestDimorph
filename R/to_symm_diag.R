<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 3ffea2d74b01d4945f77061b801bd9783c34951e
>>>>>>> 49bab1ae1b9ab0ccdc09b835509661bcaf96caeb

#' @keywords internal
=======
>>>>>>> 6dea5f88f050c15e544b721a245e69cf92b48861
# Get pooled within group correlation matrix & standard deviations
to_symm_diag <- function(x) {
  # A utility function to take an upper triangular matrix stored by rows
  # and convert to a symmetric matrix For example, the matrix: [,1] [,2]
  # [,3] [1,] 3089 1079 785 [2,] 1079 574 351 [3,] 785 351 330 Should be
  # passed as the vector: c(3089, 1079, 785, 574, 351, 350)

  x <- as.numeric(x)
  k <- length(x)
  n.dim <- (sqrt(8 * k + 1) - 1) / 2
  if (n.dim %% 1 != 0) {
    return("Vector is wrong length for a triangular matrix")
  }
  sym.mat <- diag(n.dim)
  ip <- 0
  for (i in 1:n.dim) {
    for (j in i:n.dim)
    {
      ip <- ip + 1
      sym.mat[i, j] <- sym.mat[j, i] <- x[ip]
    }
  }
  return(sym.mat)
}
