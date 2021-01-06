#' @title Visualization Of t-Greene Pairwise Comparisons
#' @description Returns a graphical or numerical correlational matrix of p
#' values for the interpopulation degree of sexual dimorphism as measured by
#' Greene t test
#' @return Graphical or numerical matrix of p values from Greene t test
#' pairwise comparisons. N.B: contrary to the usual corrplots where higher
#' values indicate stronger correlation, here lower values indicate
#' significance
#' @details Data is entered as a data frame of summary statistics where
#' the column containing population names is chosen by position (first by
#' default), other columns of summary data should have specific names (case
#' sensitive) similar to [baboon.parms_df]
#' @seealso
#'  [TestDimorph-deprecated()]
#' @name pMatrix-deprecated
#' @keywords internal
NULL
#' @rdname TestDimorph-deprecated
#' @section `pMatrix`:
#' For `pMatrix`, use [t_greene()]
#' @export
#' @importFrom stats pt
#' @importFrom corrplot corrplot
pMatrix <- function(x,
                    Pop = 1,
                    plot = FALSE,
                    padjust = p.adjust.methods,
                    CI = 0.95,
                    digits = 4,
                    alternative = c("two.sided", "less", "greater"),
                    ...) {
  if (!identical(Sys.getenv("TESTTHAT"), "true")) {
    .Deprecated("t_greene")
  }
  if (!(is.data.frame(x))) {
    stop("x should be a dataframe")
  }
  if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f") %in% names(x))) {
    stop(
      "colnames must contain:
            M.mu= Male mean
            F.mu=Female mean
            M.sdev=Male sd
            F.sdev=Female sd
            m= Male sample size
            f=Female sample size
            N.B: colnames are case sensitive"
    )
  }
  alternative <-
    match.arg(alternative, c("two.sided", "less", "greater"))
  padjust <- match.arg(padjust, p.adjust.methods)
  if (!(Pop %in% seq_along(x))) {
    stop("Pop should be number from 1 to ncol(x)")
  }
  if (!is.logical(plot)) {
    stop("pairwise should be either TRUE or FALSE")
  }
  x <- x %>%
    drop_na() %>%
    as.data.frame()
  x$Pop <- x[, Pop]
  x$Pop <- factor(x$Pop)
  levels(x$Pop) <- sort(levels(x$Pop))
  if (length(unique(x$Pop)) != length(x$Pop[which(!is.na(x$Pop))])) {
    stop("Population names are not unique")
  }

  mat <-
    matrix(
      data = NA,
      nrow = nlevels(x$Pop),
      ncol = nlevels(x$Pop),
      byrow = TRUE,
      dimnames = list(levels(x$Pop), levels(x$Pop))
    )

  for (i in seq(x$Pop)) {
    mat[i, ] <-
      t_test(
        m = x[, "m"],
        f = x[, "f"],
        M.mu = x[, "M.mu"],
        F.mu = x[, "F.mu"],
        M.sdev = x[, "M.sdev"],
        F.sdev = x[, "F.sdev"],
        m2 = x[i, "m"],
        f2 = x[i, "f"],
        M.mu2 = x[i, "M.mu"],
        F.mu2 = x[i, "F.mu"],
        M.sdev2 = x[i, "M.sdev"],
        F.sdev2 = x[i, "F.sdev"],
        N = ((nlevels(x$Pop)^2 - nlevels(x$Pop)) / 2),
        padjust = padjust,
        alternative = alternative,
        CI = CI,
        digits = digits,
        es = "none"
      )$p.value
  }


  if (isTRUE(plot)) {
    corrplot::corrplot(
      corr = mat,
      p.mat = mat,
      is.corr = FALSE,
      sig.level = 1 - CI,
      ...
    )
  } else {
    mat
  }
}
