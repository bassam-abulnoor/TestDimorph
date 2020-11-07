#' @title Univariate Analysis Of Sexual Dimorphism
#' @description Calculation and visualization of the differences in degree
#' sexual dimorphism between multiple populations using a modified one way
#' ANOVA and summary statistics as input
#' @inheritParams t_greene
#' @inheritParams t_test
#' @param lower.tail Logical; if TRUE probabilities are `P[X <= x]`,
#' otherwise, `P[X > x]`., Default: FALSE
#' @param pairwise Logical; if TRUE runs multiple pairwise comparisons on
#' different populations using [t_greene] test, Default: FALSE
#' @param ... Additional arguments that could be passed to the [t_greene]
#' function
#' @return Tibble of ANOVA results
#' @details Data is entered as a tibble/data frame of summary statistics where
#' the column containing population names is chosen by position (first by
#' default), other columns of summary data should have specific names (case
#' sensitive) similar to [baboon.parms_df]
#' @examples
#' \donttest{
#' # Comparisons of femur head diameter in four populations
#' library(TestDimorph)
#' df <-
#'   data.frame(
#'     Pop = c("Turkish", "Bulgarian", "Greek", "Portuguese "),
#'     m = c(150.00, 82.00, 36.00, 34.00),
#'     M.mu = c(49.39, 48.33, 46.99, 45.20),
#'     M.sdev = c(3.01, 2.53, 2.47, 2.00),
#'     f = c(150.00, 58.00, 34.00, 24.00),
#'     F.mu = c(42.91, 42.89, 42.44, 40.90),
#'     F.sdev = c(2.90, 2.84, 2.26, 2.90)
#'   )
#' univariate(df, pairwise = TRUE, padjust = "bonferroni")
#' }
#' @rdname univariate
#' @export
#' @importFrom stats pf
#' @importFrom tibble as_tibble
univariate <- function(x,
                       Pop = 1,
                       es = FALSE,
                       pairwise = FALSE,
                       padjust = p.adjust.methods,
                       ...,
                       lower.tail = FALSE,
                       N = NULL,
                       digits = 4) {
  # Data preparation --------------------------------------------------------

  if (!(is.data.frame(x))) {
    stop("x should be a dataframe")
  }
  if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f") %in% names(x))) {
    stop(
      "colnames should contain:
            M.mu= Male mean
            F.mu=Female mean
            M.sdev=Male sd
            F.sdev=Female sd
            m= Male sample size
            f=Female sample size
            N.B: colnames are case sensitive"
    )
  }
  if (!(Pop %in% seq_along(x))) {
    stop("Pop should be number from 1 to ncol(x)")
  }
  if (nrow(x) < 2) {
    stop("x should at least have 2 rows")
  }
  if (!is.logical(es)) {
    stop("es should be either TRUE or FALSE")
  }
  if (!is.logical(pairwise)) {
    stop("pairwise should be either TRUE or FALSE")
  }
  padjust <- match.arg(padjust, choices = p.adjust.methods)
  x <- data.frame(x)
  x$Pop <- x[, Pop]
  x$Pop <- factor(x$Pop)

  # univariate analysis -----------------------------------------------------

  r <- nrow(x)
  n <- sum(x$m, x$f)
  df1 <- (r - 1)
  df2 <- (n - (2 * r))
  x$w <- (x$m * x$f) / (x$m + x$f)
  x$d <- x$M.mu - x$F.mu
  sse <-
    sum((x$m - 1) * (x$M.sdev^2) + ((x$f - 1) * (x$F.sdev^
      2)))
  ssi <-
    sum(x$w * x$d^2) - (sum((x$w * x$d))^2 / sum(x$w))
  within <- sse / df2
  between <- ssi / df1
  f <- between / within
  if (is.null(N)) {
    p <- stats::pf(f, df1, df2, lower.tail = lower.tail)
  } else {
    p <-
      padjust_n(
        p = stats::pf(f, df1, df2, lower.tail = lower.tail),
        method = padjust,
        n = N
      )
  }
  signif <-
    case_when(
      p > 0.05 ~ "ns",
      p < 0.05 & p > 0.01 ~ "*",
      p < 0.01 & p > 0.001 ~ "**",
      p < 0.001 ~ "***"
    )
  ss <- c(ssi, sse)
  df <- c(df1, df2)
  ms <- c(between, within)
  eta <- ssi / (ssi + sse)
  omega <- (ssi - (df1 * within)) / (ssi + sse + within)
  cohen <- sqrt((eta) / (1 - eta))
  if (isTRUE(es)) {
    out <-
      cbind_fill(
        "term" = c("Sex", "Residuals"),
        "df" = round(df, 1),
        "sumsq" = round(ss, digits),
        "meansq" = round(ms, digits),
        "statistic" = round(f, digits),
        "p.value" = round(p, digits),
        "signif" = signif,
        "etasq" = round(eta, digits),
        "omegasq" = round(omega, digits),
        "cohen.f" = round(cohen, digits)
      )
  } else {
    out <-
      cbind_fill(
        "term" = c("Sex", "Residuals"),
        "df" = round(df, 1),
        "sumsq" = round(ss, digits),
        "meansq" = round(ms, digits),
        "statistic" = round(f, digits),
        "p.value" = round(p, digits),
        "signif" = signif
      )
  }

  # Pairwise comparisons ----------------------------------------------------

  if (isTRUE(pairwise)) {
    out <-
      list(
        univariate = tibble::as_tibble(out),
        pairwise = t_greene(
          x,
          Pop = Pop,
          padjust = padjust,
          es = es,
          ...
        )
      )
    out
  } else {
    tibble::as_tibble(out)
  }
}
