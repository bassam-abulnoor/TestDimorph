#' @title anova_es
#' @param x ANOVA model
#' @inheritParams univariate
#' @return OUTPUT_DESCRIPTION
#' @keywords internal
#' @importFrom tibble as_tibble
anova_es <- function(x,
                     es = es,
                     digits = digits) {
  df1 <- summary(x)[[1]][1, 1]
  df2 <- summary(x)[[1]][2, 1]
  ssi <- summary(x)[[1]][1, 2]
  sse <- summary(x)[[1]][2, 2]
  within <- sse / df2
  between <- ssi / df1
  f <- summary(x)[[1]][[4]][1]
  p <- summary(x)[[1]][[5]][1]
  ss <- summary(x)[[1]][, 2]
  df <- summary(x)[[1]][, 1]
  ms <- summary(x)[[1]][, 3]
  eta <- ssi / (ssi + sse)
  omega <- (ssi - (df1 * within)) / (ssi + sse + within)
  cohen <- sqrt((eta) / (1 - eta))
  if (es == TRUE) {
    out <-
      cbind_fill(
        c("Sex", "Residuals"),
        round(df, 1),
        round(ss, digits),
        round(ms, digits),
        round(f, digits),
        round(p, digits),
        round(eta, digits),
        round(omega, digits),
        round(cohen, digits)

      )
    colnames(out) <-
      c(
        "term",
        "df",
        "sumsq",
        "meansq" ,
        "statistic",
        "p.value",
        "etasq",
        "omegasq",
        "cohen.f"
      )
  } else{
    out <-
      cbind_fill(
        c("Sex", "Residuals"),
        round(df, 1),
        round(ss, digits),
        round(ms, digits),
        round(f, digits),
        round(p, digits)

      )
    colnames(out) <-
      c("term", "df", "sumsq", "meansq" , "statistic", "p.value")


  }
  return(tibble::as_tibble(out))
}
