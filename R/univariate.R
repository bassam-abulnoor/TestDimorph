#' @title Univariate Analysis Of Sexual Dimorphism
#' @description Calculation and visualization of the differences in degree
#' sexual dimorphism between multiple populations using a modified one way
#' ANOVA and summary statistics as input
#' @param type_anova type of ANOVA test "I","II" or "III", Default:"II".
#' @param interact_anova Logical; if TRUE calculates interaction effect,
#' Default: TRUE.
#' @param es_anova Type of effect size either "f" for f squared,"eta" for eta
#' squared or "none", Default:"none".
#' @inheritParams t_greene
#' @inheritParams t_test
#' @param lower.tail Logical; if TRUE probabilities are `P[X <= x]`,
#' otherwise, `P[X > x]`., Default: FALSE
#' @param pairwise Logical; if TRUE runs multiple pairwise comparisons on
#' different populations using \link{t_greene} Default: FALSE
#' @param ... Additional arguments that could be passed to the \link{t_greene}
#' function
#' @return  ANOVA tale.
#' @details Data is entered as a data frame of summary statistics where
#' the column containing population names is chosen by position (first by
#' default), other columns of summary data should have specific names (case
#' sensitive) similar to \link{baboon.parms_df}
#' @examples
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
#' @rdname univariate
#' @export
#' @importFrom stats pf
univariate <- function(x,
                       Pop = 1,
                       type_anova = "II",
                       interact_anova = TRUE,
                       es_anova = "none",
                       pairwise = FALSE,
                       padjust = "none",
                       ...,
                       lower.tail = FALSE,
                       CI = 0.95,
                       N = NULL,
                       digits = 4) {
  padjust <- match.arg(padjust, choices = p.adjust.methods)
  es_anova <-
    match.arg(es_anova, choices = c("none", "eta", "f"))

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
  if (!is.logical(pairwise)) {
    stop("pairwise should be either TRUE or FALSE")
  }
  if (CI < 0 ||
    CI > 1 || !is.numeric(CI)) {
    stop("CI should be a number between 0 and 1")
  }
  if (isFALSE(interact_anova) && type_anova == "III") {
    stop("main effects ANOVA is only available for types (I) and (II)")
  }

  x <- x %>%
    drop_na() %>%
    as.data.frame()
  x$Pop <- x[, Pop]
  x$Pop <- factor(x$Pop)

  if (isFALSE(interact_anova)) {
    out <- switch(type_anova,
      I = anova_main_I(
        x,
        es_anova,

        digits, CI, lower.tail
      ),
      II = anova_main_II(
        x,
        es_anova,

        digits, CI, lower.tail
      )
    )
  } else {
    out <- switch(type_anova,
      I = anova_I(
        x,
        es_anova,

        digits, CI, lower.tail
      ),
      II = anova_II(
        x,
        es_anova,
        digits, CI, lower.tail
      ),
      III = anova_III(
        x,
        es_anova,

        digits, CI, lower.tail
      )
    )
  }




  if (isTRUE(pairwise)) {
    out <-
      list(
        univariate = out,
        pairwise = TestDimorph::t_greene(
          x,
          Pop = Pop,
          padjust = padjust,
          ...
        )
      )
    out
  } else {
    out
  }
}
