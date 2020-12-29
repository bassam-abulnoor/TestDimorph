#' @title Multivariate Analysis Of Sexual Dimorphism
#' @description Multivariate extension of Greene t test \link{t_greene}
#' @param x Data frame or list containing summary statistics for
#' multiple parameters measured in both sexes in two or more populations.
#' @param R.res Pooled within correlational matrix, Default: NULL
#' @param Trait Number of the column containing names of measured parameters,
#' Default: 1
#' @param Pop Number of the column containing populations' names, Default: 2
#' @param type_manova type of MANOVA test "I","II" or "III", Default:"II".
#' @param manova_test_statistic type of test statistic used either "W" for "Wilks","P"
#' for "Pillai", "HL" for "Hotelling-Lawley" or "R" for "Roy's largest root",
#' Default: "W".
#' @param interact_manova Logical; if TRUE calculates MANOVA for main effects,
#' Default: FALSE.
#' @param es_manova effect size either ,"eta" for eta squared, or "none"for
#' not reporting an effect size, Default:"none".
#' @param univariate Logical; if TRUE conducts multiple univariate analyses on
#' different parameters separately, Default: FALSE
#' @inheritParams univariate

#' @param ... Additional arguments that could be passed to \link{univariate}
#' @return MANOVA table. When the term is followed by `(E)` an exact f-value
#' is calculated.
#' @details Data can be entered either as a data frame of summary
#' statistics as in \link{baboon.parms_df}. In that case the pooled within
#' correlational matrix `R.res` should be entered as a separate argument as in
#' \link{baboon.parms_R}. Another acceptable format is a named list of
#' matrices containing different summary statistics as well as the correlational
#' matrix as in \link{baboon.parms_list}. By setting the option `univariate`
#' to `TRUE`, multiple `ANOVA`s can be run on each parameter independently with
#' the required p.value correction using \link[stats]{p.adjust.methods}.
#' @examples
#' # x is a data frame with separate correlational matrix
#' library(TestDimorph)
#' multivariate(baboon.parms_df, R.res = baboon.parms_R)
#' # x is a list with the correlational matrix included
#' library(TestDimorph)
#' multivariate(baboon.parms_list, univariate = TRUE, padjust = "bonferroni")
#' #reproduces results from Konigsberg (1991)
#' multivariate(baboon.parms_df, R.res = baboon.parms_R)[3,]
# multivariate(baboon.parms_df, R.res = baboon.parms_R,interact_manova=T)

#' @rdname multivariate
#' @export
#' @importFrom stats pf
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble is_tibble

multivariate <- function(x,
                         R.res = NULL,
                         Trait = 1,
                         Pop = 2,
                         type_manova = "II",
                         manova_test_statistic = "W",
                         interact_manova = TRUE,
                         es_manova = "none",
                         univariate = FALSE,
                         padjust = "none",
                         ...,
                         lower.tail = FALSE,
                         CI = 0.95,
                         digits = 4) {
  padjust <- match.arg(padjust, choices = p.adjust.methods)
  type_manova <- match.arg(type_manova, c("I", "II", "III"))
  manova_test_statistic <- match.arg(manova_test_statistic, c("W", "R", "P", "HL"))
  es_manova <-
    match.arg(es_manova, choices = c("none", "eta"))
  if (!(is.list(x) || is.data.frame(x))) {
    stop("x should be a list or a dataframe")
  }
  if (!is.logical(univariate)) {
    stop("univariate should be either TRUE or FALSE")
  }
  if (isFALSE(interact_manova) && type_manova == "III") {
    stop("main effects MANOVA is only available for types (I) and (II)")
  }
  if (is.data.frame(x)) {
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
    if (!(Trait %in% seq_along(x))) {
      stop("Trait should be number from 1 to ncol(x)")
    }
    if (!(Pop %in% seq_along(x))) {
      stop("Pop should be number from 1 to ncol(x)")
    }

    x <- data.frame(x)
    if (is.null(R.res)) {
      stop("R.res should be supplied if x is a dataframe")
    }
    if (!is.matrix(R.res)) {
      stop("R.res should be a matrix")
    }

    x <- dataframe2list(
      x = x,
      R.res = R.res,
      Trait = Trait,
      Pop = Pop
    )
  }
  if (is.list(x) && !(is.data.frame(x))) {
    if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f", "R.res") %in%
      names(x))) {
      stop(
        "List should have the following named matricies:
            M.mu= Male mean
            F.mu=Female mean
            M.sdev=Male sd
            F.sdev=Female sd
            m= Male sample size
            f=Female sample size
            R.res=Pooled within correlational matrix
            N.B: names are case sensitive"
      )
    }
    R <- x$R.res
    M <- x$M.mu
    F <- x$F.mu
    nM <- x$m
    nF <- x$f
    M.sd <- x$M.sdev
    F.sd <- x$F.sdev
  }

  if (isFALSE(interact_manova)) {
    out <- switch(type_manova,
      I = manova_main_I(
        x,
        es_manova,
        test = manova_test_statistic,
        digits, CI, lower.tail
      ),
      II = manova_main_II(
        x,
        es_manova,
        test = manova_test_statistic,
        digits, CI, lower.tail
      )
    )
  } else {
    out <- switch(type_manova,
      I = manova_I(
        x,
        es_manova,
        test = manova_test_statistic,
        digits, CI, lower.tail
      ),
      II = manova_II(
        x,
        es_manova,
        test = manova_test_statistic,
        digits, CI, lower.tail
      ),
      III = manova_III(
        x,
        es_manova,
        test = manova_test_statistic,
        digits, CI, lower.tail
      )
    )
  }
  out

  if (isTRUE(univariate)) {
    univariate_pairwise(
      x = x, out = out, padjust = padjust, digits = digits,
      lower.tail = lower.tail, CI = CI, ...
    )
  } else {
    out
  }
}
