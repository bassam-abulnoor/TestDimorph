#'@title Univariate Analysis Of Sexual Dimorphism
#'@description Calculation and visualization of the differences in degree sexual
#'  dimorphism between multiple populations using a modified one-way ANOVA and
#'  summary statistics as input
#'@inheritParams Tg
#'@param lower.tail Logical; if TRUE probabilities are \code{P[X <= x]},
#'  otherwise, \code{P[X > x]}., Default: FALSE
#'@param pairwise Logical; if TRUE runs multiple pairwise comparisons on
#'  different populations using [Tg] test, Default: FALSE
#'@param ... Additional arguments that could be passed to the [Tg] function
#'@return Tibble of ANOVA results
#'@details Data is entered as a tibble/data frame of summary statistics where
#'  the column containing population names is chosen by position (first by
#'  default), other columns of summary data should have specific names (case
#'  sensitive) similar to [baboon.parms_df]
#' @examples
#' # Comparisons of femur head diameter in four populations
#' library(TestDimorph)
#' m <- c(150.00, 82.00, 36.00, 34.00)
#' f <- c(150.00, 58.00, 34.00, 24.00)
#' M.mu <- c(49.39, 48.33, 46.99, 45.20)
#' F.mu <- c(42.91, 42.89, 42.44, 40.90)
#' M.sdev <- c(3.01, 2.53, 2.47, 2.00)
#' F.sdev <- c(2.90, 2.84, 2.26, 2.90)
#' df <-
#'cbind.data.frame(
#'    Pop = c('Turkish', 'Bulgarian', 'Greek', 'Portuguese '),
#'    m,
#'    f,
#'    M.mu,
#'    F.mu,
#'    M.sdev,
#'    F.sdev,
#'    stringsAsFactors = TRUE
#')
#' univariate(df, pairwise = TRUE, padjust = 'bonferroni')
#'@rdname univariate
#'@export
#'@importFrom stats pf
#'@importFrom tibble as_tibble is_tibble
#'@importFrom rlang abort
#'@references \insertRef{konigsberg1991historical}{TestDimorph}
#'  \insertRef{timonov2014study}{TestDimorph}
#'  \insertRef{curate2017sex}{TestDimorph}
#'  \insertRef{kranioti2009sex}{TestDimorph}
#'  \insertRef{gulhan2015new}{TestDimorph}
univariate <-     function(x,
                           Pop = 1,
                           es = FALSE,
                           pairwise = FALSE,
                           padjust = "none",
                           ...,
                           lower.tail = FALSE,
                           N = NULL,
                           digits = 4) {
    if (!(is.data.frame(x) || tibble::is_tibble(x))) {
        rlang::abort("x must be a tibble or a dataframe")
    }
    if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f") %in% names(x))) {
        rlang::abort(
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
    if (!(Pop %in% seq_along(1:ncol(x))))   {
        rlang::abort("Pop should be number from 1 to ncol(x)")

    }
    if (!(es %in% c(TRUE, FALSE)))   {
        rlang::abort("es should be either TRUE or FALSE")

    }
    if (!(pairwise %in% c(TRUE, FALSE)))   {
        rlang::abort("pairwise should be either TRUE or FALSE")

    }
    x <- data.frame(x)
    x$Pop <- x[, Pop]
    x$Pop <- factor(x$Pop)
    r <- nrow(x)
    n <- sum(x[, "m"], x[, "f"])
    df1 <- (r - 1)
    df2 <- (n - (2 * r))
    x[, "w"] <- (x[, "m"] * x[, "f"]) / (x[, "m"] + x[, "f"])
    x[, "d"] <- x[, "M.mu"] - x[, "F.mu"]
    sse <-
        sum((x[, "m"] - 1) * (x[, "M.sdev"] ^ 2) + ((x[, "f"] - 1) * (x[, "F.sdev"] ^
                                                                          2)))
    ssi <-
        sum(x[, "w"] * x[, "d"] ^ 2) - (sum((x[, "w"] * x[, "d"])) ^ 2 / sum(x[, "w"]))
    within <- sse / df2
    between <- ssi / df1
    f <- between / within
    if (is.null(N)) {
        p <- stats::pf(f, df1, df2, lower.tail = lower.tail)
    } else{
        p <-
            padjust_n(
                p = stats::pf(f, df1, df2, lower.tail = lower.tail),
                method =  padjust,
                n = N
            )
    }
    ss <- c(ssi, sse)
    df <- c(df1, df2)
    ms <- c(between, within)
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
    if (pairwise == TRUE) {
        out <-
            list(tibble::as_tibble(out),
                 Tg(
                     x,
                     Pop = Pop,
                     padjust = padjust,
                     es = es,
                     ...
                 ))
        names(out) <- c("univariate", "pairwise")
        return(out)
    } else {
        return(tibble::as_tibble(out))
    }

}
