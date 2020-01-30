#'@title Greene t-test of Sexual Dimorphism
#'@description Calculation and visualization of the differences in degree sexual
#'  dimorphism between two populations using summary statistics as input.
#'@param x Tibble/data frame containing summary statistics, Default: NULL
#'@param Pop Number of the column containing populations' names, Default: 1
#'@param es Logical; if TRUE effect size is included in the output , Default:
#'  FALSE
#'@param plot Logical; if TRUE graphical matrix of p-values, Default: TRUE
#'@param ... additional arguments that can be passed to
#'  [corrplot][corrplot::corrplot] function.
#'@param alternative a character string specifying the alternative hypothesis,
#'  must be one of "two.sided", "greater" or "less", Default: 'two.sided'
#'@param padjust Method of p.value adjustment for multiple comparisons following
#'  [p.adjust.methods], Default: 'none'
#'@param letters Logical; if TRUE returns letters for pairwise comparisons where
#'  significantly different populations are given different letters, Default:
#'  FALSE'
#'@param digits Number of significant digits, Default: 4
#'@param sig.level Critical p.value, Default: 0.05
#'@param N Number of pairwise comparisons for [p.adjust.methods], if left `NULL`
#'  it will follow the formula `n(n-1)/2` where `n` is the number of populations
#'  , Default: NULL
#'@param m Number of male sample size in the first population, Default: NULL
#'@param m2 Number of male sample size in the second population, Default: NULL
#'@param f Number of female sample size in the first population, Default: NULL
#'@param f2 Number of female sample size in the second population, Default: NULL
#'@param M.mu Means for males in the first population, Default: NULL
#'@param M.mu2 Means for males in the second population, Default: NULL
#'@param F.mu Means for females in the first population, Default: NULL
#'@param F.mu2 Means for females in the second population, Default: NULL
#'@param M.sdev Standard deviation for males in the first population, Default:
#'  NULL
#'@param M.sdev2 Standard deviation for males in the second population, Default:
#'  NULL
#'@param F.sdev Standard deviation for females in the first population, Default:
#'  NULL
#'@param F.sdev2 Standard deviation for females in the second population,
#'  Default: NULL
#'@return Tibble of t.test results
#'@details Summary statistics can be entered directly as arguments in case of
#'  comparing two populations or as a tibble/data frame of summary statistics
#'  where the column containing population names is chosen by position (first by
#'  default), other columns of summary data should have specific names (case
#'  sensitive) similar to [baboon.parms_df]
#' @examples
#' #Comparisons of femur head diameter in four populations
#' library(TestDimorph)
#' m <- c(150.00, 82.00, 36.00, 34.00)
#' f <- c(150.00, 58.00, 34.00, 24.00)
#' M.mu <- c(49.39, 48.33, 46.99, 45.20)
#' F.mu <- c(42.91, 42.89, 42.44, 40.90)
#' M.sdev <- c(3.01, 2.53, 2.47, 2.00)
#' F.sdev <- c(2.90, 2.84, 2.26, 2.90)
#' df <- cbind.data.frame(
#'   Pop = c('Turkish', 'Bulgarian', 'Greek', 'Portuguese '),
#'   m,
#'   f,
#'   M.mu,
#'   F.mu,
#'   M.sdev,
#'   F.sdev,
#'   stringsAsFactors = TRUE
#' )
#' Tg(
#'    df,
#'    plot = TRUE,
#'    method = 'ellipse',
#'   type = 'lower',
#'    col = c(
#'        '#AEB6E5',
#'        '#B1A0DB',
#'        '#B788CD',
#'        '#BC6EB9',
#'        '#BC569E',
#'        '#B6407D',
#'        '#A93154'
#'    ),
#'    tl.cex = 0.8,
#'    tl.col = 'black',
#'   insig =
#'        'label_sig',
#'    tl.srt = 0.1,
#'    pch.cex = 2.5,
#'    tl.pos = 'ld',
#'    win.asp = 1,
#'    number.cex = 0.5,
#'    na.label = 'NA'
#')
#'@seealso
#'\code{\link[multcompView]{multcompLetters}}
#'\code{\link[corrplot]{corrplot}}
#'@rdname Tg
#'@export
#'@importFrom stats qt pt
#'@importFrom utils combn
#'@importFrom purrr map map_dfr
#'@importFrom tibble as_tibble is_tibble
#'@importFrom multcompView multcompLetters vec2mat
#'@importFrom corrplot corrplot
#'@importFrom rlang abort
#'@references \insertRef{greene1989comparison}{TestDimorph}
#'  \insertRef{timonov2014study}{TestDimorph}
#'  \insertRef{gulhan2015new}{TestDimorph}
Tg <-     function(x = NULL,
                   Pop = 1,
                   es = FALSE,
                   plot = FALSE,
                   ...,
                   alternative = "two.sided",
                   padjust = "none",
                   letters = FALSE,
                   digits = 4,
                   sig.level = 0.05,
                   N = NULL,
                   m = NULL,
                   m2 = NULL,
                   f = NULL,
                   f2 = NULL,
                   M.mu = NULL,
                   M.mu2 = NULL,
                   F.mu = NULL,
                   F.mu2 = NULL,
                   M.sdev = NULL,
                   M.sdev2 = NULL,
                   F.sdev = NULL,
                   F.sdev2 = NULL) {
    t <-
        function(m,
                 f,
                 m2,
                 f2,
                 M.mu,
                 F.mu,
                 M.mu2,
                 F.mu2,
                 M.sdev,
                 F.sdev,
                 M.sdev2,
                 F.sdev2,
                 padjust = padjust) {
            tg <-
                ((M.mu - F.mu) - (M.mu2 - F.mu2)) / (sqrt(((((m - 1) * M.sdev ^ 2) + ((f -
                                                                                           1) * F.sdev ^ 2) + ((m2 - 1) * M.sdev2 ^ 2) + ((f2 - 1) * F.sdev2 ^
                                                                                                                                              2)
                )) / (m + f +
                          m2 + f2 - 4)) * sqrt((1 / m) + (1 / f) + (1 / m2) + (1 / f2)))

            df <- (m + f + m2 + f2 - 4)
            sdp <-
                (sqrt(((((m - 1) * M.sdev ^ 2) + ((f - 1) * F.sdev ^ 2) + ((m2 - 1) * M.sdev2 ^
                                                                               2) + ((f2 - 1) * F.sdev2 ^ 2)
                )) / (m + f + m2 + f2 - 4)))
            mean_diff <- ((M.mu - F.mu) - (M.mu2 - F.mu2))
            d <- abs(mean_diff / sdp)

            if (!(alternative %in% c("less", "greater", "two.sided")))   {
                rlang::abort("alternative should be one of `less`,`greater` or `two.sided`")

            }
            if (alternative == "less") {
                upper <-
                    mean_diff + (abs(stats::qt(
                        p = (1 - sig.level), df = df
                    )) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
                lower <-
                    mean_diff - (abs(stats::qt(
                        p = (1 - sig.level), df = df
                    )) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
                p <-
                    (stats::pt(abs(tg), df, lower.tail = TRUE))
            }
            if (alternative == "greater") {
                upper <-
                    mean_diff + (abs(stats::qt(
                        p = (1 - sig.level), df = df
                    )) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
                lower <-
                    mean_diff - (abs(stats::qt(
                        p = (1 - sig.level), df = df
                    )) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
                p <-
                    (stats::pt(abs(tg), df, lower.tail = FALSE))
            }
            if (alternative == "two.sided") {
                upper <-
                    mean_diff + (abs(stats::qt(
                        p = (1 - (sig.level / 2)), df = df
                    )) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
                lower <-
                    mean_diff - (abs(stats::qt(
                        p = (1 - (sig.level / 2)), df = df
                    )) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
                p <-
                    (2 * stats::pt(abs(tg), df, lower.tail = FALSE))
            }
            if (!is.null(x)) {
                if (is.null(N)) {
                    p <-
                        padjust_n(p = p,
                                  method = padjust,
                                  n = ((
                                      nlevels(x$Pop) ^ 2 - nlevels(x$Pop)
                                  ) / 2))
                } else{
                    p <-
                        padjust_n(p = p,
                                  method = padjust,
                                  n = N)
                }
            } else{
                p <-
                    padjust_n(p = p,
                              method = padjust,
                              n = N)


            }
            if (!(es %in% c(TRUE, FALSE)))   {
                rlang::abort("es should be either TRUE or FALSE")

            }
            if (es == TRUE) {
                return(
                    data.frame(
                        "df" = round(df, digits),
                        "mean.diff" = round(mean_diff, digits),
                        "conf.low" = round(lower, digits),
                        "conf.high" = round(upper, digits),
                        "statistic" = round(tg, digits),
                        "p.value" = round(p, digits),
                        "cohen.d" = round(d, digits)
                    )
                )

            } else{
                return(
                    data.frame(
                        "df" = round(df, digits),
                        "mean.diff" = round(mean_diff, digits),
                        "conf.low" = round(lower, digits),
                        "conf.high" = round(upper, digits),
                        "statistic" = round(tg, digits),
                        "p.value" = round(p, digits)
                    )
                )


            }

        }
    if (is.null(x)) {
        return(
            t(
                m,
                f,
                m2,
                f2,
                M.mu,
                F.mu,
                M.mu2,
                F.mu2,
                M.sdev,
                F.sdev,
                M.sdev2,
                F.sdev2,
                padjust = padjust
            )
        )
    } else {
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
        x <- data.frame(x)
        x$Pop <- x[, Pop]
        x$Pop <- factor(x$Pop)
        if (length(unique(x$Pop)) != length(x$Pop[which(!is.na(x$Pop))])) {
            rlang::warn("Population names are not unique")
        }
        pairs <- utils::combn(x$Pop, 2, simplify = FALSE)

        names(pairs) <- purrr::map(pairs, paste, collapse = "-")

        tg <-
            purrr::map_dfr(pairs, .id = "populations", function(y) {
                t(
                    m = x[y[1], "m"],
                    f = x[y[1], "f"],
                    m2 = x[y[2], "m"],
                    f2 = x[y[2], "f"],
                    M.mu = x[y[1], "M.mu"],
                    F.mu = x[y[1], "F.mu"],
                    M.mu2 = x[y[2], "M.mu"],
                    F.mu2 = x[y[2], "F.mu"],
                    M.sdev = x[y[1], "M.sdev"],
                    F.sdev = x[y[1],
                               "F.sdev"],
                    M.sdev2 = x[y[2], "M.sdev"],
                    F.sdev2 = x[y[2], "F.sdev"],
                    padjust = padjust
                )
            })
        pval <- tg$p.value
        names(pval) <- tg$populations
        pmatrix <- multcompView::vec2mat(pval)
        if (!(letters %in% c(TRUE, FALSE)))   {
            rlang::abort("letters should be either TRUE or FALSE")

        }
        if (letters == TRUE) {
            tg <-
                list(
                    "t.test" = tibble::as_tibble(tg),
                    "pairwise letters" = tibble::rownames_to_column(
                        data.frame(
                            "letters" = multcompView::multcompLetters(pval
                                                                      ,
                                                                      threshold = sig.level)[[1]]
                        ),
                        var = "populations"
                    )
                )

        } else{
            tg <- tibble::as_tibble(tg)
        }
        if (!(plot %in% c(TRUE, FALSE)))   {
            rlang::abort("plot should be either TRUE or FALSE")

        }
        if (plot == TRUE) {
            list(
                "p.matrix" =
                    corrplot::corrplot(
                        corr = pmatrix,
                        p.mat = pmatrix,
                        sig.level = sig.level,
                        number.digits = digits,
                        ...
                    ),
                "t.greene" = tg
            )
        } else{
            return(tg)

        }
    }
}
