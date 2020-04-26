#' @title Sex-Specific One-way ANOVA From Summary statistics
#' @description Calculates sex specific one-way ANOVA from summary statistics.
#' @inheritParams Tg
#' @param pairwise Logical; if TRUE runs multiple pairwise comparisons on
#'   different populations using post hoc test of choice, Default: TRUE
#' @param method Type of post hoc test implemented by [PostHocTest], Default:
#'   'hsd'
#' @return Sex specific ANOVA tables and pairwise comparisons in tidy format.
#' @details Data is entered as a tibble/data frame of summary statistics where
#'   the column containing population names is chosen by position (first by
#'   default), other columns of summary data should have specific names (case
#'   sensitive) similar to [baboon.parms_df]
#' @examples
#'   # Comparisons of femur head diameter in four populations
#'   library(TestDimorph)
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
#' aovSS(x = df)
#' @seealso
#'  \code{\link[DescTools]{PostHocTest}}
#' @rdname aovSS
#' @export
#' @importFrom stats rnorm aov
#' @importFrom tibble rownames_to_column as_tibble is_tibble
#' @importFrom DescTools PostHocTest
#' @importFrom rlang abort
#' @importFrom multcompView multcompLetters
aovSS <-
    function(x,
             Pop = 1,
             pairwise = TRUE,
             letters = FALSE,
             es = FALSE,
             digits = 4,
             method = "hsd",
             sig.level = 0.05) {
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
        if (!(pairwise %in% c(TRUE, FALSE)))   {
            rlang::abort("pairwise should be either TRUE or FALSE")

        }
        if (!(letters %in% c(TRUE, FALSE)))   {
            rlang::abort("letters should be either TRUE or FALSE")

        }
        if (!(es %in% c(TRUE, FALSE)))   {
            rlang::abort("es should be either TRUE or FALSE")

        }
        if (sig.level < 0 ||
            sig.level > 1 || !is.numeric(sig.level))  {
            rlang::abort("sig.level should be a number between 0 and 1")

        }
        x <- data.frame(x)
        x$Pop <- x[, Pop]
        x$Pop <- factor(x$Pop)
        if (length(unique(x$Pop)) != length(x$Pop[which(!is.na(x$Pop))])) {
            rlang::abort("Populations names'must be unique")
        }
        aov_sum <- function(.mu, .sdev, n) {
            N <- length(.mu)
            Pop <- factor(rep(levels(x$Pop), n))
            df <- lapply(1:N, function(i) {
                scale(stats::rnorm(n[i])) * .sdev[i] + .mu[i]
            })
            x <- do.call(rbind, df)
            out <- data.frame(Pop, x)
            return(out)
        }
        Male <- aov_sum(x$M.mu, x$M.sdev, x$m)
        Female <- aov_sum(x$F.mu, x$F.sdev, x$f)
        av_M <- stats::aov(x ~ Pop, data = Male)
        av_F <- stats::aov(x ~ Pop, data = Female)
        M <-
            anova_es(x = av_M, es = es, digits = digits)
        F <-
            anova_es(x = av_F, es = es, digits = digits)

        M1 <-
            tibble::as_tibble(tibble::rownames_to_column(data.frame(
                DescTools::PostHocTest(av_M, method = method)[[1]]
            )))
        colnames(M1) <-
            c("populations",
              "mean.diff",
              "conf.low",
              "conf.high",
              "p.value")
        post_M1 <- M1$p.value
        names(post_M1) <- M1$populations
        post_M1 <-
            tibble::rownames_to_column(data.frame(
                "letters" = multcompView::multcompLetters(post_M1, threshold = sig.level)[[1]]
            ),
            var = "populations")
        F1 <-
            tibble::as_tibble(tibble::rownames_to_column(data.frame(
                DescTools::PostHocTest(av_F, method = method)[[1]]
            )))
        colnames(F1) <-
            c("populations",
              "mean.diff",
              "conf.low",
              "conf.high",
              "p.value")

        F1 <-
            tibble::as_tibble(tibble::rownames_to_column(data.frame(
                DescTools::PostHocTest(av_F, method = method)[[1]]
            )))
        colnames(F1) <-
            c("populations",
              "mean.diff",
              "conf.low",
              "conf.high",
              "p.value")
        post_F1 <- F1$p.value
        names(post_F1) <- F1$populations
        post_F1 <-
            tibble::rownames_to_column(data.frame(
                "letters" = multcompView::multcompLetters(post_F1, threshold = sig.level)[[1]]
            ),
            var = "populations")


        if (pairwise == TRUE) {
            if (letters == TRUE) {
                list(
                    "Male model" = M,
                    "Male posthoc" = M1,
                    "Male letters" = post_M1,
                    "Female model" = F,
                    "Female posthoc" = F1,
                    "Female letters" = post_F1
                )

            } else{
                list(
                    "Male model" = M,
                    "Male posthoc" = M1,
                    "Female model" = F,
                    "Female posthoc" = F1
                )

            }
        } else {
            out <- list(M, F)
            names(out) <- NULL
            names(out) <- c("Male model", "Female model")
            return(out)
        }
    }
