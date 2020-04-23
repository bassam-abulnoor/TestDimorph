#' @title Visualization Of t-Greene Pairwise Comparisons
#' @description Returns a graphical or numerical correlational matrix of
#'   p-values for the interpopulation degree of sexual dimorphism as measured by
#'   Greene t-test
#'@param x Tibble/data frame containing summary statistics, Default: NULL
#'@param Pop Number of the column containing populations' names, Default: 1
#'@param plot Logical; if TRUE graphical matrix of p-values, Default: TRUE
#'@param alternative a character string specifying the alternative hypothesis,
#'  must be one of "two.sided", "greater" or "less", Default: 'two.sided'
#'@param padjust Method of p.value adjustment for multiple comparisons following
#'  [p.adjust.methods], Default: 'none'
#' @param ... additional arguments that can be passed to
#'   [corrplot][corrplot::corrplot] or [Tg] functions.
#' @return Graphical or numerical matrix of p-values from Greene t-test pairwise
#'   comparisons.
#'   N.B: contrary to the usual corrplots where higher values indicate stronger
#'   correlation, here lower values indicate significance
#' @details Data is entered as a tibble/data frame of summary statistics where
#'   the column containing population names is chosen by position (first by
#'   default), other columns of summary data should have specific names (case
#'   sensitive) similar to [baboon.parms_df]
#' @seealso
#'  \code{\link[corrplot]{corrplot}}
#'  \code{\link{TestDimorph-deprecated}}
#' @name pMatrix-deprecated
#' @keywords internal
NULL
#' @rdname TestDimorph-deprecated
#' @section \code{pMatrix}:
#' For \code{pMatrix}, use \code{\link{Tg}}.
#' @export
#' @importFrom stats pt
#' @importFrom plyr adply
#' @importFrom tibble is_tibble
#' @importFrom corrplot corrplot
#' @importFrom rlang abort
#'
pMat <- function(x = NULL,
                    Pop = 1,
                    plot = FALSE,
                    padjust = "none",
                    alternative = "two.sided",
                    ...) {
    .Deprecated("Tg")
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
    if (!(plot %in% c(TRUE, FALSE)))   {
        rlang::abort("pairwise should be either TRUE or FALSE")

    }
    x <- data.frame(x)
    x$Pop <- x[, Pop]
    x$Pop <- factor(x$Pop)
    if (length(unique(x$Pop)) != length(x$Pop[which(!is.na(x$Pop))])) {
        rlang::abort("Population names must be unique")
    }

    pMat <- function(x = x, y = x) {
        Mat <- plyr::adply(
            .data = x,
            .margins = 1,
            .fun = function(x) {
                Tg(
                    m = x[1, "m"],
                    f = x[1, "f"],
                    M.mu = x[1, "M.mu"],
                    F.mu = x[1, "F.mu"],
                    M.sdev = x[1, "M.sdev"],
                    F.sdev = x[1, "F.sdev"],
                    m2 = y[, "m"],
                    f2 = y[,
                           "f"],
                    M.mu2 = y[, "M.mu"],
                    F.mu2 = y[, "F.mu"],
                    M.sdev2 = y[, "M.sdev"],
                    F.sdev2 = y[, "F.sdev"],
                    N = ((
                        nlevels(x$Pop) ^ 2 - nlevels(x$Pop)
                    ) / 2),
                    padjust = padjust,
                    alternative = alternative,
                    ...
                )$p.value
            }
        )

        Mat <- as.matrix(Mat[, -(1:ncol(x))])
        rownames(Mat) <- levels(x$Pop)
        colnames(Mat) <- levels(x$Pop)
        return(Mat)
    }
    if (plot == TRUE) {
        cor <- pMat(x, x)
        corrplot::corrplot(corr = cor,
                           p.mat = cor,
                           ...)

    } else{
        pMat(x, x)
    }
}
