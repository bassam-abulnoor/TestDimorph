#' @title Multivariate Analysis Of Sexual Dimorphism
#' @description Multivariate extension of Greene t-test [Tg]
#' @param x Tibble/Data frame or list containing summary statistics for multiple
#'   parameters measured in both sexes in two or more populations.
#' @param R.res Pooled within correlational matrix, Default: NULL
#' @param Parms Number of the column containing names of measured parameters,
#'   Default: 1
#' @param Pop Number of the column containing populations' names, Default: 2
#' @inheritParams univariate
#' @param univariate Logical; if TRUE conducts multiple univariate analyses on
#'   different parameters separately, Default: FALSE
#' @param ... Additional arguments that could be passed to the [univariate]
#'   function
#' @return Tibble of MANOVA results
#' @details Data can be entered either as a tibble/data frame of summary
#'   statistics as in [baboon.parms_df] . In that case the pooled within
#'   correlational matrix `R.res` should be entered as a separate argument as in
#'   [R]. Another acceptable format is a named list of matrices containing
#'   different summary statistics as well as the correlational matrix as in
#'   [baboon.parms_list]. By setting the option `univariate` to `TRUE`, multiple
#'   `ANOVA`s can be run on each parameter independently with the required p
#'   value correction using [p.adjust.methods].
#' @examples
#'  # x is a data frame with separate correlational matrix
#'  library(TestDimorph)
#'  multivariate(baboon.parms_df, R.res = R)
#'  # x is a list with the correlational matrix included
#'  library(TestDimorph)
#'  multivariate(baboon.parms_list, univariate = TRUE, padjust = 'bonferroni')
#' @rdname multivariate
#' @export
#' @importFrom stats pf
#' @importFrom reshape2 melt
#' @importFrom tibble as_tibble is_tibble
#' @importFrom rlang abort
#' @references \insertRef{konigsberg1991historical}{TestDimorph}
multivariate <-
    function(x,
             R.res = NULL,
             Parms = 1,
             Pop = 2,
             es = FALSE,
             univariate = FALSE,
             padjust = "none",
             ...,
             lower.tail = FALSE,
             digits = 4) {
        if (!(is.list(x) || is.data.frame(x) || tibble::is_tibble(x))) {
            rlang::abort("x must be a list a dataframe or a tibble")
        }
        if (!(es %in% c(TRUE, FALSE)))   {
            rlang::abort("es should be either TRUE or FALSE")

        }
        if (!(univariate %in% c(TRUE, FALSE)))   {
            rlang::abort("univariate should be either TRUE or FALSE")

        }
        if (is.data.frame(x) || tibble::is_tibble(x)) {
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
            if (!(Parms %in% seq_along(1:ncol(x))))   {
                rlang::abort("Parms should be number from 1 to ncol(x)")

            }
            if (!(Pop %in% seq_along(1:ncol(x))))   {
                rlang::abort("Pop should be number from 1 to ncol(x)")

            }
            x <- data.frame(x)
            if (is.null(R.res)) {
                rlang::abort("R.res must be supplied if x is a dataframe or a tibble")
            }
            if (!is.matrix(R.res)) {
                rlang::abort("R.res must be a matrix")

            }
            R <- R.res
            i <- list(levels(x[, Pop]), levels(x[, Parms]))
            M <-
                matrix(
                    data = x$M.mu,
                    nrow = nlevels(x[, Pop]),
                    ncol = nlevels(x[, Parms]),
                    dimnames = i
                )
            F <-
                matrix(
                    data = x$F.mu,
                    nrow = nlevels(x[, Pop]),
                    ncol = nlevels(x[, Parms]),
                    dimnames = i
                )
            nM <- as.vector(x$m)[1:nlevels(x[, Pop])]
            nF <- as.vector(x$f)[1:nlevels(x[, Pop])]
            nM <- nM[!is.na(nM)]
            nF <- nF[!is.na(nF)]
            M.sd <-
                matrix(
                    data = x$M.sdev,
                    nrow = nlevels(x[, Pop]),
                    ncol = nlevels(x[,
                                     Parms]),
                    dimnames = i
                )
            F.sd <-
                matrix(
                    data = x$F.sdev,
                    nrow = nlevels(x[, Pop]),
                    ncol = nlevels(x[,
                                     Parms]),
                    dimnames = i
                )
            x <-
                list(
                    R.res = R,
                    M.mu = M,
                    F.mu = F,
                    m = nM,
                    f = nF,
                    M.sdev = M.sd,
                    F.sdev = F.sd
                )
        }
        if (is.list(x) &&
            !(is.data.frame(x) || tibble::is_tibble(x))) {
            if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f", "R.res") %in% names(x))) {
                rlang::abort(
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
        p <- NROW(R)
        r <- NROW(M)
        o.p <- rep(1, p)
        o.r <- rep(1, r)
        J <- matrix(1, nr <- r, nc <- r)
        D <- M - F
        N <- sum(nM) + sum(nF)
        w <- (nM * nF) / (nM + nF)

        weighted.D <- as.numeric(t(D) %*% w)
        SSCPSex <- weighted.D %o% weighted.D / sum(w)
        SSCPi <- t(D) %*% (w %*% t(o.p) * D) - SSCPSex

        T <-
            diag(sqrt(apply((nM - 1) %o% o.p * M.sd ^ 2 + (nF - 1) %o% o.p * F.sd ^
                                2, 2, sum
            )))
        SSCPe <- T %*% R %*% T

        Xm <- nM %*% t(o.p) * M
        Xf <- nF %*% t(o.p) * F
        SSCPsamp <-
            t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) - t(Xf) %*% J %*%
            Xf / sum(nF) - SSCPi

        Lambda <- det(SSCPe) / det(SSCPi + SSCPe)
        Lambda[2] <-
            det(SSCPi + SSCPe) / det(SSCPSex + SSCPi + SSCPe)
        Lambda[3] <- det(SSCPe) / det(SSCPsamp + SSCPe)

        vh <- r - 1
        vh[2] <- 1
        vh[3] <- vh[1]

        ve <- N - 2 * r
        ve[2] <- N - r - 1
        ve[3] <- ve[1]


        m <- ve + vh - (p + vh + 1) / 2
        s <- sqrt(((p * vh) ^ 2 - 4) / (p ^ 2 + vh ^ 2 - 5))

        DF1 <- p * vh
        DF2 <- m * s - p * vh / 2 + 1

        Rao.R <-
            (1 - Lambda ^ (1 / s)) / (Lambda ^ (1 / s)) * (m * s - p * vh / 2 + 1) /
            (p * vh)
        eta <- (1 - Lambda ^ (1 / s))
        cohen <- (Lambda ^ (-1 / s) - 1)

        p <- stats::pf(Rao.R, DF1, DF2, lower.tail = lower.tail)



        if (es == TRUE) {
            out <-
                cbind.data.frame(
                    c("Sex:Population", "Sex", "Population"),
                    round(vh, 1),
                    round(Lambda, digits),
                    round(Rao.R, digits),
                    round(DF1, 1),
                    round(DF2,
                          digits),
                    round(p, digits),
                    round(eta, digits),
                    round(cohen, digits)
                )
            colnames(out) <-
                c(
                    "term",
                    "df",
                    "wilks",
                    "statistic",
                    "num.df",
                    "den.df",
                    "p.value",
                    "eta.squared",
                    "cohen.f"
                )
        } else{
            out <-
                cbind.data.frame(
                    c("Sex:Population", "Sex", "Population"),
                    round(vh, 1),
                    round(Lambda, digits),
                    round(Rao.R, digits),
                    round(DF1, 1),
                    round(DF2,
                          digits),
                    round(p, digits)
                )
            colnames(out) <-
                c("term",
                  "df",
                  "wilks",
                  "statistic",
                  "num.df",
                  "den.df",
                  "p.value")

        }
        if (univariate == TRUE) {
            M <- reshape2::melt(data = x$M.mu)
            F <- reshape2::melt(data = x$F.mu)
            m <- as.data.frame(as.numeric(x$m))
            f <- as.data.frame(as.numeric(x$f))
            Msd <- reshape2::melt(data = x$M.sdev)
            Fsd <- reshape2::melt(data = x$F.sdev)
            df <-
                cbind_fill(
                    colnames(x$M.mu),
                    rownames(x$M.mu),
                    M$value,
                    F$value,
                    Msd$value,
                    Fsd$value,
                    m,
                    f
                )
            colnames(df) <-
                c("Parms",
                  "Pop",
                  "M.mu",
                  "F.mu",
                  "M.sdev",
                  "F.sdev",
                  "f",
                  "m")
            if (is.data.frame(x)) {
                out <-
                    list(
                        tibble::as_tibble(out),
                        by(
                            df,
                            INDICES = df$Parms,
                            univariate,
                            N = nlevels(x[, Parms]),
                            padjust = padjust,
                            es = es,
                            digits = digits,
                            lower.tail = lower.tail,
                            ...
                        )
                    )
            } else{
                out <-
                    list(
                        tibble::as_tibble(out),
                        by(
                            df,
                            df$Parms,
                            TestDimorph::univariate,
                            N = ncol(x[["M.mu"]]),
                            padjust = padjust,
                            es = es,
                            digits = digits,
                            lower.tail = lower.tail,
                            ...
                        )
                    )


            }
            names(out) <- c("multivariate", "univariate")
            return(out)
        } else {
            return(tibble::as_tibble(out))
        }
    }
