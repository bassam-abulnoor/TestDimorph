#' @title Multivariate Analysis Of Sexual Dimorphism
#' @description Multivariate extension of Greene t test [t_greene]
#' @param x Tibble/Data frame or list containing summary statistics for
#' multiple parameters measured in both sexes in two or more populations.
#' @param R.res Pooled within correlational matrix, Default: NULL
#' @param Parms Number of the column containing names of measured parameters,
#' Default: 1
#' @param Pop Number of the column containing populations' names, Default: 2
#' @inheritParams univariate
#' @param univariate Logical; if TRUE conducts multiple univariate analyses on
#' different parameters separately, Default: FALSE
#' @param ... Additional arguments that could be passed to the [univariate]
#' function
#' @return Tibble of MANOVA results
#' @details Data can be entered either as a tibble/data frame of summary
#' statistics as in [baboon.parms_df] . In that case the pooled within
#' correlational matrix `R.res` should be entered as a separate argument as in
#' [R]. Another acceptable format is a named list of matrices containing
#' different summary statistics as well as the correlational matrix as in
#' [baboon.parms_list]. By setting the option `univariate` to `TRUE`, multiple
#' `ANOVA`s can be run on each parameter independently with the required p
#' value correction using [p.adjust.methods].
#' @examples
#' \donttest{
#' # x is a data frame with separate correlational matrix
#' library(TestDimorph)
#' multivariate(baboon.parms_df, R.res = R)
#' # x is a list with the correlational matrix included
#' library(TestDimorph)
#' multivariate(baboon.parms_list, univariate = TRUE, padjust = "bonferroni")
#' }
#' @rdname multivariate
#' @export
#' @importFrom stats pf
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble is_tibble

multivariate <- function(x,
                         R.res = NULL,
                         Parms = 1,
                         Pop = 2,
                         es = FALSE,
                         univariate = FALSE,
                         padjust = p.adjust.methods,
                         ...,
                         lower.tail = FALSE,
                         digits = 4) {
  # Data preparation --------------------------------------------------------

  if (!(is.list(x) || is.data.frame(x))) {
    stop("x should be a list or a dataframe")
  }
  if (!is.logical(es)) {
    stop("es should be either TRUE or FALSE")
  }
  if (!is.logical(univariate)) {
    stop("univariate should be either TRUE or FALSE")
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
    if (!(Parms %in% seq_along(x))) {
      stop("Parms should be number from 1 to ncol(x)")
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
    padjust <- match.arg(padjust, choices = p.adjust.methods)
    x <- dataframe2list(
      x = x,
      R.res = R.res,
      Parms = Parms,
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

    # MANOVA ------------------------------------------------------------------

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

  Gmean <-
    apply(nM %o% o.p * M + nF %o% o.p * F, 2, sum) / (sum(nM + nF))

  d <- M[1, ] - Gmean
  between <- d %o% d * nM[1]
  d <- F[1, ] - Gmean
  between <- between + d %o% d * nF[1]
  for (i in 2:r)
  {
    d <- M[i, ] - Gmean
    between <- between + d %o% d * nM[i]
    d <- F[i, ] - Gmean
    between <- between + d %o% d * nF[i]
  }

  J <- matrix(1, nrow = r, ncol = r)
  D <- M - F
  N <- sum(nM) + sum(nF)
  w <- (nM * nF) / (nM + nF)


  weighted.D <- as.numeric(t(D) %*% w)
  SSCPi <-
    t(D) %*% (w %*% t(o.p) * D) - weighted.D %o% weighted.D / sum(w)

  T <-
    diag(sqrt(apply((nM - 1) %o% o.p * M.sd^2 + (nF - 1) %o% o.p *
      F.sd^2, 2, sum)))
  SSCPe <- T %*% R %*% T
  TotSSCP <- SSCPe + between

  Xm <- nM %*% t(o.p) * M
  Xf <- nF %*% t(o.p) * F
  SSCPsamp <-
    t(Xm) %*% M + t(Xf) %*% F - t(Xm) %*% J %*% Xm / sum(nM) -
    t(Xf) %*% J %*% Xf / sum(nF) - SSCPi

  SSCPsex <- TotSSCP - SSCPsamp - SSCPe - SSCPi

  Lambda <- det(SSCPe) / det(SSCPi + SSCPe)
  Lambda[2] <- det(SSCPe) / det(SSCPsex + SSCPe)
  Lambda[3] <- det(SSCPe) / det(SSCPsamp + SSCPe)

  vh <- r - 1
  vh[2] <- 1
  vh[3] <- vh[1]

  ve <- N - 2 * r
  for (i in 2:3) {
    ve[i] <- ve[1]
  }

  m <- ve + vh - (p + vh + 1) / 2
  s <- sqrt(((p * vh)^2 - 4) / (p^2 + vh^2 - 5))

  DF1 <- p * vh
  DF2 <- m * s - p * vh / 2 + 1

  Rao.R <-
    (1 - Lambda^(1 / s)) / (Lambda^(1 / s)) * (m * s - p * vh / 2 + 1) /
      (p *
        vh)
  eta <- (1 - Lambda^(1 / s))
  cohen <- (Lambda^(-1 / s) - 1)

  p <- stats::pf(Rao.R, DF1, DF2, lower.tail = lower.tail)
  signif <-
    case_when(
      p > 0.05 ~ "ns",
      p < 0.05 & p > 0.01 ~ "*",
      p < 0.01 & p > 0.001 ~ "**",
      p < 0.001 ~ "***"
    )
  if (isTRUE(es)) {
    out <- data.frame(
      "term" = c("Sex:Population", "Sex", "Population"),
      "df" = round(vh, 1),
      "wilks" = round(Lambda, digits),
      "statistic" = round(Rao.R, digits),
      "num.df" = round(DF1, 1),
      "den.df" = round(DF2, digits),
      "p.value" = round(p, digits),
      "signif" = signif,
      "eta.squared" = round(
        eta,
        digits
      ),
      "cohen.f" = round(cohen, digits)
    )
  } else {
    out <- data.frame(
      "term" = c("Sex:Population", "Sex", "Population"),
      "df" = round(vh, 1),
      "wilks" = round(Lambda, digits),
      "statistic" = round(Rao.R, digits),
      "num.df" = round(DF1, 1),
      "den.df" = round(DF2, digits),
      "p.value" = round(p, digits),
      "signif" = signif
    )
  }

  # Pairwise univariate tests -----------------------------------------------

  if (isTRUE(univariate)) {
    m_mean <-
      x$M.mu %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      pivot_longer(cols = 2:(ncol(as.data.frame(x$M.mu)) +
        1))
    f_mean <-
      x$F.mu %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      pivot_longer(cols = 2:(ncol(as.data.frame(x$F.mu)) +
        1))
    Msd <-
      x$M.sdev %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      pivot_longer(cols = 2:(ncol(as.data.frame(x$M.sdev)) +
        1))
    Fsd <-
      x$F.sdev %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      pivot_longer(cols = 2:(ncol(as.data.frame(x$F.sdev)) +
        1))
    A <- (ncol(x$M.mu) * nrow(x$M.mu)) / length(x$m)
    B <- (ncol(x$F.mu) * nrow(x$F.mu)) / length(x$f)
    m <- rep((as.numeric(x$m)), A)
    f <- rep((as.numeric(x$f)), B)
    df <-
      cbind.data.frame(
        "Parms" = colnames(x$M.mu),
        "Pop" = rownames(x$M.mu),
        "M.mu" = m_mean$value,
        "F.mu" = f_mean$value,
        "M.sdev" = Msd$value,
        "F.sdev" = Fsd$value,
        "m" = m,
        "f" = f
      )
    df$Pop <- factor(df$Pop)
    df$Parms <- factor(df$Parms)
    A <- ((length(df) - nlevels(df$Parms)) / nlevels(df$Parms))
    B <- ((length(df) - nlevels(df$Pop)) / nlevels(df$Pop))
    df$Parms <- rep(levels(df$Parms), A)
    df$Pop <- rep(levels(df$Pop), B)
    if (is.data.frame(x)) {
      out <- list(
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
    } else {
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
    out
  } else {
    tibble::as_tibble(out, pillar.sigfig = 6)
  }
}
