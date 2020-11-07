# t_test ------------------------------------------------------------------
#' @title t-greene test
#' @param m Number of male sample size in the first population, Default: NULL
#' @param m2 Number of male sample size in the second population, Default:
#' NULL
#' @param f Number of female sample size in the first population, Default:
#' NULL
#' @param f2 Number of female sample size in the second population, Default:
#' NULL
#' @param M.mu Means for males in the first population, Default: NULL
#' @param M.mu2 Means for males in the second population, Default: NULL
#' @param F.mu Means for females in the first population, Default: NULL
#' @param F.mu2 Means for females in the second population, Default: NULL
#' @param M.sdev Standard deviation for males in the first population,
#' Default: NULL
#' @param M.sdev2 Standard deviation for males in the second population,
#' Default: NULL
#' @param F.sdev Standard deviation for females in the first population,
#' Default: NULL
#' @param F.sdev2 Standard deviation for females in the second population,
#' Default: NULL
#' @param N Number of pairwise comparisons for [p.adjust.methods], if left
#' `NULL` it will follow the formula `n(n 1)/2` where `n` is the number of
#' populations , Default: NULL
#' @inheritParams t_greene
#' @keywords internal
t_test <-
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
           padjust = padjust,
           N,
           digits,
           sig.level,
           alternative,
           es) {
    tg <-
      ((M.mu - F.mu) - (M.mu2 - F.mu2)) / (sqrt(((((m - 1) * M.sdev^2) + ((f -
        1) * F.sdev^2) + ((m2 - 1) * M.sdev2^2) + ((f2 - 1) * F.sdev2^
        2)
      )) / (m + f +
        m2 + f2 - 4)) * sqrt((1 / m) + (1 / f) + (1 / m2) + (1 / f2)))
    df <- (m + f + m2 + f2 - 4)
    sdp <-
      (sqrt(((((m - 1) * M.sdev^2) + ((f - 1) * F.sdev^2) + ((m2 - 1) * M.sdev2^
        2) + ((f2 - 1) * F.sdev2^2)
      )) / (m + f + m2 + f2 - 4)))
    mean_diff <- ((M.mu - F.mu) - (M.mu2 - F.mu2))
    d <- abs(mean_diff / sdp)
    n1<-m+f
    n2<-m2+f2
    sigma_d <- sqrt(((n1+n2)/(n1*n2))+(d^2/(2*(n1+n2))))
    if (sig.level < 0 ||
      sig.level > 1 || !is.numeric(sig.level)) {
      stop("sig.level should be a number between 0 and 1")
    }
    crit_d <- stats::qnorm(sig.level/2,lower.tail = FALSE)
    lower_d <- d-(crit_d*sigma_d)
    upper_d <- d+(crit_d*sigma_d)
    alternative <-
      match.arg(alternative, choices = c("two.sided", "less", "greater"))
    padjust <-
      match.arg(padjust, choices = p.adjust.methods)
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
        mean_diff + (abs(stats::qt(p = (
          1 - (sig.level / 2)
        ), df = df)) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
      lower <-
        mean_diff - (abs(stats::qt(p = (
          1 - (sig.level / 2)
        ), df = df)) * sdp * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
      p <-
        (2 * stats::pt(abs(tg), df, lower.tail = FALSE))
    }
    if (!is.null(N)) {
      p <-
        padjust_n(
          p = p,
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
    if (!is.logical(es)) {
      stop("es should be either TRUE or FALSE")
    }
    if (isTRUE(es)) {
      data.frame(
        "df" = round(df, digits),
        "mean.diff" = round(mean_diff, digits),
        "conf.low" = round(lower, digits),
        "conf.high" = round(upper, digits),
        "statistic" = round(tg, digits),
        "p.value" = round(p, digits),
        "signif" = signif,
        "cohen.d" = round(d, digits),
        "conf.low.d"=round(lower_d, digits),
        "conf.high.d" = round(upper_d, digits)
      )
    } else {
      data.frame(
        "df" = round(df, digits),
        "mean.diff" = round(mean_diff, digits),
        "conf.low" = round(lower, digits),
        "conf.high" = round(upper, digits),
        "statistic" = round(tg, digits),
        "p.value" = round(p, digits),
        "signif" = signif
      )
    }
  }

# multi_raw ---------------------------------------------------------------

#' @title multivariate data generation
#' @description multivariate data generation helper function
#' @inheritParams raw_gen
#' @keywords internal
multi_raw <- function(x,
                      format,
                      complete_cases,
                      dist,
                      lower,
                      upper) {
  if (dist == "log") {
    stop(
      "Transformation of raw summary data to logged data is only possible for univariate distribution"
    )
  }
  R <- x$R.res
  m <- x$m
  f <- x$f
  m_mean <- x$M.mu
  f_mean <- x$F.mu
  M.sdev <- x$M.sdev
  F.sdev <- x$F.sdev
  traits <- colnames(m_mean)
  n.t <- length(traits)
  pops <- row.names(m_mean)
  n.pops <- length(pops)
  Sex <-
    rep(rep(c("M", "F"), n.pops), times = as.vector(rbind(m, f)))
  Pop <- rep(pops, times = m + f)
  X <- matrix(NA, nrow = sum(m) + sum(f), ncol = n.t)
  start <- 0
  stop <- 0
  for (i in 1:n.pops) {
    start <- stop + 1
    S <- diag(M.sdev[i, ])
    V <- S %*% R %*% S
    stop <- stop + m[i]
    X[start:stop, ] <-
      tmvtnorm::rtmvnorm(
        n = m[i],
        mean = m_mean[i, ],
        sigma = V,
        lower = rep(lower, n.t),
        upper = rep(upper, n.t)
      )
    start <- stop + 1
    S <- diag(F.sdev[i, ])
    V <- S %*% R %*% S
    stop <- stop + f[i]
    X[start:stop, ] <-
      tmvtnorm::rtmvnorm(
        n = f[i],
        mean = f_mean[i, ],
        sigma = V,
        lower = rep(lower, n.t),
        upper = rep(upper, n.t)
      )
  }
  colnames(X) <- traits
  X <- data.frame(Sex, Pop, X)

  if (format == "wide") {
    if (isTRUE(complete_cases)) {
      return(as_tibble(tidyr::drop_na(X)))
    } else {
      return(tibble::as_tibble(X))
    }
  } else {
    tibble::as_tibble(
      pivot_longer(
        data = X,
        cols = -c("Sex", "Pop"),
        names_to = "Parms",
        values_drop_na = complete_cases
      )
    )
  }
}

# dataframe2list ----------------------------------------------------------

#' @title converts multivariate dataframe to list
#' @description helper function for multivariate analysis
#' @inheritParams multivariate
#' @keywords internal
dataframe2list <- function(x, R.res, Parms, Pop) {
  x <- data.frame(x)

  R <- R.res
  i <- list(levels(x[, Pop]), levels(x[, Parms]))
  m_mean <-
    matrix(
      data = x$M.mu,
      nrow = nlevels(x[, Pop]),
      ncol = nlevels(x[, Parms]),
      dimnames = i
    )
  f_mean <-
    matrix(
      data = x$F.mu,
      nrow = nlevels(x[, Pop]),
      ncol = nlevels(x[, Parms]),
      dimnames = i
    )
  nM <- as.vector(x$m)[seq(nlevels(x[, Pop]))]
  nF <- as.vector(x$f)[seq(nlevels(x[, Pop]))]
  nM <- nM[!is.na(nM)]
  nF <- nF[!is.na(nF)]
  M.sd <-
    matrix(
      data = x$M.sdev,
      nrow = nlevels(x[, Pop]),
      ncol = nlevels(x[, Parms]),
      dimnames = i
    )
  F.sd <-
    matrix(
      data = x$F.sdev,
      nrow = nlevels(x[, Pop]),
      ncol = nlevels(x[, Parms]),
      dimnames = i
    )
  x <-
    list(
      R.res = R,
      M.mu = m_mean,
      F.mu = f_mean,
      m = nM,
      f = nF,
      M.sdev = M.sd,
      F.sdev = F.sd
    )
  x
}

# cbind_fill --------------------------------------------------------------

#' @title cbind_fill
#' @description cbind columns with unequal length
#' @param ... columns with unequal length
#' @keywords internal
cbind_fill <- function(...) {
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  df <- do.call(cbind, lapply(nm, function(x) {
  rbind.data.frame(x, matrix(data = NA, n - nrow(x), ncol(x)))
  }))
  names(df) <- names(nm)
  df
}

# anov_es -----------------------------------------------------------------

#' @title anova_es
#' @param x ANOVA model
#' @inheritParams univariate
#' @keywords internal
#' @importFrom tibble as_tibble
#' @importFrom dplyr case_when
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
  signif <- dplyr::case_when(
    p > 0.05 ~ "ns",
    p < 0.05 & p > 0.01 ~ "*",
    p < 0.01 & p > 0.001 ~ "**",
    p < 0.001 ~ "***"
  )
  ss <- summary(x)[[1]][, 2]
  df <- summary(x)[[1]][, 1]
  ms <- summary(x)[[1]][, 3]
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
  tibble::as_tibble(out)
}


# padjust_n ---------------------------------------------------------------

#' @title padjust_n
#' @param p pvalue to be adjusted
#' @param method method of adjustment, Default: p.adjust.methods
#' @param n number of pairwise comparisons, Default: length(p)
#' @description doesn't return an error when n >= lp
#' @rdname padjust_n
#' @importFrom stats p.adjust.methods setNames
#' @keywords internal
padjust_n <- function(p, method = p.adjust.methods, n = length(p)) {
  method <- match.arg(method)
  if (method == "fdr") {
    method <- "BH"
  }
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) {
    nna <- TRUE
  }
  p <- p[nna]
  lp <- length(p)
  if (n <= 1) {
    return(p0)
  }
  if (n == 2 && method == "hommel") {
    method <- "hochberg"
  }
  p0[nna] <- switch(
    method,
    bonferroni = pmin(1, n * p),
    holm = {
      i <- seq_len(lp)
      o <- order(p)
      ro <- order(o)
      pmin(1, cummax((n + 1L - i) * p[o]))[ro]
    },
    hommel = {
      if (n > lp) {
        p <- c(p, rep.int(1, n - lp))
      }
      i <- seq_len(n)
      o <- order(p)
      p <- p[o]
      ro <- order(o)
      q <- pa <- rep.int(min(n * p / i), n)
      for (j in (n - 1L):2L) {
        ij <- seq_len(n - j + 1L)
        i2 <- (n - j + 2L):n
        q1 <- min(j * p[i2] / (2L:j))
        q[ij] <- pmin(j * p[ij], q1)
        q[i2] <- q[n - j + 1L]
        pa <- pmax(pa, q)
      }
      pmax(pa, p)[if (lp < n) {
        ro[1L:lp]
      } else {
        ro
      }]
    },
    hochberg = {
      i <- lp:1L
      o <- order(p, decreasing = TRUE)
      ro <- order(o)
      pmin(1, cummin((n + 1L - i) * p[o]))[ro]
    },
    BH = {
      i <- lp:1L
      o <- order(p, decreasing = TRUE)
      ro <- order(o)
      pmin(1, cummin(n / i * p[o]))[ro]
    },
    BY = {
      i <- lp:1L
      o <- order(p, decreasing = TRUE)
      ro <- order(o)
      q <- sum(1 / (1L:n))
      pmin(1, cummin(q * n / i * p[o]))[ro]
    },
    none = p
  )
  p0
}
