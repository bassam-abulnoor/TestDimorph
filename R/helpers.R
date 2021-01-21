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
           CI,
           alternative,
           es) {
    if (CI < 0 ||
      CI > 1 || !is.numeric(CI)) {
      stop("CI should be a number between 0 and 1")
    }
    CI <- 1 - CI
    es <-
      match.arg(es, choices = c("none", "d", "g"))
    alternative <-
      match.arg(alternative, choices = c("two.sided", "less", "greater"))
    padjust <-
      match.arg(padjust, choices = p.adjust.methods)
    df <- (m + f + m2 + f2 - 4)
    sd_pooled <-
      sqrt(((((m - 1) * M.sdev^2) + ((f - 1) * F.sdev^2) + ((m2 - 1) *
        M.sdev2^2) + ((f2 - 1) * F.sdev2^2)
      )) / df)
    mean_diff <- ((M.mu - F.mu) - (M.mu2 - F.mu2))
    tg <-
      mean_diff / (sd_pooled * sqrt((1 / m) + (1 / f) + (1 / m2) + (1 / f2)))
    d <- abs(2 * tg / sqrt(df))
    var_d <-
      (m + f + m2 + f2) / ((m + f) * (m2 + f2)) + (d^2) / (2 * (m + f + m2 + f2))
    j <- 1 - (3 / (4 * df - 1))
    g <- j * d
    var_g <- j^2 * var_d
    if (es == "d") {
      eff <- d
      var_eff <- var_d
    } else {
      eff <- g
      var_eff <- var_g
    }
    if (alternative == "two.sided") {
      crit <- stats::qt(p = (1 - (CI / 2)), df = df)
    } else {
      crit <- stats::qt(p = (1 - CI), df = df)
    }
    upper <-
      mean_diff + (abs(crit) * sd_pooled * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
    lower <-
      mean_diff - (abs(crit) * sd_pooled * sqrt(sum(1 / m, 1 / f, 1 / m2, 1 / f2)))
    upper_eff <- eff + crit * sqrt(var_eff)
    lower_eff <- eff - crit * sqrt(var_eff)
    p <- switch(
      alternative,
      less = stats::pt(abs(tg), df, lower.tail = TRUE),
      greater = stats::pt(abs(tg), df, lower.tail = FALSE),
      two.sided = 2 * stats::pt(abs(tg), df, lower.tail = FALSE)
    )

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

    if (es != "none") {
      out <- data.frame(
        "df" = round(df, digits),
        "mean.diff" = round(mean_diff, digits),
        "conf.low" = round(lower, digits),
        "conf.high" = round(upper, digits),
        "statistic" = round(tg, digits),
        "p.value" = round(p, digits),
        "signif" = signif,
        round(eff, digits),
        "conf.es.low" = round(lower_eff, digits),
        "conf.es.high" = round(upper_eff, digits)
      )
      names(out)[8] <- es
      return(out)
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
      "Transformation of raw summary data to logged data is only possible for
      univariate distribution"
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
      return(tidyr::drop_na(X))
    } else {
      return(X)
    }
  } else {
    pivot_longer(
      data = X,
      cols = -c("Sex", "Pop"),
      names_to = "Parms",
      values_drop_na = complete_cases
    )
  }
}

# dataframe2list ----------------------------------------------------------

#' @title converts multivariate dataframe to list
#' @description helper function for multivariate analysis
#' @inheritParams multivariate
#' @keywords internal
dataframe2list <- function(x, R.res, Trait, Pop) {
  x <- x %>%
    rename("Pop" = Pop, "Trait" = Trait) %>%
    mutate(Pop = factor(.data$Pop), Trait = factor(.data$Trait)) %>%
    as.data.frame()
  x$Pop <- droplevels(x$Pop)
  x$Trait <- droplevels(x$Trait)
  R <- R.res
  name_mat <- list(levels(x$Pop), levels(x$Trait))
  z <- arrange(x, Trait, Pop)
  m_mean <-
    matrix(
      data = x$M.mu,
      nrow = nlevels(x$Pop),
      ncol = nlevels(x$Trait),
      dimnames = name_mat
    )
  f_mean <-
    matrix(
      data = x$F.mu,
      nrow = nlevels(x$Pop),
      ncol = nlevels(x$Trait),
      dimnames = name_mat
    )
  nM <- as.vector(x$m)[seq(nlevels(x$Pop))]
  nF <- as.vector(x$f)[seq(nlevels(x$Pop))]
  nM <- nM[!is.na(nM)]
  nF <- nF[!is.na(nF)]
  M.sd <-
    matrix(
      data = x$M.sdev,
      nrow = nlevels(x$Pop),
      ncol = nlevels(x$Trait),
      dimnames = name_mat
    )
  F.sd <-
    matrix(
      data = x$F.sdev,
      nrow = nlevels(x$Pop),
      ncol = nlevels(x$Trait),
      dimnames = name_mat
    )
  x <-
    list(
      z = z,
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
#' @description cbind columns with unequal length with naming
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
#' @title cbind_fill2
#' @description cbind columns with unequal length without naming
#' @param ... columns with unequal length
#' @keywords internal
cbind_fill2 <- function(...) {
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function(x) {
    as.data.frame(rbind(x, matrix(data = NA, n - nrow(x), ncol(x))))
  }))
}

# anov_es -----------------------------------------------------------------

#' @title anova_es
#' @param x ANOVA model
#' @inheritParams univariate
#' @keywords internal
#' @importFrom dplyr case_when
anova_es <- function(x,
                     es_anova = es_anova,
                     digits = digits,
                     CI = CI) {
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
  SST <- sum(ss)
  df <- summary(x)[[1]][, 1]
  N <- sum(df)
  ms <- summary(x)[[1]][, 3]
  eta <- f * df[1] / (f * df[1] + df[2])
  omega <- (ssi - (df1 * within)) / (ssi + sse + within)
  cohen_squared <- (eta) / (1 - eta)
  if (es_anova != "none") {
    eff <- switch(es_anova,
      f = cohen_squared,
      eta = eta
    )
    eff <-
      eff_CI(
        f = f,
        CI = CI,
        eff = eff,
        df1 = df1,
        df2 = df2,
        es_type = es_anova
      )
    lower_eff <- eff[[2]]
    upper_eff <- eff[[3]]
    eff <- eff[[1]]
    out <-
      cbind_fill(
        "term" = c("Populations", "Residuals"),
        "df" = round(df, 1),
        "sumsq" = round(ss, digits),
        "meansq" = round(ms, digits),
        "statistic" = round(f, digits),
        "p.value" = round(p, digits),
        "signif" = signif,
        round(eff, digits),
        "conf.es.low" = round(lower_eff, digits),
        "conf.es.high" = round(upper_eff, digits)
      )
    names(out)[8] <- es_anova
  } else {
    out <-
      cbind_fill(
        "term" = c("Populations", "Residuals"),
        "df" = round(df, 1),
        "sumsq" = round(ss, digits),
        "meansq" = round(ms, digits),
        "statistic" = round(f, digits),
        "p.value" = round(p, digits),
        "signif" = signif
      )
  }
  as.data.frame(out)
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

# pooled_cov --------------------------------------------------------------

#' pooled_cov
#' @description Pooled covariance matrix
#' @param x A matrix with continuous data
#' @param ina A numerical vector indicating the groups
#' @importFrom Morpho covW
#' @keywords internal
pooled_cov <- function(x, ina) {
  Morpho::covW(x, as.factor(ina))
}
# Confidence interval for anova and manova effect sizes -------------------
#' eff_CI
#' @description Confidence intervals for ANOVA and MANOVA effect sizes
#' @param f critical F-value
#' @param CI confidence level, Default 0.95
#' @param eff effect size
#' @param df1 numerator degree of freedom
#' @param df2 denominator degree of freedom
#' @param SS sum of squares
#' @param SST total SS
#' @param N sum of df
#' @param es_type type of effect size
#' @keywords internal
eff_CI <- function(f, CI, eff, df1, df2, es_type = "eta") {
  get_NCP <- function(F, df.1, df.2, CI) {
    # From the FORTRAN code in:
    # Guirguis, G. H. (1990). A note on computing the noncentrality
    # parameter of the noncentral F-distribution.
    # Communications in Statistics-Simulation and Computation, 19(4), 1497-1511.

    #### ANORM function  #########
    ANORM <- function(X, DFN, DFD, FL, QUANT) {
      A <- DFN + FL
      B <- 0.22222 * (1 + FL / A) / A
      A <- exp(log(X * DFN / A) / 3)
      anorm <- (A * (1 - 0.22222 / DFD) - (1 - B)) / sqrt(B + 0.22222 / DFD * A^2) - QUANT
      return(anorm)
    }

    ### GUESS function #########
    GUESS <- function(X = 20, DFN = 2, DFD = 2, FX = 0.01) {
      ACC <- 0.01
      N <- 50
      QUANT <- qnorm(FX)
      FA <- pf(X, DFN, DFD)
      if (FA - FX <= 0) {
        return(0)
      }

      REFQ <- ANORM(X, DFN, DFD, 0, QUANT)
      FL <- 2 * log(FA / FX)
      FL <- max(c(FL, 1))
      FLO <- 1.e30

      for (I in 1:50) {
        REF <- ANORM(X, DFN, DFD, FL, QUANT)
        if (abs(FLO - FL) < ACC) {
          if (FL < 0) FL <- 0
          return(FL)
        }
        FLO <- FL
        FL <- FL * REFQ / (REFQ - REF)
      }
      if (FL < 0) FL <- 0
      return(FL)
    }

    ### FLAMDA function #####
    FLAMDA <- function(X = 20, DFN = 2, DFD = 2, FX = 0.01) {
      ACC <- 1.e-06

      FL <- GUESS(X, DFN, DFD, FX)
      if (FL == 0) {
        return(0)
      }

      B <- X * DFN / (DFN + 2)

      for (I in 1:50) {
        FA <- pf(X, DFN, DFD, FL)
        FC <- pf(B, DFN + 2, DFD, FL)
        APROX <- 2 * FA / (FC - FA) * log(FX / FA)
        FL <- FL + APROX
        if (abs(APROX) < ACC * max(c(FL, 1))) {
          return(FL)
        }
      }
    }
    ##################################################

    hi <- (1 - CI) / 2
    lo <- 1 - hi
    lo <- FLAMDA(F, df.1, df.2, lo)
    hi <- FLAMDA(F, df.1, df.2, hi)
    return(c(lo, hi))
  }

  eta_squared <- function(f, CI, eff, df1, df2) {
    ncps <- get_NCP(f, df1, df2, CI)
    N <- df1 + df2 + 1
    ci <- ncps / (ncps + N)
    lower_eff <- ci[1]
    upper_eff <- ci[2]
    cbind.data.frame(lower_eff = lower_eff, upper_eff = upper_eff)
  }
  f_squared <- function(f, CI, eff, df1, df2) {
    ncps <- get_NCP(f, df1, df2, CI)
    N <- df1 + df2 + 1
    ci <- ncps / (ncps + N)
    lower_eff <- ci[1] / (1 - ci[1])
    upper_eff <- ci[2] / (1 - ci[2])
    cbind.data.frame(lower_eff = lower_eff, upper_eff = upper_eff)
  }


  out <- switch(es_type,
    none = eta_squared(f, CI, eff, df1, df2),
    eta = eta_squared(f, CI, eff, df1, df2),
    f = f_squared(f, CI, eff, df1, df2)
  )

  cbind.data.frame(
    ES = eff,
    conf.es.low = out$lower_eff,
    conf.es.high = out$upper_eff
  )
}


# van_vark from raw data --------------------------------------------------
#' Van_vark_raw
#' @description runs van_vark function on raw data
#' @inheritParams extract_sum
#' @keywords internal
Van_vark_raw <- function(x, Sex, Pop, firstX, ...) {
  vec_sum <- function(x, i) {
    x$Pop <- factor(x$Pop, levels = unique(x$Pop))
    cbind.data.frame(
      Trait = rep(names(x)[i], nlevels(x$Pop)),
      extract_sum(x, run = FALSE, firstX = i, Sex = Sex, Pop = Pop, test = 1)
    )
  }
  vec_sum <-
    Vectorize(vec_sum, vectorize.args = "i", SIMPLIFY = FALSE)
  df <-
    do.call(rbind.data.frame, c(vec_sum(x, firstX:ncol(x))))
  df$Trait <- factor(df$Trait, levels = unique(df$Trait))
  df <- df %>% relocate(.data$Trait, .before = 1)
  x <- as.data.frame.list(x)
  x <- x %>% arrange(Pop, Sex)
  sex <- as.numeric(x$Sex) - 1
  pop <- as.numeric(x$Pop)
  pop.names <- names(table(x$Pop))
  N.pops <- length(pop.names)
  ina <- pop + N.pops * sex
  X <- x[, -(1:(firstX - 1))]
  list("W" = pooled_cov(as.matrix(X), ina), "x" = df)
}

# univariate pairwise -----------------------------------------------------
#' post hoc univariate analysis to MANOVA
#'
#' @inheritParams univariate
#' @inheritDotParams multivariate
#' @param out output of multivariate function
#' @param ... other arguments that are passed to univariate function
#' @keywords internal
univariate_pairwise <- function(x, out, padjust, digits, lower.tail, ...) {
  Parms <- NULL
  x$R.res <- NULL
  A <- (ncol(x$M.mu) * nrow(x$M.mu)) / length(x$m)
  B <- (ncol(x$F.mu) * nrow(x$F.mu)) / length(x$f)
  m_dat <- cbind.data.frame(m = rep(x$m, A), Pop = factor(rep(rownames(x$M.mu), A),
    levels = rownames(x$M.mu)
  ))
  f_dat <- cbind.data.frame(f = rep(x$f, A), Pop = factor(rep(rownames(x$M.mu), A),
    levels = rownames(x$M.mu)
  ))
  m_dat <- m_dat[order(m_dat$Pop), ]
  f_dat <- f_dat[order(f_dat$Pop), ]
  x$m <- NULL
  x$f <- NULL
  rownames(x$M.mu) -> r
  colnames(x$M.mu) -> cl
  x <- lapply(x[which(names(x) != "z")], function(y) {
    matrix(y, nrow = nrow(y), ncol = ncol(y), dimnames = list(r, cl))
  })
  m_mean <-
    x$M.mu %>%
    as.data.frame() %>%
    rown_col(var = "Pop") %>%
    pivot_longer(cols = 2:(ncol(as.data.frame(x$M.mu)) +
      1), names_to = "Parms", values_to = "M.mu") %>%
    mutate(Pop = factor(.data$Pop, levels = rownames(x$M.mu)), Parms = factor(
      .data$Parms,
      levels = colnames(x$M.mu)
    ))

  m_mean <- with(m_mean, m_mean[order(Parms, Pop, M.mu), ])
  f_mean <-
    x$F.mu %>%
    as.data.frame() %>%
    rown_col(var = "Pop") %>%
    pivot_longer(
      cols = 2:(ncol(as.data.frame(x$F.mu)) + 1), names_to = "Parms",
      values_to = "F.mu"
    ) %>%
    mutate(Pop = factor(.data$Pop,
      levels =
        rownames(x$F.mu)
    ), Parms = factor(.data$Parms, levels = colnames(x$F.mu)))
  f_mean <- with(f_mean, f_mean[order(Parms, Pop, F.mu), ])
  Msd <-
    x$M.sdev %>%
    as.data.frame() %>%
    rown_col(var = "Pop") %>%
    pivot_longer(cols = 2:(ncol(as.data.frame(x$M.sdev)) + 1), names_to = "Parms", values_to = "M.sdev") %>%
    mutate(Pop = factor(.data$Pop, levels = rownames(x$M.sdev)), Parms = factor(
      .data$Parms,
      levels = colnames(x$M.sdev)
    ))
  Msd <- with(Msd, Msd[order(Parms, Pop, M.sdev), ])
  Fsd <-
    x$F.sdev %>%
    as.data.frame() %>%
    rown_col(var = "Pop") %>%
    pivot_longer(cols = 2:(ncol(as.data.frame(x$F.sdev)) +
      1), names_to = "Parms", values_to = "F.sdev") %>%
    mutate(Pop = factor(.data$Pop, levels = rownames(x$F.sdev)), Parms = factor(
      .data$Parms,
      levels = colnames(x$F.sdev)
    ))
  Fsd <- with(Fsd, Fsd[order(Parms, Pop, F.sdev), ])
  my_merge <- function(df1, df2) {
    merge(df1, df2, by = c("Pop", "Parms"), all.x = TRUE, all.y = TRUE)
  }
  list_frames <- Reduce(my_merge, list(m_mean, f_mean, Msd, Fsd)) %>%
    mutate(Pop = factor(.data$Pop, levels = rownames(x$M.mu))) %>%
    relocate(.data$Parms,
      .before = 1
    ) %>%
    arrange(.data$Pop, .data$Parms)
  df <- cbind.data.frame(list_frames, m = m_dat$m, f = f_dat$f)
  if (is.data.frame(x)) {
    out <- list(
      out,
      by(
        df,
        INDICES = df$Parms,
        univariate,
        N = nlevels(x[, Parms]),
        padjust = padjust,
        digits = digits,
        lower.tail = lower.tail,
        ...
      )
    )
  } else {
    out <-
      list(
        out,
        by(
          df,
          df$Parms,
          TestDimorph::univariate,
          N = ncol(x[["M.mu"]]),
          padjust = padjust,
          digits = digits,
          lower.tail = lower.tail,
          ...
        )
      )
  }
  names(out) <- c("multivariate", "univariate")
  out
}

# add_sig -----------------------------------------------------------------

#' add_sig
#'
#' @param x output of ANOVA and MANOVA tests
#'
#' @return markers for significance values
#' @keywords internal

add_sig <- function(x) {
  x %>%
    as.data.frame() %>%
    mutate(p.value = as.numeric(.data$p.value), signif = case_when(
      .data$p.value > 0.05 ~ "ns",
      .data$p.value < 0.05 & p.value > 0.01 ~ "*",
      .data$p.value < 0.01 & p.value > 0.001 ~ "**",
      .data$p.value < 0.001 ~ "***"
    )) %>%
    relocate(.data$signif, .after = .data$p.value) -> x
}

# rown_col -----------------------------------------------------------------

#' rown_col
#'
#' @param x dataframe with rownames
#' @param var new column name
#' @details convert rownames to column
#' @keywords internal
rown_col <- function(x, var) {
  x <- as.data.frame(x)
  rownames(x) -> var2
  x[, var] <- var2
  relocate(x, var, .before = 1) -> x
  rownames(x) <- NULL
  as.data.frame(x)
}
