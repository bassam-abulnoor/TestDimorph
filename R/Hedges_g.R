#' @title Hedges' g
#' @description quantifies the size of difference between sexes in measured traits.
#' @inheritParams D_index
#' @return a table of Hedge's g values with confidence interval for different traits.
#' @details Calculates Hedges' (1981) g and its confidence intervals using
#' the pooled standard deviation and correcting for bias.  See Goulet-Pelletier
#' and Cousineau (2018) for details of the calculations and \link{D_index} for
#' description of the bootstrap.
#' @examples
#' library(TestDimorph)
#' data("Cremains_measurements")
#' # Confidence intervals with non-central t distribution
#' Hedges_g(Cremains_measurements[1, ])
#' \dontrun{
#' # confidence interval with bootstrapping
#' Hedges_g(Cremains_measurements[1, ], rand = FALSE, B = 1000)
#' }
#' @rdname Hedges_g
#' @references Hedges, L. V. (1981). Distribution theory for Glass's estimator
#' of effect size and related estimators. Journal of Educational Statistics,
#' 6(2), 107-128.
#'
#'
#' Goulet-Pelletier, J.-C., & Cousineau, D. (2018). A review of effect sizes and
#' their confidence intervals, part I: The Cohen's d family. The Quantitative
#' Methods for Psychology, 14(4), 242-265.
#' @export

Hedges_g <- function(x, Trait = 1, CI = 0.95, B = NULL, rand = TRUE, digits = 4) {
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
  if (!(Trait %in% seq_along(x))) {
    stop("Trait should be number from 1 to ncol(x)")
  }
  if (length(unique(x$Trait)) != length(which(!is.na(x$Trait)))) {
    stop("Each trait should be represented by a single raw with a unique name")
  }
  if (!is.numeric(CI) && !is.null(CI)) {
    stop("confidence level should be a number from 0 to 1")
  }
  if (is.numeric(CI) && any(CI < 0, CI > 1)) {
    stop("confidence level should be a number from 0 to 1")
  }
  x <- x %>%
    drop_na() %>%
    as.data.frame()
  x$Trait <- x[, Trait]
  hedge <- function(x) {
    m <- x$m[1]
    M.mu <- x$M.mu[1]
    M.sdev <- x$M.sdev[1]
    f <- x$f[1]
    F.mu <- x$F.mu[1]
    F.sdev <- x$F.sdev[1]
    Trait <- Trait[1]

    # From: Hedges, L. V. (1981) Distribution theory for Glass's
    #       estimator of effect size and related estimators.
    #       Journal of Educational Statistics, 6(2), 107-128.

    # Get Hedge's g using pooled standard deviation and bias correction
    n.tilde <- m * f / (m + f)
    dif <- abs(M.mu - F.mu)
    nu <- m + f - 2
    sd.pool <- sqrt((M.sdev^2 * (m - 1) + F.sdev^2 * (f - 1)) / nu)
    c.m <- sqrt(2) * gamma(nu / 2) / (sqrt(nu) * gamma((nu - 1) / 2))
    g <- dif / sd.pool * c.m
    if (is.null(B)) {
      # Find confidence interval using non-central t-distribution
      CI <- 0.5 * c(1 - CI, 1 + CI)
      CI <- qt(CI, nu, g * sqrt(n.tilde)) / sqrt(n.tilde)
      cbind.data.frame(g, lower = CI[1], upper = CI[2])
    } else {
      if (isFALSE(rand)) {
        set.seed(42)
      }
      sto.boot <- rep(NA, B)
      for (i in 1:B) {
        males <- rnorm(m, M.mu, M.sdev)
        females <- rnorm(f, F.mu, F.sdev)
        M.mu.boot <- mean(males)
        M.sdev.boot <- sd(males)
        F.mu.boot <- mean(females)
        F.sdev.boot <- sd(females)
        dif_boot <- abs(M.mu.boot - F.mu.boot)
        sd.pool_boot <- sqrt((M.sdev.boot^2 * (m - 1) + F.sdev.boot^2 * (f - 1)) / nu)
        boot.D <- dif_boot / sd.pool_boot * c.m
        obs.D <- g
        sto.boot[i] <- boot.D
        cat(paste("\r", i, " of ", B, "bootstraps"))
        flush.console()
      }
      half.CI <- CI / 2
      bot <- 0.5 - half.CI
      top <- 0.5 + half.CI
      b <- -qnorm(sum(sto.boot < obs.D) / B)
      alpha.1 <- pnorm(qnorm(bot) - 2 * b)
      alpha.2 <- pnorm(qnorm(top) - 2 * b)
      bounds <- quantile(sto.boot, c(alpha.1, alpha.2))
      bounds <- as.numeric(bounds)
      cbind.data.frame(g, lower = bounds[1], upper = bounds[2])
    }
  }
  hedge_list <- lapply(split.data.frame(x, x$Trait), hedge)
  name_hedge_list <- names(hedge_list)
  hedge_df <- do.call(rbind.data.frame, hedge_list)
  rownames(hedge_df) <- name_hedge_list
  rown_col(hedge_df, "Trait") %>%
    relocate(.data$lower,
      .before =
        .data$g
    ) %>%
    mutate(across(-1, round, digits)) %>%
    as.data.frame() %>%
    drop_na()
}
