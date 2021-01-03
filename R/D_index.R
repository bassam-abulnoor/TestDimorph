#' @title Dissimilarity index
#' @description Visual and statistical computation of the area of non-overlap in
#' the trait distribution between two sex groups.
#' @inheritParams univariate
#' @inheritParams multivariate
#' @param plot logical; if true a plot of densities for both sexes is returned,
#' Default: FALSE
#' @param fill Specify which sex's density to be filled with color in the plot;
#' either "male" in blue color, "female" in pink color or "both", Default: 'female'
#' @param B number of bootstrap samples for generating confidence intervals. Higher
#' number means greater accuracy but slower execution. If NULL bootstrap confidence
#' intervals are not produced, Default:NULL
#' @param rand logical; if TRUE, uses random seed.  If FALSE, then set.seed(42)
#' for repeatability, Default: TRUE
#' @details Chakraborty and Majumder's (1982) D index.  The calculations are done
#' using Inman and Bradley's (1989) equations, and the relationship that
#' D =  1 - OVL where OVL is the overlap coefficient described in Inman and Bradley.
#' A parametric bootstrap was used assuming normal distributions. The method is
#' known as the "bias-corrected percentile method" (Efron, 1981) or
#' the "bias-corrected percentile interval" (Tibshirani, 1984)
#' @return a table and a graphical representation of the selected traits and
#' their corresponding dissimilarity indices, confidence intervals and
#' significance tests.
#' @import dplyr
#' @import ggplot2
#' @importFrom stats dnorm integrate na.omit qnorm contr.sum pnorm quantile
#' @importFrom tidyr drop_na
#' @importFrom utils flush.console
#' @export
#' @references Chakraborty, Ranajit, and Partha P. Majumder.(1982) "On Bennett's
#' measure of sex dimorphism." American Journal of Physical Anthropology
#' 59.3 : 295-298.
#'
#'
#' Inman, Henry F., and Edwin L. Bradley Jr.(1989) "The overlapping coefficient as a
#' measure of agreement between probability distributions and point estimation
#' of the overlap of two normal densities." Communications in Statistics-Theory
#' and Methods 18.10:3851-3874.
#'
#' Efron, B. (1981). Nonparametric standard errors and confidence intervals.
#' Canadian Journal of Statistics, 9(2), 139-158.
#'
#' Tibshirani, R. J. (1984). Bootstrap confidence intervals. Technical Report
#' No. 3, Laboratory for Computational Statistics, Department of Statistics,
#' Stanford University.
#'
#' @examples
#' library(TestDimorph)
#' data("Cremains_measurements")
#' # plot and test of significance
#' D_index(Cremains_measurements[1, ], plot = TRUE)
#' \dontrun{
#' # confidence interval with bootstrapping
#' D_index(Cremains_measurements[1, ], rand = FALSE, B = 1000)
#' }
#'
D_index <- function(x, plot = FALSE, fill = "female", Trait = 1, B = NULL, CI = 0.95,
                    rand = TRUE,digits=4) {
  fill <- match.arg(fill, choices = c("female", "male", "both"))
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
  if (!is.logical(plot)) {
    stop("plot should be either TRUE or FALSE")
  }
  if (!is.logical(rand)) {
    stop("rand should be either TRUE or FALSE")
  }
  if (!is.numeric(B) && !is.null(B)) {
    stop("number of bootstraps should be numeric")
  }
  if (!is.numeric(CI) && !is.null(CI)) {
    stop("confidence level should be a number from 0 to 1")
  }
  if (is.numeric(CI) && any(CI < 0, CI > 1)) {
    stop("confidence level should be a number from 0 to 1")
  }
  if (is.null(CI) && !is.null(B)) {
    warning("confidence level is not spicified")
  }
  if (isFALSE(rand)) {
    set.seed(42)
  }
  x <- x %>%
    drop_na() %>%
    as.data.frame()
  x$Trait <- x[, Trait]
  x$Trait <- factor(x$Trait, levels = x$Trait)
  m <- x$m
  M.mu <- x$M.mu
  M.sdev <- x$M.sdev
  f <- x$f
  F.mu <- x$F.mu
  F.sdev <- x$F.sdev
  InBr.D <- function(Trait, m, M.mu, M.sdev, f, F.mu, F.sdev, i) {
    m <- x$m[i]
    M.mu <- x$M.mu[i]
    M.sdev <- x$M.sdev[i]
    f <- x$f[i]
    F.mu <- x$F.mu[i]
    F.sdev <- x$F.sdev[i]
    Trait <- x$Trait[i]

    from <- min(qnorm(1e-04,
      mean = c(M.mu, F.mu), sd = c(M.sdev, F.sdev),
      lower.tail = TRUE
    ))
    to <-
      max(qnorm(1e-04,
        mean = c(M.mu, F.mu), sd = c(M.sdev, F.sdev), lower.tail =
          FALSE
      ))
    abs.dens <- function(x) {
      abs(dnorm(x, c(M.mu, F.mu)[1], c(M.sdev, F.sdev)[1]) - dnorm(x, c(
        M.mu, F.mu
      )[2], c(M.sdev, F.sdev)[2]))
    }
    D <-
      integrate(Vectorize(abs.dens), lower = from, upper = to)$val / 2
    z <- seq(from, to, (to - from) / 1000)
    trait <- rep(Trait, length(z))
    dn_male <- dnorm(z, as.numeric(na.omit(M.mu)), as.numeric(na.omit(M.sdev))) %>% as.data.frame()
    dn_male <- cbind.data.frame(z = z, dn = dn_male, sex = rep("M", nrow(
      dn_male
    )))
    dn_female <- dnorm(z, as.numeric(na.omit(F.mu)), as.numeric(na.omit(F.sdev))) %>% as.data.frame()
    dn_female <- cbind.data.frame(z = z, dn = dn_female, sex = rep("F", nrow(dn_female)))
    df <- rbind.data.frame(dn_male, dn_female)
    names(df) <- c("z", "dn", "sex")
    df <- cbind.data.frame(trait = trait, df)
    if (!is.null(CI) && !is.null(B)) {
      sto_boot <- rep(NA, B)
      for (i in 1:B) {
        males <- rnorm(m, M.mu, M.sdev)
        females <- rnorm(f, F.mu, F.sdev)
        M.mu_boot <- mean(males)
        M.sdev_boot <- sd(males)
        F.mu_boot <- mean(females)
        F.sdev_boot <- sd(females)
        from_boot <-
          min(qnorm(1e-04, mean = c(M.mu_boot, F.mu_boot), sd = c(
            M.sdev_boot,
            F.sdev_boot
          ), lower.tail = TRUE))
        to_boot <-
          max(qnorm(1e-04, mean = c(M.mu_boot, F.mu_boot), sd = c(
            M.sdev_boot,
            F.sdev_boot
          ), lower.tail = FALSE))
        abs.dens_boot <- function(x) {
          abs(dnorm(x, c(M.mu_boot, F.mu_boot)[1], c(M.sdev_boot, F.sdev_boot)[1]) -
            dnorm(x, c(M.mu_boot, F.mu_boot)[2], c(M.sdev_boot, F.sdev_boot)[2]))
        }
        D_boot <-
          integrate(Vectorize(abs.dens_boot),
            lower = from_boot, upper =
              to_boot
          )$val / 2
        sto_boot[i] <- D_boot
        cat(paste("\r", i, " of ", B, "bootstraps"))
        flush.console()
      }
      half_CI <- CI / 2
      bot <- 0.5 - half_CI
      top <- 0.5 + half_CI
      b <- -qnorm(sum(sto_boot < D) / B)
      alpha.1 <- pnorm(qnorm(bot) - 2 * b)
      alpha.2 <- pnorm(qnorm(top) - 2 * b)
      bounds <- quantile(sto_boot, c(alpha.1, alpha.2))
      bounds <- as.numeric(bounds)
      upper <- bounds[2]
      lower <- bounds[1]

      D_list <-
        list(D = cbind.data.frame(lower = lower, D = D, upper = upper), df = df)
    } else {
      D_list <-
        list(D = data.frame(D = D), df = df)
    }
    D_list
  }
  vec_inBr <- Vectorize(InBr.D, vectorize.args = "i", SIMPLIFY = F)
  vec_inBr(i = seq_len(nrow(x))) -> all_list
  D_list <- lapply(all_list, function(x) x[[1]])
  names(D_list) <- levels(x$Trait)
  df <- lapply(all_list, function(x) x[[2]])
  df <- do.call(rbind.data.frame, df)
  do.call(rbind, D_list) -> D_df
  fill_list <- switch(fill, female = list(
    col = c("pink1", "white"),
    col2 = c("pink1", "light blue"),
    sex_levels = c("F", "M")
  ),
  male = list(
    col = c("light blue", "white"),
    col2 = c("light blue", "pink1"),
    sex_levels = c("M", "F")
  ),
  both = list(
    col = c("pink1", "light blue"),
    col2 = c("pink1", "light blue"),
    sex_levels = c("F", "M")
  )
  )
  col <- fill_list$col
  col2 <- fill_list$col2
  df$sex <- factor(df$sex, levels = fill_list$sex_levels)
  names(col) <- levels(df$sex)
  names(col2) <- levels(df$sex)
  scale <- scale_fill_manual(name = "sex", values = col)
  scale2 <- scale_color_manual(name = "sex", values = col2)
  p <- ggplot(data = df, aes(
    x = .data$z,
    y = .data$dn,
    color = .data$sex
  )) +
    geom_polygon(aes(
      fill =
        .data$sex
    )) +
    scale +
    scale2 +
    geom_density(stat = "identity") +
    facet_wrap(~trait, scales = "free") +
    ylab("Density") +
    xlab("x") +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  D_df <- rown_col(as.data.frame(D_df), var = "Trait") %>% mutate(
    across(-1,round,digits)) %>% as.data.frame()
  if(isTRUE(plot)){
  list(D = D_df, plot = p)
  }else{

    return(D_df)

}
}
