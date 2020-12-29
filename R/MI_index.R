#' @title Mixture Index ("MI")
#' @description Ipina and Durand's (2010) mixture intersection (MI) measure of
#' sexual dimorphism.  This measure is an overlap coefficient where the sum of
#'  the frequency of males and the frequency of females equals 1.0.  Ipina and
#'  Durand (2010) also define a normal intersection (NI) measure which is the
#'  overlap coefficient of two normal distributions, equivalent to Inman and
#'  Bradley's (1989) overlap coefficient
#' @details see \link{D_index} for bootstrap method.
#' @inheritParams D_index
#' @param p.f proportion of sample that is female (if p.f>0 then
#' p.m=1-p.f, where p.m is the proportion of males and bootstrap won't
#' be available) , Default: 0
#' @param index_type type of coefficient (if "MI" it fits the mixture index.
#' If = "NI" it fits the overlap coefficient for two normal distributions,
#' which is equal to 1 – D_index, Default: 'MI'
#' @return returns a table of Ipina and Durand's (2010) mixture index ("MI")
#' for different traits with graphical representation.
#' @examples
#' library(TestDimorph)
#' data("Cremains_measurements")
#' # plot and test of significance
#' MI_index(Cremains_measurements[1, ], plot = TRUE)
#' #' #NI index
#' MI_index(Cremains_measurements[1, ], index_type = "NI")
#' \dontrun{
#' # confidence interval was bootstrapping
#' MI_index(Cremains_measurements[1, ], rand = FALSE, B = 1000)
#' }
#' @references Inman, H. F., & Bradley Jr, E. L. (1989). The overlapping
#' coefficient as a measure of agreement between probability distributions and
#' point estimation of the overlap of two normal densities. Communications in
#' Statistics-Theory and Methods, 18(10), 3851-3874.
#' @rdname MI_index
#' @import dplyr
#' @import ggplot2
#' @importFrom stats dnorm integrate na.omit qnorm contr.sum pnorm quantile
#' @importFrom tidyr drop_na
#' @importFrom utils flush.console
#' @export

MI_index <- function(x, plot = FALSE, fill = "both", Trait = 1, B = NULL, CI = 0.95,
                     p.f = 0, index_type = "MI", rand = TRUE) {
  fill <- match.arg(fill, choices = c("female", "male", "both"))
  index_type <- match.arg(index_type, choices = c("MI", "NI"))
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
  if (!is.numeric(p.f)) {
    stop("p.f should be a number")
  }

  if (is.null(CI) && !is.null(B)) {
    warning("confidence level is not spicified")
  }
  if (isFALSE(rand)) {
    set.seed(42)
  }
  x %>% drop_na() %>% as.data.frame()
  m <- x$m
  M.mu <- x$M.mu
  M.sdev <- x$M.sdev
  f <- x$f
  F.mu <- x$F.mu
  F.sdev <- x$F.sdev
  Trait <- x[, Trait]
  if (p.f == 0) {
    p.f <- f / (f + m)
    p.m <- 1 - p.f
  }

  else {
    p.m <- 1 - p.f
  }

  if (index_type == "NI") {
    p.f <- 1
    p.m <- 1
  }

  Ipina_Durand <- function(Trait, m, M.mu, M.sdev, f, F.mu, F.sdev, i) {
    m <- x$m[i]
    M.mu <- x$M.mu[i]
    M.sdev <- x$M.sdev[i]
    f <- x$f[i]
    F.mu <- x$F.mu[i]
    F.sdev <- x$F.sdev[i]
    Trait <- x$Trait[i]
    left <- min(c(qnorm(0.000001, F.mu, F.sdev), qnorm(0.000001, M.mu, M.sdev)))
    right <- max(c(qnorm(0.999999, F.mu, F.sdev), qnorm(0.999999, M.mu, M.sdev)))

    ID <- function(x) min(c(p.f * dnorm(x, F.mu, F.sdev), p.m * dnorm(x, M.mu, M.sdev)))

    MI <- round(integrate(Vectorize(ID), left, right)$val, 4)
    z <- seq(left, right, 0.01)
    trait <- rep(Trait, length(z))
    dn_male <- p.m * dnorm(z, as.numeric(na.omit(M.mu)), as.numeric(na.omit(
      M.sdev
    ))) %>% as.data.frame()
    dn_male <- cbind.data.frame(z = z, dn = dn_male, sex = rep("M", nrow(dn_male)))
    dn_female <- p.f * dnorm(z, as.numeric(na.omit(F.mu)), as.numeric(na.omit(
      F.sdev
    ))) %>% as.data.frame()
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
        left_boot <- min(c(qnorm(0.000001, F.mu_boot, F.sdev_boot), qnorm(
          0.000001, M.mu_boot,
          M.sdev_boot
        )))
        right_boot <- max(c(qnorm(0.999999, F.mu_boot, F.sdev_boot), qnorm(
          0.999999, M.mu_boot,
          M.sdev_boot
        )))
        ID_boot <- function(x) {
          min(c(p.f * dnorm(x, F.mu_boot, F.sdev_boot), p.m * dnorm(
            x,
            M.mu_boot, M.sdev_boot
          )))
        }
        MI_boot <- round(integrate(Vectorize(ID_boot), left_boot, right_boot)$val, 4)
        sto_boot[i] <- MI_boot
        cat(paste("\r", i, " of ", B, "bootstraps"))
        flush.console()
      }
      half_CI <- CI / 2
      bot <- 0.5 - half_CI
      top <- 0.5 + half_CI
      b <- -qnorm(sum(sto_boot < MI) / B)
      alpha.1 <- pnorm(qnorm(bot) - 2 * b)
      alpha.2 <- pnorm(qnorm(top) - 2 * b)
      bounds <- quantile(sto_boot, c(alpha.1, alpha.2))
      bounds <- as.numeric(bounds)
      upper <- bounds[2]
      lower <- bounds[1]

      MI <- cbind.data.frame(lower = lower, MI, upper = upper)
      names(MI)[2] <- index_type
    } else {
      MI <- data.frame(MI)
      names(MI)[1] <- index_type
    }
    IM_list <- list(MI, df = df)
    IM_list
  }
  vec_Ipina_Durand <- Vectorize(Ipina_Durand, vectorize.args = "i", SIMPLIFY = F)
  all_list <- vec_Ipina_Durand(i = seq_len(nrow(x)))
  IM_list <- lapply(all_list, function(x) x[[1]])
  x$Trait <- factor(x$Trait, levels = x$Trait)
  names(IM_list) <- levels(x$Trait)
  df <- lapply(all_list, function(x) x[[2]])
  df <- do.call(rbind.data.frame, df)
  IM_df <- do.call(rbind, IM_list)
  fill_list <- switch(fill,
    female = list(
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
  IM_df <- rownames_to_column(as.data.frame(IM_df), var = "Trait")
  if (isTRUE(plot)) {
    list(index = as_tibble(IM_df), plot = p)
  } else {
    as_tibble(IM_df)
  }
}
