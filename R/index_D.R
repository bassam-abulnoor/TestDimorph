#' @title Dissimilarity index
#' @description FUNCTION_DESCRIPTION
#' @inheritParams t_greene
#' @param plot Logical; if true a plot of densities is returned, Default:
#' FALSE
#' @param fill Specify which sex to be filled with color in the plot; either
#' "male" in blue color, "female" in pink color or "both", Default: 'female'
#' @param r Number of row for a desired trait, Default: NULL
#' @param Trait Number of column containing trait names as factors, Default: 1
#' @return A tibble and a graphical representation of the chosen traits and
#' their corresponding dissimilarity indices
#' @details DETAILS
#' @rdname index_D
#' @import dplyr
#' @import ggplot2
#' @importFrom stats dnorm integrate na.omit qnorm
#' @importFrom tidyr drop_na pivot_longer
<<<<<<< HEAD
#' @keywords internal
=======
>>>>>>> 6dea5f88f050c15e544b721a245e69cf92b48861
index_D <- function(x,
                    plot = FALSE,
                    fill = c("female", "male", "both"),
                    r = NULL,
                    Trait = 1) {
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
  fill <- match.arg(fill, choices = c("female", "male", "both"))
  InBr.D <- function(Trait, m, M.mu, M.sdev, f, F.mu, F.sdev) {
    from <-
      min(qnorm(
        1e-04,
        mean = c(M.mu, F.mu),
        sd = c(M.sdev, F.sdev),
        lower.tail = TRUE
      ))
    to <-
      max(qnorm(
        1e-04,
        mean = c(M.mu, F.mu),
        sd = c(M.sdev, F.sdev),
        lower.tail = FALSE
      ))
    abs.dens <- function(x) {
      abs(dnorm(x, c(M.mu, F.mu)[1], c(M.sdev, F.sdev)[1]) - dnorm(x, c(M.mu, F.mu)[2], c(M.sdev, F.sdev)[2]))
    }
    D <-
      integrate(Vectorize(abs.dens), lower = from, upper = to)$val /
        2
    z <- seq(from, to, (to - from) / 1000)
    trait <- rep(Trait, length(z))
    list_D <-
      list(z, D, F.mu, M.mu, F.sdev, M.sdev, trait)

    names(list_D) <-
      c(
        "z",
        "D",
        "F.mu",
        "M.mu",
        "F.sdev",
        "M.sdev",
        "trait"
      )
    list_D
  }
  D1_fun <- function(x) {
    m <- x$m[1]
    M.mu <- x$M.mu[1]
    M.sdev <- x$M.sdev[1]
    f <- x$f[1]
    F.mu <- x$F.mu[1]
    F.sdev <- x$F.sdev[1]
    Trait <- x$Trait[1]
    InBr.D(Trait, m, M.mu, M.sdev, f, F.mu, F.sdev)
  }
  D2_fun <- function(x) {
    df <-
      cbind_fill(x[["z"]], x[["trait"]], x[["M.mu"]], x[["M.sdev"]], x[["F.mu"]], x[["F.sdev"]], x[["D"]])
    names(df) <-
      c(
        "z",
        "trait",
        "M.mu",
        "M.sdev",
        "F.mu",
        "F.sdev",
        "D"
      )
    df
  }
  m_fun <- function(x) {
    df <-
      dnorm(x[["z"]], as.numeric(na.omit(x[["M.mu"]])), as.numeric(na.omit(x[["M.sdev"]])))
    df <-
      cbind.data.frame(x[["z"]], df, x[["trait"]], "M")
    names(df) <- c("z", "dn", "trait", "sex")
    df
  }
  f_fun <- function(x) {
    df <-
      dnorm(x[["z"]], as.numeric(na.omit(x[["F.mu"]])), as.numeric(na.omit(x[["F.sdev"]])))
    df <-
      cbind.data.frame(x[["z"]], df, x[["trait"]], "F")
    names(df) <- c("z", "dn", "trait", "sex")
    df
  }
  if (!is.null(r)) {
    x <- x[r, ]
  }
  x <- data.frame(x)
  x$Trait <- x[, Trait]
  x$Trait <- factor(x$Trait)
  x$Trait <- droplevels(x$Trait)
  x <- dplyr::arrange(x, Trait)
  D1_list <- lapply(split.data.frame(x, x$Trait), D1_fun)
  D2_list <- lapply(D1_list, D2_fun)
  male <- do.call(rbind.data.frame, lapply(D2_list, m_fun))
  female <- do.call(rbind.data.frame, lapply(D2_list, f_fun))
  df <- rbind.data.frame(male, female)

  if (fill == "female") {
    col <- c("pink1", "white")
    col2 <- c("pink1", "light blue")
    df$sex <- factor(df$sex, levels = c("F", "M"))
  }
  if (fill == "male") {
    col <- c("light blue", "white")
    col2 <- c("light blue", "pink1")
    df$sex <- factor(df$sex, levels = c("M", "F"))
  }
  if (fill == "both") {
    col <- c("pink1", "light blue")
    col2 <- c("pink1", "light blue")
    df$sex <- factor(df$sex, levels = c("F", "M"))
  }
  names(col) <- levels(df$sex)
  names(col2) <- levels(df$sex)
  scale <- scale_fill_manual(name = "sex", values = col)
  scale2 <- scale_color_manual(name = "sex", values = col2)
  p <-
    ggplot(data = df, aes(
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

  extract_D <- function(x) {
    as.vector(round(x[["D"]], 4))
  }
  out <-
    cbind.data.frame(levels(x$Trait), do.call(rbind.data.frame, lapply(D1_list, extract_D)))
  names(out) <- c("Trait", "D")
  out$Trait <- factor(out$Trait)
  out <- as_tibble(out)
  if (isTRUE(plot)) {
    list(out, p)
  } else {
    return(out)
  }
}
