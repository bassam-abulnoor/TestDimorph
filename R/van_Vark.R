#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @inheritParams extract_sum
#' @param N Number of column containing sample size, Default: 3
#' @param W Pooled within-group variance-covariance matrix
#' @param q Number of canonical variates to retain for chi square test,
#' Default: 2
#' @param plot PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname van_Vark
#' @import ggplot2
#' @import dplyr
#' @importFrom stats pchisq
#' @importFrom tibble is_tibble rownames_to_column
#' @importFrom utils combn
#' @keywords internal
van_Vark <- function(x, Sex = 1, Pop = 2, N = 3, W, q = 2, plot = TRUE) {
  if (!(is.data.frame(x) || tibble::is_tibble(x))) {
    stop("x  should be a tibble or a dataframe")
  }
  if (!(Sex %in% seq_along(x))) {
    stop("Sex should be number from 1 to ncol(x)")
  }
  if (!(Pop %in% seq_along(x))) {
    stop("Pop should be number from 1 to ncol(x)")
  }
  if (!(N %in% seq_along(x))) {
    stop("N should be number from 1 to ncol(x)")
  }
  if (!is.logical(plot)) {
    stop("plot should be either TRUE or FALSE")
  }
  if (q == 1 && isTRUE(plot)) {
    warning("plot can't be produced if q=1")
  }
  x <- data.frame(x)
  x$N <- x[, N]
  x$Pop <- x[, Pop]
  x$Pop <- factor(x$Pop)
  x$Pop <- droplevels(x$Pop)
  levels(x$Pop) <- sort(levels(x$Pop))
  x$Sex <- x[, Sex]
  x$Sex <- factor(x$Sex, levels = c("F", "M"))
  x <- arrange(x, x$Pop, x$Sex)
  p <- NCOL(x) - 3
  g <- NROW(x)
  Rank <- min(c(g - 1, p))
  if (q > Rank) {
    q <- Rank
    warning("q has been decreased to match the matrix rank")
  }

  calc.B0 <- function() {
    y <- x[, -c(Pop, Sex, N)]
    G.mean <- apply(y, 2, mean)
    B0 <- matrix(0, ncol = p, nrow = p)
    for (i in 1:g) {
      d <- as.numeric(y[i, ]) - G.mean
      B0 <- B0 + (d %*% t(d))
    }
    B0 <- B0 / (g - 1)
    return(B0)
  }

  B0 <- calc.B0()
  eig.struc <- eigen(solve(W) %*% B0)
  e.vecs <- Re(eig.struc$vec[, 1:Rank])
  e.vecs <- e.vecs / (rep(1, p) %o% sqrt(diag(t(e.vecs) %*% W %*% e.vecs)))

  Center <- diag(g) - 1 / g
  y <- x[, -c(Pop, Sex, N)]
  X <- Center %*% as.matrix(y)
  CVs <- X %*% e.vecs[, 1:min(c(q, Rank))]
  CVs <- data.frame(CVs)
  names(CVs) <- paste0("x", seq_along(CVs))
  CVs$Pop <- factor(x$Pop)
  CVs$Sex <- factor(x$Sex)

  levels(CVs$Pop) <- sort(levels(CVs$Pop))
  CVs <- arrange(CVs, CVs$Pop, CVs$Sex)
  line1 <- vector()
  line2 <- vector()
  point1 <- vector()
  point2 <- vector()

  for (i in seq(1, (g - 1), 2)) {
    line1[i] <- CVs[i:(i + 1), 1][1]
    line2[i] <- CVs[i:(i + 1), 2][1]
  }
  for (i in seq(2, g, 2)) {
    point1[i] <- CVs[i, 1]
    point2[i] <- CVs[i, 2]
  }
  line1 <- line1 %>%
    as.data.frame() %>%
    mutate(n = row_number()) %>%
    filter(n %% 2 != 0) %>%
    select(-n) %>%
    drop_na()
  line2 <- line2 %>%
    as.data.frame() %>%
    mutate(n = row_number()) %>%
    filter(n %% 2 != 0) %>%
    select(-n) %>%
    drop_na()
  point1 <- point1 %>%
    as.data.frame() %>%
    mutate(n = row_number()) %>%
    filter(n %% 2 == 0) %>%
    select(-n) %>%
    drop_na()
  point2 <- point2 %>%
    as.data.frame() %>%
    mutate(n = row_number()) %>%
    filter(n %% 2 == 0) %>%
    select(-n) %>%
    drop_na()
  male <- seq(2, g, 2)
  female <- seq(1, (g - 1), 2)

  CVs <- cbind_fill(CVs, line1, line2, point1, point2, female, male)

  names(CVs)[(q + 3):ncol(CVs)] <- c(
    "line1", "line2", "point1", "point2",
    "female", "male"
  )
  p <- ggplot(CVs, aes(x = .data$x1, y = .data$x2, color = Sex, fill = Sex)) +
    geom_segment(aes(
      xend = point1,
      yend = point2, x = line1, y = line2
    ), color = "black", na.rm = TRUE) +
    geom_point(aes(x = point1, y = point2), color = "blue", na.rm = TRUE) +
    geom_point(aes(x = line1, y = line2), color = "red", na.rm = TRUE) +
    geom_text(
      na.rm = TRUE, aes(x = point1, y = point2, label = male),
      nudge_y = 0.2, color = "blue", check_overlap = T
    ) +
    geom_text(aes(
      x = line1,
      y = line2, label = female
    ),
    check_overlap = T, color = "red", na.rm = TRUE,
    nudge_y = 0.2
    ) +
    scale_color_manual(values = c("red", "blue")) +
    scale_fill_manual(values = c("red", "blue")) +
    guides(fill = guide_legend(override.aes = list(color = c(
      "red",
      "blue"
    )))) +
    xlab("First Canonical Axis") +
    ylab("Second Canonical Axis") +
    theme_minimal() +
    theme(legend.title = element_blank())
  pairs <- utils::combn(levels(x$Pop), 2, simplify = FALSE)

  names(pairs) <- lapply(pairs, paste, collapse = "-")
  chi <- function(y) {
    dif <- as.numeric(CVs[CVs$Pop == y[1] & CVs$Sex == "M", 1:q]) +
      as.numeric(CVs[CVs$Pop == y[2] & CVs$Sex == "F", 1:q]) - as.numeric(CVs[CVs$Pop ==
        y[1] & CVs$Sex == "F", 1:q]) - as.numeric(CVs[CVs$Pop == y[2] &
        CVs$Sex == "M", 1:q])
    X2 <- as.numeric(t(dif) %*% dif) / (1 / x[x$Pop == y[1] & x$Sex ==
      "F", N] + 1 / x[x$Pop == y[1] & x$Sex == "M", N] + 1 / x[x$Pop ==
      y[2] & x$Sex == "F", N] + 1 / x[
      x$Pop == y[2] & x$Sex == "M",
      N
    ])
    df <- min(c(q, Rank))
    p <- round(pchisq(X2, df, lower.tail = FALSE), 4)
    out <- tibble(X2, df, p)
    out$signif <- ifelse(out$p < 0.05, yes = "*", no = "ns")
    names(out) <- c("statistic", "df", "p.value", "signif")
    return(out)
  }
  pairs_list <- lapply(pairs, chi)
  pairs_df <- do.call(rbind.data.frame, pairs_list)
  pairs_df <- tibble::rownames_to_column(pairs_df, var = "populations")
  if (isTRUE(plot) && q > 1) {
    list(pairs_df, p)
  } else {
    pairs_df
  }
}
