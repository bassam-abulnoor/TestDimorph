#' @title Raw Data Generation By Log-normal Or Truncated Distribution
#' @description  Generates raw data from summary statistics using
#'   uni/multivariate log/truncated normal distribution
#' @inheritParams multivariate
#' @param dist univariate distribution used for data generation either `log` for
#'   log-normal or `trunc` for truncated, Default: 'trunc'
#' @param lower vector of lower bounds, Default: -Inf
#' @param upper vector of upper bounds, Default: Inf
#' @param format form of the resultant tibble either 'long' or 'wide', Default:
#'   'wide'
#' @param complete_cases Logical; if TRUE rows with missing values will be
#'   removed, Default: FALSE
#' @return tibble of raw data
#' @details If data generation is desired using multivariate distribution data
#'   is entered in the form of a list of summary statistics and pooled within
#'   correlational matrix as in [baboon.parms_list], or the summary statistics
#'   are entered separately in the form of a data frame/tibble as in
#'   [baboon.parms_df] with a separate correlational matrix as in [R]. If data
#'   frame/tibble is entered without a correlational matrix, data generation is
#'   carried out using univariate distribution.
#'   N.B: Transformation of raw summary data to logged data is only possible for
#'   univariate distribution and if multivariate log-normal distribution is desired
#'   logged values should be entered directly with `dist` set to `trunc`.
#' @examples
#'  # Data generation using univariate distribution
#'  library(TestDimorph)
#'  RawGen(baboon.parms_df)
#'  # Data generation using multivariate distribution
#'  library(TestDimorph)
#'  RawGen(baboon.parms_list)
#' @rdname RawGen
#' @export
#' @importFrom purrr map
#' @importFrom truncnorm rtruncnorm
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble is_tibble
#' @importFrom rlang abort
#' @importFrom tidyr drop_na
#' @importFrom reshape2 melt
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom stats rlnorm
#' @references \insertRef{HUSSEIN2019}{TestDimorph}
RawGen <- function(x,
                   Parms = 1,
                   Pop = 2,
                   R.res = NULL,
                   dist = "trunc",
                   lower = -Inf,
                   upper = Inf,
                   format = "wide",
                   complete_cases = FALSE)
{
  if (!(is.list(x) || is.data.frame(x) || tibble::is_tibble(x))) {
    rlang::abort("x must be a list a dataframe or a tibble")
  }
  if (!(dist %in% c("log", "trunc")))   {
    rlang::abort("distribution should be one of `log` or `trunc`")

  }
  if (!(format %in% c("wide", "long")))   {
    rlang::abort("format should be one of `wide` or `long`")

  }
  if (!(complete_cases %in% c(TRUE, FALSE)))   {
    rlang::abort("complete_cases should be either TRUE or FALSE")

  }
  mult <- function(y,
                   format,
                   complete_cases,
                   dist = dist,
                   lower = lower,
                   upper = upper) {
    if (dist == "log")   {
      rlang::abort(
        "Transformation of raw summary data to logged data is only possible for univariate distribution"
      )
    }
    R <- y$R.res
    m <- y$m
    f <- y$f
    M <- y$M.mu
    F <- y$F.mu
    M.sdev <- y$M.sdev
    F.sdev <- y$F.sdev
    traits <- colnames(M)
    n.t <- length(traits)
    pops <- row.names(M)
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
      X[start:stop,] <-
        tmvtnorm::rtmvnorm(
          n = m[i],
          mean = M[i,],
          sigma = V,
          lower = rep(lower, n.t),
          upper = rep(upper, n.t)
        )
      start <- stop + 1
      S <- diag(F.sdev[i, ])
      V <- S %*% R %*% S
      stop <- stop + f[i]
      X[start:stop,] <-
        tmvtnorm::rtmvnorm(
          n = f[i],
          mean = F[i,],
          sigma = V,
          lower = rep(lower, n.t),
          upper = rep(upper, n.t)
        )
    }
    colnames(X) <- traits
    X <- data.frame(Sex, Pop, X)

    if (format == "wide") {
      if (complete_cases == TRUE) {
        return(tidyr::drop_na(X))
      } else{
        return(tibble::as_tibble(X))
      }
    } else{
      tibble::as_tibble(reshape2::melt(
        X,
        id = c("Sex", "Pop"),
        variable = "Parms",
        na.rm = complete_cases
      ))


    }
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
    if (is.null(R.res)) {
      x <- data.frame(x)
      x$Pop <- x[, Pop]
      x$Pop <- factor(x$Pop)
      x$Parms <- x[, Parms]
      x$Parms <- factor(x$Parms)
      M <- function(x) {
        df <- by(x, list(x$Parms), list)
        if (dist == "log") {
          df <- do.call(cbind_fill, c(lapply(
            purrr::map(df, function(x) {
              rlnorm(
                n = x$m[1],
                meanlog = log(x$M.mu ^ 2 / sqrt(x$M.sdev ^ 2 + x$M.mu ^ 2)),
                sdlog = sqrt(log(1 + (
                  x$M.sdev ^ 2 / x$M.mu ^ 2
                )))
              )
            }), as.data.frame
          )))
        } else{
          df <- do.call(cbind_fill, c(lapply(
            purrr::map(df, function(x) {
              truncnorm::rtruncnorm(
                n = x$m[1],
                a = lower,
                b = upper,
                mean = x$M.mu[1],
                sd = x$M.sdev[1]
              )
            }), as.data.frame
          )))
        }
        colnames(df) <- levels(x$Parms)
        return(df)
      }
      male <- do.call(rbind, by(x, x$Pop, M))
      male$Pop <-
        as.factor(stringr::str_split(
          rownames(male),
          pattern = "\\.",
          n = 2,
          simplify = T
        )[, 1])
      male$Sex <- as.factor(rep("M", nrow(male)))
      male <-
        male[, c(ncol(male), ncol(male) - 1, 1:nlevels(x$Parms))]

      F <- function(x) {
        df <- by(x, list(x$Parms), list)
        if (dist == "log") {
          df <-
            do.call(cbind_fill, c(lapply(
              purrr::map(df, function(x) {
                rlnorm(
                  n = x$f[1],
                  meanlog = log(x$F.mu ^ 2 / sqrt(x$F.sdev ^ 2 + x$F.mu ^ 2)),
                  sdlog = sqrt(log(1 + (
                    x$F.sdev ^ 2 / x$F.mu ^ 2
                  )))
                )
              }), as.data.frame
            )))
        } else{
          df <-
            do.call(cbind_fill, c(lapply(
              purrr::map(df, function(x) {
                truncnorm::rtruncnorm(
                  n = x$f[1],
                  a = lower,
                  b = upper,
                  mean = x$F.mu[1],
                  sd = x$F.sdev[1]
                )
              }), as.data.frame
            )))
        }
        colnames(df) <- levels(x$Parms)
        return(df)
      }
      female <- do.call(rbind, by(x, x$Pop, F))
      female$Pop <-
        as.factor(stringr::str_split(
          rownames(female),
          pattern = "\\.",
          n = 2,
          simplify = TRUE
        )[, 1])
      female$Sex <- as.factor(rep("F", nrow(female)))
      female <-
        female[, c(ncol(female) , ncol(female) - 1, 1:nlevels(x$Parms))]
      wide <- tibble::as_tibble(rbind.data.frame(male, female))
      if (format == "wide") {
        if (complete_cases == TRUE) {
          return(tidyr::drop_na(wide))
        } else{
          return(wide)
        }
      }
      if (format == "long") {
        long <- tibble::as_tibble(reshape2::melt(
          wide,
          id = c("Sex", "Pop"),
          variable = "Parms",
          na.rm = complete_cases
        ))

        return(long)
      }
    }
    if (!is.null(R.res)) {
      if (!is.matrix(R.res)) {
        rlang::abort("R.res must be a matrix")

      }
      x <- data.frame(x)
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
      mult(
        y = x,
        format = format,
        complete_cases = complete_cases,
        dist = dist,
        lower = lower,
        upper = upper
      )
    }
  }
  if (!(is.data.frame(x) || tibble::is_tibble(x))) {
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
    mult(
      y = x,
      format = format,
      complete_cases = complete_cases,
      dist = dist,
      lower = lower,
      upper = upper
    )

  }
}
