#' @title RawGen
#' @seealso
#'  [TestDimorph-deprecated()]
#' @name RawGen-deprecated
#' @keywords internal
NULL
#' @rdname TestDimorph-deprecated
#' @section `RawGen`:
#' For `RawGen`, use [raw_gen()].
#' @export
RawGen <- function(x,
                   Trait = 1,
                   Pop = 2,
                   R.res = NULL,
                   dist = c("truncated", "log"),
                   lower = -Inf,
                   upper = Inf,
                   format = c("wide", "long"),
                   complete_cases = FALSE) {
  if (!identical(Sys.getenv("TESTTHAT"), "true"))
  .Deprecated("raw_gen")
  if (!(is.list(x) || is.data.frame(x))) {
    stop("x should be a list or a dataframe")
  }
  dist <- match.arg(dist, choices = c("truncated", "log"))
  format <- match.arg(format, choices = c("wide", "long"))
  if (!is.logical(complete_cases)) {
    stop("complete_cases should be either TRUE or FALSE")
  }
  # univariate data generation ----------------------------------------------

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
    if (!(Trait %in% seq_along(x))) {
      stop("Trait should be number from 1 to ncol(x)")
    }
    if (!(Pop %in% seq_along(x))) {
      stop("Pop should be number from 1 to ncol(x)")
    }
    if (is.null(R.res)) {
      x <- x %>%
        drop_na() %>%
        as.data.frame()
      x$Pop <- x[, Pop]
      x$Pop <- factor(x$Pop)
      x$Trait <- x[, Trait]
      x$Trait <- factor(x$Trait)

      # Data generation --------------------------------------------------

      if (dist == "log") {
        message("Data generation was done using univariate log distribution")
        gen_m <- function(x) {
          rlnorm(
            n = x$m[1],
            meanlog = log(x$M.mu^2 / sqrt(x$M.sdev^2 + x$M.mu^2)),
            sdlog = sqrt(log(1 + (
              x$M.sdev^2 / x$M.mu^2
            )))
          )
        }
        gen_f <- function(x) {
          rlnorm(
            n = x$f[1],
            meanlog = log(x$F.mu^2 / sqrt(x$F.sdev^2 + x$F.mu^2)),
            sdlog = sqrt(log(1 + (
              x$F.sdev^2 / x$F.mu^2
            )))
          )
        }
      } else {
        message("Data generation was done using univariate truncated distribution")
        gen_m <- function(x) {
          truncnorm::rtruncnorm(
            n = x$m[1],
            a = lower,
            b = upper,
            mean = x$M.mu[1],
            sd = x$M.sdev[1]
          )
        }
        gen_f <- function(x) {
          truncnorm::rtruncnorm(
            n = x$f[1],
            a = lower,
            b = upper,
            mean = x$F.mu[1],
            sd = x$F.sdev[1]
          )
        }
      }
      m_function <- function(x) {
        df <- by(x, list(x$Trait), list)
        df <- lapply(df, gen_m)
        df <- lapply(df, as.data.frame)
        df <- do.call(cbind_fill2, df)
        colnames(df) <- levels(x$Trait)
        df
      }
      f_function <- function(x) {
        df <- by(x, list(x$Trait), list)
        df <- lapply(df, gen_f)
        df <- lapply(df, as.data.frame)
        df <- do.call(cbind_fill2, df)
        colnames(df) <- levels(x$Trait)
        df
      }
      pops <- split.data.frame(x, x$Pop)
      male <- lapply(pops, m_function)
      male <- do.call(rbind.data.frame, male)
      female <- lapply(pops, f_function)
      female <- do.call(rbind.data.frame, female)
      males <- strsplit(rownames(male), split = "\\.")
      females <- strsplit(rownames(female), split = "\\.")
      male$Pop <- as.factor(do.call(rbind.data.frame, males)[, 1])
      male$Sex <- as.factor(rep("M", nrow(male)))
      female$Pop <-
        as.factor(do.call(rbind.data.frame, females)[, 1])
      female$Sex <- as.factor(rep("F", nrow(female)))
      male <-
        male[, c(ncol(male), ncol(male) - 1, seq(nlevels(x$Trait)))]


      female <-
        female[, c(ncol(female), ncol(female) - 1, seq(nlevels(x$Trait)))]

      # Joining both datasets ---------------------------------------------------

      wide <- rbind.data.frame(male, female)
      if (format == "wide") {
        if (isTRUE(complete_cases)) {
          return(tidyr::drop_na(wide))
        } else {
          return(wide)
        }
      }
      if (format == "long") {
        long <-
          pivot_longer(
            data = wide,
            cols = -c("Sex", "Pop"),
            names_to = "Trait",
            values_drop_na = complete_cases
          )

        return(long)
      }
    }

    # multivariate generation with data.frame and correlation matrix -----------

    if (!is.null(R.res)) {
      if (!is.matrix(R.res)) {
        stop("R.res should be a matrix")
      }
      x <- dataframe2list(
        x = x,
        R.res = R.res,
        Trait = Trait,
        Pop = Pop
      )
    }
  }

  # multivariate data generation with list input -------------------------------------

  if (!(is.data.frame(x))) {
    if (!all(c("M.mu", "F.mu", "M.sdev", "F.sdev", "m", "f", "R.res") %in% names(x))) {
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
    message("Data generation was done using multivariate truncated distribution")
    multi_raw(
      x = x,
      format = format,
      complete_cases = complete_cases,
      dist = dist,
      lower = lower,
      upper = upper
    )
  }
}
