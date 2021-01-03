#' @title Tg
#' @seealso
#'  [TestDimorph-deprecated()]
#' @name Tg-deprecated
#' @keywords internal
NULL
#' @rdname TestDimorph-deprecated
#' @section `Tg`:
#' For `Tg`, use [t_greene()].
#' @export
Tg <- function(x = NULL,
               Pop = 1,
               es = "none",
               plot = FALSE,
               ...,
               alternative = c("two.sided", "less", "greater"),
               padjust = "none",
               letters = FALSE,
               digits = 4,
               CI = 0.95) {
  if (!identical(Sys.getenv("TESTTHAT"), "true"))
  .Deprecated("t_greene")
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
  if (!(Pop %in% seq_along(x))) {
    stop("Pop should be number from 1 to ncol(x)")
  }
  if (nrow(x) < 2) {
    stop("x should at least have 2 rows")
  }
  padjust <- match.arg(padjust, choices = p.adjust.methods)
  x <- x %>%
    drop_na() %>%
    as.data.frame()
  x$Pop <- x[, Pop]
  if (length(dplyr::contains("-", vars = x$Pop)) != 0) {
    x$Pop <-
      gsub(
        x = x$Pop,
        pattern = "[^[:alnum:]]",
        replacement = "_"
      )
  }
  x$Pop <- factor(x$Pop,levels = x$Pop)
  x$Pop <- droplevels(x$Pop)
  if (length(unique(x$Pop)) != length(which(!is.na(x$Pop)))) {
    warning("Population names are not unique")
  }
  pairs <- utils::combn(x$Pop, 2, simplify = FALSE)

  names(pairs) <- lapply(pairs, paste, collapse = "-")

  tg <-
    lapply(pairs, function(y) {
      t_test(
        m = x[y[1], "m"],
        f = x[y[1], "f"],
        m2 = x[y[2], "m"],
        f2 = x[y[2], "f"],
        M.mu = x[y[1], "M.mu"],
        F.mu = x[y[1], "F.mu"],
        M.mu2 = x[y[2], "M.mu"],
        F.mu2 = x[y[2], "F.mu"],
        M.sdev = x[y[1], "M.sdev"],
        F.sdev = x[
          y[1],
          "F.sdev"
        ],
        M.sdev2 = x[y[2], "M.sdev"],
        F.sdev2 = x[y[2], "F.sdev"],
        padjust = padjust,
        N = ((
          nlevels(x$Pop)^2 - nlevels(x$Pop)
        ) / 2),
        alternative = alternative,
        CI = CI,
        digits = digits,
        es = es
      )
    })
  tg <- do.call(rbind.data.frame, tg)
  tg <- rown_col(tg, var = "populations")

  # Pairwise comparisons and corrplot ---------------------------------------

  pval <- tg$p.value
  names(pval) <- tg$populations
  pmatrix <- multcompView::vec2mat(pval)
  if (!is.logical(letters)) {
    stop("letters should be either TRUE or FALSE")
  }
  if (isTRUE(letters)) {
    tg <-
      list(
        "t.test" = tg,
        "pairwise letters" = rown_col(
          data.frame(
            "letters" = multcompView::multcompLetters(pval,
                                                      threshold = CI
            )[[1]]
          ),
          var = "populations"
        )
      )
  } else {
    tg <- tg
  }
  if (!is.logical(plot)) {
    stop("plot should be either TRUE or FALSE")
  }
  if (isTRUE(plot)) {
    plot_list <- structure(list(
      "t.greene" = tg,
      corrplot::corrplot(
        corr = pmatrix,
        p.mat = pmatrix,
        sig.level = 1-CI,
        number.digits = digits,
        is.corr = FALSE,
        ...
      ),
      class = "tg"
    ))
    print_tg <- function(x) {
      x[[1]]
    }

    print_tg(plot_list)
  } else {
    tg
  }
}
