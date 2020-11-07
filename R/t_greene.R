#' @title Greene t test of Sexual Dimorphism
#' @description Calculation and visualization of the differences in degree
#' sexual dimorphism between two populations using summary statistics as
#' input.
#' @param x A data frame containing summary statistics.
#' @param Pop Number of the column containing populations' names, Default: 1
#' @param es Logical; if TRUE effect size is included in the output , Default:
#' FALSE
#' @param plot Logical; if TRUE graphical matrix of p values, Default: FALSE
#' @param ... additional arguments that can be passed to
#' [corrplot][corrplot::corrplot] function.
#' @param alternative a character string specifying the alternative
#' hypothesis, must be one of "two.sided", "greater" or "less".
#' @param padjust Method of p.value adjustment for multiple comparisons
#' following [p.adjust.methods]
#' @param letters Logical; if TRUE returns letters for pairwise comparisons
#' where significantly different populations are given different letters,
#' Default: FALSE'
#' @param digits Number of significant digits, Default: 4
#' @param sig.level Critical p.value, Default: 0.05
#' @return Tibble of t.test results
#' @details The input is a data frame of summary statistics where the column
#' containing population names is chosen by position (first by default), other
#' columns of summary data should have specific names (case sensitive) similar
#' to [baboon.parms_df]
#' @examples
#' \donttest{
#' # Comparisons of femur head diameter in four populations
#' library(TestDimorph)
#' df <- data.frame(
#'   Pop = c("Turkish", "Bulgarian", "Greek", "Portuguese "),
#'   m = c(150.00, 82.00, 36.00, 34.00),
#'   f = c(150.00, 58.00, 34.00, 24.00),
#'   M.mu = c(49.39, 48.33, 46.99, 45.20),
#'   F.mu = c(42.91, 42.89, 42.44, 40.90),
#'   M.sdev = c(3.01, 2.53, 2.47, 2.00),
#'   F.sdev = c(2.90, 2.84, 2.26, 2.90)
#' )
#' t_greene(
#'   df,
#'   plot = TRUE,
#'   method = "ellipse",
#'   padjust = "none",
#'   type = "lower",
#'   col = c(
#'     "#AEB6E5",
#'     "#B1A0DB",
#'     "#B788CD",
#'     "#BC6EB9",
#'     "#BC569E",
#'     "#B6407D",
#'     "#A93154"
#'   ),
#'   tl.cex = 0.8,
#'   tl.col = "black",
#'   insig =
#'     "label_sig",
#'   tl.srt = 0.1,
#'   pch.cex = 2.5,
#'   tl.pos = "ld",
#'   win.asp = 1,
#'   number.cex = 0.5,
#'   na.label = "NA"
#' )
#' }
#' @seealso
#' [multcompView::multcompLetters()]
#' [corrplot::corrplot()]
#' @rdname t_greene
#' @export
#' @importFrom stats qt pt
#' @importFrom utils combn
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom multcompView multcompLetters vec2mat
#' @importFrom corrplot corrplot
#' @importFrom dplyr contains

t_greene <- function(x,
                     Pop = 1,
                     es = FALSE,
                     plot = FALSE,
                     ...,
                     alternative = c("two.sided", "less", "greater"),
                     padjust = p.adjust.methods,
                     letters = FALSE,
                     digits = 4,
                     sig.level = 0.05) {
  # t-test for a data.frame -------------------------------------------------

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
  x <- data.frame(x)
  x$Pop <- x[, Pop]
  if (length(dplyr::contains("-", vars = x$Pop)) != 0) {
    x$Pop <-
      gsub(
        x = x$Pop,
        pattern = "-",
        replacement = "_"
      )
  }
  x$Pop <- factor(x$Pop)
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
        sig.level = sig.level,
        digits = digits,
        es = es
      )
    })
  tg <- do.call(rbind.data.frame, tg)
  tg <- tibble::rownames_to_column(tg, var = "populations")

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
        "t.test" = tibble::as_tibble(tg),
        "pairwise letters" = tibble::rownames_to_column(
          data.frame(
            "letters" = multcompView::multcompLetters(pval,
              threshold = sig.level
            )[[1]]
          ),
          var = "populations"
        )
      )
  } else {
    tg <- tibble::as_tibble(tg)
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
        sig.level = sig.level,
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
