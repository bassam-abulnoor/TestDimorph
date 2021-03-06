#' @title aovSS
#' @seealso
#'  [TestDimorph-deprecated()]
#' @name aovSS-deprecated
#' @keywords internal
NULL
#' @rdname TestDimorph-deprecated
#' @section `aovSS`:
#' For `aovSS`, use [aov_ss()]
#' @export
aovSS <-
  function(x,
           Pop = 1,
           pairwise = TRUE,
           letters = FALSE,
           es_anova = "none",
           digits = 4,
           CI = 0.95) {
    if (!identical(Sys.getenv("TESTTHAT"), "true")) {
      .Deprecated("aov_ss")
    }
    # Data preparation --------------------------------------------------------

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
    if (!is.logical(pairwise)) {
      stop("pairwise should be either TRUE or FALSE")
    }
    if (!is.logical(letters)) {
      stop("letters should be either TRUE or FALSE")
    }
    es_anova <-
      match.arg(es_anova, choices = c("none", "eta", "f"))
    if (CI < 0 ||
      CI > 1 || !is.numeric(CI)) {
      stop("CI should be a number between 0 and 1")
    }
    x <- data.frame(x)
    x$Pop <- x[, Pop]
    x$Pop <-
      gsub(
        x = x$Pop,
        pattern = "[^[:alnum:]]",
        replacement = "_"
      )
    x$Pop <- factor(x$Pop)
    levels(x$Pop) <- droplevels(x$Pop)
    if (length(unique(x$Pop)) != length(x$Pop[which(!is.na(x$Pop))])) {
      stop("Populations names'should be unique")
    }

    # ANOVA from summary data -------------------------------------------------

    aov_sum <- function(.mu, .sdev, n) {
      N <- length(.mu)
      Pop <- factor(rep(levels(x$Pop), n))
      df <- lapply(seq_len(N), function(i) {
        scale(stats::rnorm(n[i])) * .sdev[i] + .mu[i]
      })
      x <- do.call(rbind, df)
      out <- data.frame(Pop, x)
      return(out)
    }
    av_M <- aov_sum(x$M.mu, x$M.sdev, x$m)
    av_F <- aov_sum(x$F.mu, x$F.sdev, x$f)
    av_M <- stats::aov(x ~ Pop, data = av_M)
    av_F <- stats::aov(x ~ Pop, data = av_F)
    M_model <-
      anova_es(x = av_M, es_anova = es_anova, digits = digits, CI = CI)
    F_model <-
      anova_es(x = av_F, es_anova = es_anova, digits = digits, CI = CI)

    # Pairwise comparisons ----------------------------------------------------

    M_post <-
      data.frame(
        TukeyHSD(av_M, conf.level = CI)[[1]]
      )
    rownames(M_post) -> populations
    rownames(M_post) <- NULL
    M_post <- apply(M_post, 2, round, digits) %>% as.data.frame()
    M_post$populations <- populations
    M_post <- relocate(M_post, populations, .before = 1)
    colnames(M_post) <-
      c(
        "populations",
        "mean.diff",
        "conf.low",
        "conf.high",
        "p.value"
      )
    p <- M_post$p.value
    M_post$signif <-
      case_when(
        p > 0.05 ~ "ns",
        p < 0.05 & p > 0.01 ~ "*",
        p < 0.01 & p > 0.001 ~ "**",
        p < 0.001 ~ "***"
      )

    M_letters <- M_post$p.value
    names(M_letters) <- M_post$populations
    M_letters <-
      rown_col(data.frame(
        "letters" = multcompView::multcompLetters(M_letters, threshold = CI)[[1]]
      ),
      var = "populations"
      )
    F_post <-
      data.frame(
        TukeyHSD(av_F, conf.level = CI)[[1]]
      )
    rownames(F_post) -> populations
    rownames(F_post) <- NULL
    F_post <- apply(F_post, 2, round, digits) %>% as.data.frame()
    F_post$populations <- populations
    F_post <- relocate(F_post, populations, .before = 1)
    colnames(F_post) <-
      c(
        "populations",
        "mean.diff",
        "conf.low",
        "conf.high",
        "p.value"
      )
    p <- F_post$p.value
    F_post$signif <-
      case_when(
        p > 0.05 ~ "ns",
        p < 0.05 & p > 0.01 ~ "*",
        p < 0.01 & p > 0.001 ~ "**",
        p < 0.001 ~ "***"
      )

    F_letters <- F_post$p.value
    names(F_letters) <- F_post$populations
    F_letters <-
      rown_col(data.frame(
        "letters" = multcompView::multcompLetters(F_letters, threshold = CI)[[1]]
      ),
      var = "populations"
      )


    if (isTRUE(pairwise)) {
      if (isTRUE(letters)) {
        list(
          "Male model" = M_model,
          "Male posthoc" = M_post,
          "Male letters" = M_letters,
          "Female model" = F_model,
          "Female posthoc" = F_post,
          "Female letters" = F_letters
        )
      } else {
        list(
          "Male model" = M_model,
          "Male posthoc" = M_post,
          "Female model" = F_model,
          "Female posthoc" = F_post
        )
      }
    } else {
      list("Male model" = M_model, "Female model" = F_model)
    }
  }
