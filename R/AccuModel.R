#' @title AccuModel
#' @seealso
#'  [TestDimorph-deprecated()]
#' @name AccuModel-deprecated
#' @keywords internal
NULL
#' @rdname TestDimorph-deprecated
#' @section `AccuModel`:
#' For `AccuModel`, use [accu_model()].
#' @export
AccuModel <-
  function(f,
           x,
           y = NULL,
           method = "lda",
           res_method = "repeatedcv",
           prior = NULL,
           p = 0.75,
           nf = 10,
           nr = 3,
           plot = FALSE,
           Sex = 1,
           Pop = NULL,
           byPop = FALSE,
           ref. = "F",
           post. = "M",
           ...) {
    if (!identical(Sys.getenv("TESTTHAT"), "true")) {
      .Deprecated("accu_model")
    }
    prob <- NULL
    if (!(is.data.frame(x))) {
      stop("x and y should be dataframes")
    }
    if (!is.logical(byPop)) {
      stop("byPop should be either TRUE or FALSE")
    }
    if (isTRUE(byPop) && is.null(Pop)) {
      warning("Pop column number should be specified if byPop is TRUE")
    }
    if (!is.logical(plot)) {
      stop("plot should be either TRUE or FALSE")
    }
    if (!(Sex %in% seq_along(x))) {
      stop("Sex should be a number from 1 to ncol(x)")
    }
    x <- x %>%
      drop_na() %>%
      as.data.frame()
    x$Sex <- x[, Sex]
    x$Sex <- factor(x$Sex, levels = c(ref., post.))
    if (is.null(Pop)) {
      x$Pop <- as.factor(rep("pop_1", nrow(x)))
    } else {
      if (!(Pop %in% seq_along(x))) {
        stop("Pop should be a number from 1 or ncol(x)")
      }
      x$Pop <- x[, Pop]
      x$Pop <- factor(x$Pop)
      x <- dplyr::arrange(x, x$Pop, x$Sex)
    }
    if (is.null(y)) {
      x <- x %>% mutate(id = row_number())
      z <- x
      z$Pop <- NULL
      smp_size <- floor(p * nrow(z))
      train_ind <- sample(seq_len(nrow(z)), size = smp_size)
      train.data <- z[train_ind, ]
      test.data <- z[-train_ind, ]


      train.control <- caret::trainControl(
        classProbs = TRUE,
        method = res_method,
        number = nf,
        repeats = nr
      )

      model <- caret::train(f,
        data = train.data,
        method = method,
        trControl = train.control
      )
      preds <-
        data.frame(
          "id" = test.data$id, "class" = predict(model, test.data),
          "prob" = predict(model, test.data, type = "prob")[, 2]
        )
      df <-
        dplyr::full_join(data.frame("id" = test.data$id, "Sex" = test.data$Sex),
          preds,
          by = "id"
        )
      df <-
        dplyr::right_join(data.frame("id" = x$id, "Pop" = x$Pop), df,
          by =
            "id"
        )
    } else {
      if (!(is.data.frame(y))) {
        stop("x and y should be dataframes")
      }

      if (!(Sex %in% seq_along(y))) {
        stop("Sex should be number from 1 or ncol(y)")
      }

      y <- y %>%
        drop_na() %>%
        as.data.frame()
      y$Sex <- y[, Sex]
      y$Sex <- factor(y$Sex)
      if (is.null(Pop)) {
        y$Pop <- as.factor(rep("pop_1", nrow(y)))
      } else {
        if (!(Pop %in% seq_along(y))) {
          stop("Pop should be number from 1 or ncol(y)")
        }
        y$Pop <- y[, Pop]
        y$Pop <- factor(y$Pop)
        y <- dplyr::arrange(y, y$Pop, y$Sex)
      }
      if (!(ref. %in% c("M", "F"))) {
        stop("ref. should be one of `M` or `F`")
      }
      if (!(post. %in% c("M", "F"))) {
        stop("post. should be one of `M` or `F`")
      }
      x$Sex <- stats::relevel(x$Sex, ref = ref.)
      y$Sex <- stats::relevel(y$Sex, ref = ref.)
      if (length(unique(x$Sex)) != 2 &&
        (!(levels(x$Sex) %in% c("M", "F")))) {
        stop("Sex column should be a factor with only 2 levels `M` and `F`")
      }
      if (length(unique(y$Sex)) != 2 &&
        (!(levels(y$Sex) %in% c("M", "F")))) {
        stop("Sex column should be a factor with only 2 levels `M` and `F`")
      }
      train.control <- caret::trainControl(
        classProbs = TRUE,
        method = res_method,
        number = nf,
        repeats = nr
      )
      model <- caret::train(f,
        data = x,
        method = method,
        trControl = train.control
      )
      preds <- cbind.data.frame(
        class = predict(model, newdata = y),
        prob = predict(model, newdata = y, type = "prob")[, 2]
      )

      df <- data.frame(
        "Sex" = y$Sex,
        "class" = preds$class,
        "prob" = preds$prob,
        "Pop" = y$Pop,
        stringsAsFactors = TRUE
      )
    }

    if (isTRUE(byPop)) {
      list <- by(
        df,
        df$Pop,
        FUN = function(x) {
          table(x$class, x$Sex, dnn = c("Prediction", "Reference"))
        }
      )

      df$Sex <- as.numeric(df$Sex)
      df$class <- as.numeric(df$class)


      cutpoint <- cutpointr::cutpointr(
        data = df,
        x = prob,
        class = Sex,
        subgroup = Pop,
        pos_class = 2,
        neg_class = 1,
        silent = TRUE, method = cutpointr::maximize_metric,
        metric = cutpointr::sum_sens_spec
      )


      roc <-
        cutpointr::plot_roc(cutpoint) +
        theme(legend.title = ggplot2::element_blank()) + labs(title = NULL, subtitle = NULL)
      conf <-
        lapply(list,
          caret::confusionMatrix,
          positive = post.,
          reference = ref.,

          ...
        )
      cutpoint <- cbind.data.frame(cutpoint[, 1], cutpoint[, 3])
      names(cutpoint) <- c("populations", "cutpoint")
      if (isTRUE(plot)) {
        list(cutpoint = cutpoint, conf, roc)
      } else {
        list(cutpoint = cutpoint, conf)
      }
    } else {
      xtab <- table(df$class, df$Sex, dnn = c("Prediction", "Reference"))
      df$Sex <- factor(df$Sex, levels = c(ref., post.))
      df$Sex <- as.numeric(df$Sex)
      df$class <- as.numeric(df$class)

      cutpoint <- cutpointr::cutpointr(
        data = df,
        x = prob,
        class = Sex,
        pos_class = 2,
        neg_class = 1,
        silent = TRUE, method = cutpointr::maximize_metric,
        metric = cutpointr::sum_sens_spec
      )
      roc <-
        cutpointr::plot_roc(cutpoint) +
        theme(legend.title = ggplot2::element_blank()) + labs(title = NULL, subtitle = NULL)
      conf <-
        caret::confusionMatrix(
          xtab,
          positive = post.,
          reference = ref.,
          dnn = c("Prediction", "Reference"),
          ...
        )
      cutpoint <- pull(cutpoint, 2)
      if (isTRUE(plot)) {
        list(cutpoint = round(cutpoint, 4), conf, roc)
      } else {
        list(cutpoint = round(cutpoint, 4), conf)
      }
    }
  }
