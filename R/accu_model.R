#' @title Evaluation Of Sex prediction Accuracy
#' @description Testing, cross validation and visualization of the accuracy of
#' different sex prediction models using the
#' [confusionMatrix][caret::confusionMatrix] and roc curves
#' @param f Formula in the form `groups ~ x1 + x2 + ...`. The grouping factor
#' is placed to the left hand side while the numerical measurements are placed
#' to the right hand side
#' @param x Data frame to be fitted to the model
#' @param y New data frame to be tested, if `NULL` `x` is splitted to test and
#' training data seta, Default: NULL
#' @param method A string specifying which classification or regression model
#' to use,
#' @param res_method The resampling method used by
#' [trainControl][caret::trainControl], Default: 'repeatedcv'
#' @param p Percentage of `x` for testing the model in case `y` is NULL,
#' Default: 0.75
#' @param nf number of folds or of resampling iterations, Default: 10
#' @param nr Number of repeats for repeated k fold cross validation, Default:
#' 3
#' @param plot Logical; if TRUE returns an roc curve for model accuracy,
#' Default:
#'  FALSE
#' @param Pop Number of the column containing populations' names, Default:
#' NULL
#' @inheritParams extract_sum
#' @param byPop Logical; if TRUE returns the accuracy in different populations
#' of the new data frame, Default: FALSE.
#' @param ref. reference category in the grouping factor, Default: 'F'
#' @param post. positive category in the grouping factor, Default: 'M'
#' @param ... additional arguments that can passed to modeling,
#' [confusionMatrix][caret::confusionMatrix] function and roc curve generated
#' by [plot_roc][cutpointr::plot_roc]
#' @return Visual and numerical accuracy parameters for the tested model
#' @details Data frames to be entered as input need to be arranged in a
#' similar manner to [Howells] dataset.
#' @examples
#' \donttest{
#' # Using a single dataset
#' library(TestDimorph)
#' accu_model(
#'   Sex ~ GOL + NOL + BNL,
#'   x = Howells,
#'   method = "lda",
#'   plot = FALSE
#' )
#' }
#' @seealso [cutpointr::plot_roc()]
#' [caret::confusionMatrix()]
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom stats relevel predict
#' @importFrom cutpointr plot_roc cutpointr
#' @importFrom caret confusionMatrix trainControl train
accu_model <-
  function(f,
           x,
           y = NULL,
           method = "lda",
           res_method = "repeatedcv",
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
    # First data.frame preparation --------------------------------------------------------

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
    x <- data.frame(x)
    x$Sex <- x[, Sex]
    x$Sex <- as.factor(x$Sex)
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

    # Cross validation --------------------------------------------------------
    if (is.null(y)) {
      x <- x %>% mutate(id = row_number())
      z <- x
      z$Pop <- NULL
      smp_size <- floor(p * nrow(z))
      train_ind <- sample(seq_len(nrow(z)), size = smp_size)
      train.data <- z[train_ind, ]
      test.data <- z[-train_ind, ]


      train.control <- caret::trainControl(
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
        data.frame("id" = test.data$id, "class" = predict(model, test.data))
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

      # Second data.frame preparation -------------------------------------------
    } else {
      if (!(is.data.frame(y))) {
        stop("x and y should be dataframes")
      }

      if (!(Sex %in% seq_along(y))) {
        stop("Sex should be number from 1 or ncol(y)")
      }

      y <- data.frame(y)
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

      # Modeling ----------------------------------------------------------------

      model <- caret::train(f,
        data = x,
        method = method
      )
      preds <- stats::predict(model, newdata = y)

      df <- data.frame(
        "Sex" = y$Sex,
        "class" = preds,
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

      # ROC curve and confusion matrix ------------------------------------------
      df$Sex <- as.numeric(df$Sex)
      df$class <- as.numeric(df$class)

      roc <-
        cutpointr::plot_roc(
          cutpointr::cutpointr(
            data = df,
            x = class,
            class = Sex,
            subgroup = Pop,
            pos_class = 2,
            neg_class = 1,
            silent = TRUE
          )
        ) +
        theme(legend.title = ggplot2::element_blank()) + labs(title = NULL, subtitle = NULL)
      conf <-
        lapply(list,
          caret::confusionMatrix,
          positive = post.,
          reference = ref.,

          ...
        )
      if (isTRUE(plot)) {
        list(roc, conf)
      } else {
        conf
      }
    } else {
      xtab <- table(df$class, df$Sex, dnn = c("Prediction", "Reference"))
      df$Sex <- as.numeric(df$Sex)
      df$class <- as.numeric(df$class)
      roc <-
        cutpointr::plot_roc(
          cutpointr::cutpointr(
            data = df,
            x = class,
            class = Sex,
            pos_class = 2,
            neg_class = 1,
            silent = TRUE
          )
        ) +
        theme(legend.title = ggplot2::element_blank()) + labs(title = NULL, subtitle = NULL)
      conf <-
        caret::confusionMatrix(
          xtab,
          positive = post.,
          reference = ref.,
          dnn = c("Prediction", "Reference"),
          ...
        )

      if (isTRUE(plot)) {
        list(roc, conf)
      } else {
        conf
      }
    }
  }
