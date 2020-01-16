#'@title Evaluation Of Sex-prediction Accuracy
#'@description Testing and visualization of the accuracy of different sex
#'  prediction models using the [confusionMatrix][caret::confusionMatrix] and
#'  roc curves
#'@param f Formula in the form `groups ~ x1 + x2 + ...`. The grouping factor is
#'  placed to the left hand side while the numerical measurements are placed to
#'  the right hand side
#'@param x Data frame to be fitted to the model
#'@param y New data frame to be tested
#'@inheritParams extract_sum
#'@param byPop Logical; if TRUE returns the accuracy in different populations of
#'  the new data frame, Default: TRUE.
#'@param method Different methods of modeling see `details` , Default:'lda'
#'@param plot Logical; if TRUE returns an roc curve for model accuracy, Default:
#'  FALSE
#'@param cutoff cutoff value when using logistic regression, Default: 0.5
#'@param ref. reference category in the grouping factor, Default: 'F'
#'@param post. positive category in the grouping factor, Default: 'M'
#'@param ... additional arguments that can passed to modeling,
#'  [confusionMatrix][caret::confusionMatrix] function and roc curve generated
#'  by [geom_roc][plotROC::geom_roc]
#'@return Visual and numerical accuracy parameters for the tested model
#'@details Tibble/data frames to be entered as input need to be arranged in a
#'  similar manner to [Howells] dataset. Methods used for modeling are:
#'  \describe{ \item{[lda][MASS::lda]}{linear discriminant analysis}
#'  \item{[qda][MASS::qda]}{quadratic discriminant analysis}
#'  \item{[mda][mda::mda]}{mixture discriminant analysis}
#'  \item{[fda][mda::fda]}{flexible discriminant analysis}
#'  \item{[rda][klaR::rda]}{regularized discriminant analysis}
#'  \item{[glm][stats::glm]}{binomial logistic regression}
#'  \item{[raf][randomForest::randomForest]}{random forest}}
#' @examples
#' #Splitting Howells dataset into training and test datasets
#' smp_size <- floor(0.5 * nrow(Howells))
#' set.seed(123)
#' train_ind <- sample(seq_len(nrow(Howells)), size = smp_size)
#' train <- Howells[train_ind, ]
#' test <- Howells[-train_ind, ]
#' library(TestDimorph)
#' AccuModel(
#'Sex ~ GOL + NOL + BNL,
#'x = train,
#'y = test,
#'byPop = FALSE,
#'method = "lda",
#'plot = FALSE
#')
#'@seealso \code{\link[MASS]{lda}},\code{\link[MASS]{qda}}
#'\code{\link[mda]{mda}},\code{\link[mda]{fda}} \code{\link[klaR]{rda}}
#'\code{\link[randomForest]{randomForest}} \code{\link[plotROC]{GeomRoc}}
#'\code{\link[caret]{confusionMatrix}}
#'@export
#'@importFrom stats relevel predict glm binomial
#'@importFrom MASS lda qda
#'@importFrom mda mda fda
#'@importFrom klaR rda
#'@importFrom ggplot2 ggplot geom_abline xlab ylab facet_wrap aes theme
#'  element_blank
#'@importFrom plotROC geom_roc
#'@importFrom purrr map
#'@importFrom tibble is_tibble
#'@importFrom rlang abort
#'@importFrom caret confusionMatrix
#'@importFrom randomForest randomForest
AccuModel <-
  function(f,
           x,
           y,
           Sex = 1,
           Pop = 2,
           byPop = TRUE,
           method = "lda",
           plot = FALSE,
           cutoff = 0.5,
           ref. = "F",
           post. = "M",
           ...) {
    if (!(is.data.frame(x) || tibble::is_tibble(x))) {
      rlang::abort("x and y should be tibbles or dataframes")
    }
    if (!(is.data.frame(y) || tibble::is_tibble(y))) {
      rlang::abort("x and y should be tibbles or dataframes")
    }
    if (!(byPop %in% c(TRUE, FALSE)))   {
      rlang::abort("byPop should be either TRUE or FALSE")

    }
    if (!(plot %in% c(TRUE, FALSE)))   {
      rlang::abort("plot should be either TRUE or FALSE")

    }
    if (!(Sex %in% seq_along(1:ncol(x))))   {
      rlang::abort("Sex should be number from 1 to ncol(x)")

    }
    if (!(Pop %in% seq_along(1:ncol(x))))   {
      rlang::abort("Pop should be number from 1 to ncol(x)")

    }
    if (!(Sex %in% seq_along(1:ncol(y))))   {
      rlang::abort("Sex should be number from 1 or ncol(y)")

    }
    if (!(Pop %in% seq_along(1:ncol(y))))   {
      rlang::abort("Pop should be number from 1 or ncol(y)")

    }
    x <- data.frame(x)
    y <- data.frame(y)
    x$Pop <- x[, Pop]
    x$Sex <- x[, Sex]
    y$Pop <- y[, Pop]
    y$Sex <- y[, Sex]
    x$Pop <- as.factor(x$Pop)
    y$Pop <- as.factor(y$Pop)
    x$Sex <- as.factor(x$Sex)
    y$Sex <- as.factor(y$Sex)
    if (!(ref. %in% c("M", "F")))   {
      rlang::abort("ref. should be one of `M` or `F`")

    }
    if (!(post. %in% c("M", "F")))   {
      rlang::abort("post. should be one of `M` or `F`")

    }
    x$Sex <- stats::relevel(x$Sex, ref = ref.)
    y$Sex <- stats::relevel(y$Sex, ref = ref.)
    if (length(unique(x$Sex)) != 2 &&
        (!(levels(x$Sex) %in% c("M", "F")))) {
      rlang::abort("Sex column should be a factor with only 2 levels `M` and `F`")
    }
    if (length(unique(y$Sex)) != 2 &&
        (!(levels(y$Sex) %in% c("M", "F")))) {
      rlang::abort("Sex column should be a factor with only 2 levels `M` and `F`")
    }
    if (!(method %in% c("lda", "qda", "mda", "fda", "rda", "glm", "raf")))   {
      rlang::abort("method should be one of `lda`, `qda`,`mda`,`fda`,`rda`,`glm`,`raf`")

    }
    if (method == "raf") {
      model <- randomForest::randomForest(f, data = x, ...)
      preds <- stats::predict(model, newdata = y)
      preds <- as.data.frame(preds)
      colnames(preds)[1] <- "class"
    }
    if (method == "lda") {
      model <- MASS::lda(f, data = x, ...)
      preds <- stats::predict(model, newdata = y)
    }
    if (method == "qda") {
      model <- MASS::qda(f, data = x, ...)
      preds <- stats::predict(model, newdata = y)
    }
    if (method == "mda") {
      model <- mda::mda(f, data = x, ...)
      preds <- stats::predict(model, newdata = y)
      preds <- as.data.frame(preds)
      colnames(preds)[1] <- "class"
    }
    if (method == "fda") {
      model <- mda::fda(f, data = x, ...)
      preds <- stats::predict(model, newdata = y)
      preds <- as.data.frame(preds)
      colnames(preds)[1] <- "class"
    }
    if (method == "rda") {
      model <- klaR::rda(f, data = x, ...)
      preds <- stats::predict(model, newdata = y)
    }
    if (method == "glm") {
      model <-
        stats::glm(
          f,
          family = stats::binomial(link = 'logit'),
          data = x,
          maxit = 100,
          ...
        )
      preds <-
        stats::predict(model, newdata = y, type = 'response')
      class <- ifelse(test = preds > cutoff,
                      yes = "M",
                      no = "F")
      preds <- cbind.data.frame(preds, class)
    }
    df <- data.frame("Sex" = y$Sex,
                     "class" = preds$class,
                     "Pop" = y$Pop)
    if (byPop == TRUE) {
      list <- by(
        df,
        df$Pop,
        FUN = function(x) {
          table(x$class, x$Sex)
        }
      )
      roc <-
        ggplot2::ggplot(df, aes(
          m = as.numeric(df$Sex),
          d = as.numeric(df$class),
          color = df$Pop,
          ...
        )) + plotROC::geom_roc(n.cuts = 0, ...) + ylab("Sensitivity") + xlab("1-Specificity") +
        geom_abline(...) + facet_wrap( ~ Pop) + ggplot2::theme(legend.title = ggplot2::element_blank())
      conf <-
        purrr::map(list,
                   caret::confusionMatrix,
                   positive = post.,
                   reference = ref.,
                   ...)
      if (plot == TRUE) {
        list(roc, conf)
      } else{
        return(conf)
      }
    } else{
      xtab <- table(df$class, df$Sex)
      roc <-
        ggplot2::ggplot(df, aes(m = as.numeric(df$Sex), d = as.numeric(df$class), ...)) + plotROC::geom_roc(n.cuts = 0, ...) +
        ylab("Sensitivity") + xlab("1-Specificity") + geom_abline(...)
      conf <-
        caret::confusionMatrix(xtab, positive = post., reference = ref., ...)

      if (plot == TRUE) {
        list(roc, conf)
      } else{
        return(conf)
      }


    }

  }
