# scMiscellaneous
#
# Set of tools for easy score for credit scorecard build.
#
# Author : Tri Le <lmtleminh@gmail.com>
#
#' Miscellaneous tools for Credit Scoring
#'
#' These are set of tools for easy credit scorecard build. The set includes \code{\link{sc.score}},
#' \code{\link{sc.point}} and \code{\link{sc.corplot}}.
#'
#' @param data A data frame which contains target varible as well as predictor variables.
#' @param model A model built.
#' @param pdo Point double odds.
#' @param score Score at which the desired odd is.
#' @param odd The desired odd.
#'
#' @return \code{\link{sc.score}}. The output is a data frame with chosen variables in the model. All variables are in score point.
#'   There is also a final score.
#'
#' @examples
#' \dontrun{
#' # Load library
#' library(easysc)
#'
#' # Generate final score data based on the model built with data.
#' scoredata <- sc.score(data = df, model = model.Lasso, pdo = 100, score = 800, odd = 5)
#'
#' # Generate final scorecard based on the model built with data.
#' scorecard <- sc.point(data = df, woe = woe, model = model.Lasso, pdo = 100, score = 800, odd = 5)
#' }
#' @name sc.Misc
NULL
#' @rdname sc.Misc
#' @export
#score calculation
sc.score <- function(data, model, pdo, score, odd) {
  if (any(class(model) %in% c('glm', 'lm', 'glmboost', 'mboost'))) {
    var <- coef(model)
  } else if (any(class(model) == 'cv.glmnet')) {
    var <- as.vector(coef(model, s = 'lambda.1se'))
    names(var) <- rownames(coef(model, s = 'lambda.1se'))
    var <- var[var!=0]
  }
  data %<>%
    dplyr::select(names(var)[-1])
  factor_ <- pdo / log(2)
  offset_ <- score - (factor_ * log(odd))
  m <- ncol(data)
  for (i in 1:m) {
    data[[i]] <- round((var[[names(data)[i]]] * data[[i]] +
                          var[['(Intercept)']] / m) * factor_ +
                         offset_ / m, 0)
  }
  data$SCORE <- rowSums(data)
  return(data)
}

#' @rdname sc.Misc
#' @param woe A woe result which is calculated by \code{\link[klaR:woe]{woe}}.
#' @return \code{\link{sc.point}}. The output is a list with chosen variables in the model. Variables with their values and score point.
#extract score point
#' @export
sc.point <- function(data, woe, model, pdo, score, odd) {
  data_0 <- finalize(data)
  data_1 <- predict(woe, as.data.frame(data_0))
  if (any(class(model) %in% c('glm', 'lm', 'glmboost', 'mboost'))) {
    var <- coef(model)
  } else if (any(class(model) == 'cv.glmnet')) {
    var <- as.vector(coef(model, s = 'lambda.1se'))
    names(var) <- rownames(coef(model, s = 'lambda.1se'))
    var <- var[var!=0]
  }
  data_0 %<>%
    dplyr::select(names(var)[-1])
  data_1 %<>%
    dplyr::select(paste0('woe_', names(var)[-1]))
  factor_ <- pdo / log(2)
  offset_ <- score - (factor_ * log(odd))
  m <- ncol(data_1)
  for (i in 1:m) {
    data_1[[i]] <- round((var[[str_replace(names(data_1)[i], 'woe_', '')]] * data_1[[i]] +
                            var[['(Intercept)']] / m) * factor_ +
                           offset_ / m, 0)
  }
  data_2 <- cbind(data_0, data_1)
  x <- list()
  for (i in 1:m) {
    nam <- paste0('n', i)
    assign(nam, unique(data_2[,c(i, i+m)]))
    x[[names(data_0)[[i]]]] <- get(nam)
  }
  return(x)
}

#' @rdname sc.Misc
#' @param s In case of using \code{\link[glmnet:cv.glmnet]{Lasso}} model. Lambda min or
#'   lambda at 1se is chosen.
#' @return \code{\link{sc.corplot}}. The output plot is generated showing correlated variables.
#correlation plot
#' @export
sc.corplot <- function(model, data, s = NULL) {
  if (is.null(s))
    s = 'lambda.1se'
  var <- data.frame(COEF = as.vector(coef(model, s)),
                    VAR = rownames(coef(model, s)))
  var <- var[var$COEF != 0,]
  na <- var %>%
    filter(VAR != '(Intercept)') %>%
    pull(VAR) %>%
    as.vector()
  na <- na[order(na)]
  corr <- data[,na] %>% cor()

  corr[upper.tri(corr)] <- NA

  corr %>%
    broom::tidy() %>%
    tidyr::gather('.colnames', 'value', c(2:ncol(.))) %>%
    dplyr::filter(!is.na(value)) %>%
    ggplot2::ggplot(mapping = aes(x = .rownames, y = .colnames)) +
    geom_point(mapping = aes(size = value, color = value)) +
    geom_text(mapping = aes(label =round(value, 1))) +
    scale_colour_gradient(low = 'blue', high = 'red')
}

#' @export
#performance calculation
sc.performance <- function(model, X, y, s = NULL, plot = F) {
  if (any(class(model) %in% c('glm', 'lm', 'glmboost', 'mboost'))) {
    prob <- predict(model, X, type = 'response')
  } else if (any(class(model) == 'cv.glmnet')) {
    if (is.null(s))
      s = 'lambda.1se'
    Xt <- model.matrix(BAD ~., data = cbind(X %>%
                                              dplyr::select(rownames(coef(model, s = s))[rownames(coef(model, s = s)) != '(Intercept)']),
                                            BAD = y))
    prob <- predict(model, Xt, s = s, type = 'response')
  } else if (any(class(model) %in% c('train', 'train.formula'))) {
    prob <- predict(model, X, type = 'prob')[1]
  }
  pred <- ROCR::prediction(1-prob, y)
  perf <- ROCR::performance(pred, 'tpr', 'fpr')
  auc <- attr(ROCR::performance(pred, 'auc'), 'y.values')[[1]]
  ks <- max(attr(perf, 'y.values')[[1]] - attr(perf, 'x.values')[[1]])
  gini <- 2 * auc - 1
  if (plot) {
    plot(perf)
    abline(0, 1)
  }
  return(list(GINI = gini, KS = ks))
}
