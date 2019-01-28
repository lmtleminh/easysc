# scFeatureSelection
#
# Features Selection for easy credit scorecard build.
#
# Author : Tri Le <lmtleminh@gmail.com>
#
#' Feature Selection for Credit Scoring
#'
#' Using multiple techniques for features selection.
#'
#' @param data A data frame which contains target varible as well as predictor variables.
#' @param target Target variable name.
#' @param median A logical scalar.
#' @param p_cor A logical scalar.
#' @param s_cor A logical scalar.
#' @param logreg A logical scalar.
#' @param gini_rf A logical scalar.
#' @param parallel A logical scalar.
#'
#' @return The output is a data frame of variables with important measurement. All variables are arranged in point order.
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
#' @name sc.featureSelection
NULL
#' @rdname sc.featureSelection
#' @export
sc.featureSelection <- function(data, target, median = T, p_cor = T, s_cor = T, logreg = T, gini_rf = F, parallel = F) {
  target <- deparse(substitute(target))
  y <- data[[target]]
  X <- data %>% dplyr::select(-target)
  classes <- unique(y)
  if (median) {
    print('Start median')
    w_test <- base::sapply(colnames(X),
                           function(x) wilcox.test(data[data[,target] == classes[1],x],
                                                   data[data[,target] == classes[2],x])$p.value)
    imp1 <- cbind(1:length(w_test), 1 - w_test)
    imp1[,2] <- imp1[,2]/max(imp1[,2])
    print('Done!')
  }
  else {
    imp1 <- cbind(1:ncol(X), rep(0, ncol(X)))
  }
  cor_func <- function(method) {
    #http://www.jmlr.org/papers/volume5/yu04a/yu04a.pdf
    cor <- cor(data, method = method)
    cor_c <- as.matrix(cbind(1:(nrow(cor) - 1), abs(cor[-nrow(cor),ncol(cor)])))
    cor_c <- cor_c[order(cor_c[,2], decreasing = T),]
    for (i in 1:(nrow(cor_c)-1)) {
      for (j in 1:(nrow(cor_c)-i)) {
        if (cor_c[row.names(cor_c)[i],2] != 0 &
            cor[row.names(cor_c)[i], row.names(cor_c)[i+j]] >= cor_c[row.names(cor_c)[i],2]) {
          cor_c[row.names(cor_c)[i+j],2] <- 0
        }
      }
    }
    return(cor_c)
  }
  if (p_cor) {
    print('Start pearson')
    imp2 <- cor_func('pearson')
    imp2[,2] <- imp2[,2]/max(imp2[,2])
    imp2 <- imp2[order(imp2[,1], decreasing = F),]
    print('Done!')
  }
  else {
    imp2 <- cbind(1:ncol(X), rep(0, ncol(X)))
  }
  if (s_cor) {
    print('Start spearman')
    imp3 <- cor_func('spearman')
    imp3[,2] <- imp3[,2]/max(imp3[,2])
    imp3 <- imp3[order(imp3[,1], decreasing = F),]
    print('Done!')
  }
  else {
    imp3 <- cbind(1:ncol(X), rep(0, ncol(X)))
  }
  if (logreg) {
    print("Start LogReg")
    X_std <- as.data.frame(sapply(colnames(X), function(x) (X[,x] - mean(X[,x]))/var(X[,x])))
    #X_std <- cbind(X, y)
    model_lg <- glm(as.factor(y) ~ ., data = X_std, family = binomial(link = 'logit'))
    corr <- abs(model_lg$coefficients)[-1]
    imp4 <- cbind(1:length(corr), corr)
    imp4[is.na(imp4[,2]),2] <- 0
    imp4[,2] <- imp4[,2]/max(imp4[,2])
    print('Done!')
  }
  else {
    imp4 <- cbind(1:ncol(X), rep(0, ncol(X)))
  }
  if (gini_rf) {
    print("Start RF")
    start.run <- Sys.time()
    if (parallel) {
      mc <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(mc)
      doParallel::registerDoParallel(cl)
      rf <- foreach(ntree = rep(1000, 3), .combine = randomForest::combine,
                    .multicombine = TRUE, .packages = 'randomForest') %dopar%
        randomForest::randomForest(x = X, y = as.factor(y), importance = T, replace = F,
                                   ntree = ntree)


      parallel::stopCluster(cl)
    } else {
      rf <- foreach(ntree = rep(1000, 3), .combine = randomForest::combine,
                    .multicombine = TRUE, .packages = 'randomForest') %do%
        randomForest::randomForest(x = X, y = as.factor(y), importance = T, replace = F,
                                   ntree = ntree)
    }
    end.run <- Sys.time()
    diff <- end.run - start.run
    print(diff)
    imp5 <- cbind(1:ncol(X), rf$importance[,3])
    imp5[,2] <- imp5[,2]/max(imp5[,2])
    imp6 <- cbind(1:ncol(X), rf$importance[,4])
    imp6[,2] <- imp6[,2]/max(imp6[,2])
    print('Done!')
  }
  else {
    imp5 <- cbind(1:ncol(X), rep(0, ncol(X)))
    imp6 <- cbind(1:ncol(X), rep(0, ncol(X)))
  }
  print('Returing result...')
  var_imp <- cbind(imp1[,2], imp2[,2], imp3[,2], imp4[,2], imp5[,2], imp6[,2])
  colnames(var_imp) <- c('Median', 'Pearson', 'Spearman', 'LogReg', 'RF_ACC', 'RF_GINI')
  structure(var_imp, class = 'varimp')
}
#' @method plot varimp
#' @export
plot.varimp <- function(fs) {
  fs %>%
    as.data.frame() %>%
    tibble::rownames_to_column('Variables') %>%
    tidyr::gather(Method, VarImp, -Variables) %>%
    ggplot() +
    geom_bar(mapping = aes(x = reorder(Variables, VarImp), y = VarImp, fill = reorder(Method,VarImp)),
             stat = 'identity') +
    labs(title = 'Feature Selection', x = 'Variables') +
    scale_fill_discrete(name = 'Method') +
    coord_flip() +
    theme_classic()
}
