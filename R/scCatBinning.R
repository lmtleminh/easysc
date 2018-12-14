# scCatBinning
#
# This is an easy binning solution for credit scorecard build.
#
# Author : Tri Le <lmtleminh@gmail.com>
#
#' Categorical binning for Credit Scoring
#'
#' This is an easy categorical binning solution for credit scorecard build. It is designed to
#' group the optimal categories by utilizing the \href{https://cran.r-project.org/package=partykit}{Recursive Partitioning}
#' which is only applied on factor variables.
#'
#' @param data A data frame which contains target varible as well as predictor variables.
#' @param target Target variable name.
#' @param n Number of bootstrap iterations. Default 10 times.
#' @param uni The minimum number of unique values within predictor variables. Default 4.
#' @param parallel A logical scalar. Use parallel backend. Default FALSE.
#'
#' @return The output is a list of group plan which can be applied to the orginal data frame via
#'   the \code{\link{predict}} function.
#'   The user can also update the cut plan via the \code{\link{update}} function.
#'
#' @examples
#' \dontrun{
#' # Load library
#' library(easysc)
#'
#' # Generate a grouping plan which maximize IV via 500 bootstrap resampling
#' group.plan <- sc.cat.binning(data = df, target = BAD, n = 500, parallel = TRUE)
#' # Update the grouping plan
#' update(group.plan, MARRIAGE = list(c('SINGLE', 'DIVORCE'),
#'                                    c('MARRIED', 'WIDOW'))
#' # Apply to the data frame
#' predict(group.plan, df, keepTarget = TRUE)
#' }

#for categorical binning
#' @export
sc.cat.binning <- function(data, target, n = 10, uni = 4, parallel = FALSE) {
  start_time <- Sys.time()
  auto_seg <- function(X, y, n, uni, na, parallel = FALSE) {
    print(paste0(na, '...'))
    if (class(X) != 'factor' | length(unique(X)) <= uni) return(NULL)
    start_time <- Sys.time()
    dataX <- data.frame(X, y)
    seeds <- c(0, round(runif(n) * as.numeric(paste('1e', ceiling(log10(n)) + 2, sep = '')), 0))
    seg <- function(seed, data) {
      b <- sample(nrow(data), nrow(data), replace = T)
      tr <- party::ctree(factor(y) ~ X, data=data[b,])
      m <- tr@tree[['psplit']]$splitpoint
      m1 <- tr@tree[['left']]$psplit$splitpoint
      m2 <- tr@tree[['right']]$psplit$splitpoint
      check.null <- function(x) {
        if (is.null(x)) rep(1, length(levels(data$X)))
        else x
      }
      if (!is.null(m))
        if(any(levels(data$X) != attr(m, 'levels')))
          stop('Levels are not matched!')
      M <- cbind(PRO = levels(data$X),
                 M1 = check.null(m),
                 M2 = check.null(m1),
                 M3 = check.null(m2))
      p <- c()
      for (i in 1:(nrow(M) - 1)) {
        for (j in (i+1):nrow(M)) {
          if (M[[i,'M1']] == M[[j,'M1']] &
              M[[i,'M2']] == M[[j,'M2']] &
              M[[i,'M3']] == M[[j,'M3']]) {
            o <- cbind(PROS = paste0(M[[i,1]], ':', M[[j,1]]), VAL = 1)
          } else {
            o <- cbind(PROS = paste0(M[[i,1]], ':', M[[j,1]]), VAL = 0)
          }
          p <- rbind(p, o)
        }
      }
      q <- as.matrix(as.numeric(p[,2]))
      rownames(q) <- p[,1]
      return(q)
    }
    if (parallel) {
      mc <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(mc)
      doParallel::registerDoParallel(cl)
      foreach::foreach(s = seeds, .combine = '+') %dopar% {
        seg(s, dataX)
      } -> q
      parallel::stopCluster(cl)
    } else {
      foreach::foreach(s = seeds, .combine = '+') %do% {
        seg(s, dataX)
      } -> q
    }
    g1 <- c()
    g2 <- c()
    g3 <- c()
    g4 <- c()
    g5 <- c()
    if (is.null(nrow(q))) return(NULL)
    q %<>% as.data.frame() %>% tibble::rownames_to_column('PROS')
    q %<>%
      dplyr::arrange(desc(V1)) %>%
      dplyr::mutate(LEFT_PRO = sapply(stringr::str_split(PROS, pattern = ':'), function(x) x[1]),
                    RIGHT_PRO = sapply(stringr::str_split(PROS, pattern = ':'), function(x) x[2])) %>%
      dplyr::select(-PROS)
    for (i in 1:nrow(q)) {
      if (q[i, 'LEFT_PRO'] %in% c(g1, g2, g3, g4, g5) |
          q[i, 'RIGHT_PRO'] %in% c(g1, g2, g3, g4, g5)) {
        for (j in 1:5) {
          g <- paste('g', j, sep = '')
          if (q[i, 'LEFT_PRO'] %in% get(g) & !(q[i, 'RIGHT_PRO'] %in% c(g1, g2, g3, g4, g5))) {
            assign(g, c(get(g), q[i, 'RIGHT_PRO']))
          } else if (!(q[i, 'LEFT_PRO'] %in% c(g1, g2, g3, g4, g5)) & (q[i, 'RIGHT_PRO'] %in% get(g))) {
            assign(g, c(get(g), q[i, 'LEFT_PRO']))
          }
        }
      } else if (is.null(c(g1, g2, g3, g4, g5))) {
        g <- 'g1'
        assign(g, c(get(g), q[i, 'RIGHT_PRO']))
        assign(g, c(get(g), q[i, 'LEFT_PRO']))
      } else if(is.null(c(g2, g3, g4, g5))) {
        g <- 'g2'
        assign(g, c(get(g), q[i, 'RIGHT_PRO']))
        assign(g, c(get(g), q[i, 'LEFT_PRO']))
      } else if(is.null(c(g3, g4, g5))) {
        g <- 'g3'
        assign(g, c(get(g), q[i, 'RIGHT_PRO']))
        assign(g, c(get(g), q[i, 'LEFT_PRO']))
      } else if(is.null(c(g4, g5))) {
        g <- 'g4'
        assign(g, c(get(g), q[i, 'RIGHT_PRO']))
        assign(g, c(get(g), q[i, 'LEFT_PRO']))
      } else if(is.null(c(g5))) {
        g <- 'g5'
        assign(g, c(get(g), q[i, 'RIGHT_PRO']))
        assign(g, c(get(g), q[i, 'LEFT_PRO']))
      }
    }
    end_time <- Sys.time()
    dif <- end_time - start_time
    print(dif)
    g_final <- plyr::compact(list(g1, g2, g3, g4, g5))
    return(g_final)
  }
  target <- deparse(substitute(target))
  if (!(target %in% names(data)))
    stop(paste0(target, ' is not exist!'))
  y <- data[[target]]
  data[,target] <- NULL
  segmt <- lapply(names(data), function(x)
    auto_seg(data[[x]], y, n, uni, x, parallel)
  )
  names(segmt) <- names(data)
  end_time <- Sys.time()
  dif <- end_time - start_time
  print(dif)
  attr(segmt, 'target') <- target
  structure(segmt, class = 'group.plan')
}

#' @export update
update <- function(...){
  UseMethod('update')
}
#' @method update group.plan
update.group.plan <- function(segmt, ...){
  un_arg <- list(...)
  if (class(segmt) != 'group.plan') {
    stop('Not a group plan!')
  } else if (length(un_arg) != 0) {
    for (i in 1:length(un_arg)) {
      if (!is.na(names(segmt[names(un_arg)[i]])))
        segmt[names(un_arg)[i]] <- un_arg[i]
      else {
        segmt[names(un_arg)[i]] <- un_arg[i]
        print(paste(names(un_arg)[i], 'is just added!'))
      }
    }
  } else {
    print('Nothing happens!')
  }
  return(segmt)
}

#' @export predict
predict <- function(...){
  UseMethod('predict')
}
#' @method predict group.plan
predict.group.plan <- function(segmt, data, keepTarget = FALSE) {
  if (class(segmt) != 'group.plan') {
    stop('Not a group plan!')
  } else {
    NewData <- data.frame(A = rep('', nrow(data)))
    for (i in 1:length(segmt)) {
      if (!is.null(segmt[[i]])) {
        if (class(data[[names(segmt)[i]]]) %in% c('factor')) {
          col <- data[[names(segmt)[i]]]
          for (j in 1:length(segmt[[i]])) {
            levels(col) <- c(levels(col), paste0('G', j))
            col[col %in% segmt[[i]][[j]]] <- paste0('G', j)
            col <- droplevels(col)
          }
        } else {
          stop(paste0(names(segmt)[i], ' is not factor'))
        }
      } else {
        col <- data[,names(segmt)[i]]
      }
      NewData <- cbind(NewData, col)
    }
    NewData['A'] <- NULL
    names(NewData) <- names(segmt)
    if (keepTarget & attr(segmt, 'target') %in% names(data))
      NewData[attr(segmt, 'target')] <- data[attr(segmt, 'target')]
    if (!(attr(segmt, 'target') %in% names(data)))
      warning(paste0(attr(segmt, 'target'), ' is not exist!'))
    return(NewData)
  }
}
