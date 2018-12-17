# scBinning
#
# Author : Tri Le <lmtleminh@gmail.com>
#
#' @importFrom magrittr "%>%"
#' @importFrom foreach "%do%" "%dopar%"
NULL

#' Numeric binning for Credit Scoring
#'
#' This is an easy numerical binning solution for credit scorecard build. It is designed to
#' choose the optimal binning solution by utilizing the \href{https://cran.r-project.org/package=partykit}{Recursive Partitioning}.
#' This will only bin numeric or integer variables and ignore factor or character variables.
#'
#' @param data A data frame which contains target varible as well as predictor variables.
#' @param target Target variable name.
#' @param n Number of bootstrap iterations. Default 10 times.
#' @param p The minimum percentage of observation per bin. Default 3\%.
#' @param thres Threshold differences of target between bins. Default 0.5\%.
#' @param uni The minimum number of unique values within predictor variables. Default 4.
#' @param best A logical scalar. Use different methods which maximize IV. Default TRUE.
#' @param parallel A logical scalar. Use parallel backend. Default FALSE.
#'
#' @return The output is a list of cut plan which can be applied to the orginal data frame via
#'   the \code{\link{predict}} function.
#'   The user can also update the cut plan via the \code{\link{update}} function.
#'
#' @examples
#' \dontrun{
#' # Load library
#' library(easysc)
#'
#' # Generate a cut plan which maximize IV via 500 bootstrap resampling
#' cut.plan <- sc.binning(data = df, target = BAD, n = 500, p = 5, best = TRUE, parallel = TRUE)
#' # Update the cut plan
#' update(cut.plan, AGE = c(20, 30, 40))
#' # Apply to the data frame
#' predict(cut.plan, df, keepTarget = TRUE)
#' }
#for numeric binning
#' @export
sc.binning <- function(data, target, n = 10, p = 3, thres = .5, uni = 4, best = TRUE, parallel = FALSE) {
  start_time <- Sys.time()
  target <- deparse(substitute(target))
  if (!(target %in% names(data)))
    stop(paste0(target, ' is not exist!'))
  y <- data[[target]]
  data[,target] <- NULL
  woeZ <- function(X, y) {
    df <- tibble::tibble(X, yZ = y) %>%
      dplyr::group_by(X) %>%
      dplyr::summarise(
        pct_bad = sum(yZ)/sum(y),
        pct_good = sum(1-yZ)/sum(1-y),
        woe = log(pct_good/pct_bad),
        iv = (pct_good - pct_bad)*woe
      )
    return(sum(df$iv))
  }
  superbin <- function(X, y, n = 10, p = 3, thres = .5, parallel = FALSE) {
    data <- data.frame(X, y)
    tree_bin <- function(seed, data) {
      if (seed == 0) b <- 1:nrow(data)
      else {
        set.seed(seed)
        b <- sample(nrow(data), nrow(data), replace = T)
      }
      t <- party::ctree(factor(y) ~ ., data = data[b,],
                        controls = party::ctree_control(maxdepth = 2))
      c(t@tree$psplit$splitpoint,
        t@tree$left$psplit$splitpoint,
        t@tree$right$psplit$splitpoint)
    }
    seeds <- c(0, round(runif(n) * as.numeric(paste('1e', ceiling(log10(n)) + 2, sep = '')), 0))
    x <- c()
    if (parallel) {
      mc <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(mc)
      doParallel::registerDoParallel(cl)
      foreach::foreach(s = seeds, .combine = 'c') %dopar% {
        tree_bin(s, data)
      } -> x
      parallel::stopCluster(cl)
    } else {
      foreach::foreach(s = seeds, .combine = 'c') %do% {
        tree_bin(s, data)
      } -> x
    }
    m <- table(x) %>% as.matrix()
    if (!is.null(m) & length(m) > 1) {
      m <- data.frame(x = as.numeric(rownames(m)), Freq = m)
      m$Check <- 1
      for (i in 1:(nrow(m)-1)) {
        if (m[i + 1,'x'] == m[i,'x'] + 1 & m[i + 1, 'Freq'] >= m[i, 'Freq'])
          m[i,'Check'] <- 0
        if (i >= 2) {
          if (m[i - 1,'x'] == m[i,'x'] - 1 & m[i - 1, 'Freq'] >= m[i, 'Freq'])
            m[i,'Check'] <- 0
        }
      }
      m <- m[m$Check == 1,]
      m$w <- (m$Freq - min(m$Freq)) / (max(m$Freq) - min(m$Freq))
      m <- m[m$w >=.2, 'x']
      woeC <- woeZ(cut(X, breaks = c(-Inf, m, Inf)), y)
      len <- length(m)
      findMonotonic <- function(m) {
        l1 <- length(m)
        tbl <- data %>%
          dplyr::mutate(X = cut(X, breaks = c(-Inf, m, Inf))) %>%
          dplyr::group_by(X) %>%
          dplyr::summarise(
            pct_bad = mean(y) * 100,
            pct_n = n()/nrow(data) * 100
          )
        tbl$dif <- c(diff(tbl$pct_bad), 0)
        tbl$sign <- sign(tbl$dif)
        tbl$m <- c(m, m[length(m)])
        tbl$rk <- sign(sum(tbl$sign)) * 1:nrow(tbl)
        reg <- isoreg(tbl$rk, tbl$pct_bad)
        cut <- knots(as.stepfun(reg))
        m <- unique(tbl$m[tbl$rk %in% cut])
        l2 <- length(m)
        if (l1 == l2) {
          m <- unique(tbl$m[abs(tbl$dif) >= thres])
          l2 <- length(m)
        }
        if (l1 == l2) {
          m <- unique(tbl$m[tbl$pct_n >= p])
        }
        return(m)
      }
      while(length(m) > length(findMonotonic(m))) {
        m <- findMonotonic(m)
        woeC <- c(woeC, woeZ(cut(X, breaks = c(-Inf, m, Inf)), y))
        len <- c(len, length(m))
      }
    } else {
      m <- rownames(m)
    }
    return(m)
  }
  #based on this post https://statcompute.wordpress.com/2018/11/25/improving-binning-by-bootstrap-bumping/
  bump_bin <- function(X, y, n, p, parallel = FALSE) {
    n1 <- round(p * length(y), 0)
    n2 <- 10
    #set.seed(2019)
    seeds <- c(0, round(runif(n) * as.numeric(paste('1e', ceiling(log10(n)) + 2, sep = '')), 0))
    df1 <- data.frame(X, y)
    df2 <- df1[!is.na(df1[, 'X']), c('X', 'y')]
    cor <- cor(df2[, 2], df2[, 1], method = "spearman", use = "complete.obs")
    ### MONOTONIC BINNING WITH BOOTSTRAP SAMPLES ###
    bin <- function(seed, df2) {
      if (seed == 0) b <- 1:nrow(df2)
      else {
        set.seed(seed)
        b <- sample(nrow(df2), nrow(df2), replace = T)
      }
      smp <- df2
      reg <- isoreg(smp[, 1], cor / abs(cor) * smp[, 2])
      cut <- knots(as.stepfun(reg))
      df2$cut <- cut(df2[['X']], breaks = unique(cut), include.lowest = T)
      df3 <- Reduce(rbind,
                    lapply(split(df2, df2$cut),
                           function(x) data.frame(n = nrow(x), b = sum(x[['y']]), g = sum(1 - x[['y']]),
                                                  maxx = max(x[['X']]), minx = min(x[['X']]))))
      df4 <- df3[which(df3[["n"]] > n1 & df3[["b"]] > n2 & df3[["g"]] > n2), ]
      df2$good <- 1 -  df2[['y']]
      out <- smbinning::smbinning.custom(df2, "good", 'X', cuts = df4[['maxx']][-nrow(df4)])
      if (out == "No Bins") return(NULL)
      out <- out$ivtable
      return(data.frame(iv = out[['IV']][length(out[['IV']])], nbin = nrow(out) - 2,
                        cuts = I(list(df4[['maxx']][-nrow(df4)])),
                        abs_cor = abs(cor(as.numeric(row.names(out)[1:(nrow(out) - 2)]),
                                          out[['WoE']][1:(nrow(out) - 2)], method = "spearman"))))
    }
    if (parallel) {
      mc <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(mc)
      doParallel::registerDoParallel(cl)
      foreach::foreach(s = seeds, .combine = 'rbind') %dopar% {
        bin(s, df2)
      } -> bump_out
      parallel::stopCluster(cl)
    } else {
      foreach::foreach(s = seeds, .combine = 'rbind') %do% {
        bin(s, df2)
      } -> bump_out
    }
    #bump_out <- Reduce(rbind, lapply(seeds, bin))
    ### FIND THE CUT MAXIMIZING THE INFORMATION VALUE ###
    if (is.null(bump_out)) return(NULL)
    cut2 <- bump_out[order(-bump_out["abs_cor"], -bump_out["iv"], bump_out["nbin"]), ]$cuts[[1]]
    return(cut2)
  }
  bestbin <- function(X, y, n = 10, p = 3, thres = .5, uni = 4, name = NULL, best = T, parallel = FALSE) {
    print(paste0(name, '...'))
    u <- unique(X)
    if (!(class(X) %in% c('numeric', 'integer')) | length(u) <= uni) {
      return(NULL)
    }
    spb <- superbin(X, y, n, p, thres, parallel)
    if (!is.null(spb)) attr(spb, 'method') <- 'spb'
    finalBin <- spb
    if (best) {
      smb <- smbinning::smbinning(data.frame(X, y), 'y', 'X', p = p / 100)
      #fixing smb give no split
      if (smb == "No significant splits") {
        smb.cuts <- 0
        smb.iv <- -1
      } else {
        smb.cuts <- smb$cuts
        smb.iv <- smb$iv
        attr(smb.cuts, 'method') <- 'smb'
      }
      bpb <- bump_bin(X, y, n, p = p / 100, parallel)
      if (is.null(bpb)) {
        bpb.iv <- -1
      } else {
        bpb.iv <- woeZ(cut(X, breaks = c(-Inf, bpb, Inf)), y)
        attr(bpb, 'method') <- 'bpb'
      }
      iv <- c(spb = woeZ(cut(X, breaks = c(-Inf, spb, Inf)), y),
              smb = smb.iv,
              bpb = bpb.iv)
      finalBin <- switch(names(which.max(iv)),
                         spb = spb,
                         smb = smb.cuts,
                         bpb = bpb)
    }
    return(finalBin)
  }
  cut_plan <- lapply(names(data), function(x)
    bestbin(data[[x]], y, n, p, thres, uni, x, best, parallel)
  )
  names(cut_plan) <- names(data)
  end_time <- Sys.time()
  diff = end_time - start_time
  print(diff)
  attr(cut_plan, 'target') <- target
  structure(cut_plan, class = 'cut.plan')
}

#' @method update cut.plan
#' @export
#for manual updating
update.cut.plan <- function(cut_plan, ...) {
  un_arg <- list(...)
  if (class(cut_plan) != 'cut.plan') {
    stop('Not a cut plan!')
  } else if (length(un_arg) != 0) {
    for (i in 1:length(un_arg)) {
      if (!is.na(names(cut_plan[names(un_arg)[i]]))) {
        cut_plan[[names(un_arg)[i]]] <- un_arg[[i]]
        attr(cut_plan[names(un_arg)[i]], 'method') <- 'manual'
      } else
        stop('Column names are not matched!')
    }
  } else {
    print('Nothing happens!')
  }
  return(cut_plan)
}

#' @method predict cut.plan
#' @export
#for applying to data frame
predict.cut.plan <- function (cut_plan, data, keepTarget = FALSE) {
  if (class(cut_plan) != 'cut.plan') {
    stop('Not a cut plan!')
  } else {
    NewData <- data.frame(A = rep(0, nrow(data)))
    for (i in 1:length(cut_plan)) {
      if (!is.null(cut_plan[[i]])) {
        if (class(data[[names(cut_plan)[i]]]) %in% c('numeric', 'integer')) {
          col <- cut(data[[names(cut_plan)[i]]], breaks = c(-Inf, cut_plan[[i]], Inf))
        } else {
          stop(paste0(names(cut_plan)[i], ' is not numeric'))
        }
      } else if (names(cut_plan)[i] %in% colnames(data)){
        col <- data[,names(cut_plan)[i]]
      } else {
        col <- rep(NA, nrow(data))
      }
      NewData <- cbind(NewData, col)
    }
    NewData['A'] <- NULL
    names(NewData) <- names(cut_plan)
    if (keepTarget & attr(cut_plan, 'target') %in% names(data))
      NewData[attr(cut_plan, 'target')] <- data[attr(cut_plan, 'target')]
    if (!(attr(cut_plan, 'target') %in% names(data)))
      warning(paste0(attr(cut_plan, 'target'), ' is not exist!'))
    return(NewData)
  }
}
