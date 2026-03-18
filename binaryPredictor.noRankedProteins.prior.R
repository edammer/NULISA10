#library(caret)  #ADDI Windows03 VM, R v4.1.3
#library(glmnet)
#library(xgboost)
#library(ranger)
#library(progressr)
#library(dplyr)
#library(future)
#library(doFuture)
#library(doRNG)

## Prevent RStudio-specific error at very end of function ("error in evaluating the argument 'x' in selecting a method for function 'print': object 'spaces' not found")
#if (!exists("spaces", envir = asNamespace("cli")))
#    assign("spaces", c("", vapply(1:20, strrep, "", x = " ")),
#           envir = asNamespace("cli"))
## Just don't run in RStudio.

binaryPredictor <- function(expr, APOE_gt,
                                  nfold = 5, nrep = 5,
                                  ncores = parallel::detectCores() - 1,
                                  target=c("Control"),   # run one genotype vs all others
                                  target_ppv=0.95,
                                  seed   = 1) {

  set.seed(seed)

  if (!"package:cli" %in% search()) suppressPackageStartupMessages(library(cli))

  memLimit=4*1024^3
  options(future.globals.maxSize= memLimit)  #4GB Total size of all global objects that need to be exported - up from 500MB
  Sys.setenv(R_FUTURE_GLOBALS_MAXSIZE=memLimit) #inherited by workers

  ## Helper function - threshold minimum
  pick_thr <- function(prob, truth, target_ppv = 0.95,
                       min_tp = 30, floor = 0.80) {
    ok  <- !is.na(prob)                       # drop rows with NA prob
    prob <- prob[ok]
    truth <- truth[ok]

    if (sum(truth) == 0)  return(1)           # no positives
    if (all(prob == 0))   return(1)           # degenerate

    o   <- order(prob, decreasing = TRUE)
    tp  <- cumsum(truth[o] == 1)
    fp  <- cumsum(truth[o] == 0)
    ppv <- tp / (tp + fp)
    ok1  <- which(ppv >= target_ppv & tp >= min_tp)
    thr <- if (length(ok1)) prob[o[max(ok1)]] else quantile(prob[truth==1], .90, na.rm=TRUE)
    max(thr, floor)
  }

  ## ------------------------------------------------------------------------
  ## 0.  create resampling indices once
  ## ------------------------------------------------------------------------
  cvIndex <- createMultiFolds(APOE_gt, k = nfold, times = nrep)
  nTasks  <- length(cvIndex)
  ## ------------------------------------------------------------------------
  ## 1.  start the cluster
  ## ------------------------------------------------------------------------
  handlers("progress")
  handlers(global=TRUE)  # set handler for progress bar before the cluster

  plan(multisession, workers = ncores)  #, globals.maxSize=memLimit)  # or multicore on linux/macOS
  registerDoFuture()

  ## ------------------------------------------------------------------------
  ## 2.  run each fold in a worker  -----------------------------------------
  ## ------------------------------------------------------------------------

  with_progress({                              # << all progress lives here
    n_sub <- 3   # Progress bar increments 3x per fold
    p <- progressor(steps=length(cvIndex) * n_sub )     #along = cvIndex)          # one step per fold

    p(sprintf("Initializing %d workers...", ncores), amount=0)


  metrics <- foreach(fold = seq_along(cvIndex),
                     .combine   = rbind,
                     .options.future = list(  #expr not exported explicitly; captured only once.
                        globals = list(APOE_gt  = APOE_gt,
                                       pick_thr = pick_thr,
                                       target   = target,
                                       seed     = seed)
                     ),
                     .export = "p",            # let workers see 'p'
                     .packages  = c("progressr","glmnet","xgboost","ranger","dplyr","cli")) %dorng% {

    ## ----------------- announce fold start -----------------------------
    p(sprintf("fold %d/%d  -  started", fold, length(cvIndex)), amount=0)

    set.seed(seed + fold)                        # reproducible inside worker
    
    tr <- cvIndex[[fold]]
    te <- setdiff(seq_len(nrow(expr)), tr)
    
    tr_idx <- cvIndex[[fold]]
    te_idx <- setdiff(seq_len(nrow(expr)), tr_idx)

    ## --------- preprocessing  (Z scale, no missing data in input) -------
    prep <- function(m) scale(m[, colMeans(is.na(m)) <= .20, drop = FALSE])
    X_tr <- prep(expr[tr, ])
    X_te <- scale(expr[te, colnames(X_tr)],
                  center = attr(X_tr,"scaled:center"),
                  scale  = attr(X_tr,"scaled:scale"))
    X_te[is.na(X_te)] <- 0

    y_bin_tr <- factor(ifelse(APOE_gt[tr] == target, "pos", "neg"))
    y_bin_te <- factor(ifelse(APOE_gt[te] == target, "pos", "neg"))
    w_bin    <- ifelse(y_bin_tr == "pos", 8, 1)   # keep your weight idea

    # detect folds without positives OR without negatives
    if (length(unique(y_bin_tr)) < 2) {
      p(sprintf("fold %d/%d  *  skipped (only one class)", fold,
                length(cvIndex)), amount = 3)             # step the bar
      return(data.frame(Precision = NA, Recall = NA, Fold = fold))
    }
    
    # ----------   learner 1   glmnet  -------------------------------------
    glm_cv <- cv.glmnet(X_tr, y_bin_tr, family = "binomial",
                        weights = w_bin, type.measure = "class")
    p_glm <-    drop(predict(glm_cv, X_te, s = "lambda.min", type = "response"))
    p(sprintf("fold %d/%d  *  glmnet done", fold, length(cvIndex)), amount=1)
    
    # ----------   learner 2   xgboost  ------------------------------------
    dtr <- xgb.DMatrix(X_tr, label = as.numeric(y_bin_tr) - 1, weight = w_bin)
    dte <- xgb.DMatrix(X_te, label = as.numeric(y_bin_te) - 1)
    xpar <- list(objective = "binary:logistic", eta = 0.1,
                 max_depth = 6, subsample = 0.8, colsample_bytree = 0.8,
                 nthread = 1, eval_metric = "logloss")
    xgb <- xgb.train(xpar, dtr, watchlist=list(train=dtr, eval=dte), nrounds = 200,
                     verbose = 0, early_stopping_rounds = 20)
    p_xgb <- drop(predict(xgb, dte))
    p(sprintf("fold %d/%d  *  xgboost done", fold, length(cvIndex)), amount=1)

    # ----------   learner 3   ranger   ------------------------------------
    rf  <- ranger(y_bin_tr ~ ., data = data.frame(y_bin_tr, X_tr),
                  probability   = TRUE,
                  num.trees     = 500, num.threads = 1,
                  class.weights = c(neg = 1, pos = 8))
    p_rf <- predict(rf, data.frame(X_te))$predictions[,"pos"]
    p(sprintf("fold %d/%d  *  ranger done", fold, length(cvIndex)), amount=1)

    # ----------   averaged prob & hard threshold --------------------------
    # probabilities on training set 
    p_glm_tr <- drop(predict(glm_cv, X_tr, s = "lambda.min", type = "response"))
    p_xgb_tr <- drop(predict(xgb, dtr))
    p_rf_tr  <-        predict(rf, data.frame(X_tr))$predictions[,"pos"]
    p_avg_tr <- (p_glm_tr + p_xgb_tr + p_rf_tr) / 3    # length == length(y_bin_tr)

    thr   <- pick_thr(p_avg_tr, y_bin_tr == "pos", target_ppv)   # use helper
    
    # probabilities on test set 
    p_avg_te <- (p_glm + p_xgb + p_rf) / 3  # already computed on X_te
    pred  <- factor(ifelse(p_avg_te >= thr, target, NA_character_),
                    levels = c(target))
    
    tp <- sum(pred == target & y_bin_te == "pos", na.rm = TRUE)
    fp <- sum(pred == target & y_bin_te == "neg", na.rm = TRUE)
    fn <- sum(pred != target & y_bin_te == "pos", na.rm = TRUE)
    
    prec <- if ((tp+fp) > 0) tp/(tp+fp) else NA_real_
    rec  <- if ((tp+fn) > 0) tp/(tp+fn) else NA_real_
    
    data.frame(Precision = prec, Recall = rec, Fold = fold)
  } # foreach
  }) # with_progress

#  stopCluster(cl)                               # tidy up
#  registerDoSEQ()                               # back to sequential

  cat(sprintf("CV %s - Precision %.3f +/- %.3f | Recall %.3f +/- %.3f\n\n",
              target,
              mean(metrics$Precision, na.rm=TRUE), sd(metrics$Precision, na.rm=TRUE),
              mean(metrics$Recall,    na.rm=TRUE), sd(metrics$Recall,    na.rm=TRUE)))

  if (!exists("binaryPredictionMetrics", envir=.GlobalEnv)) assign("binaryPredictionMetrics",list(), envir=.GlobalEnv)
  binaryPredictionMetrics[[target]] <<- data.frame(Precision=metrics$Precision, Recall=metrics$Recall)

  ## ------------------------------------------------------------------------
  ## 3.  fit final ensemble on all data  (sequential) -----------------------
  ## ------------------------------------------------------------------------
  prep_expr <- function(mat) scale(mat[, colMeans(is.na(mat)) <= .20, drop = FALSE])
  X_all <- prep_expr(expr)
  y_all <- APOE_gt                                   # keep as character

  with_progress({
    p2 <- progressor(steps=length(target)*3)
    p2(sprintf("Starting binary fit of final ensemble on all data (serial/non-parallel)..."), amount=0)

    ## ----------  each genotype specified in fallback-targets, vs rest  -------
    ovr <- list()
    
    for (tg in target) {
    
      y_bin  <- factor(ifelse(y_all == tg, "pos", "neg"))
      w_pos  <- if (tg %in% c("MED13L","TCF4","CACNA1A")) 12 else 8  # higher weight for rarer genotypes
      w_bin  <- ifelse(y_bin == "pos", w_pos, 1)
    
      ## ---- (i) glmnet --------------------------------------------------------
      glm_tg <- cv.glmnet(X_all, y_bin, family = "binomial",
                          weights = w_bin, type.measure = "class")
      p2(sprintf("GLM fit on full data  *  finished"), amount=1)
    
      ## ---- (ii) xgboost ------------------------------------------------------
      d_bin  <- xgb.DMatrix(X_all, label = as.numeric(y_bin) - 1, weight = w_bin)
      xpar_b <- list(objective = "binary:logistic", eta = 0.1,
                     max_depth = 6, subsample = 0.8,
                     colsample_bytree = 0.8, nthread = ncores, eval_metric = "logloss")
      xgb_tg <- xgb.train(xpar_b, d_bin, watchlist=list(train=d_bin), nrounds = 200,
                          verbose = 0, early_stopping_rounds = 20)
      p2(sprintf("XGboost fit on full data  *  finished"), amount=1)

      ## ---- (iii) ranger ------------------------------------------------------
      rf_tg  <- ranger(y_bin ~ ., data = data.frame(y_bin, X_all),
                       probability = TRUE, num.trees = 500,
                       class.weights = c("neg" = 1, "pos" = w_pos),
                       num.threads  = ncores)
      p2(sprintf("Random Forest fit on full data  *  finished"), amount=1)
    
      ## ---- averaged prob on the *training* set ------------------------------
      p_glm <-      drop(predict(glm_tg, X_all, s = "lambda.min", type="response"))
      p_xgb <-      drop(predict(xgb_tg, d_bin))
      p_rf  <- rf_tg$predictions[,"pos"]
    
      p_avg <- (p_glm + p_xgb + p_rf) / 3
    
      ## learn a threshold that guarantees >= target_ppv on the 15 k labelled rows
      thr   <- pick_thr(p_avg, y_bin == "pos", target_ppv)
    
      ovr[[tg]] <- list(glm = glm_tg, xgb = xgb_tg, rf = rf_tg, thr = thr)
    }

    p2(sprintf("Binary fit on full data  *  finished"), amount=0)


## ---------------- feature importance -------------------------------- ###
    imp_glm <- abs(as.matrix(coef(ovr[[target]]$glm, s="lambda.min"))[-1,1])
    imp_xgb <- {
         g <- xgb.importance(model = ovr[[target]]$xgb)
         setNames(g$Gain, g$Feature)
    }
    imp_rf  <- ovr[[target]]$rf$variable.importance
  
    ## put everything on the same scale and average
    all_feats <- union(names(imp_glm), union(names(imp_xgb), names(imp_rf)))
    top_k=length(all_feats)

    imp_mat   <- cbind(
       glm = imp_glm[all_feats], xgb = imp_xgb[all_feats], rf = imp_rf[all_feats])
    imp_mat[is.na(imp_mat)] <- 0
    imp_scaled <- scale(imp_mat)
    imp_mean   <- rowMeans(imp_scaled)
    top_feats  <- sort(imp_mean, decreasing = TRUE) #[1:top_k]

  }) # with_progress
  
  ## -------- export top_k predictive proteins --------------------------
    if (!exists("rankedProteins", envir=.GlobalEnv)) assign("rankedProteins",list(), envir=.GlobalEnv)
    rankedProteins[[target]] <<- data.frame(feature=names(sort(top_feats, decreasing = TRUE)), importance=sort(top_feats, decreasing = TRUE))
    names(top_feats)<-gsub("\\|","_",names(top_feats))
#    print(knitr::kable(data.frame(Protein = names(top_feats)[c(1:10,(top_k-9):top_k)],
#                                  Importance = round(top_feats,3)[c(1:10,(top_k-9):top_k)]),
#                       caption = sprintf("Top 10 and bottom 10 (of top %d) predictive proteins - (all exported to list rankedProteins)", top_k)))
  ### ------------------------------------------------------------------- ###


  ## ------------------------------------------------------------------------
  ## 4.  prediction wrapper -------------------------------------------------
  ## ------------------------------------------------------------------------
    ## helper - probability only ---------------------------------------------
     prob_fun <- function(new_expr){
       # keep exactly the columns the model was trained with ------------------
       new_expr <- as.matrix(new_expr)[ , colnames(X_all), drop = FALSE]
       # use the same centering/scaling that is stored inside X_all ----------
       new_expr <- scale(new_expr,
                         center = attr(X_all, "scaled:center"),
                         scale  = attr(X_all, "scaled:scale"))
       new_expr[is.na(new_expr)] <- 0

       p_g <-      drop(predict(ovr[[target]]$glm, new_expr,
                                s = "lambda.min", type = "response"))
      p_x <-      drop(predict(ovr[[target]]$xgb, new_expr))
       p_r <- predict(ovr[[target]]$rf , data.frame(new_expr))$predictions[,"pos"]

       (p_g + p_x + p_r) / 3                # numeric vector of probabilities
     }

    ## main wrapper - returns factor but *carries* the probability -----------
     predict_wrapper <- function(new_expr){
       p_bin <- prob_fun(new_expr)
       pred  <- factor(ifelse(p_bin >= thr, target, NA_character_),
                       levels = target)
       attr(pred, "prob") <- p_bin          # <-- hand the prob back to caller
       pred
     }
 
     ## expose helper & threshold on the function object itself ---------------
     attr(predict_wrapper, "prob") <- prob_fun
     attr(predict_wrapper, "thr")  <- thr
     ## ---- remember the 50 features that were used to train this model ----
     attr(predict_wrapper, "features") <- colnames(X_all) #rankedProteins.prior[[target]][1:50, "feature"]

     predict_wrapper
}

# The above function is used to generate the 6 genotype-specific feature importance lists with all 7334 assays as input, and target ppv 0.95
#  (see code after genotype accounting below):

## (Don't have these vars until after running): binaryPredictionMetrics, rankedProteins






# not in memory:  rankedProteins.prior <- rankedProteins.prior[genotypes]


################################################################################
## 1.  Binary learner 2----------------------------------------------------------
################################################################################
binaryPredictor2 <- function(expr, APOE_gt,
                            target,                    # e.g. "e4/e4"
                            target_ppv = 0.95,         # PPV for hard-threshold
                            nfold = 5, nrep = 5,
                            ncores = parallel::detectCores() - 1,
                            seed   = 42)
{
  set.seed(seed)

  ## ---------- helpers --------------------------------------------------------
  prep_expr <- function(mat, ref = NULL) {   #Z scaling (run intra-fold) and handling of missing data (there is none in our input)
    keep <- which(colMeans(is.na(mat)) <= 0.20)
    if (is.null(ref)) {
      X <- scale(mat[, keep, drop = FALSE])
      list(x = X,
           center = attr(X, "scaled:center"),
           scale  = attr(X, "scaled:scale"),
           vars   = colnames(X))
    } else {
      X <- scale(mat[, ref$vars, drop = FALSE],
                 center = ref$center,
                 scale  = ref$scale)
      X[is.na(X)] <- 0
      list(x = X)
    }
  }

  pick_thr <- function(prob, truth, target_ppv = 0.95,
                       min_tp = 30, floor = 0.80) {
    ok <- !is.na(prob)
    prob <- prob[ok]; truth <- truth[ok]
    if (sum(truth) == 0 || all(prob == 0)) return(1)
    ord <- order(prob, decreasing = TRUE)
    tp <- cumsum(truth[ord] == 1); fp <- cumsum(truth[ord] == 0)
    ppv <- tp / (tp + fp)
    idx <- which(ppv >= target_ppv & tp >= min_tp)
    thr <- if (length(idx)) prob[ord[max(idx)]] else
                         quantile(prob[truth == 1], .90, na.rm = TRUE)
    max(thr, floor)
  }

  ## ---------- inner CV  ------------------------------------------------------
  cvIndex <- caret::createMultiFolds(APOE_gt, k = nfold, times = nrep)

  ## -- future backend (no nested progress objects) ------------
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(multisession, workers = ncores)
  doFuture::registerDoFuture()
  doRNG::registerDoRNG(seed)                # reproducible but no %dorng

    cv_stats <- foreach(fold = seq_along(cvIndex),
            .packages = c("glmnet", "xgboost", "ranger", "caret", "progressr", "dplyr", "cli"),
            .combine  = dplyr::bind_rows,
            .export   = c("pick_thr", "APOE_gt", "target", "seed"),
            .options.RNG = seed) %dopar% {

      expr_loc <- expr
      tr <- cvIndex[[fold]]
      te <- setdiff(seq_len(nrow(expr_loc)), tr)

      ## ---------- scaling ----------------------------------------------------
      X_tr <- prep_expr(expr[tr, ])$x
      ref  <- list(center = attr(X_tr, "scaled:center"),
                   scale  = attr(X_tr, "scaled:scale"),
                   vars   = colnames(X_tr))
      X_te <- prep_expr(expr[te, ], ref)$x

      ## ---------- binarise labels -------------------------------------------
      y_tr <- factor(ifelse(APOE_gt[tr] == target, "pos", "neg"))
      y_te <- factor(ifelse(APOE_gt[te] == target, "pos", "neg"))
      w_tr <- ifelse(y_tr == "pos", 8, 1)

      if (length(unique(y_tr)) < 2) {
#        p(message = sprintf("fold %d – skipped (single class)", fold), amount = 3)
        return(data.frame(Precision = NA, Recall = NA))
      }

      ## ---------- glmnet -----------------------------------------------------
      glm_cv <- glmnet::cv.glmnet(X_tr, y_tr, family = "binomial",
                                  weights = w_tr, type.measure = "class")
      p_glm <- drop(predict(glm_cv, X_te, s = "lambda.min",
                            type = "response"))
#      p(amount = 1)

      ## ---------- xgboost ----------------------------------------------------
      dtr <- xgboost::xgb.DMatrix(X_tr, label = as.numeric(y_tr) - 1,
                                  weight = w_tr)
      dte <- xgboost::xgb.DMatrix(X_te, label = as.numeric(y_te) - 1)
      xpar <- list(objective = "binary:logistic",
                   eta = 0.1, max_depth = 6,
                   subsample = 0.8, colsample_bytree = 0.8,
                   nthread = 1, eval_metric = "logloss")
      bst <- xgboost::xgb.train(xpar, dtr,
                                watchlist = list(train = dtr, eval = dte),
                                nrounds = 200, verbose = 0,
                                early_stopping_rounds = 20)
      p_xgb <- drop(predict(bst, dte))
#      p(amount = 1)

      ## ---------- ranger -----------------------------------------------------
      rf <- ranger::ranger(y_tr ~ ., data = data.frame(y_tr, X_tr),
                           probability   = TRUE,
                           num.trees     = 500,
                           num.threads   = 1,
                           class.weights = c(neg = 1, pos = 8))
      p_rf <- predict(rf, data.frame(X_te))$predictions[, "pos"]
#      p(amount = 1)

      ## ---------- evaluate ---------------------------------------------------
      p_avg_tr <- rowMeans(cbind(
        drop(predict(glm_cv, X_tr, s = "lambda.min", type = "response")),
        drop(predict(bst, dtr)),
        predict(rf, data.frame(X_tr))$predictions[, "pos"]))

      thr <- pick_thr(p_avg_tr, APOE_gt[tr] == target, target_ppv)

      p_avg_te <- rowMeans(cbind(p_glm, p_xgb, p_rf))
      pred <- ifelse(p_avg_te >= thr, "pos", "neg")

      tp <- sum(pred == "pos" & y_te == "pos")
      fp <- sum(pred == "pos" & y_te == "neg")
      fn <- sum(pred == "neg" & y_te == "pos")

      data.frame(
        Precision = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
        Recall    = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_)
    } #}) no more with_progress

  message(sprintf("[%s]  CV Precision %.3f +/- %.3f | Recall %.3f +/- %.3f",
                  target,
                  mean(cv_stats$Precision, na.rm = TRUE),
                  sd  (cv_stats$Precision, na.rm = TRUE),
                  mean(cv_stats$Recall,    na.rm = TRUE),
                  sd  (cv_stats$Recall,    na.rm = TRUE)))

  ## ---------- fit once on all data ------------------------------------------
  prep_all <- prep_expr(expr)
  X_all <- prep_all$x
  y_all <- factor(ifelse(APOE_gt == target, "pos", "neg"))
  w_all <- ifelse(y_all == "pos",
                  if (target %in% c("MED13L","TCF4","CACNA1A")) 12 else 8,
                  1)

  glm_fit <- glmnet::cv.glmnet(X_all, y_all, family = "binomial",
                               weights = w_all, type.measure = "class")

  d_all  <- xgboost::xgb.DMatrix(X_all, label = as.numeric(y_all) - 1,
                                 weight = w_all)
  xpar_b <- list(objective = "binary:logistic",
                 eta = 0.1, max_depth = 6,
                 subsample = 0.8, colsample_bytree = 0.8,
                 nthread = ncores, eval_metric = "logloss")
  xgb_fit <- xgboost::xgb.train(xpar_b, d_all,
                                watchlist = list(train = d_all),
                                nrounds   = 200,
                                verbose   = 0,
                                early_stopping_rounds = 20)

  rf_fit <- ranger::ranger(y_all ~ ., data = data.frame(y_all, X_all),
                           probability   = TRUE,
                           num.trees     = 500,
                           class.weights = c(neg = 1,
                                             pos = if (target %in%
                                                       c("MED13L","TCF4","CACNA1A"))
                                                        12 else 8),
                           num.threads   = ncores)

  p_all <- rowMeans(cbind(
    drop(predict(glm_fit, X_all, s = "lambda.min", type = "response")),
    drop(predict(xgb_fit, d_all)),
    rf_fit$predictions[, "pos"]))

  thr <- pick_thr(p_all, APOE_gt == target, target_ppv)

  ## ---------- wrapper with probability accessor -----------------------------
  prob_fun <- function(new_expr) {
    new_expr <- prep_expr(new_expr, prep_all)$x
    rowMeans(cbind(
      drop(predict(glm_fit, new_expr, s = "lambda.min", type = "response")),
      drop(predict(xgb_fit, new_expr)),
      predict(rf_fit, data.frame(new_expr))$predictions[, "pos"]))
  }

  predict_fun <- function(new_expr) {
    p_bin <- prob_fun(new_expr)
    res   <- factor(ifelse(p_bin >= thr, target, NA_character_),
                    levels = target)
    attr(res, "prob") <- p_bin
    res
  }

  attr(predict_fun, "prob")     <- prob_fun
  attr(predict_fun, "thr")      <- thr
  attr(predict_fun, "features") <- colnames(expr)

  predict_fun
}

################################################################################
## 2.  Helper – train the 10 binary models & build a 10-class predictor ---------
################################################################################
train_NULISA10_full <- function(expr, APOE_gt,
                             rankedProteins.prior,        # list of top-50 feats
                             target_ppv = 0.95,
                             ncores = parallel::detectCores() - 1,
                             topX = 10,
                             seed = 42)
{
  genotypes <- names(rankedProteins.prior)

  ## one binary fit per genotype – each with *its own* topX columns ------------
  fit_fns <- lapply(genotypes, function(gt) {
               feats <- rankedProteins.prior[[gt]][1:topX, "feature"]
               binaryPredictor2(expr   = expr[, feats, drop = FALSE],
                               APOE_gt = APOE_gt,
                               target  = gt,
                               target_ppv = target_ppv,
                               ncores  = ncores,
                               seed    = seed)
             })
  names(fit_fns) <- genotypes

  ## ----------- helper to query *all* probabilities for one sample ----------
  score_all <- function(new_expr_row) {
    sapply(genotypes, function(gt) {
      fn    <- fit_fns[[gt]]
      feats <- attr(fn, "features")
      res   <- fn(new_expr_row[, feats, drop = FALSE])
      attr(res, "prob")
    })
  }

  ## ----------- final 6-class wrapper ---------------------------------------
  predict_NULISA10 <- function(new_expr,
                               ncores = parallel::detectCores()) {

    if (is.data.frame(new_expr)) new_expr <- as.matrix(new_expr)

    fit_fns   <- attr(predict_NULISA10, "fit_fns")
    genotypes <- attr(predict_NULISA10, "genotypes")

    # helper: safe threshold extractor
    get_thr <- function(gt) {
      fn <- fit_fns[[gt]]
      thr <- attr(fn, "thr")

      if (is.null(fn)) {
        warning(sprintf("Missing model for genotype '%s'", gt))
        return(NA_real_)
      }
      if (is.null(thr) || length(thr) != 1 || !is.numeric(thr)) {
        warning(sprintf("Invalid or missing threshold for genotype '%s'", gt))
        return(NA_real_)
      }
      thr
    }

    thr_vec <- vapply(genotypes, get_thr, numeric(1))

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(multisession, workers = ncores)

    prob_mat <- t(future.apply::future_apply(
      new_expr, 1L, future.seed = TRUE,
      future.packages = c("glmnet", "xgboost", "ranger"),
      function(row_vec) {
        mat <- matrix(row_vec, nrow = 1,
                      dimnames = list(NULL, colnames(new_expr)))
        score_all(mat)
      }))

    res <- character(nrow(prob_mat))

    for (k in seq_len(nrow(prob_mat))) {
      p_row <- prob_mat[k, ]
      names(p_row) <- genotypes

      above_thr <- p_row >= thr_vec & !is.na(thr_vec)

      if (any(above_thr)) {
        margin <- p_row[above_thr] - thr_vec[above_thr]
        res[k] <- names(margin)[which.max(margin)]
      } else {
        res[k] <- names(p_row)[which.max(p_row)]
      }
    }

    factor(res, levels = genotypes)
  }

  ## ship out: the 6-class wrapper plus access to internals ------------------
  attr(predict_NULISA10, "fit_fns")   <- fit_fns
  attr(predict_NULISA10, "genotypes") <- genotypes
  predict_NULISA10
}



