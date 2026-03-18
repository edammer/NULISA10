#####################################################################################
## WGCNA Coexpression Network Build - Seyfried Systems Biology Pipeline
##
## input cleanDat from:  plasma transposed 2 batch variable regressed NULISA NPQ data
##
## Eric Dammer -- Seyfried Lab Proteomics (2026)  5b run 03/18/2026
##             -- Epilepsy prediction confusion matrix 100 fold crossvalidation 80/20
#####################################################################################

rootdir2="F:/OneDrive - Emory/Faundez_NULISA10/5.prediction/"
setwd(rootdir2)

load("../4.subtyping/4b.saved.image.NULISA10.plasma-eigensampleWGCNA.Rdata")
#contains: cleanDat,numericMeta
rootdir<-rootdir2


dim(cleanDat)
#[1] 465 130
#^transposed (465 samples as rows)


source("binaryPredictor.noRankedProteins.prior.R")

library(caret)
library(glmnet)
library(xgboost)
library(ranger)
library(progressr)
library(dplyr)
library(future)
library(doFuture)
library(doRNG)


### Replacement train_NULISA10_full() hardened to:
# 1. resolve requested features safely against available columns
# 2. tolerate fewer than topX available features
# 3. pad missing prediction-time columns with NA
# 4. force the subset matrices eagerly so future does not later re-evaluate a fragile promise
# 5. optionally keep the outer loop going even if a fold still fails for another reason

train_NULISA10_full <- function(expr, APOE_gt,
                                rankedProteins.prior,
                                target_ppv = 0.95,
                                ncores = parallel::detectCores() - 1,
                                topX = 10,
                                seed = 42,
                                verbose = TRUE)
{
  expr <- as.matrix(expr)
  if (is.null(colnames(expr))) stop("'expr' must have column names.")
  genotypes <- names(rankedProteins.prior)

  ## ---------- helpers ----------------------------------------------------
  norm_name <- function(x) gsub("\\|", "_", as.character(x))

  resolve_features <- function(feats, available_cols) {
    feats <- unique(stats::na.omit(as.character(feats)))
    if (!length(feats)) return(character(0))

    ## exact matches first
    exact <- feats[feats %in% available_cols]

    ## normalized fallback: treat "|" and "_" as equivalent
    miss <- setdiff(feats, exact)
    if (!length(miss)) return(unique(exact))

    avail_norm <- norm_name(available_cols)
    miss_norm  <- norm_name(miss)
    idx <- match(miss_norm, avail_norm)
    mapped <- available_cols[idx[!is.na(idx)]]

    unique(c(exact, mapped))
  }

  align_newdata_to_features <- function(new_expr_row, feats_needed) {
    new_expr_row <- as.matrix(new_expr_row)
    if (is.null(colnames(new_expr_row))) {
      stop("Prediction matrix/row must have column names.")
    }

    feats_needed <- unique(as.character(feats_needed))
    present <- intersect(feats_needed, colnames(new_expr_row))
    missing <- setdiff(feats_needed, present)

    if (length(missing)) {
      pad <- matrix(NA_real_, nrow = nrow(new_expr_row), ncol = length(missing),
                    dimnames = list(rownames(new_expr_row), missing))
      new_expr_row <- cbind(new_expr_row, pad)
    }

    new_expr_row[, feats_needed, drop = FALSE]
  }

  ## ---------- fit one binary model per class ------------------------------
  fit_fns <- vector("list", length(genotypes))
  names(fit_fns) <- genotypes

  for (j in seq_along(genotypes)) {
    gt <- genotypes[j]

    if (!all(c("feature") %in% colnames(rankedProteins.prior[[gt]]))) {
      stop(sprintf("rankedProteins.prior[['%s']] must contain a 'feature' column.", gt))
    }

    raw_feats <- rankedProteins.prior[[gt]][seq_len(min(topX, nrow(rankedProteins.prior[[gt]]))), "feature"]
    feats_use <- resolve_features(raw_feats, colnames(expr))

    if (length(feats_use) == 0) {
      stop(sprintf("No usable features found for target '%s'. Requested: %s",
                   gt, paste(raw_feats, collapse = ", ")))
    }

    if (verbose && length(feats_use) < min(topX, length(stats::na.omit(raw_feats)))) {
      warning(sprintf("[%s] only %d/%d requested features were found in expr for this fold.",
                      gt, length(feats_use), min(topX, length(stats::na.omit(raw_feats)))))
    }

    ## force the subset now so future does not later try to realize a fragile promise
    expr_sub <- expr[, feats_use, drop = FALSE]
    expr_sub <- as.matrix(expr_sub)
    force(expr_sub)
    force(feats_use)

    fit_fns[[gt]] <- local({
      gt_local   <- gt
      expr_local <- expr_sub
      feats_local <- feats_use
      seed_local <- seed + j

      fn <- binaryPredictor2(
        expr       = expr_local,
        APOE_gt    = APOE_gt,
        target     = gt_local,
        target_ppv = target_ppv,
        ncores     = ncores,
        seed       = seed_local
      )

      attr(fn, "features") <- colnames(expr_local)
      attr(fn, "resolved_features") <- feats_local
      fn
    })
  }

  ## ---------- score all classes for one sample ----------------------------
  score_all <- function(new_expr_row) {
    sapply(genotypes, function(gt) {
      fn <- fit_fns[[gt]]
      feats_needed <- attr(fn, "features")
      row_aligned <- align_newdata_to_features(new_expr_row, feats_needed)
      res <- fn(row_aligned)
      attr(res, "prob")
    })
  }

  ## ---------- final multiclass predictor ---------------------------------
  predict_NULISA10 <- function(new_expr,
                               ncores = parallel::detectCores()) {

    if (is.data.frame(new_expr)) new_expr <- as.matrix(new_expr)
    if (is.null(colnames(new_expr))) stop("'new_expr' must have column names.")

    fit_fns   <- attr(predict_NULISA10, "fit_fns")
    genotypes <- attr(predict_NULISA10, "genotypes")

    thr_vec <- vapply(genotypes, function(gt) {
      fn <- fit_fns[[gt]]
      thr <- attr(fn, "thr")
      if (is.null(thr) || length(thr) != 1 || !is.numeric(thr)) NA_real_ else thr
    }, numeric(1))

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = ncores)

    prob_mat <- t(future.apply::future_apply(
      new_expr, 1L, future.seed = TRUE,
      future.packages = c("glmnet", "xgboost", "ranger"),
      function(row_vec) {
        mat <- matrix(row_vec, nrow = 1, dimnames = list(NULL, colnames(new_expr)))
        score_all(mat)
      }
    ))

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

  attr(predict_NULISA10, "fit_fns")   <- fit_fns
  attr(predict_NULISA10, "genotypes") <- genotypes
  predict_NULISA10
} 

## Previously, all known status samples used for training, and predicted those same knowns.
## Here, we do an 80/20 hold-out evaluation, so reported performance is out-of-sample without resubstitution.

## rankedProteins.epi is not computed from all 370 labeled subjects to avoid information leakage for feature selection

## Wrapper:

## ------------------------------------------------------------------
## Helper: run stage-1 epilepsy feature ranking on one training split
##         and return only the fold-specific rankedProteins list
## ------------------------------------------------------------------
get_epilepsy_rankings_fold <- function(expr_train, y_train,
                                       stage1_ppv = 0.95,
                                       ncores = 8,
                                       seed = 1) {
  targets <- c("Epilepsy", "No")

  ## preserve any existing globals that binaryPredictor() mutates
  had_rp  <- exists("rankedProteins", envir = .GlobalEnv, inherits = FALSE)
  had_bpm <- exists("binaryPredictionMetrics", envir = .GlobalEnv, inherits = FALSE)

  old_rp  <- if (had_rp)  get("rankedProteins", envir = .GlobalEnv) else NULL
  old_bpm <- if (had_bpm) get("binaryPredictionMetrics", envir = .GlobalEnv) else NULL

  on.exit({
    if (had_rp) {
      assign("rankedProteins", old_rp, envir = .GlobalEnv)
    } else if (exists("rankedProteins", envir = .GlobalEnv, inherits = FALSE)) {
      rm("rankedProteins", envir = .GlobalEnv)
    }

    if (had_bpm) {
      assign("binaryPredictionMetrics", old_bpm, envir = .GlobalEnv)
    } else if (exists("binaryPredictionMetrics", envir = .GlobalEnv, inherits = FALSE)) {
      rm("binaryPredictionMetrics", envir = .GlobalEnv)
    }
  }, add = TRUE)

  invisible(lapply(seq_along(targets), function(j) {
    gt <- targets[j]
    binaryPredictor(
      expr       = expr_train,
      APOE_gt    = y_train,
      target     = gt,
      ncores     = ncores,
      target_ppv = stage1_ppv,
      seed       = seed + j
    )
  }))

  rp <- get("rankedProteins", envir = .GlobalEnv)
  rp[targets]
}


## ============================================================
## Helper: safely get top epilepsy features from one dataset
##        using stage-1 binaryPredictor(target = "Epilepsy")
## ============================================================
get_epilepsy_top_features <- function(expr, y,
                                      topX = 12,
                                      stage1_ppv = 0.95,
                                      ncores = 8,
                                      seed = 1) {
  expr <- as.matrix(expr)
  y <- as.character(y)
  y[y == "Yes"] <- "Epilepsy"
  y <- factor(y, levels = c("Epilepsy", "No"))

  ## preserve globals that binaryPredictor() mutates
  had_rp  <- exists("rankedProteins", envir = .GlobalEnv, inherits = FALSE)
  had_bpm <- exists("binaryPredictionMetrics", envir = .GlobalEnv, inherits = FALSE)

  old_rp  <- if (had_rp)  get("rankedProteins", envir = .GlobalEnv) else NULL
  old_bpm <- if (had_bpm) get("binaryPredictionMetrics", envir = .GlobalEnv) else NULL

  on.exit({
    if (had_rp) {
      assign("rankedProteins", old_rp, envir = .GlobalEnv)
    } else if (exists("rankedProteins", envir = .GlobalEnv, inherits = FALSE)) {
      rm("rankedProteins", envir = .GlobalEnv)
    }

    if (had_bpm) {
      assign("binaryPredictionMetrics", old_bpm, envir = .GlobalEnv)
    } else if (exists("binaryPredictionMetrics", envir = .GlobalEnv, inherits = FALSE)) {
      rm("binaryPredictionMetrics", envir = .GlobalEnv)
    }
  }, add = TRUE)

  invisible(
    binaryPredictor(
      expr       = expr,
      APOE_gt    = y,
      target     = "Epilepsy",
      ncores     = ncores,
      target_ppv = stage1_ppv,
      seed       = seed
    )
  )

  rp <- get("rankedProteins", envir = .GlobalEnv)
  if (!"Epilepsy" %in% names(rp)) {
    stop("Stage-1 ranking did not produce rankedProteins[['Epilepsy']].")
  }
  if (!"feature" %in% colnames(rp[["Epilepsy"]])) {
    stop("rankedProteins[['Epilepsy']] lacks a 'feature' column.")
  }

  unique(as.character(head(rp[["Epilepsy"]][, "feature"], topX)))
}


## ============================================================
## Main: repeated 80/20 hold-out CV for one-model epilepsy call
##
## Mode A (recommended): fixed_features = NULL
##   -> re-rank top features inside each outer training split
##
## Mode B (leaky / sensitivity analysis): fixed_features = ...
##   -> use the same feature list in every fold
## ============================================================
epilepsy_binary_holdout_cv <- function(cleanDat,
                                       epilepsy,
                                       n_outer = 100,
                                       hold_frac = 0.20,
                                       topX = 12,
                                       stage1_ppv = 0.95,
                                       stage2_ppv = 0.80,
                                       ncores = 8,
                                       seed = 1,
                                       fixed_features = NULL,
                                       continue_on_error = TRUE) {
  targets <- c("Epilepsy", "No")

  ## align known-label samples to expression matrix
  y <- as.character(epilepsy)
  names(y) <- names(epilepsy)
  y[y == "Yes"] <- "Epilepsy"
  y <- factor(y, levels = targets)

  expr_all <- t(na.omit(t(cleanDat)))
  keep_ids <- intersect(names(y), rownames(expr_all))
  expr_all <- as.matrix(expr_all[keep_ids, , drop = FALSE])
  y <- factor(as.character(y[keep_ids]), levels = targets)

  if (nrow(expr_all) != length(y)) stop("Expression matrix and labels are misaligned.")
  if (anyNA(y)) stop("Labels contain NA after alignment.")
  if (length(unique(y)) < 2) stop("Need both Epilepsy and No classes.")

  ## optional fixed features: normalize once
  if (!is.null(fixed_features)) {
    fixed_features <- unique(as.character(fixed_features))
    fixed_features <- intersect(fixed_features, colnames(expr_all))
    if (!length(fixed_features)) stop("None of fixed_features were found in expr_all.")
  }

  outer_res <- vector("list", n_outer)

  for (i in seq_len(n_outer)) {
    fit_one_fold <- function() {
      set.seed(seed + i)

      train_idx <- caret::createDataPartition(y, p = 1 - hold_frac, list = FALSE)[, 1]
      test_idx  <- setdiff(seq_len(nrow(expr_all)), train_idx)

      expr_tr <- expr_all[train_idx, , drop = FALSE]
      expr_te <- expr_all[test_idx,  , drop = FALSE]
      y_tr    <- factor(as.character(y[train_idx]), levels = targets)
      y_te    <- factor(as.character(y[test_idx]),  levels = targets)

      ## choose features
      feats <- if (is.null(fixed_features)) {
        get_epilepsy_top_features(
          expr       = expr_tr,
          y          = y_tr,
          topX       = topX,
          stage1_ppv = stage1_ppv,
          ncores     = ncores,
          seed       = seed + 1000 + i
        )
      } else {
        fixed_features
      }

      feats <- intersect(feats, colnames(expr_tr))
      if (!length(feats)) stop(sprintf("Fold %d: no usable features after intersection.", i))

      ## fit one binary Epilepsy-vs-rest model on training split
      epi_model <- binaryPredictor2(
        expr       = expr_tr[, feats, drop = FALSE],
        APOE_gt    = y_tr,
        target     = "Epilepsy",
        target_ppv = stage2_ppv,
        ncores     = ncores,
        seed       = seed + 2000 + i
      )

      ## predict hold-out
      pred_pos <- epi_model(expr_te[, feats, drop = FALSE])

      ## binaryPredictor2 returns target or NA; map NA -> No
      y_hat <- ifelse(is.na(as.character(pred_pos)), "No", "Epilepsy")
      y_hat <- factor(y_hat, levels = targets)

      ## rows = Predicted, cols = Truth
      cm <- table(Predicted = y_hat, Truth = y_te)
      cm_full <- matrix(0L, nrow = 2, ncol = 2,
                        dimnames = list(Predicted = targets, Truth = targets))
      cm_full[rownames(cm), colnames(cm)] <- cm

      TP <- unname(cm_full["Epilepsy", "Epilepsy"])
      FP <- unname(cm_full["Epilepsy", "No"])
      FN <- unname(cm_full["No",       "Epilepsy"])
      TN <- unname(cm_full["No",       "No"])

      acc  <- (TP + TN) / sum(cm_full)
      sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
      spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
      ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
      npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
      f1   <- if ((2 * TP + FP + FN) > 0) 2 * TP / (2 * TP + FP + FN) else NA_real_
      bal  <- mean(c(sens, spec), na.rm = TRUE)

      list(
        ok = TRUE,
        Fold = i,
        n_train = length(train_idx),
        n_test = length(test_idx),
        features = feats,
        n_features = length(feats),
        ConfMat = cm_full,
        TP = TP, FP = FP, FN = FN, TN = TN,
        Accuracy = acc,
        Sensitivity = sens,
        Specificity = spec,
        PPV = ppv,
        NPV = npv,
        F1 = f1,
        BalancedAccuracy = bal,
        error = ""
      )
    }

    if (continue_on_error) {
      outer_res[[i]] <- tryCatch(
        fit_one_fold(),
        error = function(e) {
          message(sprintf("Outer fold %d failed: %s", i, conditionMessage(e)))
          list(
            ok = FALSE,
            Fold = i,
            n_train = NA_integer_,
            n_test = NA_integer_,
            features = character(0),
            n_features = NA_integer_,
            ConfMat = matrix(NA_integer_, nrow = 2, ncol = 2,
                             dimnames = list(Predicted = targets, Truth = targets)),
            TP = NA_integer_, FP = NA_integer_, FN = NA_integer_, TN = NA_integer_,
            Accuracy = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
            PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, BalancedAccuracy = NA_real_,
            error = conditionMessage(e)
          )
        }
      )
    } else {
      outer_res[[i]] <- fit_one_fold()
    }
  }

  ## per-fold summary table
  per_fold <- do.call(rbind, lapply(outer_res, function(x) {
    data.frame(
      Fold = x$Fold,
      ok = x$ok,
      n_train = x$n_train,
      n_test = x$n_test,
      n_features = x$n_features,
      TP = x$TP, FP = x$FP, FN = x$FN, TN = x$TN,
      Accuracy = x$Accuracy,
      Sensitivity = x$Sensitivity,
      Specificity = x$Specificity,
      PPV = x$PPV,
      NPV = x$NPV,
      F1 = x$F1,
      BalancedAccuracy = x$BalancedAccuracy,
      error = x$error,
      stringsAsFactors = FALSE
    )
  }))

  good <- which(vapply(outer_res, function(x) isTRUE(x$ok), logical(1)))

  cm_arr <- array(
    0L,
    dim = c(2, 2, length(good)),
    dimnames = list(Predicted = targets, Truth = targets, Fold = paste0("fold", good))
  )
  for (j in seq_along(good)) {
    cm_arr[, , j] <- outer_res[[good[j]]]$ConfMat
  }

  pooled_confmat <- if (length(good)) apply(cm_arr, c(1, 2), sum) else NULL
  mean_confmat   <- if (length(good)) apply(cm_arr, c(1, 2), mean) else NULL
  sd_confmat     <- if (length(good)) apply(cm_arr, c(1, 2), sd) else NULL

  pooled_counts <- if (length(good)) {
    c(
      TP = unname(pooled_confmat["Epilepsy", "Epilepsy"]),
      FP = unname(pooled_confmat["Epilepsy", "No"]),
      FN = unname(pooled_confmat["No",       "Epilepsy"]),
      TN = unname(pooled_confmat["No",       "No"])
    )
  } else {
    c(TP = NA, FP = NA, FN = NA, TN = NA)
  }

  metric_summary <- data.frame(
    Metric = c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "F1", "BalancedAccuracy"),
    Mean = c(
      mean(per_fold$Accuracy, na.rm = TRUE),
      mean(per_fold$Sensitivity, na.rm = TRUE),
      mean(per_fold$Specificity, na.rm = TRUE),
      mean(per_fold$PPV, na.rm = TRUE),
      mean(per_fold$NPV, na.rm = TRUE),
      mean(per_fold$F1, na.rm = TRUE),
      mean(per_fold$BalancedAccuracy, na.rm = TRUE)
    ),
    SD = c(
      sd(per_fold$Accuracy, na.rm = TRUE),
      sd(per_fold$Sensitivity, na.rm = TRUE),
      sd(per_fold$Specificity, na.rm = TRUE),
      sd(per_fold$PPV, na.rm = TRUE),
      sd(per_fold$NPV, na.rm = TRUE),
      sd(per_fold$F1, na.rm = TRUE),
      sd(per_fold$BalancedAccuracy, na.rm = TRUE)
    )
  )

  list(
    per_fold = per_fold,
    metric_summary = metric_summary,
    pooled_confmat = pooled_confmat,   # rows = Predicted, cols = Truth
    mean_confmat = mean_confmat,
    sd_confmat = sd_confmat,
    pooled_counts = pooled_counts,     # TP, FP, FN, TN
    confmat_array = cm_arr,
    features_by_fold = lapply(outer_res, `[[`, "features"),
    fixed_features = fixed_features,
    n_successful_folds = length(good),
    n_failed_folds = n_outer - length(good),
    outer_res = outer_res
  )
}

######################################
### Stage 1 + 2 - predict epilepsy


## known epilepsy labels only
epilepsy <- numericMeta$Epilepsy[!is.na(numericMeta$Epilepsy)]
names(epilepsy) <- rownames(numericMeta)[!is.na(numericMeta$Epilepsy)]
epilepsy[epilepsy == "Yes"] <- "Epilepsy"

table(epilepsy)
#Epilepsy       No 
#     149      222

length(epilepsy)
# 371


## Top features of all known label data
fixed_epi_features <- get_epilepsy_top_features(
  expr       = t(na.omit(t(cleanDat)))[names(epilepsy), , drop = FALSE],
  y          = epilepsy,
  topX       = 12,
  stage1_ppv = 0.95,
  ncores     = 25,
  seed       = 1
)

fixed_epi_features
#[1] "HBA1"   "IGFBP7" "CRH"    "NEFL"   "FLT1"   "ANXA5"  "DDC"    "IL10"   "VCAM1"  "CCL4"   "CRP"    "MAPT"


## Run CV with fixed list for sensitivity analysis
epi_cv_fixed <- epilepsy_binary_holdout_cv(
  cleanDat    = t(na.omit(t(cleanDat)))[names(epilepsy), , drop = FALSE],
  epilepsy    = epilepsy,
  n_outer     = 100,
  hold_frac   = 0.20,
  topX        = 12,
  stage1_ppv  = 0.95,
  stage2_ppv  = 0.80,
  ncores      = 25,
  seed        = 1,
  fixed_features = fixed_epi_features
)

epi_cv_fixed$metric_summary
#            Metric      Mean         SD
#1         Accuracy 0.6909589 0.03721442
#2      Sensitivity 0.3389655 0.08202025
#3      Specificity 0.9229545 0.03765706
#4              PPV 0.7477038 0.10290657
#5              NPV 0.6801999 0.02744837
#6               F1 0.4611859 0.08730187
#7 BalancedAccuracy 0.6309600 0.04284795

## how to read the 2x2 confusion matrix returned:
## rows = Predicted, cols = Truth
epi_cv_fixed$pooled_confmat
#          Truth
#Predicted  Epilepsy   No
#  Epilepsy      983  339
#  No           1917 4061
  
## therefore:
TP <- epi_cv_fixed$pooled_confmat["Epilepsy", "Epilepsy"]
FP <- epi_cv_fixed$pooled_confmat["Epilepsy", "No"]
FN <- epi_cv_fixed$pooled_confmat["No",       "Epilepsy"]
TN <- epi_cv_fixed$pooled_confmat["No",       "No"]

epi_cv_fixed$pooled_counts
#  TP   FP   FN   TN 
# 983  339 1917 4061
epi_cv_fixed$mean_confmat
#          Truth
#Predicted  Epilepsy    No
#  Epilepsy     9.83  3.39
#  No          19.17 40.61
epi_cv_fixed$sd_confmat
#          Truth
#Predicted  Epilepsy       No
#  Epilepsy 2.378587 1.656911
#  No       2.378587 1.656911

head(epi_cv_fixed$per_fold)





## Confirm top features for positive prediction of epilepsy

#known epilepsy sample data, no missing values:
#cleanDat.epilepsy=t(na.omit(t(cleanDat)))[names(epilepsy),]

epilepsy.fit_fn <- lapply(c("Epilepsy","No"), function(gt)
    binaryPredictor(expr =t(na.omit(t(cleanDat)))[names(epilepsy),] ,
    #                    = t(training.cleanDat.noNA[rankedProteins.prior[[gt]][1:50,"feature"],]),     # your 15 000 × p matrix
                          APOE_gt= epilepsy,   # training.gt.APOE,  # factor of length 15 000
                          target = gt,
                          ncores=8,  # parallel::detectCores() - 1,
                          target_ppv=0.95,  # HERE WE USE DEFAULT, 0.95; in stage 2, with only 50 features, we set this to 0.50 (but the floor is 0.80 in the helper)
                          seed   = 1))

#(previously) CV Yes - Precision 0.887 +/- 0.148 | Recall 1.000 +/- 0.000
#CV Epilepsy - Precision 0.931 +/- 0.149 | Recall 1.000 +/- 0.000
#CV No - Precision 0.876 +/- 0.141 | Recall 1.000 +/- 0.000

names(epilepsy.fit_fn) <- c("positive","negative")
## (new vars in environment after running): binaryPredictionMetrics, rankedProteins

lastColNum=length(names(rankedProteins))
names(rankedProteins)[c((lastColNum-1):lastColNum)] <- c("Epilepsy","No")

head(rankedProteins[["Epilepsy"]],20)

predFeats.epi<-matrix(NA,nrow=20,ncol=0)
for (gt in c("Epilepsy","No")) if(dim(head(rankedProteins[[gt]],20))[1]==20) {
  this.set=data.frame(head(rankedProteins[[gt]],20))
  colnames(this.set)<-c(gt,paste0(gt,"_importance"))
  predFeats.epi<-cbind(predFeats.epi,this.set)
}

predFeats.epi
#       Epilepsy Epilepsy_importance     No No_importance
#HBA1       HBA1           4.5929784   HBA1     4.1053412
#IGFBP7   IGFBP7           3.8952989   KLK6     3.2597948
#CRH         CRH           3.0479254    DDC     3.2357337
#NEFL       NEFL           3.0331622   FLT1     3.1115648
#FLT1       FLT1           2.7012729 IGFBP7     2.7293426
#ANXA5     ANXA5           2.1313775   CCL4     2.6181717
#DDC         DDC           2.0796786   NEFL     2.4537298
#IL10       IL10           1.7124297   IL10     2.0205561
#VCAM1     VCAM1           1.6837726  ANXA5     1.6069774
#CCL4       CCL4           1.6501709    CRH     1.5085497
#CRP         CRP           1.3544448    MME     1.0655302
#MAPT       MAPT           1.1239118   CST3     1.0509662
#CST3       CST3           0.9071813  NPTX1     0.8762381
#GDF15     GDF15           0.7421629   MAPT     0.7716189
#UCHL1     UCHL1           0.7088304  VCAM1     0.7699722
#MME         MME           0.6836618    TNF     0.6800827
#CHIT1     CHIT1           0.6279624   IL13     0.6631550
#TREM1     TREM1           0.5653328  CCL13     0.5302359
#IL33       IL33           0.5333263 PDLIM5     0.5065959
#ACHE       ACHE           0.4561012    IL9     0.4416678

# same top 12 positive prediction features as:
fixed_epi_features
