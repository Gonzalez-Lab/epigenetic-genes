library(limma)
library(lightgbm)
library(caret)
library(Metrics)
library(dplyr)

############################################################
# Datos
############################################################

y <- covars$age

X <- t(EpiGenes_expr_corrected)
X <- t(log_cpm_values_corrected)

stopifnot(nrow(X) == length(y))

cat("Número de muestras:", nrow(X), "\n")
cat("Número de genes:", ncol(X), "\n")

############################################################
# Función: corregir batch dentro de cada fold
# Input: matrices muestras x genes
############################################################

correct_batch_train_test <- function(X_train, X_test, covars_train, covars_test) {
  
  expr_train <- t(X_train)  # genes x muestras
  expr_test  <- t(X_test)
  
  covars_train$batch <- factor(covars_train$batch)
  covars_test$batch  <- factor(covars_test$batch, levels = levels(covars_train$batch))
  
  design_train <- model.matrix(
    ~ age + batch,
    data = covars_train
  )
  
  fit <- lmFit(expr_train, design_train)
  
  batch_cols <- grep("^batch", colnames(design_train))
  
  batch_effect_train <- fit$coefficients[, batch_cols, drop = FALSE] %*%
    t(design_train[, batch_cols, drop = FALSE])
  
  expr_train_corrected <- expr_train - batch_effect_train
  
  design_test <- model.matrix(
    ~ age + batch,
    data = covars_test
  )
  
  design_test <- design_test[, colnames(design_train), drop = FALSE]
  
  batch_effect_test <- fit$coefficients[, batch_cols, drop = FALSE] %*%
    t(design_test[, batch_cols, drop = FALSE])
  
  expr_test_corrected <- expr_test - batch_effect_test
  
  list(
    train = t(expr_train_corrected), # vuelve a muestras x genes
    test = t(expr_test_corrected)
  )
}

############################################################
# Cross Validation
############################################################

set.seed(111)

folds <- createFolds(y, k = 5)

results <- data.frame(
  Fold = integer(),
  MAE = numeric(),
  RMSE = numeric(),
  R2 = numeric(),
  Baseline_MAE = numeric()
)

importance_list <- list()

############################################################
# CV
############################################################

for(i in seq_along(folds)){
  
  cat("Fold", i, "\n")
  
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)
  
  X_train_raw <- X[train_idx, , drop = FALSE]
  X_test_raw  <- X[test_idx, , drop = FALSE]
  
  covars_train <- covars[train_idx, ]
  covars_test  <- covars[test_idx, ]
  
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  corrected <- correct_batch_train_test(
    X_train = X_train_raw,
    X_test = X_test_raw,
    covars_train = covars_train,
    covars_test = covars_test
  )
  
  X_train <- corrected$train
  X_test  <- corrected$test
  
  cat("Genes usados:", ncol(X_train), "\n")
  
  dtrain <- lgb.Dataset(
    data = as.matrix(X_train),
    label = y_train
  )
  
  params <- list(
    objective = "regression",
    metric = "rmse",
    learning_rate = 0.01,
    num_leaves = 31,
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 5,
    min_data_in_leaf = 10,
    verbosity = -1
  )
  
  model <- lgb.train(
    params = params,
    data = dtrain,
    nrounds = 1000
  )
  
  imp <- lgb.importance(model)
  imp$Fold <- i
  importance_list[[i]] <- imp
  
  pred <- predict(model, as.matrix(X_test))
  
  mae_fold <- mae(y_test, pred)
  rmse_fold <- rmse(y_test, pred)
  r2_fold <- cor(y_test, pred)^2
  
  baseline_pred <- rep(mean(y_train), length(y_test))
  baseline_mae <- mae(y_test, baseline_pred)
  
  results <- rbind(
    results,
    data.frame(
      Fold = i,
      MAE = mae_fold,
      RMSE = rmse_fold,
      R2 = r2_fold,
      Baseline_MAE = baseline_mae
    )
  )
}

############################################################
# Resultados
############################################################

results

cat("\n")
cat("MODEL MAE:", mean(results$MAE), "\n")
cat("BASELINE MAE:", mean(results$Baseline_MAE), "\n")
cat("MODEL RMSE:", mean(results$RMSE), "\n")
cat("MODEL R2:", mean(results$R2), "\n")

improvement <- 100 * (
  mean(results$Baseline_MAE) - mean(results$MAE)
) / mean(results$Baseline_MAE)

cat("Improvement over baseline:", round(improvement, 2), "%\n")

wilcox.test(results$Baseline_MAE, results$MAE, paired = TRUE, alternative = "greater")
t.test(results$Baseline_MAE, results$MAE, paired = TRUE, alternative = "greater")

############################################################
# Feature importance
############################################################

importance_all <- bind_rows(importance_list)

importance_summary <- importance_all %>%
  group_by(Feature) %>%
  summarise(
    Mean_Gain = mean(Gain, na.rm = TRUE),
    SD_Gain = sd(Gain, na.rm = TRUE),
    Mean_Cover = mean(Cover, na.rm = TRUE),
    Mean_Frequency = mean(Frequency, na.rm = TRUE),
    N_Folds = n_distinct(Fold),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Gain))

head(importance_summary, 30)

robust_importance <- importance_summary %>%
  filter(N_Folds == 5) %>%
  arrange(desc(Mean_Gain))

head(robust_importance, 30)

plot(y_test, pred)
abline(0,1,col="red")
