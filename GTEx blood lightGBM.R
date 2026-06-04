#library("RNAAgeCalc")

library(lightgbm)
library(caret)
library(Metrics)
library(dplyr)

############################################################
# Edad como valor medio de década
############################################################

age <- ifelse(blood.covars$age.decade=="20-29",25,
              ifelse(blood.covars$age.decade=="30-39",35,
                     ifelse(blood.covars$age.decade=="40-49",45,
                            ifelse(blood.covars$age.decade=="50-59",55,
                                   ifelse(blood.covars$age.decade=="60-69",65,75)))))

y <- age

############################################################
# Datos
############################################################

X <- t(EpiGenes_expr_corrected)
# X <- t(log_cpm_values_corrected)

############################################################
# Cross Validation
############################################################

set.seed(111)

folds <- createFolds(y, k=5)

results <- data.frame(
  Fold = integer(),
  MAE = numeric(),
  RMSE = numeric(),
  R2 = numeric(),
  Baseline_MAE = numeric()
)

importance_list <- list()

############################################################
# CV regresión
############################################################

for(i in seq_along(folds)){
  
  cat("Fold", i, "\n")
  
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)
  
  X_train <- X[train_idx,]
  X_test  <- X[test_idx,]
  
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  ##########################################################
  # Feature selection SOLO en training
  ##########################################################
  
  cors <- apply(
    X_train,
    2,
    function(x)
      cor(as.numeric(x), y_train, method="spearman")
  )
  
  cors_sorted <- sort(abs(cors), decreasing=TRUE)
  
  n_features <- min(240, length(cors_sorted))
  age.genes <- names(cors_sorted)[1:n_features]
  
  cat("Genes seleccionados:", length(age.genes), "\n")
  
  X_train <- X_train[, age.genes, drop=FALSE]
  X_test  <- X_test[, age.genes, drop=FALSE]
  
  ##########################################################
  # LightGBM regresión
  ##########################################################
  
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
  
  ##########################################################
  # Feature importance por fold
  ##########################################################
  
  imp <- lgb.importance(model)
  imp$Fold <- i
  
  importance_list[[i]] <- imp
  
  ##########################################################
  # Predicción
  ##########################################################
  
  pred <- predict(model, as.matrix(X_test))
  
  mae_fold <- mae(y_test, pred)
  rmse_fold <- rmse(y_test, pred)
  r2_fold <- cor(y_test, pred)^2
  
  ##########################################################
  # Baseline
  ##########################################################
  
  baseline_pred <- rep(mean(y_train), length(y_test))
  baseline_mae <- mae(y_test, baseline_pred)
  
  results <- rbind(
    results,
    data.frame(
      Fold=i,
      MAE=mae_fold,
      RMSE=rmse_fold,
      R2=r2_fold,
      Baseline_MAE=baseline_mae
    )
  )
}

############################################################
# Resultados regresión
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

cat("Improvement over baseline:", round(improvement,2), "%\n")

wilcox.test(
  results$Baseline_MAE,
  results$MAE,
  paired = TRUE,
  alternative = "greater"
)

t.test(
  results$Baseline_MAE,
  results$MAE,
  paired = TRUE,
  alternative = "greater"
)

############################################################
# Feature importance global
############################################################

importance_all <- bind_rows(importance_list)

importance_summary <- importance_all %>%
  group_by(Feature) %>%
  summarise(
    Mean_Gain = mean(Gain, na.rm=TRUE),
    Mean_Cover = mean(Cover, na.rm=TRUE),
    Mean_Frequency = mean(Frequency, na.rm=TRUE),
    N_Folds = n_distinct(Fold),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Gain))

head(importance_summary, 30)

write.csv(
  importance_summary,
  "lightgbm_regression_feature_importance.csv",
  row.names = FALSE
)

lgb.plot.importance(
  importance_all,
  top_n = 30,
  measure = "Gain"
)