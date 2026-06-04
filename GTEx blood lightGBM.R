#library("RNAAgeCalc")

library(lightgbm)
library(caret)
library(Metrics)
library(dplyr)

############################################################
# Edad como valor medio de década
############################################################

age <- ifelse(blood.covars$age.decade == "20-29", 25,
              ifelse(blood.covars$age.decade == "30-39", 35,
                     ifelse(blood.covars$age.decade == "40-49", 45,
                            ifelse(blood.covars$age.decade == "50-59", 55,
                                   ifelse(blood.covars$age.decade == "60-69", 65, 75)))))

y <- age

############################################################
# Datos: todos los genes epigenéticos
############################################################

X <- t(EpiGenes_expr_corrected)

X <- t(log_cpm_values_corrected) #transcriptoma completo

# Verificación básica
stopifnot(nrow(X) == length(y))

cat("Número de muestras:", nrow(X), "\n")
cat("Número de genes epigenéticos:", ncol(X), "\n")

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
# CV regresión
############################################################

for(i in seq_along(folds)){
  
  cat("Fold", i, "\n")
  
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)
  
  ##########################################################
  # Split
  ##########################################################
  
  X_train <- X[train_idx, , drop = FALSE]
  X_test  <- X[test_idx, , drop = FALSE]
  
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  cat("Genes usados:", ncol(X_train), "\n")
  
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
  
  pred <- predict(
    model,
    as.matrix(X_test)
  )
  
  mae_fold <- mae(y_test, pred)
  rmse_fold <- rmse(y_test, pred)
  r2_fold <- cor(y_test, pred)^2
  
  ##########################################################
  # Baseline: edad promedio del training
  ##########################################################
  
  baseline_pred <- rep(mean(y_train), length(y_test))
  baseline_mae <- mae(y_test, baseline_pred)
  
  ##########################################################
  # Guardar resultados
  ##########################################################
  
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
# Resultados regresión
############################################################

results

cat("\n")
cat("MODEL MAE:", mean(results$MAE), "\n")
cat("BASELINE MAE:", mean(results$Baseline_MAE), "\n")
cat("MODEL RMSE:", mean(results$RMSE), "\n")
cat("MODEL R2:", mean(results$R2), "\n")

############################################################
# Mejora porcentual sobre baseline
############################################################

improvement <- 100 * (
  mean(results$Baseline_MAE) - mean(results$MAE)
) / mean(results$Baseline_MAE)

cat("Improvement over baseline:", round(improvement, 2), "%\n")

############################################################
# Tests estadísticos contra baseline
############################################################

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
    Mean_Gain = mean(Gain, na.rm = TRUE),
    SD_Gain = sd(Gain, na.rm = TRUE),
    Mean_Cover = mean(Cover, na.rm = TRUE),
    Mean_Frequency = mean(Frequency, na.rm = TRUE),
    N_Folds = n_distinct(Fold),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Gain))

head(importance_summary, 30)

############################################################
# Genes importantes robustos
############################################################

robust_importance <- importance_summary %>%
  filter(N_Folds == 5) %>%
  arrange(desc(Mean_Gain))

head(robust_importance, 30)

############################################################
# Guardar resultados
############################################################

write.csv(
  results,
  "lightgbm_epigenetic_genes_regression_CV_results.csv",
  row.names = FALSE
)

write.csv(
  importance_summary,
  "lightgbm_epigenetic_genes_feature_importance.csv",
  row.names = FALSE
)

write.csv(
  robust_importance,
  "lightgbm_epigenetic_genes_robust_feature_importance.csv",
  row.names = FALSE
)

############################################################
# Plot importancia
############################################################

lgb.plot.importance(
  importance_all,
  top_n = 30,
  measure = "Gain"
)

