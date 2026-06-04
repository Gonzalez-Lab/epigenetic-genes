#library("RNAAgeCalc")

library(lightgbm)
library(caret)
library(Metrics)

############################################################
# Edad
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

# Opción 1
#X <- EpiGenes_expr

# Opción 2
X <- t(log_cpm_values)

############################################################
# Cross Validation
############################################################

set.seed(123)

folds <- createFolds(y, k=5)

results <- data.frame(
  Fold = integer(),
  MAE = numeric(),
  RMSE = numeric(),
  R2 = numeric(),
  Baseline_MAE = numeric()
)

############################################################
# CV
############################################################

for(i in seq_along(folds)){
  
  cat("Fold",i,"\n")
  
  test_idx <- folds[[i]]
  
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)
  
  ##########################################################
  # Split
  ##########################################################
  
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
      cor(as.numeric(x),
          y_train,
          method="spearman")
  )
  
  cors_sorted <- sort(abs(cors),
                      decreasing=TRUE)
  
  age.genes <- names(
    which(cors_sorted > 0.2)
  )
  
  X_train <- X_train[, age.genes, drop=FALSE]
  X_test  <- X_test[, age.genes, drop=FALSE]
  
  ##########################################################
  # LightGBM
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
  # Predicción modelo
  ##########################################################
  
  pred <- predict(
    model,
    as.matrix(X_test)
  )
  
  mae_fold <- mae(y_test,pred)
  
  rmse_fold <- rmse(y_test,pred)
  
  r2_fold <- cor(y_test,pred)^2
  
  ##########################################################
  # Baseline
  ##########################################################
  
  baseline_pred <- rep(
    mean(y_train),
    length(y_test)
  )
  
  baseline_mae <- mae(
    y_test,
    baseline_pred
  )
  
  ##########################################################
  # Guardar
  ##########################################################
  
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
# Resultados
############################################################

results

cat("\n")
cat("MODEL MAE :",mean(results$MAE),"\n")
cat("BASELINE MAE :",mean(results$Baseline_MAE),"\n")

cat("MODEL RMSE :",mean(results$RMSE),"\n")
cat("MODEL R2 :",mean(results$R2),"\n")

############################################################
# Mejora porcentual
############################################################

improvement <- 100 * (
  mean(results$Baseline_MAE) -
    mean(results$MAE)
) / mean(results$Baseline_MAE)

cat(
  "Improvement over baseline:",
  round(improvement,2),
  "%\n"
)

############################################################
# Test estadístico
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
