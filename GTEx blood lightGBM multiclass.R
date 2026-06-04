library(lightgbm)
library(caret)
library(dplyr)

############################################################
# Clases
############################################################

y_factor <- factor(
  blood.covars$age.decade,
  levels = c(
    "20-29",
    "30-39",
    "40-49",
    "50-59",
    "60-69",
    "70-79"
  ),
  ordered = TRUE
)

y <- as.numeric(y_factor) - 1

class_names <- levels(y_factor)

############################################################
# Datos
############################################################

X <- t(EpiGenes_expr_corrected)

# X <- t(log_cpm_values_corrected)

############################################################
# Cross validation
############################################################

set.seed(111)

folds <- createFolds(y_factor, k = 5)

results <- data.frame(
  Fold = integer(),
  Accuracy = numeric(),
  Accuracy_pm1 = numeric(),
  Baseline_Accuracy = numeric()
)

importance_list <- list()

############################################################
# CV
############################################################

for(i in seq_along(folds)){
  
  cat("Fold", i, "\n")
  
  test_idx <- folds[[i]]
  
  train_idx <- setdiff(
    seq_len(nrow(X)),
    test_idx
  )
  
  ##########################################################
  # Split
  ##########################################################
  
  X_train <- X[train_idx, ]
  X_test  <- X[test_idx, ]
  
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  ##########################################################
  # Feature selection SOLO TRAIN
  ##########################################################
  
  cors <- apply(
    X_train,
    2,
    function(x)
      cor(
        as.numeric(x),
        y_train,
        method = "spearman"
      )
  )
  
  cors_sorted <- sort(
    abs(cors),
    decreasing = TRUE
  )
  
  n_features <- min(
    240,
    length(cors_sorted)
  )
  
  age.genes <- names(cors_sorted)[1:n_features]
  
  cat(
    "Genes seleccionados:",
    length(age.genes),
    "\n"
  )
  
  X_train <- X_train[, age.genes, drop = FALSE]
  X_test  <- X_test[, age.genes, drop = FALSE]
  
  ##########################################################
  # Dataset
  ##########################################################
  
  dtrain <- lgb.Dataset(
    data = as.matrix(X_train),
    label = y_train
  )
  
  ##########################################################
  # LightGBM multiclass
  ##########################################################
  
  params <- list(
    objective = "multiclass",
    
    num_class = 6,
    
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
  # Feature importance
  ##########################################################
  
  imp <- lgb.importance(model)
  
  imp$Fold <- i
  
  importance_list[[i]] <- imp
  
  ##########################################################
  # Predicción
  ##########################################################
  
  pred_prob <- predict(
    model,
    as.matrix(X_test)
  )
  
  pred_matrix <- matrix(
    pred_prob,
    ncol = 6,
    byrow = TRUE
  )
  
  pred_class <- max.col(
    pred_matrix
  ) - 1
  
  ##########################################################
  # Accuracy exacta
  ##########################################################
  
  acc <- mean(
    pred_class == y_test
  )
  
  ##########################################################
  # Accuracy ±1 década
  ##########################################################
  
  acc_pm1 <- mean(
    abs(pred_class - y_test) <= 1
  )
  
  ##########################################################
  # Baseline
  ##########################################################
  
  majority_class <- as.numeric(
    names(
      which.max(
        table(y_train)
      )
    )
  )
  
  baseline_pred <- rep(
    majority_class,
    length(y_test)
  )
  
  baseline_acc <- mean(
    baseline_pred == y_test
  )
  
  ##########################################################
  # Guardar
  ##########################################################
  
  results <- rbind(
    results,
    data.frame(
      Fold = i,
      Accuracy = acc,
      Accuracy_pm1 = acc_pm1,
      Baseline_Accuracy = baseline_acc
    )
  )
}

############################################################
# Resultados
############################################################

results

cat("\n")

cat(
  "Mean Accuracy:",
  mean(results$Accuracy),
  "\n"
)

cat(
  "Mean Accuracy ±1 decade:",
  mean(results$Accuracy_pm1),
  "\n"
)

cat(
  "Baseline Accuracy:",
  mean(results$Baseline_Accuracy),
  "\n"
)

############################################################
# Test estadístico
############################################################

wilcox.test(
  results$Accuracy,
  results$Baseline_Accuracy,
  paired = TRUE,
  alternative = "greater"
)

t.test(
  results$Accuracy,
  results$Baseline_Accuracy,
  paired = TRUE,
  alternative = "greater"
)

############################################################
# Importancia global
############################################################

importance_all <- bind_rows(
  importance_list
)

importance_summary <- importance_all %>%
  group_by(Feature) %>%
  summarise(
    Mean_Gain = mean(Gain),
    Mean_Cover = mean(Cover),
    Mean_Frequency = mean(Frequency),
    N_Folds = n_distinct(Fold),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Gain))

head(
  importance_summary,
  30
)

write.csv(
  importance_summary,
  "LightGBM_multiclass_importance.csv",
  row.names = FALSE
)