#library("RNAAgeCalc")

#### Upload epigenetic genes data
samples.df <- read.csv2("data frame para tfi.csv")

rownames(samples.df) = samples.df$X
samples.df = samples.df[,-1]

table(samples.df$age.decade)

head(samples.df)

age = ifelse(samples.df$age.decade=="20-29",25,
                    ifelse(samples.df$age.decade=="30-39",35,
                           ifelse(samples.df$age.decade=="40-49",45,
                                  ifelse(samples.df$age.decade=="50-59",55,
                                         ifelse(samples.df$age.decade=="60-69",65,75)))))

age <- as.numeric(age)

#######################################################################
library(lightgbm)
library(caret)
library(Metrics)

set.seed(123)

#-----------------------------------
# Datos
#-----------------------------------

y <- age

#epigenetic genes
X <- samples.df[,-c(1:12)]

#all genes
X <- t(blood_expr)
colnames(X) <- make.unique(colnames(X))

X <- as.matrix(X)

cors <- apply(X, 2, function(x)
  cor(as.numeric(x), y, method="spearman"))

summary(abs(cors))

cors_sorted <- sort(abs(cors), decreasing = TRUE)

age.genes <- names(which(cors_sorted>0.15))

X <- X[,age.genes]

#-----------------------------------
# Folds
#-----------------------------------

folds <- createFolds(y, k = 5)

results <- data.frame(
  Fold = integer(),
  MAE = numeric(),
  RMSE = numeric(),
  R2 = numeric()
)

#-----------------------------------
# Cross-validation
#-----------------------------------

for(i in seq_along(folds)){
  
  test_idx <- folds[[i]]
  
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)
  
  X_train <- X[train_idx, ]
  X_test  <- X[test_idx, ]
  
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  dtrain <- lgb.Dataset(
    data = X_train,
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
  
  pred <- predict(model, X_test)
  
  mae_fold <- mae(y_test, pred)
  
  rmse_fold <- rmse(y_test, pred)
  
  r2_fold <- cor(y_test, pred)^2
  
  results <- rbind(
    results,
    data.frame(
      Fold = i,
      MAE = mae_fold,
      RMSE = rmse_fold,
      R2 = r2_fold
    )
  )
}

results

cat("Mean MAE :", mean(results$MAE), "\n")
cat("Mean RMSE:", mean(results$RMSE), "\n")
cat("Mean R2  :", mean(results$R2), "\n")


#####################################################################
