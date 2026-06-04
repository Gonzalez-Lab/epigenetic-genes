library(nnet)
library(caret)

#samples.df <- read.csv2("data frame para tfi.csv")
#rownames(samples.df) = samples.df$X
#samples.df = samples.df[,-1]

table(samples.df$age)

head(samples.df)


age = ifelse(blood.covars$age.decade=="20-29",25,
             ifelse(blood.covars$age.decade=="30-39",35,
                    ifelse(blood.covars$age.decade=="40-49",45,
                           ifelse(blood.covars$age.decade=="50-59",55,
                                  ifelse(blood.covars$age.decade=="60-69",65,75)))))

age = ifelse(samples.df$age=="20-29","20s",
             ifelse(samples.df$age=="30-39","30s",
                    ifelse(samples.df$age=="40-49","40s",
                           ifelse(samples.df$age=="50-59","50s",
                                  ifelse(samples.df$age=="60-69","60s","70s")))))

#samples.df$age <- as.factor(age)

#crear particion
set.seed(111); particion=createDataPartition (y=samples.df$age, p=0.8, list=FALSE)

entreno=samples.df[particion, ]

testeo=samples.df[-particion, ]

##############################################
#paso las variables numericas a class numeric
entreno[-which(names(entreno) %in% c("age","sex"))] <- 
  lapply(entreno[-which(names(entreno) %in% c("age","sex"))], as.numeric)

testeo[-which(names(testeo) %in% c("age","sex"))] <- 
  lapply(testeo[-which(names(testeo) %in% c("age","sex"))], as.numeric)
#############################################

set.seed(111); red=nnet(age~., 
                        entreno, 
                        size=35,
                        maxit=20000, 
                        MaxNWts=20000)

pred=predict(red, testeo, type="class")

pred <- factor(pred, levels = levels(testeo$age))

confusionMatrix(pred, testeo$age)


library(MASS)
#este modelo respeta ordinalidad

modelo <- polr(age ~ ., data = entreno, Hess=TRUE) 

pred=predict(modelo, testeo, type="class")

pred <- factor(pred, levels = levels(testeo$age))

confusionMatrix(pred, testeo$age)


library(glmnet)

x_train <- model.matrix(age ~ ., data = entreno)[,-1]
y_train <- entreno$age.decade

cv <- cv.glmnet(x_train, y_train, family="multinomial")

x_test <- model.matrix(age.decade ~ ., data = testeo)[,-1]
y_test <- testeo$age.decade

x_test <- x_test[, colnames(x_train)]

plot(cv)

pred <- predict(cv, 
                newx = x_test,
                s = "lambda.min",
                type = "class")

pred <- as.vector(pred)
pred <- factor(pred, levels = levels(y_train))

confusionMatrix(pred, y_test)


boxplot(as.numeric(samples.df$DNMT1)~samples.df$age.decade)
