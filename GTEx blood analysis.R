
#PCA analysis

age = ifelse(blood.covars$age.decade=="20-29",25,
             ifelse(blood.covars$age.decade=="30-39",35,
                    ifelse(blood.covars$age.decade=="40-49",45,
                           ifelse(blood.covars$age.decade=="50-59",55,
                                  ifelse(blood.covars$age.decade=="60-69",65,75)))))

table(colMeans(EpiGenes_expr > 1) > 0.2)
keep <- colMeans(EpiGenes_expr > 1) > 0.2
EpiGenes_expr2 <- EpiGenes_expr[, keep]

samples.matrix <- cbind(EpiGenes_expr2,age)

PC <- prcomp(EpiGenes_expr2, scale = TRUE)
PC$rotation
summary(PC)

cor.test(PC$x[,1], age, method="spearman")
cor.test(PC$x[,2], age, method="spearman") #captura senial temporal
cor.test(PC$x[,3], age, method="spearman") #captura senial temporal
cor.test(PC$x[,4], age, method="spearman")

summary(lm(age ~ PC$x[,2] + PC$x[,3]))

library(limma)

covars <- samples.df[,c(1,3:12)]

covars$sex <- as.factor(covars$sex)

cell_cols <- c("T.cells","CD8.T.cells","Cytotoxic.lymphocytes",
               "B.lineage","NK.cells","Monocytic.lineage",
               "Myeloid.dendritic.cells","Neutrophils",
               "Endothelial.cells","Fibroblasts")

covars[cell_cols] <- lapply(covars[cell_cols], function(x)
  as.numeric(as.character(x)))

design <- model.matrix(~ sex + T.cells + 
                         CD8.T.cells +
                         Cytotoxic.lymphocytes +
                         B.lineage + NK.cells +
                         Monocytic.lineage +
                         Myeloid.dendritic.cells +
                         Neutrophils +
                         Endothelial.cells +
                         Fibroblasts,
                       data = covars)

EpiGenes_expr_log.corrected <- removeBatchEffect(EpiGenes_expr_log2, covariates = design[,-1])

PCnew <- prcomp(t(EpiGenes_expr_log.corrected), scale = TRUE)
PCnew$rotation
summary(PCnew)

cor.test(PCnew$x[,1], age.decade, method="spearman")
cor.test(PCnew$x[,2], age.decade, method="spearman") 
cor.test(PCnew$x[,3], age.decade, method="spearman") #captura senial temporal
cor.test(PCnew$x[,4], age.decade, method="spearman")

sort(PCnew$rotation[,3], decreasing=TRUE)[1:25]
sort(PCnew$rotation[,3])[1:25]

cor.test(PCnew$x[,3], covars$Neutrophils)
cor.test(PCnew$x[,3], covars$T.cells)
cor.test(PCnew$x[,3], covars$B.lineage)

plot(age.decade, PCnew$x[,3])
smooth.spline(age.decade, PCnew$x[,3])

summary(lm(PCnew$x[,3] ~ age.decade + sex, data=covars))

plot(age.decade, PCnew$x[,3])
abline(lm(PCnew$x[,3] ~ age.decade))

cor.test(EpiGenes_expr_log.corrected["EZH2",], age.decade, method="spearman")

loadings <- PCnew$rotation[,3]
top_genes <- sort(abs(loadings), decreasing=TRUE)
top_genes[which(top_genes>0.1)]

#Seleccionar por percentil: top 10%:
threshold <- quantile(abs(loadings), 0.9)
selected <- loadings[abs(loadings) >= threshold]

anova(lm(PCnew$x[,3] ~ age.decade))


contrib <- loadings^2
ordered <- sort(contrib, decreasing=TRUE)
cumvar <- cumsum(ordered) / sum(ordered)
selected <- names(ordered[cumvar <= 0.5])
loadings[selected]

summary(lm(PCnew$x[,3] ~ age.decade))

###############################################################################
library(nnet)
library(caret)

#samples.df <- read.csv2("data frame para tfi.csv")
#rownames(samples.df) = samples.df$X
#samples.df = samples.df[,-1]

table(samples.df$age.decade)

head(samples.df)

age.decade = ifelse(samples.df$age.decade=="20-29","20s",
                    ifelse(samples.df$age.decade=="30-39","30s",
                           ifelse(samples.df$age.decade=="40-49","40s",
                                  ifelse(samples.df$age.decade=="50-59","50s",
                                         ifelse(samples.df$age.decade=="60-69","60s","70s")))))

#samples.df$age.decade <- as.factor(age.decade)

#crear particion
set.seed(111); particion=createDataPartition (y=samples.df$age.decade, p=0.8, list=FALSE)

entreno=samples.df[particion, ]

testeo=samples.df[-particion, ]

##############################################
#paso las variables numericas a class numeric
entreno[-which(names(entreno) %in% c("age.decade","sex"))] <- 
  lapply(entreno[-which(names(entreno) %in% c("age.decade","sex"))], as.numeric)

testeo[-which(names(testeo) %in% c("age.decade","sex"))] <- 
  lapply(testeo[-which(names(testeo) %in% c("age.decade","sex"))], as.numeric)
#############################################

set.seed(111); red=nnet(age.decade~., 
                        entreno, 
                        size=35,
                        maxit=20000, 
                        MaxNWts=20000)

pred=predict(red, testeo, type="class")

pred <- factor(pred, levels = levels(testeo$age.decade))

confusionMatrix(pred, testeo$age.decade)


library(MASS)
#este modelo respeta ordinalidad

modelo <- polr(age.decade ~ ., data = entreno, Hess=TRUE) 

pred=predict(modelo, testeo, type="class")

pred <- factor(pred, levels = levels(testeo$age.decade))

confusionMatrix(pred, testeo$age.decade)


library(glmnet)

x_train <- model.matrix(age.decade ~ ., data = entreno)[,-1]
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