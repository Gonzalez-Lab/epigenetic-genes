# Acces sample data 
blood_data <- read.csv("expression Allen.csv")

# Access sample metadata
blood_metadata <- read.csv("metadata Allen.csv")

table(blood_metadata$subject.biologicalSex)
table(blood_metadata$subject.ageAtFirstDraw)

# 1. Leer sin modificar nombres
blood_data <- read.csv("expression Allen.csv", check.names = FALSE)

# 2. Ver estructura
dim(blood_data)
head(blood_data[,1:5])
str(blood_data[,1:5])

# 3. Armar matriz de counts
counts <- blood_data[,-1]
rownames(counts) <- blood_data[,1]

# 4. Forzar todo a numérico
counts <- as.data.frame(lapply(counts, function(x) as.numeric(as.character(x))))
rownames(counts) <- blood_data[,1]

# 5. Detectar columnas problema
lib_sizes <- colSums(counts, na.rm = TRUE)

summary(lib_sizes)
sum(lib_sizes == 0)
which(lib_sizes == 0)[1:20]
colnames(counts)[lib_sizes == 0][1:20]

counts <- counts[, lib_sizes > 0]

# eliminar genes sin counts en todas las muestras
counts <- counts[rowSums(counts, na.rm = TRUE) > 0, ]
counts <- t(counts)

library(edgeR)

dge <- DGEList(counts = counts)
dge <- normLibSizes(dge)

log_cpm_values <- cpm(dge, log = TRUE, prior.count = 1)


# View first few rows
head(log_cpm_values)

############################################################################
#calculo de poblaciones celulares con MCPcounter
library(MCPcounter)

CellEstimates = MCPcounter.estimate(log_cpm_values,
                                    featuresType="HUGO_symbols")

CellEstimates <- t(CellEstimates)

#############################################################################
#armado del dataset de genes epigeneticos en sangre

library(readxl)

EpiGenes <- read_xlsx("epigenetic code genes HG19.xlsx")

EpiGenes_expr <- log_cpm_values[which(rownames(log_cpm_values) %in% EpiGenes$Symbol),]

EpiGenes_expr <- as.data.frame(EpiGenes_expr)

EpiGenes_expr <- t(EpiGenes_expr)

#write.csv2(EpiGenes_expr, "Epigenetic genes blood expression.csv")

#matriz de covariables para modelos
blood.covars <- cbind(age = blood_metadata$subject.ageAtFirstDraw,
                      sex = blood_metadata$subject.biologicalSex,
                      CellEstimates)

blood.covars <- as.data.frame(blood.covars)

rownames(EpiGenes_expr) == rownames(blood.covars)

#write.csv(blood.covars,"blood metadata and covars.csv")

############################################################################
#Cell composition correction
library(limma)

covars <- blood.covars[, -1]

covars$sex <- as.factor(covars$sex)

cell_cols <- c(
  "T cells",
  "CD8 T cells",
  "Cytotoxic lymphocytes",
  "B lineage",
  "NK cells",
  "Monocytic lineage",
  "Myeloid dendritic cells",
  "Neutrophils",
  "Endothelial cells",
  "Fibroblasts"
)

covars[cell_cols] <- lapply(covars[cell_cols], function(x) {
  as.numeric(as.character(x))
})

design <- model.matrix(
  ~ sex +
    `T cells` +
    `CD8 T cells` +
    `Cytotoxic lymphocytes` +
    `B lineage` +
    `NK cells` +
    `Monocytic lineage` +
    `Myeloid dendritic cells` +
    Neutrophils +
    `Endothelial cells` +
    Fibroblasts,
  data = covars
)

head(design)

EpiGenes_expr_corrected <- removeBatchEffect(t(EpiGenes_expr), 
                                             covariates = design[,-1])

log_cpm_values_corrected <- removeBatchEffect(log_cpm_values, 
                                              covariates = design[,-1])

# eliminar rownames NA
log_cpm_values_corrected <- log_cpm_values_corrected[
  !is.na(rownames(log_cpm_values_corrected)),
  ,
  drop = FALSE
]

# verificar duplicados
sum(duplicated(rownames(log_cpm_values_corrected)))


#PCA analysis

PC <- prcomp(t(EpiGenes_expr_corrected), scale = TRUE)
#PC$rotation
summary(PC)

age <- as.numeric(blood.covars$age)

cor.test(PC$x[,1], age, method="spearman") #captura señal temporal
cor.test(PC$x[,2], age, method="spearman") #captura señal temporal
cor.test(PC$x[,3], age, method="spearman") 
cor.test(PC$x[,4], age, method="spearman")

sort(PC$rotation[,1], decreasing=TRUE)[1:30]
sort(PC$rotation[,1])[1:25]

sort(PC$rotation[,2], decreasing=TRUE)[1:30]
sort(PC$rotation[,2])[1:25]

plot(age, PC$x[,2])
smooth.spline(age, PC$x[,2])

summary(lm(PC$x[,2] ~ age + sex, data=covars))

plot(age, PC$x[,3])
abline(lm(PC$x[,3] ~ age))

summary(lm(PC$x[,3] ~ age + sex, data=covars))

#para ver correlacion de genes individuales
#cor.test(EpiGenes_expr_corrected["EZH2",], age, method="spearman")

loadings <- PC$rotation[,2]

#Seleccionar por percentil: top 10%:
threshold <- quantile(abs(loadings), 0.9)
selected <- loadings[abs(loadings) >= threshold]

#analisis estadistico
anova(lm(PC$x[,2] ~ age))

#los que explican el 50% del peso en el PC
contrib <- loadings^2
ordered <- sort(contrib, decreasing=TRUE)
cumvar <- cumsum(ordered) / sum(ordered)
selected <- names(ordered[cumvar <= 0.5])
loadings[selected]