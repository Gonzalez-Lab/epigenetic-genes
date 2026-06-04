library(recount3)

# List available GTEx projects
proj_info <- available_projects()
print(proj_info)

head(subset(proj_info, file_source == "gtex")) 
head(subset(proj_info, file_source == "tcga"))

# Load GTEx project (v8 human RNA-seq)
gtex_data <- create_rse_manual(
  project = "BLOOD",
  project_home = "data_sources/gtex",
  organism = "human"
)

# View the dataset
gtex_data

# Convert raw counts to counts per million (CPM)
library(edgeR)

dge <- DGEList(counts = assays(gtex_data)$raw_counts)
dge <- calcNormFactors(dge)
log_cpm_values <- cpm(dge,log = TRUE, prior.count = 1)

# View first few rows
head(log_cpm_values)

#write.csv(log_cpm_values,"GTEx blood dataset.csv")

rownames(log_cpm_values) <- substr(rownames(log_cpm_values),1,15)

# Access sample metadata
blood_metadata <- colData(gtex_data)

table(blood_metadata$gtex.sex)
table(blood_metadata$gtex.age)

############################################################################
#calculo de poblaciones celulares con MCPcounter
library(MCPcounter)

CellEstimates = MCPcounter.estimate(log_cpm_values,
                                    featuresType="ENSEMBL_ID")

CellEstimates <- t(CellEstimates)

#############################################################################

#armado del dataset de genes epigeneticos en sangre

library(readxl)
library(org.Hs.eg.db)

EpiGenes <- read_xlsx("epigenetic code genes HG19.xlsx")

EpiGenes$ensemble <- mapIds(org.Hs.eg.db, keys = EpiGenes$Symbol, 
                            keytype = "SYMBOL", column="ENSEMBL")

EpiGenes_expr <- log_cpm_values[which(rownames(log_cpm_values) %in% EpiGenes$ensemble),]

EpiGenes_expr <- as.data.frame(EpiGenes_expr)

EpiGenes_expr <- t(EpiGenes_expr)

colnames(EpiGenes_expr) <- mapIds(org.Hs.eg.db, keys = colnames(EpiGenes_expr), 
                                  keytype = "ENSEMBL", column="SYMBOL")

#write.csv2(EpiGenes_expr, "Epigenetic genes blood expression.csv")

#matriz de covariables para modelos
blood.covars <- cbind(age.decade = blood_metadata$gtex.age,
                      sex = blood_metadata$gtex.sex,
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

#PCA analysis

PC <- prcomp(t(EpiGenes_expr_corrected), scale = TRUE)
#PC$rotation
summary(PC)

cor.test(PC$x[,1], age, method="spearman")
cor.test(PC$x[,2], age, method="spearman") #captura señal temporal
cor.test(PC$x[,3], age, method="spearman") #captura señal temporal
cor.test(PC$x[,4], age, method="spearman")

sort(PC$rotation[,2], decreasing=TRUE)[1:30]
sort(PC$rotation[,2])[1:25]

sort(PC$rotation[,3], decreasing=TRUE)[1:30]
sort(PC$rotation[,3])[1:25]

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