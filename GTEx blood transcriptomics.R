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

ensembl <- substr(rownames(log_cpm_values),1,15)
rownames(log_cpm_values) <- ensembl

# Access sample metadata
gtex_metadata <- colData(gtex_data)

# Filter for whole blood samples
blood_samples <- gtex_metadata[gtex_metadata$study == "BLOOD", ]

# Subset expression data for blood samples
blood_expr <- log_cpm_values[, colnames(log_cpm_values) %in% rownames(blood_samples)]

# View dimensions of filtered data
dim(blood_expr)

table(blood_samples$gtex.sex)
table(blood_samples$gtex.age)

write.csv(blood_expr,"GTEx blood dataset.csv")

library(readxl)
library(org.Hs.eg.db)

EpiGenes <- read_xlsx("epigenetic code genes HG19 TFI.xlsx")

EpiGenes$Symbol <- toupper(EpiGenes$Symbol)


EpiGenes$ensemble <- mapIds(org.Hs.eg.db, keys = EpiGenes$Symbol, 
                            keytype = "SYMBOL", column="ENSEMBL")

EpiGenes.in.blood <- rownames(blood_expr[which(substr(rownames(blood_expr),1,15) %in% EpiGenes$ensemble),])

EpiGenes_expr <- blood_expr[EpiGenes.in.blood,]

EpiGenes_expr <- as.data.frame(EpiGenes_expr)

EpiGenes_expr <- t(EpiGenes_expr)

colnames(EpiGenes_expr) <- substr(colnames(EpiGenes_expr),1,15)

colnames(EpiGenes_expr) <- mapIds(org.Hs.eg.db, keys = colnames(EpiGenes_expr), 
                                  keytype = "ENSEMBL", column="SYMBOL")

#EpiGenes_expr$Symbol <- mapIds(org.Hs.eg.db, keys = substr(rownames(EpiGenes_expr),1,15), 
#                               keytype = "ENSEMBL", column="SYMBOL")

#write.csv2(EpiGenes_expr, "Epigenetic genes blood expression.csv")

#calculo de poblaciones celulares con MCPcounter
library(MCPcounter)

CellEstimates = MCPcounter.estimate(log_cpm_values,
                                    featuresType="ENSEMBL_ID")


#matriz de datos para modelos

blood.df <- cbind(sex = blood_samples$gtex.sex,
                          age.decade = blood_samples$gtex.age,
                          t(CellEstimates),
                          EpiGenes_expr)

samples.df <- as.data.frame(blood.df)

write.csv(samples.df,"data frame para tfi.csv")
