library(Biobase)
library(GEOquery)
library(dplyr)
library(readxl)
library(org.Hs.eg.db)

#############################################################################
#Deconvolution of bulk blood eQTL effects into immune cell subpopulations. 
#BMC Bioinformatics. 2020 Jun 12;21(1):243.
blood1 <- getGEO("GSE134080")[[1]]

head(exprs(blood1))
dim(exprs(blood1))

head(pData(blood1))
names(pData(blood1))
pData(blood1)$'age:ch1' #edad exacta
pData(blood1)$'Sex:ch1'
table(pData(blood1)$'age:ch1') #100 muestras

df1 <- read.csv("GSE134080_500FG_RNASeq_counts.csv")
rownames(df1) <- df1$X
df1 <- df1[,-1]

covar1 <- data.frame(
  id=colnames(df1),
  age= pData(blood1)$'age:ch1',
  sex= pData(blood1)$'Sex:ch1'
)



#############################################################################

blood2 <- getGEO("GSE190125")[[1]]

head(exprs(blood2))
dim(exprs(blood2))

head(pData(blood2))
names(pData(blood2))
pData(blood2)$'age:ch1' #edad exacta
pData(blood2)$'Sex:ch1'
table(pData(blood2)$'karyotype:ch1') #96 controles
control.idx2 <- which(pData(blood2)$'karyotype:ch1'=="Control")

df2 <- read.delim("GSE190125_raw_counts_GRCh38.p13_NCBI.tsv")
rownames(df2) <- df2$GeneID
df2 <- df2[,-1]
df2 <- df2[,control.idx2]

covar2 <- data.frame(
  id=colnames(df2),
  age= round(as.numeric(pData(blood2)$'age:ch1'[control.idx2])),
  sex= pData(blood2)$'Sex:ch1'[control.idx2]
)

df2$ensembl <- mapIds(org.Hs.eg.db,
                     keys = rownames(df2),
                     column = "ENSEMBL",
                     keytype = "ENTREZID",
                     multiVals = "first")

df2 <- df2[!is.na(df2$ensembl), ]
df2 <- df2[!duplicated(df2$ensembl), ]

rownames(df2) <- df2$ensembl
df2$ensembl <- NULL

#############################################################################

blood3 <- getGEO("GSE183701")[[1]]

head(exprs(blood3))
dim(exprs(blood3))

head(pData(blood3))
names(pData(blood3))
pData(blood3)$'age at draw:ch1' #edad exacta
pData(blood3)$'Sex:ch1'
table(pData(blood3)$'characteristics_ch1') #26 controles

control.idx3 <- which(pData(blood3)$'characteristics_ch1'=="registry: Control")

df3 <- read.csv("GSE183701_combinedRawCounts.csv")
rownames(df3) <- df3$X
df3 <- df3[,-1]
df3 <- df3[,control.idx3]

covar3 <- data.frame(
  id=colnames(df3),
  age= round(as.numeric(pData(blood3)$'age at draw:ch1'[control.idx3])),
  sex= pData(blood3)$'Sex:ch1'[control.idx3]
)


############################################################################
############################################################################
#armado del ds GEO

covars <- data.frame(
  batch=c( rep(1,ncol(df1)),rep(2,ncol(df2)),rep(3,ncol(df3)) ),
  id= c(covar1$id,covar2$id,covar3$id),
  age=as.numeric(c(covar1$age,covar2$age,covar3$age)),
  sex=c(covar1$sex,covar2$sex,covar3$sex)
)

#covars$sex <- tolower(covars$sex)
#hist(covars$age)
#summary(covars$age)
#table(cut(covars$age,breaks = c(0,10,20,30,40,50,60,70,80)))
#boxplot(age ~ batch, data=covars)

#matriz de counts
common_genes <- Reduce(
  intersect,
  list(
    rownames(df1),
    rownames(df2),
    rownames(df3) #,
    #rownames(df4)
  )
)

length(common_genes)

counts <- cbind(
  df1[common_genes, ],
  df2[common_genes, ],
  df3[common_genes, ]#,
  #df4[common_genes, ]
)

##################  MATRIZ DE COUNTS
counts <- as.matrix(counts) #muestras en columnas genes en filas
storage.mode(counts) <- "numeric" # asegurar numérico
#######################################

# Convert raw counts to counts per million (CPM)
library(edgeR)

dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
log_cpm_values <- cpm(dge,log = TRUE, prior.count = 1)

# View first few rows
head(log_cpm_values)

#write.csv(log_cpm_values,"GEO blood dataset.csv")

############################################################################
#calculo de poblaciones celulares con MCPcounter

library(org.Hs.eg.db)
library(AnnotationDbi)
library(MCPcounter)

symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(log_cpm_values),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

log_cpm_symbols <- log_cpm_values[!is.na(symbols), ]
rownames(log_cpm_symbols) <- symbols[!is.na(symbols)]

log_cpm_symbols <- log_cpm_symbols[!duplicated(rownames(log_cpm_symbols)), ]

CellEstimates <- MCPcounter.estimate(
  log_cpm_symbols,
  featuresType = "HUGO_symbols"
)

CellEstimates <- t(CellEstimates)


#############################################################################

#armado del dataset de genes epigeneticos en sangre

EpiGenes <- read.csv2("epigenetic code genes HG19.csv")

EpiGenes_expr <- log_cpm_values[which(rownames(log_cpm_values) %in% EpiGenes$Gene_ID),]

EpiGenes_expr <- as.data.frame(EpiGenes_expr)

EpiGenes_expr <- t(EpiGenes_expr)

#write.csv2(EpiGenes_expr, "Epigenetic genes blood expression.csv")

#matriz de covariables para modelos
blood.covars <- cbind(covars,
                      CellEstimates)

blood.covars <- as.data.frame(blood.covars)

rownames(EpiGenes_expr) == rownames(blood.covars)

#write.csv(blood.covars,"blood metadata and covars.csv")

############################################################################
#Cell composition correction
library(limma)

blood.covars <- blood.covars[,-c(1:3)]

covars$sex <- as.factor(covars$sex)

design <- model.matrix(
  ~ sex + `T cells` +
    `CD8 T cells` +
    `Cytotoxic lymphocytes` +
    `B lineage` +
    `NK cells` +
    `Monocytic lineage` +
    `Myeloid dendritic cells` +
    `Neutrophils` +
    `Endothelial cells` +
    `Fibroblasts`,
  data = blood.covars
)

head(design)

EpiGenes_expr_corrected <- removeBatchEffect(t(EpiGenes_expr), 
                                             covariates = design[,-1])

log_cpm_values_corrected <- removeBatchEffect(log_cpm_values, 
                                              covariates = design[,-1])



"%notin%" <- Negate("%in%")

absent.genes = EpiGenes[which(EpiGenes$Gene_ID %notin% rownames(EpiGenes_expr_corrected)),]$Symbol

length(absent.genes)

absent.genes


###############################################################################

blood4 <- getGEO("GSE75337")[[1]]

head(exprs(blood4))
dim(exprs(blood4))

head(pData(blood4))
names(pData(blood4))
pData(blood4)$'age:ch1' #NO tiene edad exacta
pData(blood4)$'characteristics_ch1.2' #no encuentro el sexo
table(pData(blood4)$'characteristics_ch1')

blood5 <- getGEO("GSE103232")[[1]]

head(exprs(blood5))
dim(exprs(blood5))

head(pData(blood5))
names(pData(blood5))
pData(blood5)$'age:ch1' #NO tiene edad exacta
#pData(blood5)$'Sex:ch1' #no encuentro el sexo
table(pData(blood5)$'characteristics_ch1')