#iSee
#load libraries
library(tidyverse)
library(DESeq2)
library(edgeR)
library(limma)
library(biomaRt)
library(ggfortify)
library(pheatmap)
library(cowplot)
library(DT)
library(gt)
library(plotly)
library(msigdbr)
library(clusterProfiler)
library(ggrepel)
library(topGO)
library("org.Hs.eg.db")
library(lattice)
library(magrittr)
library(iSEE)
library(ComplexHeatmap)
library(enrichplot)
library(ggplot2)
library(ggridges)

#import raw counts csv
Exp1_countsimport<-read.csv(file = "rawcounts2021.csv", header = TRUE)
Exp2_countsimport<-read.csv(file = "rawcounts2022.csv", header = TRUE)
Hep1_countsimport<-read.csv(file = "HepRawCounts_2021.csv", header = TRUE)
Hep2_countsimport<-read.csv(file = "HEPRawCounts_2022.csv", header = TRUE)
HEK1_countsimport<-read.csv(file = "HEKRawCounts_2021.csv", header = TRUE)
HEK2_countsimport<-read.csv(file = "HEKRawCounts_2022.csv", header = TRUE)
HeLa1_countsimport<- read.csv(file = "HELArawcounts2021.csv", header = TRUE)
HeLa2_countsimport<- read.csv(file = "HeLaRawCounts_2022.csv", header = TRUE)

#Merge Raw Counts from Batch1 and Batch2 and assign ENSGIDs to rownames
RawReadCountsE1E2<-inner_join(Exp2_countsimport,Exp1_countsimport, "gene_id") %>% column_to_rownames(var='gene_id') 
RawReadCountsE1<-Exp1_countsimport%>% column_to_rownames(var='gene_id')
RawReadCountsE2<-Exp2_countsimport%>% column_to_rownames(var='gene_id')
Hep1cts<-Hep1_countsimport%>% column_to_rownames("gene_id")
Hep2cts<-Hep2_countsimport%>% column_to_rownames("gene_id")
HEK1cts<-HEK1_countsimport%>% column_to_rownames("gene_id")
HEK2cts<-HEK2_countsimport%>% column_to_rownames("gene_id")
HeLa1cts<-HeLa1_countsimport%>% column_to_rownames("gene_id")
HeLa2cts<-HeLa2_countsimport%>% column_to_rownames("gene_id")

#Prepare cts and coldata objects
coldata<-read.csv(file = 'B1B2coldata.csv',header = TRUE) %>% column_to_rownames(var='X')
coldataE1<-coldata %>%filter(Batch=='E1')
coldataE2<-coldata %>%filter(Batch=='E2')
Hep1coldata<-read.csv("Hepcoldata_2021.csv", header  = TRUE) %>% column_to_rownames(var='X')
Hep2coldata<-read.csv("HEP2022_coldata.csv", header = TRUE) %>% column_to_rownames(var='X')
HEK1coldata<-read.csv("HEK2021_coldata.csv", header = TRUE) %>% column_to_rownames(var='X')
HEK2coldata<-read.csv("HEK2022_coldata.csv", header = TRUE) %>% column_to_rownames(var='X')
HeLa1coldata<-read.csv("HELAcoldata2021.csv", header = TRUE) %>% column_to_rownames(var='X')
HeLa2coldata<-read.csv("HeLa2022_coldata.csv", header = TRUE) %>% column_to_rownames(var='X')

#verify
all(rownames(coldata) %in% colnames(RawReadCountsE1E2))
all(rownames(coldataE1) %in% colnames(RawReadCountsE1))
all(rownames(coldata$Batch== 'E2') %in% colnames(RawReadCountsE2))
all(rownames(Hep1coldata) %in% colnames(Hep1cts))
all(rownames(Hep2coldata) %in% colnames(Hep2cts))
all(rownames(HEK1coldata) %in% colnames(HEK1cts))
all(rownames(HEK2coldata) %in% colnames(HEK2cts))
all(rownames(HeLa1coldata) %in% colnames(HeLa1cts))
all(rownames(HeLa2coldata) %in% colnames(HeLa2cts))


all(rownames(coldata) == colnames(RawReadCountsE1E2))
all(rownames(coldataE1) == colnames(RawReadCountsE1))
all(rownames(coldata$Batch== 'E2') == colnames(RawReadCountsE2))
all(rownames(Hep1coldata) == colnames(Hep1cts))
all(rownames(Hep2coldata) == colnames(Hep2cts))
all(rownames(HEK1coldata) == colnames(HEK1cts))
all(rownames(HEK2coldata) == colnames(HEK2cts))
all(rownames(HeLa1coldata) == colnames(HeLa1cts))
all(rownames(HeLa2coldata) == colnames(HeLa2cts))


RawReadCountsE1E2 <- RawReadCountsE1E2[, rownames(coldata)]
RawReadCountsE1 <- RawReadCountsE1[,rownames(coldataE1)]
RawReadCountsE2 <- RawReadCountsE2[,rownames(coldataE2)]

all(rownames(coldata) == colnames(RawReadCountsE1E2))
all(rownames(coldataE1) == colnames(RawReadCountsE1))
all(rownames(coldata$Batch== 'E2') == colnames(RawReadCountsE2))

colannotationsdf <- as.data.frame(coldata$Condition)
colnames(colannotationsdf) <- c("Condition")
row.names(colannotationsdf) <-colnames(RawReadCountsE1E2)

#Make  local "mart"
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

RawReadCountsE1E2_ENSG<-row.names(RawReadCountsE1E2)
RawReadCountsE1_ENSG<-row.names(RawReadCountsE1)
RawReadCountsE2_ENSG<-row.names(RawReadCountsE2)
Hep1cts_ENSG<-row.names(Hep1cts)
Hep2cts_ENSG<-row.names(Hep2cts)
HEK1cts_ENSG<-row.names(HEK1cts)
HEK2cts_ENSG<-row.names(HEK2cts)
HeLa1cts_ENSG<-row.names(HeLa1cts)
HeLa2cts_ENSG<-row.names(HeLa2cts)

#Make an ENSG lookup table to convert known ENSG IDs to HGNC IDs
ENSGlookup <- getBM(
  mart = mart,
  attributes = c('ensembl_gene_id',
                 'gene_biotype','hgnc_symbol'),
  filter = 'ensembl_gene_id',
  values = RawReadCountsE1E2_ENSG,
  uniqueRows = TRUE)

#prep a column to recieve ENSG IDs 
RawReadCountsE1E2_ENSG_colprep<-rownames_to_column(RawReadCountsE1E2, "ensembl_gene_id")
RawReadCountsE1_ENSG_colprep<-rownames_to_column(RawReadCountsE1, "ensembl_gene_id")
RawReadCountsE2_ENSG_colprep<-rownames_to_column(RawReadCountsE2, "ensembl_gene_id")
Hep1cts_ENSG_colprep<-rownames_to_column(Hep1cts, "ensembl_gene_id")
Hep2cts_ENSG_colprep<-rownames_to_column(Hep2cts, "ensembl_gene_id")
HEK1cts_ENSG_colprep<-rownames_to_column(HEK1cts, "ensembl_gene_id")
HEK2cts_ENSG_colprep<-rownames_to_column(HEK2cts, "ensembl_gene_id")
HeLa1cts_ENSG_colprep<-rownames_to_column(HeLa1cts, "ensembl_gene_id")
HeLa2cts_ENSG_colprep<-rownames_to_column(HeLa2cts, "ensembl_gene_id")

#inner_join with lookup table w/HGNC + biotypes
RawReadCountsE1E2_ENSG_HGNC_BT<-inner_join(RawReadCountsE1E2_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
RawReadCountsE1_ENSG_HGNC_BT<-inner_join(RawReadCountsE1_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
RawReadCountsE2_ENSG_HGNC_BT<-inner_join(RawReadCountsE2_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
Hep1cts_ENSG_HGNC_BT<-inner_join(Hep1cts_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
Hep2cts_ENSG_HGNC_BT<-inner_join(Hep2cts_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
HEK1cts_ENSG_HGNC_BT<-inner_join(HEK1cts_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
HEK2cts_ENSG_HGNC_BT<-inner_join(HEK2cts_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
HeLa1cts_ENSG_HGNC_BT<-inner_join(HeLa1cts_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)
HeLa2cts_ENSG_HGNC_BT<-inner_join(HeLa2cts_ENSG_colprep,ENSGlookup,by ="ensembl_gene_id", keep=FALSE)

#remove duplicate HGNC IDs and name rows by HGNC
RawReadCountsE1E2_ENSG_HGNC_BT<-RawReadCountsE1E2_ENSG_HGNC_BT[!duplicated(RawReadCountsE1E2_ENSG_HGNC_BT$hgnc_symbol), ]
RawReadCountsE1_ENSG_HGNC_BT<-RawReadCountsE1_ENSG_HGNC_BT[!duplicated(RawReadCountsE1_ENSG_HGNC_BT$hgnc_symbol), ]
RawReadCountsE2_ENSG_HGNC_BT<-RawReadCountsE2_ENSG_HGNC_BT[!duplicated(RawReadCountsE2_ENSG_HGNC_BT$hgnc_symbol), ]
Hep1cts_ENSG_HGNC_BT<-Hep1cts_ENSG_HGNC_BT[!duplicated(Hep1cts_ENSG_HGNC_BT$hgnc_symbol),]
Hep2cts_ENSG_HGNC_BT<-Hep2cts_ENSG_HGNC_BT[!duplicated(Hep2cts_ENSG_HGNC_BT$hgnc_symbol),]
HEK1cts_ENSG_HGNC_BT<-HEK1cts_ENSG_HGNC_BT[!duplicated(HEK1cts_ENSG_HGNC_BT$hgnc_symbol),]
HEK2cts_ENSG_HGNC_BT<-HEK2cts_ENSG_HGNC_BT[!duplicated(HEK2cts_ENSG_HGNC_BT$hgnc_symbol),]
HeLa1cts_ENSG_HGNC_BT<-HeLa1cts_ENSG_HGNC_BT[!duplicated(HeLa1cts_ENSG_HGNC_BT$hgnc_symbol),]
HeLa2cts_ENSG_HGNC_BT<-HeLa2cts_ENSG_HGNC_BT[!duplicated(HeLa2cts_ENSG_HGNC_BT$hgnc_symbol),]

row.names(RawReadCountsE1E2_ENSG_HGNC_BT)<-RawReadCountsE1E2_ENSG_HGNC_BT$hgnc_symbol
row.names(RawReadCountsE1_ENSG_HGNC_BT)<-RawReadCountsE1_ENSG_HGNC_BT$hgnc_symbol
row.names(RawReadCountsE2_ENSG_HGNC_BT)<-RawReadCountsE2_ENSG_HGNC_BT$hgnc_symbol
row.names(Hep1cts_ENSG_HGNC_BT)<-Hep1cts_ENSG_HGNC_BT$hgnc_symbol
row.names(Hep2cts_ENSG_HGNC_BT)<-Hep2cts_ENSG_HGNC_BT$hgnc_symbol
row.names(HEK1cts_ENSG_HGNC_BT)<-HEK1cts_ENSG_HGNC_BT$hgnc_symbol
row.names(HEK2cts_ENSG_HGNC_BT)<-HEK2cts_ENSG_HGNC_BT$hgnc_symbol
row.names(HeLa1cts_ENSG_HGNC_BT)<-HeLa1cts_ENSG_HGNC_BT$hgnc_symbol
row.names(HeLa2cts_ENSG_HGNC_BT)<-HeLa2cts_ENSG_HGNC_BT$hgnc_symbol

RawReadCountsE1E2_ENSG_HGNC_BT <- RawReadCountsE1E2_ENSG_HGNC_BT[,-1]
RawReadCountsE1E2_ENSG_HGNC_BT <- RawReadCountsE1E2_ENSG_HGNC_BT[,-c(80:81)]
RawReadCountsE1_ENSG_HGNC_BT <- RawReadCountsE1_ENSG_HGNC_BT[,-1]
RawReadCountsE1_ENSG_HGNC_BT <- RawReadCountsE1_ENSG_HGNC_BT[,-c(28:29)]
RawReadCountsE2_ENSG_HGNC_BT <- RawReadCountsE2_ENSG_HGNC_BT[,-1]
RawReadCountsE2_ENSG_HGNC_BT <- RawReadCountsE2_ENSG_HGNC_BT[,-c(53:54)]
Hep1cts_ENSG_HGNC_BT <- Hep1cts_ENSG_HGNC_BT[,-1]
Hep1cts_ENSG_HGNC_BT <- Hep1cts_ENSG_HGNC_BT[,-c(10,11)]
Hep2cts_ENSG_HGNC_BT <- Hep2cts_ENSG_HGNC_BT[,-1]
Hep2cts_ENSG_HGNC_BT <- Hep2cts_ENSG_HGNC_BT[,-c(19,20)]
HEK1cts_ENSG_HGNC_BT <- HEK1cts_ENSG_HGNC_BT[,-1]
HEK1cts_ENSG_HGNC_BT <- HEK1cts_ENSG_HGNC_BT[,-c(10,11)]
HEK2cts_ENSG_HGNC_BT <- HEK2cts_ENSG_HGNC_BT[,-1]
HEK2cts_ENSG_HGNC_BT <- HEK2cts_ENSG_HGNC_BT[,-c(19,20)]
HeLa1cts_ENSG_HGNC_BT <- HeLa1cts_ENSG_HGNC_BT[,-1]
HeLa1cts_ENSG_HGNC_BT <- HeLa1cts_ENSG_HGNC_BT[,-c(10,11)]
HeLa2cts_ENSG_HGNC_BT <- HeLa2cts_ENSG_HGNC_BT[,-1]
HeLa2cts_ENSG_HGNC_BT <- HeLa2cts_ENSG_HGNC_BT[,-c(17,18)]


#create S4 Object with DESeq2
# sce_E1E2 <- DESeqDataSetFromMatrix(countData = RawReadCountsE1E2_ENSG_HGNC_BT,
#                                            colData = coldata,
#                                            design = ~ Condition)
# sce_E1 <- DESeqDataSetFromMatrix(countData = RawReadCountsE1_ENSG_HGNC_BT,
#                                   colData = coldataE1,
#                                   design = ~ Condition)
# sce_E2 <- DESeqDataSetFromMatrix(countData = RawReadCountsE2_ENSG_HGNC_BT,
#                                  colData = coldataE2,
#                                  design = ~ Condition)
# sce_Hep1 <- DESeqDataSetFromMatrix(countData = Hep1cts_ENSG_HGNC_BT,
#                                    colData = Hep1coldata,
#                                    design = ~ Condition_NBE)
# sce_Hep2 <- DESeqDataSetFromMatrix(countData = Hep2cts_ENSG_HGNC_BT,
#                        colData = Hep2coldata,
#                        design = ~ Condition_NBE)
# sce_HEK1 <- DESeqDataSetFromMatrix(countData = HEK1cts_ENSG_HGNC_BT,
#                                    colData = HEK1coldata,
#                                    design = ~ Condition_NBE)
# sce_HEK2 <- DESeqDataSetFromMatrix(countData = HEK2cts_ENSG_HGNC_BT,
#                                    colData = HEK2coldata,
#                                    design = ~ Condition_NBE)
# sce_HeLa1 <- DESeqDataSetFromMatrix(countData = HeLa1cts_ENSG_HGNC_BT,
#                                    colData = HeLa1coldata,
#                                    design = ~ Condition_NBE)
# sce_HeLa2 <- DESeqDataSetFromMatrix(countData = HeLa2cts_ENSG_HGNC_BT,
#                                    colData = HeLa2coldata,
#                                    design = ~ Condition_NBE)
# #iSEE Raw Counts
# iSEE(sce_E1E2)
# iSEE(sce_E1)
# iSEE(sce_E2)
# iSEE(sce_Hep1)
# iSEE(sce_Hep2)
# iSEE(sce_HEK1)
# iSEE(sce_HEK2)
# iSEE(sce_HeLa1)
# iSEE(sce_HeLa2)
# 
# #Normalization
# E1E2_norm <- estimateSizeFactors(sce_E1E2)
# normalized_counts_E1E2 <- counts(sce_E1E2_norm, normalized = TRUE)
# sce_E1E2 <- DESeqDataSetFromMatrix(countData = round(normalized_counts_E1E2),
#                                    colData = coldata,
#                                    design = ~ Condition)
# 
# E1_norm <- estimateSizeFactors(sce_E1)
# normalized_counts_E1 <- counts(E1_norm, normalized = TRUE)
# sce_E1 <- DESeqDataSetFromMatrix(countData = round(normalized_counts_E1),
#                                    colData = coldataE1,
#                                    design = ~ Condition)
# 
# E2_norm <- estimateSizeFactors(sce_E2)
# normalized_counts_E2 <- counts(E2_norm, normalized = TRUE)
# sce_E2 <- DESeqDataSetFromMatrix(countData = round(normalized_counts_E2),
#                                  colData = coldataE2,
#                                  design = ~ Condition)
# 
# iSEE(sce_E1E2)
# iSEE(sce_E1)
# iSEE(sce_E2)
#Differential Expression
ddsE1E2 <- DESeqDataSetFromMatrix(countData = RawReadCountsE1E2_ENSG_HGNC_BT,
                                  colData = coldata,
                                  design = ~ Condition)

ddsE1 <- DESeqDataSetFromMatrix(countData = RawReadCountsE1_ENSG_HGNC_BT,
                                colData = coldataE1 ,
                                design = ~ Condition)

ddsE2 <- DESeqDataSetFromMatrix(countData = RawReadCountsE2_ENSG_HGNC_BT,
                                colData = coldataE2 ,
                                design = ~ Condition)
ddsHep1 <- DESeqDataSetFromMatrix(countData = Hep1cts_ENSG_HGNC_BT,
                                  colData = Hep1coldata,
                                  design = ~ Condition_NBE)
ddsHep2 <- DESeqDataSetFromMatrix(countData = Hep2cts_ENSG_HGNC_BT,
                                  colData = Hep2coldata,
                                  design = ~ Condition_NBE)
ddsHEK1 <- DESeqDataSetFromMatrix(countData = HEK1cts_ENSG_HGNC_BT,
                                  colData = HEK1coldata,
                                  design = ~ Condition_NBE)
ddsHEK2 <- DESeqDataSetFromMatrix(countData = HEK2cts_ENSG_HGNC_BT,
                                  colData = HEK2coldata,
                                  design = ~ Condition_NBE)
ddsHeLa1 <- DESeqDataSetFromMatrix(countData = HeLa1cts_ENSG_HGNC_BT,
                                   colData = HeLa1coldata,
                                   design = ~ Condition_NBE)
ddsHeLa2 <- DESeqDataSetFromMatrix(countData = HeLa2cts_ENSG_HGNC_BT,
                                   colData = HeLa2coldata,
                                   design = ~ Condition_NBE)

#factor condition and set levels
ddsE1E2$Condition <- factor(ddsE1E2$Condition, levels = c("HeLa_Control","HeLa_0R","HeLa_8R",
                                                          "HEK_Control",
                                                          "HEK_0R","HEK_8R",
                                                          "Hep_Control","Hep_0R",
                                                          "Hep_8R"))

ddsE1$Condition <- factor(ddsE1$Condition, levels = c("HeLa_Control","HeLa_0R","HeLa_8R",
                                                      "HEK_Control",
                                                      "HEK_0R","HEK_8R",
                                                      "Hep_Control","Hep_0R",
                                                      "Hep_8R"))

ddsE2$Condition <- factor(ddsE2$Condition, levels = c("HeLa_Control","HeLa_0R","HeLa_8R",
                                                      "HEK_Control",
                                                      "HEK_0R","HEK_8R",
                                                      "Hep_Control","Hep_0R",
                                                      "Hep_8R"))
ddsHep1$Condition_NBE <- factor(ddsHep1$Condition_NBE, levels = c("Hep_Control", "Hep_0R", "Hep_8R"))
ddsHep2$Condition_NBE <- factor(ddsHep2$Condition_NBE, levels = c("Hep_Control", "Hep_0R", "Hep_8R"))
ddsHEK1$Condition_NBE <- factor(ddsHEK1$Condition_NBE, levels = c("HEK_Control", "HEK_0R", "HEK_8R"))
ddsHEK2$Condition_NBE <- factor(ddsHEK2$Condition_NBE, levels = c("HEK_Control", "HEK_0R", "HEK_8R"))
ddsHeLa1$Condition_NBE <- factor(ddsHeLa1$Condition_NBE, levels = c("HeLa_Control", "HeLa_0R", "HeLa_8R"))
ddsHeLa2$Condition_NBE <- factor(ddsHeLa2$Condition_NBE, levels = c("HeLa_Control", "HeLa_0R", "HeLa_8R"))

#Run DESEQ function on dds object, and save to the same object
ddsE1E2 <- DESeq(ddsE1E2)
ddsE1 <- DESeq(ddsE1)
ddsE2 <- DESeq(ddsE2)
ddsHep1<-DESeq(ddsHep1)
ddsHep2<-DESeq(ddsHep2)
ddsHEK1<-DESeq(ddsHEK1)
ddsHEK2<-DESeq(ddsHEK2)
ddsHeLa1<-DESeq(ddsHeLa1)
ddsHeLa2<-DESeq(ddsHeLa2)

#variance stabilized transformation
vsdE1E2 <- vst(ddsE1E2)
vsdE1 <- vst(ddsE1)
vsdE2 <- vst(ddsE2)
vsdHep1 <- vst(ddsHep1)
vsdHep2 <- vst(ddsHep2)
vsdHEK1 <- vst(ddsHEK1)
vsdHEK2 <- vst(ddsHEK2)
vsdHeLa1 <-vst(ddsHeLa1)
vsdHeLa2 <- vst(ddsHeLa2)

#PCA
#Cell Line Variance
Hep_B1_PCA <- plotPCA(vsdHep1, intgroup=c("HeatShock"), returnData=TRUE)
Hep_B1_percentVar <- round(100 * attr(Hep_B1_PCA, "percentVar"))


Hep_B1_PCA_Plot<-ggplot(Hep_B1_PCA, aes(PC1, PC2, colour = HeatShock)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",Hep_B1_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",Hep_B1_percentVar[2],"% variance")) +
  ggtitle("Hep Batch 1 PCA")+
  theme(text = element_text(size = 20))+  
  coord_fixed()

Hep_B1_PCA_Plot

Hep_B2_PCA <- plotPCA(vsdHep2, intgroup=c("HeatShock"), returnData=TRUE)
Hep_B2_percentVar <- round(100 * attr(Hep_B2_PCA, "percentVar"))


Hep_B2_PCA_Plot<-ggplot(Hep_B2_PCA, aes(PC1, PC2, colour = HeatShock)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",Hep_B2_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",Hep_B2_percentVar[2],"% variance")) +
  ggtitle("Hep Batch 2 PCA")+
  theme(text = element_text(size = 20))+  
  coord_fixed()

Hep_B2_PCA_Plot

HEK_B1_PCA <- plotPCA(vsdHEK1, intgroup=c("HeatShock"), returnData=TRUE)
HEK_B1_percentVar <- round(100 * attr(HEK_B1_PCA, "percentVar"))


HEK_B1_PCA_Plot<-ggplot(HEK_B1_PCA, aes(PC1, PC2, colour = HeatShock)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",HEK_B1_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",HEK_B1_percentVar[2],"% variance")) +
  ggtitle("HEK Batch 1 PCA")+
  theme(text = element_text(size = 20))+  
  coord_fixed()

HEK_B1_PCA_Plot

HEK_B2_PCA <- plotPCA(vsdHEK2, intgroup=c("HeatShock"), returnData=TRUE)
HEK_B2_percentVar <- round(100 * attr(HEK_B2_PCA, "percentVar"))


HEK_B2_PCA_Plot<-ggplot(HEK_B2_PCA, aes(PC1, PC2, colour = HeatShock)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",HEK_B2_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",HEK_B2_percentVar[2],"% variance")) +
  ggtitle("HEK Batch 2 PCA")+
  theme(text = element_text(size = 20))+  
  coord_fixed()

HEK_B2_PCA_Plot

HeLa_B1_PCA <- plotPCA(vsdHeLa1, intgroup=c("HeatShockStatus"), returnData=TRUE)
HeLa_B1_percentVar <- round(100 * attr(HeLa_B1_PCA, "percentVar"))


HeLa_B1_PCA_Plot<-ggplot(HeLa_B1_PCA, aes(PC1, PC2, colour = HeatShockStatus)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",HeLa_B1_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",HeLa_B1_percentVar[2],"% variance")) +
  ggtitle("HeLa Batch 1 PCA")+
  theme(text = element_text(size = 20))+  
  coord_fixed()

HeLa_B1_PCA_Plot

HeLa_B2_PCA <- plotPCA(vsdHeLa2, intgroup=c("HeatShock"), returnData=TRUE)
HeLa_B2_percentVar <- round(100 * attr(HeLa_B2_PCA, "percentVar"))


HeLa_B2_PCA_Plot<-ggplot(HeLa_B2_PCA, aes(PC1, PC2, colour = HeatShock)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",HeLa_B2_percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",HeLa_B2_percentVar[2],"% variance")) +
  ggtitle("HeLa Batch 2 PCA")+
  theme(text = element_text(size = 20))+  
  coord_fixed()

HeLa_B2_PCA_Plot

png("HeLaB1_PCA.png", units = "in", width = 10, height = 10, res = 700)
print(HeLa_B1_PCA_Plot)
dev.off()

png("HeLaB2_PCA.png", units = "in", width = 10, height = 10, res = 700)
print(HeLa_B2_PCA_Plot)
dev.off()

png("HepB1_PCA.png", units = "in", width = 10, height = 10, res = 700)
print(Hep_B1_PCA_Plot)
dev.off()

png("HepB2_PCA.png", units = "in", width = 10, height = 10, res = 700)
print(Hep_B2_PCA_Plot)
dev.off()

png("HEKB1_PCA.png", units = "in", width = 10, height = 10, res = 700)
print(HEK_B1_PCA_Plot)
dev.off()

png("HEKB2_PCA.png", units = "in", width = 10, height = 10, res = 700)
print(HEK_B2_PCA_Plot)
dev.off()

#iSEE on vsds
# iSEE(vsdE1E2)
# iSEE(vsdE1)
# iSEE(vsdE2)
# iSEE(vsdHep1)
# iSEE(vsdHep2)
# iSEE(vsdHEK1)
# iSEE(vsdHEK2)
# iSEE(vsdHeLa1)
# iSEE(vsdHeLa2)

#heatmap
counts(ddsHEK1, normalized = T)
counts(ddsHep1, normalized = T)
counts(ddsHeLa1, normalized = T)

counts(ddsHEK2, normalized = T)
counts(ddsHep2, normalized = T)
counts(ddsHeLa2, normalized = T)

matHEK1 <- counts(ddsHEK1, normalized = T)
matHep1 <- counts(ddsHep1, normalized = T)
matHeLa1 <- counts(ddsHeLa1, normalized = T)

matHEK2 <- counts(ddsHEK2, normalized = T)
matHep2 <- counts(ddsHep2, normalized = T)
matHeLa2 <- counts(ddsHeLa2, normalized = T)

matHEK1.z <- t(apply(matHEK1, 1, scale))
matHep1.z <- t(apply(matHep1, 1, scale))
matHeLa1.z <- t(apply(matHeLa1, 1, scale))

matHEK2.z <- t(apply(matHEK2, 1, scale))
matHep2.z <- t(apply(matHep2, 1, scale))
matHeLa2.z <- t(apply(matHeLa2, 1, scale))

colnames(matHEK1.z) <- rownames(HEK1coldata)
colnames(matHep1.z) <- rownames(Hep1coldata)
colnames(matHeLa1.z) <- rownames(HeLa1coldata)

colnames(matHEK2.z) <- rownames(HEK2coldata)
colnames(matHep2.z) <- rownames(Hep2coldata)
colnames(matHeLa2.z) <- rownames(HeLa2coldata)

SRLgenes <- c("TNF", "GNRH2", "PSPN", "GNRH1", "MIA", "SEMA4D", "HBEGF", "LTA", "CNTF", "PGF", "FGF18", "SEMA7A", "JAG1")
HAgenes <- c("HSBP1", "RBBP7", "HSPA1A", "HSPA1B", "HSPA6", "HSBP1L1")

matHEK1.z_SRL <- matHEK1.z[SRLgenes,]
matHEK2.z_SRL <- matHEK2.z[SRLgenes,]
matHep1.z_SRL <- matHep1.z[SRLgenes,]
matHep2.z_SRL <- matHep2.z[SRLgenes,]
matHeLa1.z_SRL <- matHeLa1.z[SRLgenes,]
matHeLa2.z_SRL <- matHeLa2.z[SRLgenes,]

matHEK1.z_HA <- matHEK1.z[HAgenes,]
matHEK2.z_HA <- matHEK2.z[HAgenes,]
matHep1.z_HA <- matHep1.z[HAgenes,]
matHep2.z_HA <- matHep2.z[HAgenes,]
matHeLa1.z_HA <- matHeLa1.z[HAgenes,]
matHeLa2.z_HA <- matHeLa2.z[HAgenes,]

setwd("D:/Plots/Plots")
png("haHEK1.png", units = "in", width = 6, height = 5, res = 700)
haHEK1 <- Heatmap(matHEK1.z_HA, cluster_rows = T, cluster_columns = F, column_labels = colnames(matHEK1.z_HA), name = "Z-Score", row_labels = rownames(matHEK1.z_HA))
draw(haHEK1)
dev.off()

png("haHEK2.png", units = "in", width = 6, height = 5, res = 700)
haHEK2 <- Heatmap(matHEK2.z_HA, cluster_rows = T, cluster_columns = F, column_labels = colnames(matHEK2.z_HA), name = "Z-Score", row_labels = rownames(matHEK2.z_HA))
draw(haHEK2)
dev.off()

png("haHep1.png", units = "in", width = 6, height = 5, res = 700)
haHep1 <- Heatmap(matHep1.z_HA, cluster_rows = T, cluster_columns = F, column_labels = colnames(matHep1.z_HA), name = "Z-Score", row_labels = rownames(matHep1.z_HA))
draw(haHep1)
dev.off()

png("haHep2.png", units = "in", width = 6, height = 5, res = 700)
haHep2 <- Heatmap(matHep2.z_HA, cluster_rows = T, cluster_columns = F, column_labels = colnames(matHep2.z_HA), name = "Z-Score", row_labels = rownames(matHep2.z_HA))
draw(haHep2)
dev.off()

png("haHeLa1.png", units = "in", width = 6, height = 5, res = 700)
haHeLa1 <- Heatmap(matHeLa1.z_HA, cluster_rows = T, cluster_columns = F, column_labels = colnames(matHeLa1.z_HA), name = "Z-Score", row_labels = rownames(matHeLa1.z_HA))
draw(haHeLa1)
dev.off()

png("haHeLa2.png", units = "in", width = 6, height = 5, res = 700)
haHeLa2 <- Heatmap(matHeLa2.z_HA, cluster_rows = T, cluster_columns = F, column_labels = colnames(matHeLa2.z_HA), name = "Z-Score", row_labels = rownames(matHeLa2.z_HA))
draw(haHeLa2)
dev.off()

#GSEA
#Set Desired Organism
setwd("~/DEGtoGSEA/DEG Results")
organism = "org.Hs.eg.db"
#Upload Batch 1 DEG results
dfHEK0RvsCnt_B2 <- read.csv("HEK0RvsHEKControl_DEG_Results_2022.csv", header = TRUE)
dfHEK8RvsCnt_B2 <- read.csv("HEK8RvsHEKControl_DEG_Results_2022.csv", header = TRUE)
dfHEK8Rvs0R_B2 <- read.csv("HEK8RvsHEK0R_DEG_Results_2022.csv", header = TRUE)

dfHep0RvsCnt_B2 <- read.csv("Hep0RvsHepControl_DEG_Results_2022.csv", header = TRUE)
dfHep8RvsCnt_B2 <- read.csv("Hep8RvsHepControl_DEG_Results_2022.csv", header = TRUE)
dfHep8Rvs0R_B2 <- read.csv("Hep8RvsHep0R_DEG_Results_2022.csv", header = TRUE)

dfHeLa0RvsCnt_B2 <- read.csv("HeLa0RvsHeLaControl_DEG_Results_2022.csv", header = TRUE)
dfHeLa8RvsCnt_B2 <- read.csv("HeLa8RvsHeLaControl_DEG_Results_2022.csv", header = TRUE)
dfHeLa8Rvs0R_B2 <- read.csv("HeLa8RvsHeLa0R_DEG_Results_2022.csv", header = TRUE)

og_genelist_HEK0RvsCnt_B2 <- dfHEK0RvsCnt_B2$log2FoldChange
og_genelist_HEK8RvsCnt_B2 <- dfHEK8RvsCnt_B2$log2FoldChange
og_genelist_HEK8Rvs0R_B2 <- dfHEK8Rvs0R_B2$log2FoldChange

og_genelist_Hep0RvsCnt_B2 <- dfHep0RvsCnt_B2$log2FoldChange
og_genelist_Hep8RvsCnt_B2 <- dfHep8RvsCnt_B2$log2FoldChange
og_genelist_Hep8Rvs0R_B2 <- dfHep8Rvs0R_B2$log2FoldChange

og_genelist_HeLa0RvsCnt_B2 <- dfHeLa0RvsCnt_B2$log2FoldChange
og_genelist_HeLa8RvsCnt_B2 <- dfHeLa8RvsCnt_B2$log2FoldChange
og_genelist_HeLa8Rvs0R_B2 <- dfHeLa8Rvs0R_B2$log2FoldChange

names(og_genelist_HEK0RvsCnt_B2) <-dfHEK0RvsCnt_B2$ensembl_gene_id
names(og_genelist_HEK8RvsCnt_B2) <-dfHEK8RvsCnt_B2$ensembl_gene_id
names(og_genelist_HEK8Rvs0R_B2) <-dfHEK8Rvs0R_B2$ensembl_gene_id

names(og_genelist_Hep0RvsCnt_B2) <-dfHep0RvsCnt_B2$ensembl_gene_id
names(og_genelist_Hep8RvsCnt_B2) <-dfHep8RvsCnt_B2$ensembl_gene_id
names(og_genelist_Hep8Rvs0R_B2) <-dfHep8Rvs0R_B2$ensembl_gene_id

names(og_genelist_HeLa0RvsCnt_B2) <-dfHeLa0RvsCnt_B2$ensembl_gene_id
names(og_genelist_HeLa8RvsCnt_B2) <-dfHeLa8RvsCnt_B2$ensembl_gene_id
names(og_genelist_HeLa8Rvs0R_B2) <-dfHeLa8Rvs0R_B2$ensembl_gene_id

gene_list_HEK0RvsCnt_B2 <- na.omit(og_genelist_HEK0RvsCnt_B2)
gene_list_HEK8RvsCnt_B2 <- na.omit(og_genelist_HEK8RvsCnt_B2)
gene_list_HEK8Rvs0R_B2 <- na.omit(og_genelist_HEK8Rvs0R_B2)

gene_list_Hep0RvsCnt_B2 <- na.omit(og_genelist_Hep0RvsCnt_B2)
gene_list_Hep8RvsCnt_B2 <- na.omit(og_genelist_Hep8RvsCnt_B2)
gene_list_Hep8Rvs0R_B2 <- na.omit(og_genelist_Hep8Rvs0R_B2)

gene_list_HeLa0RvsCnt_B2 <- na.omit(og_genelist_HeLa0RvsCnt_B2)
gene_list_HeLa8RvsCnt_B2 <- na.omit(og_genelist_HeLa8RvsCnt_B2)
gene_list_HeLa8Rvs0R_B2 <- na.omit(og_genelist_HeLa8Rvs0R_B2)

gene_list_HEK0RvsCnt_B2 <- sort(gene_list_HEK0RvsCnt_B2, decreasing = TRUE)
gene_list_HEK8RvsCnt_B2 <- sort(gene_list_HEK8RvsCnt_B2, decreasing = TRUE)
gene_list_HEK8Rvs0R_B2 <- sort(gene_list_HEK8Rvs0R_B2, decreasing = TRUE)

gene_list_Hep0RvsCnt_B2 <- sort(gene_list_Hep0RvsCnt_B2, decreasing = TRUE)
gene_list_Hep8RvsCnt_B2 <- sort(gene_list_Hep8RvsCnt_B2, decreasing = TRUE)
gene_list_Hep8Rvs0R_B2 <- sort(gene_list_Hep8Rvs0R_B2, decreasing = TRUE)

gene_list_HeLa0RvsCnt_B2 <- sort(gene_list_HeLa0RvsCnt_B2, decreasing = TRUE)
gene_list_HeLa8RvsCnt_B2 <- sort(gene_list_HeLa8RvsCnt_B2, decreasing = TRUE)
gene_list_HeLa8Rvs0R_B2 <- sort(gene_list_HeLa8Rvs0R_B2, decreasing = TRUE)

gse_HEK8Rvs0R_B2 <- gseGO(geneList=gene_list_HEK8Rvs0R_B2, 
                           ont ="ALL", 
                           keyType = "ENSEMBL",
                           nPerm = 10000, 
                           minGSSize = 3, 
                           maxGSSize = 800, 
                           pvalueCutoff = 0.05, 
                           verbose = TRUE, 
                           OrgDb = organism, 
                           pAdjustMethod = "none")

require(DOSE)

setwd("D:/Plots/Plots/GSEA")
png("HEK8Rvs0R_B2_Dotplot.png", units = "in", width = 10, height = 12, res = 700)
HEK8Rvs0R_B2_Dotplot <- enrichplot::dotplot(gse_HEK8Rvs0R_B2, orderBy = "NES", showCategory = 10, split=".sign", x = "GeneRatio") + 
  facet_grid(.~.sign) +
  geom_point() +
  theme(text = element_text(size = 20))
print(HEK8Rvs0R_B2_Dotplot)
dev.off()

x2 <- pairwise_termsim(gse_HEK8Rvs0R_B2)

png("HEK8Rvs0R_B2_emap.png", units = "in", width = 8, height = 8, res = 700)
HEK8Rvs0R_B2_emap <- emapplot(x2, showCategory = 10, orderBy = "NES")
print(HEK8Rvs0R_B2_emap)
dev.off()

png("HEK8Rvs0R_B2_Ridgeplot.png", units = "in", width = 10, height = 10, res = 700)
HEK8Rvs0R_B2_Ridgeplot <- ridgeplot(gse_HEK8Rvs0R_B2, showCategory = 10, orderBy = "NES") + 
  labs(x = "enrichment distribution") +
  theme(text = element_text(size = 20)) 
print(HEK8Rvs0R_B2_Ridgeplot)
dev.off()

png("HeLa0RvsCnt_B1_GSEAPlot.png", units = "in", width = 10, height = 10, res = 700)
HeLa0RvsCnt_B1_GSEAPlot <- gseaplot(gse_HEK0RvsCnt_B1, by = "all", title = gse_HeLa0RvsCnt_B1$Description[2], geneSetID = 1)
print(HEK0RvsCnt_B1_GSEAPlot)
dev.off()

