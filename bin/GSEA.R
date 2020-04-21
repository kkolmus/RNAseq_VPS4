# R version: 3.4.4

# cleaning, loading packages and setting working directory

rm(list = ls())

# install.packages("readxl")
library(readxl)
# install.packages("WriteXLS")
library(WriteXLS)
# install.packages("clusterProfiler)
library(clusterProfiler)
# install.packages("org.Hs.eg.db")
library(org.Hs.eg.db)
# install.packages("tidyr")
library(tidyr)
# install.packages("dplyr")
library(dplyr)
#install.packages("imputeTS")
library(imputeTS)
#install.packages("ggplot2")
library(ggplot2)

set.seed(123)

setwd("~/Desktop/RNAseq Vps4/2018.11.09 R 3.4.4 complete reanalysis")

df <- read_excel("ES_Vps4_project_full.xlsx")
df <- df[!grepl("NR_", df$RefSeq),] # remove non-coding transcripts
df <- mutate_at(df, vars(-c(RefSeq, Symbol)), funs(as.numeric(.)))

# filtering data

# remove NA values

VPS4AB <- drop_na(df, 
                  siVps4AB_68_72.vs.NT_adj.pval, siVps4AB_66_73.vs.NT_adj.pval,
                  siVps4AB_68_72.vs.siCtrl_1_adj.pval, siVps4AB_66_73.vs.siCtrl_1_adj.pval)

# Converting RefSeq to EntrezID for significanlt upregulated genes in each condition

VPS4AB_EntrezID <- bitr(VPS4AB$RefSeq, fromType = "REFSEQ", toType = "ENTREZID",  OrgDb = "org.Hs.eg.db")

VPS4AB_input <- as.character(VPS4AB_EntrezID[,2])


# GO GSEA siVPS4AB#1 vs. NT

a <- VPS4AB[,c(1,7)]

## feature 1: numeric vector
geneList <- unname(unlist(a[,2]))
## feature 2: named vector
#names(geneList) <- as.character(d[,1])
names(geneList) <- VPS4AB_input
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

VPS4AB_1_NT_gseGO <- gseGO(geneList = geneList,
                           ont = "BP",
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           nPerm = 1000, # Permutations are used to calculate enrichment p value
                           minGSSize = 10,
                           maxGSSize = 500,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH", 
                           verbose = TRUE, 
                           seed = TRUE, 
                           by = "fgsea")

View(as.data.frame(VPS4AB_1_NT_gseGO))

write.table(as.data.frame(VPS4AB_1_NT_gseGO),
            file = "VPS4AB_1_NT_gseGO.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

gseaplot(VPS4AB_1_NT_gseGO, geneSetID = "GO:0006954", title = "Inflammatory response")
gseaplot(VPS4AB_1_NT_gseGO, geneSetID = "GO:2001239", title = "Regulation of extrinsic apoptotic signaling pathway in absence of ligand")
gseaplot(VPS4AB_1_NT_gseGO, geneSetID = "GO:0043122", title = "Regulation of I-kappaB kinase/NF-kappaB signaling")
gseaplot(VPS4AB_1_NT_gseGO, geneSetID = "GO:1901222", title = "Regulation of NIK/NF-kappaB signaling")


# GO GSEA siVPS4AB#2 vs. NT

b <- VPS4AB[,c(1,9)]

## feature 1: numeric vector
geneList_b <- unname(unlist(b[,2]))
## feature 2: named vector
names(geneList_b) <- VPS4AB_input
## feature 3: decreasing order
geneList_b <- sort(geneList_b, decreasing = TRUE)

VPS4AB_2_NT_gseGO <- gseGO(geneList = geneList_b,
                           ont = "BP",
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           nPerm = 1000, # Permutations are used to calculate enrichment p value
                           minGSSize = 10,
                           maxGSSize = 500,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH", 
                           verbose = TRUE, 
                           seed = TRUE, 
                           by = "fgsea")

View(as.data.frame(VPS4AB_2_NT_gseGO))

write.table(as.data.frame(VPS4AB_2_NT_gseGO),
            file = "VPS4AB_2_NT_gseGO.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

gseaplot(VPS4AB_2_NT_gseGO, geneSetID = "GO:0006954", title = "Inflammatory response")
gseaplot(VPS4AB_2_NT_gseGO, geneSetID = "GO:2001239", title = "Regulation of extrinsic apoptotic signaling pathway in absence of ligand")
gseaplot(VPS4AB_2_NT_gseGO, geneSetID = "GO:0043122", title = "Regulation of I-kappaB kinase/NF-kappaB signaling")
gseaplot(VPS4AB_2_NT_gseGO, geneSetID = "GO:1901222", title = "Regulation of NIK/NF-kappaB signaling")


# GO GSEA siVPS4AB#1 vs. Ctrl

c <- VPS4AB[,c(1,23)]

## feature 1: numeric vector
geneList_c <- unname(unlist(c[,2]))
## feature 2: named vector
names(geneList_c) <- VPS4AB_input
## feature 3: decreasing order
geneList_c <- sort(geneList_c, decreasing = TRUE)

VPS4AB_1_Ctrl_gseGO <- gseGO(geneList = geneList_c,
                             ont = "BP",
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             nPerm = 1000, # Permutations are used to calculate enrichment p value
                             minGSSize = 10,
                             maxGSSize = 500,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH", 
                             verbose = TRUE, 
                             seed = TRUE, 
                             by = "fgsea")

View(as.data.frame(VPS4AB_1_Ctrl_gseGO))

write.table(as.data.frame(VPS4AB_1_Ctrl_gseGO),
            file = "VPS4AB_1_Ctrl_gseGO.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

gseaplot(VPS4AB_1_Ctrl_gseGO, geneSetID = "GO:0006954", title = "Inflammatory response")
gseaplot(VPS4AB_1_Ctrl_gseGO, geneSetID = "GO:2001239", title = "Regulation of extrinsic apoptotic signaling pathway in absence of ligand")
gseaplot(VPS4AB_1_Ctrl_gseGO, geneSetID = "GO:0043122", title = "Regulation of I-kappaB kinase/NF-kappaB signaling")
gseaplot(VPS4AB_1_Ctrl_gseGO, geneSetID = "GO:1901222", title = "Regulation of NIK/NF-kappaB signaling")


# GO GSEA siVPS4AB#2 vs. Ctrl

d <- VPS4AB[,c(1,25)]

## feature 1: numeric vector
geneList_d <- unname(unlist(d[,2]))
## feature 2: named vector
names(geneList_d) <- VPS4AB_input
## feature 3: decreasing order
geneList_d <- sort(geneList_d, decreasing = TRUE)

VPS4AB_2_Ctrl_gseGO <- gseGO(geneList = geneList_d,
                             ont = "BP",
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             nPerm = 1000, # Permutations are used to calculate enrichment p value
                             minGSSize = 10,
                             maxGSSize = 500,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH", 
                             verbose = TRUE, 
                             seed = TRUE, 
                             by = "fgsea")

View(as.data.frame(VPS4AB_2_Ctrl_gseGO))

write.table(as.data.frame(VPS4AB_2_Ctrl_gseGO),
            file = "VPS4AB_2_Ctrl_gseGO.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

gseaplot(VPS4AB_2_Ctrl_gseGO, geneSetID = "GO:0006954", title = "Inflammatory response")
gseaplot(VPS4AB_2_Ctrl_gseGO, geneSetID = "GO:2001239", title = "Regulation of extrinsic apoptotic signaling pathway in absence of ligand")
gseaplot(VPS4AB_2_Ctrl_gseGO, geneSetID = "GO:0043122", title = "Regulation of I-kappaB kinase/NF-kappaB signaling")
gseaplot(VPS4AB_2_Ctrl_gseGO, geneSetID = "GO:1901222", title = "regulation of NIK/NF-kappaB signaling")




