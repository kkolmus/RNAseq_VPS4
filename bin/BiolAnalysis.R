# R 3.4.4

# cleaning
rm(list = ls())

# install.packages("readxl")
library(readxl)
# install.packages("#writeXLS")
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
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)
# install.packages("ComplexHeatmap")
library(ComplexHeatmap)
fontsize <- 0.6
# install.packages("VennDiagram")
library(VennDiagram)
# install.packages("gridExtra")
library(gridExtra)

############################################################
## read data frames normalised against control conditions ##
############################################################

setwd("~/Desktop/RNAseq Vps4/2018.11.09 R 3.4.4 complete reanalysis")

df <- read_excel("ES_Vps4_project_full.xlsx")
df <- df[!grepl("NR_", df$RefSeq),] # remove non-coding transcripts
df <- mutate_at(df, vars(-c(RefSeq, Symbol)), funs(as.numeric(.)))

# removing NA

df_Vps4AB <- drop_na(df, 
                     siVps4AB_68_72.vs.NT_adj.pval, siVps4AB_66_73.vs.NT_adj.pval,
                     siVps4AB_68_72.vs.siCtrl_1_adj.pval, siVps4AB_66_73.vs.siCtrl_1_adj.pval)


df_Vps4A <- drop_na(df, 
                    siVps4A_68.vs.NT_adj.pval, siVps4A_66.vs.NT_adj.pval,
                    siVps4A_68.vs.siCtrl_1_adj.pval, siVps4A_66.vs.siCtrl_1_adj.pval)

df_Vps4B <- drop_na(df, 
                    siVps4B_72.vs.NT_adj.pval, siVps4B_73.vs.NT_adj.pval,
                    siVps4B_72.vs.siCtrl_1_adj.pval, siVps4B_73.vs.siCtrl_1_adj.pval)

#####################################
## Filtering criteria for analysis ##
#####################################

UP = 1.5         # minimal fold change to consider the gene to be upregulated
#DOWN = 0.5       # minimal fold change to consider the gene to be downregulated
pval = 0.05      # p-value to consider change to be statistically significant

#################################
## filtering data for analysis ##
#################################

# Statistically upregulated genes after concurrent Vps4A silencing

UP_Vps4AB <- filter(df_Vps4AB, 
                    siVps4AB_68_72.vs.NT_FC >= UP, siVps4AB_68_72.vs.NT_adj.pval < pval,
                    siVps4AB_68_72.vs.siCtrl_1_FC >= UP, siVps4AB_68_72.vs.siCtrl_1_adj.pval < pval,
                    siVps4AB_66_73.vs.NT_FC >= UP, siVps4AB_66_73.vs.NT_adj.pval < pval,
                    siVps4AB_66_73.vs.siCtrl_1_FC >= UP, siVps4AB_66_73.vs.siCtrl_1_adj.pval < pval)

# write.table(as.data.frame(UP_Vps4AB),
#             file = "Genes_UP_Vps4AB.txt",
#             quote = FALSE,
#             sep = "\t",
#             row.names = TRUE)

# Statistically downregulated genes after concurrent Vps4B silencing

# DOWN_Vps4AB <- filter(df_Vps4AB, 
#                    siVps4AB_68_72.vs.NT_FC <= DOWN, siVps4AB_68_72.vs.NT_adj.pval < pval,
#                    siVps4AB_68_72.vs.siCtrl_1_FC <= DOWN, siVps4AB_68_72.vs.siCtrl_1_adj.pval < pval,
#                    siVps4AB_66_73.vs.NT_FC <= DOWN, siVps4AB_66_73.vs.NT_adj.pval < pval,
#                    siVps4AB_66_73.vs.siCtrl_1_FC <= DOWN, siVps4AB_66_73.vs.siCtrl_1_adj.pval < pval)

# write.table(as.data.frame(DOWN_Vps4AB),
#             file = "Genes_DOWN_Vps4AB.txt",
#             quote = FALSE,
#             sep = "\t",
#             row.names = TRUE)

##################
## GO analysis  ##
##################

# Converting names for upregulated genes

Vps4AB_EntrezID_UP <- bitr(UP_Vps4AB$RefSeq,
                           fromType = "REFSEQ", toType = "ENTREZID", 
                           OrgDb = "org.Hs.eg.db")

# Converting names for downregulated genes

# Vps4AB_EntrezID_DOWN <- bitr(DOWN_Vps4AB$RefSeq,
#                              fromType = "REFSEQ", toType = "ENTREZID", 
#                              OrgDb = "org.Hs.eg.db")

# selecting background genes for analysis

Vps4AB_EntrezID_project <- bitr(df_Vps4AB$RefSeq,
                                fromType = "REFSEQ", toType = "ENTREZID", 
                                OrgDb = "org.Hs.eg.db")

######################################
## Upregulated biological processes ##
######################################

Vps4AB_ego_BP_UP <- enrichGO(gene = Vps4AB_EntrezID_UP$ENTREZID, 
                             universe = Vps4AB_EntrezID_project$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff  = 0.05,
                             minGSSize = 10,
                             maxGSSize = 500,
                             readable = TRUE)

# dotplot(Vps4AB_ego_BP_UP, showCategory = 30, title = "Upregulated BP after siVps4AB")

# write.table(as.data.frame(Vps4AB_ego_BP_UP),
#             file = "BP_Vps4AB_ego_UP.txt",
#             quote = FALSE,
#             sep = "\t",
#             row.names = TRUE)

Vps4AB_ego_BP_UP_sim <- simplify(Vps4AB_ego_BP_UP, cutoff = 0.6, by = "p.adjust", select_fun = min)

# dotplot(Vps4AB_ego_BP_UP_sim, showCategory = 25, title = "Upregulated BP after siVps4AB")

BP_UP <- as.data.frame(Vps4AB_ego_BP_UP_sim) %>% mutate(GeneRatio = Count/510) %>%
  filter(Count >= 10) %>% top_n(25)

p <- ggplot(BP_UP, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) + #fct_reorder
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 12) + # size of description
  scale_colour_gradient(limits=c(10^(-10), 0.05), low="red", high = "blue") +
  ylab(NULL) + xlab("Gene ratio") + 
  labs(colour = "p-value", size = "Gene count") + 
  ggtitle("Gene Ontology Enrichment\nUpregulated Biological Processes") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", color = "black")) # size of title

p

# element_text(family = NULL, face = NULL, colour = NULL, size = NULL,
#              hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
#              color = NULL)

# theme_bw(base_size = 11, base_family = "",
#          base_line_size = base_size/22, base_rect_size = base_size/22)

# write.table(as.data.frame(Vps4AB_ego_BP_UP_sim),
#             file = "BP_Vps4AB_ego_UP_sim.txt",
#             quote = FALSE,
#             sep = "\t",
#             row.names = TRUE)

########################################
## Downregulated biological processes ##
########################################

# Vps4AB_ego_BP_DOWN <- enrichGO(gene = Vps4AB_EntrezID_DOWN$ENTREZID,
#                              universe = Vps4AB_EntrezID_project$ENTREZID,
#                              OrgDb = org.Hs.eg.db,
#                              ont = "BP",
#                              pAdjustMethod = "BH",
#                              pvalueCutoff  = 0.05,
#                              qvalueCutoff  = 0.05,
#                              minGSSize = 10,
#                              maxGSSize = 500,
#                              readable = TRUE)

# dotplot(Vps4AB_ego_BP_DOWN, showCategory = 30, title = "Downregulated BP after siVps4AB")

# write.table(as.data.frame(Vps4AB_ego_BP_DOWN),
#             file = "BP_Vps4AB_ego_DOWN.txt",
#             quote = FALSE,
#             sep = "\t",
#             row.names = TRUE)

#Vps4AB_ego_BP_DOWN_sim <- simplify(Vps4AB_ego_BP_DOWN, cutoff = 0.6, by = "p.adjust", select_fun = min)

# dotplot(Vps4AB_ego_BP_DOWN_sim, showCategory = 25, title = "Downregulated BP after siVps4AB")

# BP_DOWN <- as.data.frame(Vps4AB_ego_BP_DOWN_sim) %>% mutate(GeneRatio = Count/199) %>%
#   filter(Count >= 10) %>% top_n(25)

# p <- ggplot(BP_DOWN, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) + #fct_reorder
#   geom_point(aes(size = Count, color = p.adjust)) +
#   theme_bw(base_size = 12) +
#   scale_colour_gradient(limits=c(10^(-25), 0.05), low="red", high = "blue") +
#   ylab(NULL) + xlab("Gene ratio") + 
#   labs(colour = "p-value", size = "Gene count") + 
#   ggtitle("Gene Ontology Enrichment\nDownregulated Biological Processes") +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "black"))

# p

# write.table(as.data.frame(Vps4AB_ego_BP_DOWN_sim),
#             file = "BP_Vps4AB_ego_DOWN_sim.txt",
#             quote = FALSE,
#             sep = "\t",
#             row.names = TRUE)


#############
## Heatmap ##
#############

# TPM

samples <- df[,c(1,2,51:59)]
samples <- dplyr::rename(samples, NT = NT_reads, siCtrl = siCtrl1_reads, siCtrl_2 = siCtrl2_reads,
                         siVps4A_1 = siVps4A_66_reads, siVps4A_2 = siVps4A_68_reads,
                         siVps4B_1 = siVps4B_72_reads, siVps4B_2 = siVps4B_73_reads,
                         siVps4AB_1 = siVps4AB_68_72_reads, siVps4AB_2 = siVps4AB_66_73_reads)
samples <- samples[!grepl("NR_", samples$RefSeq),]
sum(is.na(samples$RefSeq))
samples[duplicated(samples$RefSeq),]

listMarts()

ensembl <- useMart("ensembl")

ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

RefSeqIDs <- samples$RefSeq
RefSeqIDs[duplicated(RefSeqIDs)]

length <- getBM(attributes = c("refseq_mrna", "transcript_length"), 
                filters = "refseq_mrna", values = RefSeqIDs, mart = ensembl)
length <- length[! duplicated(length$refseq_mrna),]
length <- dplyr::select(length, RefSeq = refseq_mrna, transcript_length)

samples <- merge(samples, length, by = c("RefSeq"))

# samples <- samples[,c(1,2,12,3:11)] # all samples 
samples <- samples[,c(1,2,12,3,4,6:11)] # all samples 
head(samples)

# normalisation for gene length, TPM1 = reads per kilobase = RPK
TPM1 <- mutate_at(samples, vars(-c(RefSeq, Symbol)), funs(./samples$transcript_length))
head(TPM1)
TPM1 <- TPM1[, -3]
head(TPM1)

# normalisation for sequencing depth, TPM2 = scaling factors
TPM2 <- apply(TPM1[, c(3:10)], 2, sum)
TPM2 <- TPM2/10^(6)

TPM1 <- mutate_at(TPM1, vars(-c(RefSeq, Symbol)), funs(as.numeric(.)))
apply(TPM1, 2, class)

TPM2 <- unname(unlist(TPM2))
class(TPM2)

# transcripts per million

TPM_temp <- sweep(TPM1[,-c(1,2)], MARGIN = 2, TPM2, FUN = "/")
TPM_names <- TPM1[, c(1,2)]

TPM <- cbind(TPM_names, TPM_temp)

# Heatmap for HTC116 after differential VPS4 silencing

Av <- apply(TPM[, c(3:10)], 1, mean)
StDev <- apply(TPM[, c(3:10)], 1, sd)
TPM_HM <- mutate_at(TPM, vars(-c(RefSeq, Symbol)), funs(((.)-Av)/StDev))

HeatmapData <- filter(TPM_HM, RefSeq %in% UP_Vps4AB[[1]])

HeatmapData_matrix <- as.matrix(HeatmapData[ , c(3:10)])
#HeatmapData_matrix <- t(HeatmapData_matrix)


Heatmap(HeatmapData_matrix,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_dend_side = "left",
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        column_title = "Heatmap of upregulated genes",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 12),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 12, fontface = "bold"))


####################################
### Inflammatory response heatmap

BP_UP2 <- tibble::column_to_rownames(BP_UP, var = "ID")
inflammatory_response <- BP_UP2["GO:0006954", 7]
inflammatory_response <- gsub("/", "', '", inflammatory_response)
inflammatory_response
inflammatory_response_genes <- c('LDLR', 'IL6R', 'CXCL8', 'SERPINE1', 'BDKRB2', 'TGFB1', 'CSF1', 'MGLL', 'PLAA', 'NFKB2', 'KDM6B', 'THEMIS2', 'MVK', 'CYP4F11', 'BIRC2', 'HYAL3', 'B4GALT1', 'CXCL1', 'TNFRSF9', 'IL18', 'IRAK2', 'CFB', 'F3', 'CXCL2', 'HMOX1', 'NT5E', 'SEMA7A', 'PLA2G4C', 'TNFRSF10D', 'TNFRSF10A', 'ECM1', 'ELF3', 'PTGES', 'HSPG2', 'SMAD3', 'TNIP1', 'TNFAIP3', 'RELB', 'DUSP10', 'GGT1', 'SBNO2', 'NFKBIA', 'ZC3H12A', 'HYAL1', 'IL31RA', 'BIRC3')

inflammatory_response_input <- filter(TPM_HM, Symbol %in% inflammatory_response_genes)
inflammatory_response_input <- tibble::column_to_rownames(inflammatory_response_input, var = "Symbol")

TPM_input <- as.matrix(inflammatory_response_input[ , c(2:9)])
#TPM_input <- t(TPM_input)

Heatmap(TPM_input,
        cluster_columns = FALSE,
        column_names_side = "top",
        column_dend_side = "top",
        column_names_gp = gpar(cex = fontsize),
        row_names_side = "left",
        row_dend_side = "left",
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        column_title = "Inflammatory response",
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"))


###################################################
### Positive regulation of programmed cell death

BP_UP2 <- tibble::column_to_rownames(BP_UP, var = "ID")
PCD <- BP_UP2["GO:0043068", 7]
PCD <- gsub("/", "', '", PCD)
PCD
PCD_genes <- c('SMPD1', 'TGFB1', 'HMGCR', 'NQO1', 'CYLD', 'MAPT', 'BCAP31', 'SKIL', 'BAK1', 'BIK', 'ATF3', 'B4GALT1', 'CYR61', 'F3', 'HMOX1', 'MSX1', 'PLAUR', 'KLF11', 'AKR1C3', 'PEA15', 'TNFRSF10A', 'SQSTM1', 'RHOB', 'DUSP1', 'FOSL1', 'IFI27', 'SMAD3', 'SFN', 'ATG7', 'NR4A3', 'TRIO', 'LATS2', 'PMAIP1', 'ZC3H12A', 'PINK1', 'AIFM2', 'CASP7', 'CDKN1A', 'BCL2L11', 'SIK1', 'OSGIN1')

PCD_input <- filter(TPM_HM, Symbol %in% PCD_genes)
PCD_input <- tibble::column_to_rownames(PCD_input, var = "Symbol")

TPM_input_2 <- as.matrix(PCD_input[ , c(2:9)])
#TPM_input_2 <- t(TPM_input_2)

Heatmap(TPM_input_2,
        cluster_columns = FALSE,
        column_names_side = "top",
        column_dend_side = "top",
        column_names_gp = gpar(cex = fontsize),
        row_names_side = "left",
        row_dend_side = "left",
        row_names_gp = gpar(cex = fontsize),
        row_dend_width = unit(3, "cm"),
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "ward.D",
        name = "z-score",
        column_title = "Positive regulation \nof programmed cell death",
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"))



################### 
## Venn diagrams ##
###################

Vps4AB_1_NT <- filter(df_Vps4AB, siVps4AB_68_72.vs.NT_FC >= UP, siVps4AB_68_72.vs.NT_adj.pval < pval)
Vps4AB_1_Ctrl1 <- filter(df_Vps4AB, siVps4AB_68_72.vs.siCtrl_1_FC >= UP, siVps4AB_68_72.vs.siCtrl_1_adj.pval < pval)
Vps4AB_2_NT <- filter(df_Vps4AB, siVps4AB_66_73.vs.NT_FC >= UP, siVps4AB_66_73.vs.NT_adj.pval < pval)
Vps4AB_2_Ctrl1 <- filter(df_Vps4AB, siVps4AB_66_73.vs.siCtrl_1_FC >= UP, siVps4AB_66_73.vs.siCtrl_1_adj.pval < pval)

n1 <- Vps4AB_1_NT[ , 1]
n2 <- Vps4AB_1_Ctrl1[ , 1]
n3 <- Vps4AB_2_NT[ , 1]
n4 <- Vps4AB_2_Ctrl1[ , 1]

n12 <- nrow(intersect(n1, n2))
n34 <- nrow(intersect(n3, n4))
n13 <- nrow(intersect(n1, n3))
n14 <- nrow(intersect(n1, n4))
n23 <- nrow(intersect(n2, n3))
n24 <- nrow(intersect(n2, n4))

n123 <- nrow(Reduce(intersect, list(n1, n2, n3)))
n124 <- nrow(Reduce(intersect, list(n1, n2, n4)))
n134 <- nrow(Reduce(intersect, list(n1, n3, n4)))
n234 <- nrow(Reduce(intersect, list(n2, n3, n4)))

n1234 <- nrow(intersect(intersect(n1, n2), intersect(n3, n4)))

grid.newpage()

Vps4AB <- draw.quad.venn(area1 = nrow(n1), area2 = nrow(n2), area3 = nrow(n3), area4 = nrow(n4),
                         n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234,
                         category = c("siVPS4AB#1/NT", "siVPS4AB#1/siCtrl",
                                      "siVps4AB#2/NT", "siVPS4AB#2/siCtrl"), 
                         fill = c("dodgerblue", "deeppink", "goldenrod1", "chartreuse4"),
                         fontfamily = "sans", cat.fontfamily = "sans",
                         cex = 2, cat.cex = 2,
                         cat.dist = c(0.25, 0.25, 0.15, 0.15),
                         cat.just = list(c(0.35,0.35) , c(0.70,0.35) , c(0.35,1.0) , c(0.60,1.0)))

Vps4AB



### Vps4A

Vps4A_1_NT <- filter(df_Vps4A, siVps4A_68.vs.NT_FC >= UP, siVps4A_68.vs.NT_adj.pval < pval)
Vps4A_1_Ctrl1 <- filter(df_Vps4A, siVps4A_68.vs.siCtrl_1_FC >= UP, siVps4A_68.vs.siCtrl_1_adj.pval < pval)
Vps4A_2_NT <- filter(df_Vps4A, siVps4A_66.vs.NT_FC >= UP, siVps4A_66.vs.NT_adj.pval < pval)
Vps4A_2_Ctrl1 <- filter(df_Vps4A, siVps4A_66.vs.siCtrl_1_FC >= UP, siVps4A_66.vs.siCtrl_1_adj.pval < pval)

n1 <- Vps4A_1_NT[ , 1]
n2 <- Vps4A_1_Ctrl1[ , 1]
n3 <- Vps4A_2_NT[ , 1]
n4 <- Vps4A_2_Ctrl1[ , 1]

n12 <- nrow(intersect(n1, n2))
n34 <- nrow(intersect(n3, n4))
n13 <- nrow(intersect(n1, n3))
n14 <- nrow(intersect(n1, n4))
n23 <- nrow(intersect(n2, n3))
n24 <- nrow(intersect(n2, n4))

n123 <- nrow(Reduce(intersect, list(n1, n2, n3)))
n124 <- nrow(Reduce(intersect, list(n1, n2, n4)))
n134 <- nrow(Reduce(intersect, list(n1, n3, n4)))
n234 <- nrow(Reduce(intersect, list(n2, n3, n4)))

n1234 <- nrow(intersect(intersect(n1, n2), intersect(n3, n4)))

grid.newpage()

Vps4A <- draw.quad.venn(area1 = nrow(n1), area2 = nrow(n2), area3 = nrow(n3), area4 = nrow(n4),
                        n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234,
                        category = c("siVPS4A#1/NT", "siVPS4A#1/siCtrl",
                                     "siVps4A#2/NT", "siVPS4A#2/siCtrl"), 
                        fill = c("dodgerblue", "deeppink", "goldenrod1", "chartreuse4"),
                        fontfamily = "sans", cat.fontfamily = "sans",
                        cex = 2, cat.cex = 2,
                        cat.dist = c(0.25, 0.25, 0.15, 0.15),
                        cat.just = list(c(0.35,0.35) , c(0.70,0.35) , c(0.35,1.0) , c(0.60,1.0)))

Vps4A


### Vps4B

Vps4B_1_NT <- filter(df_Vps4B, siVps4B_72.vs.NT_FC >= UP, siVps4B_72.vs.NT_adj.pval < pval)
Vps4B_1_Ctrl1 <- filter(df_Vps4B, siVps4B_72.vs.siCtrl_1_FC >= UP, siVps4B_72.vs.siCtrl_1_adj.pval < pval)
Vps4B_2_NT <- filter(df_Vps4B, siVps4B_73.vs.NT_FC >= UP, siVps4B_73.vs.NT_adj.pval < pval)
Vps4B_2_Ctrl1 <- filter(df_Vps4B, siVps4B_73.vs.siCtrl_1_FC >= UP, siVps4B_73.vs.siCtrl_1_adj.pval < pval)

n1 <- Vps4B_1_NT[ , 1]
n2 <- Vps4B_1_Ctrl1[ , 1]
n3 <- Vps4B_2_NT[ , 1]
n4 <- Vps4B_2_Ctrl1[ , 1]

n12 <- nrow(intersect(n1, n2))
n34 <- nrow(intersect(n3, n4))
n13 <- nrow(intersect(n1, n3))
n14 <- nrow(intersect(n1, n4))
n23 <- nrow(intersect(n2, n3))
n24 <- nrow(intersect(n2, n4))

n123 <- nrow(Reduce(intersect, list(n1, n2, n3)))
n124 <- nrow(Reduce(intersect, list(n1, n2, n4)))
n134 <- nrow(Reduce(intersect, list(n1, n3, n4)))
n234 <- nrow(Reduce(intersect, list(n2, n3, n4)))

n1234 <- nrow(intersect(intersect(n1, n2), intersect(n3, n4)))

grid.newpage()

Vps4B <- draw.quad.venn(area1 = nrow(n1), area2 = nrow(n2), area3 = nrow(n3), area4 = nrow(n4),
                        n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234,
                        category = c("siVPS4A#1/NT", "siVPS4A#1/siCtrl",
                                     "siVps4A#2/NT", "siVPS4A#2/siCtrl"), 
                        fill = c("dodgerblue", "deeppink", "goldenrod1", "chartreuse4"),
                        fontfamily = "sans", cat.fontfamily = "sans",
                        cex = 2, cat.cex = 2,
                        cat.dist = c(0.25, 0.25, 0.15, 0.15),
                        cat.just = list(c(0.35,0.35) , c(0.70,0.35) , c(0.35,1.0) , c(0.60,1.0)))

Vps4B


# Reactome Pathway

# source("https://bioconductor.org/biocLite.R")
# biocLite("ReactomePA")

library(ReactomePA)

Vps4AB_pathway <- enrichPathway(gene = Vps4AB_EntrezID_UP$ENTREZID, 
                                universe = Vps4AB_EntrezID_project$ENTREZID,
                                organism = "human", 
                                pAdjustMethod = "BH", 
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05, 
                                minGSSize = 10,
                                maxGSSize = 500, 
                                readable = TRUE)



Pathway_UP <- as.data.frame(Vps4AB_pathway) %>% mutate(GeneRatio = Count/342)

path <- ggplot(Pathway_UP, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) + #fct_reorder
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 12) + # size of description
  scale_colour_gradient(limits=c(10^(-10), 0.05), low="red", high = "blue") +
  ylab(NULL) + xlab("Gene ratio") + 
  labs(colour = "p-value", size = "Gene count") + 
  ggtitle("Pathways Enrichment Amongst Gene Clusters") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold", color = "black")) # size of title

path

# element_text(family = NULL, face = NULL, colour = NULL, size = NULL,
#              hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
#              color = NULL)

# theme_bw(base_size = 11, base_family = "",
#          base_line_size = base_size/22, base_rect_size = base_size/22)


# write.table(as.data.frame(Vps4AB_pathway),
#             file = "Vps4AB_pathway.txt",
#             quote = FALSE,
#             sep = "\t",
#             row.names = TRUE)
