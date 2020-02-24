# cleaning
rm(list = ls())

# Normalisation vs. NT
setwd("~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis/NT")

files <- list.files(path="~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis/NT")
genes <- read.table(files[1], row.names = 1, header=TRUE, sep="\t")[,0] # gene names
df <- do.call(cbind,
              lapply(files, function(fn) 
                read.table(fn, header=TRUE, sep="\t")[,c(2,6)]))
prefix <- gsub(".txt", "", files)
suffix <- c("FC", "adj.pval")
new_col_names <- paste(rep(prefix, each = 2), suffix, sep = "_")
colnames(df) <- new_col_names
DF_NT <- cbind(genes,df)

setwd("~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis")
write.table(DF_NT, "DF_NT.txt", quote = F, sep = "\t", row.names = T)

# cleaning
rm(list = ls())

# Normalisation vs. siCtrl1
setwd("~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis/siCtrl_1")

files <- list.files(path="~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis/siCtrl_1")
genes <- read.table(files[1], row.names = 1, header=TRUE, sep="\t")[,0] # gene names
df <- do.call(cbind,
              lapply(files, function(fn) 
                read.table(fn, header=TRUE, sep="\t")[,c(2,6)]))
prefix <- gsub(".txt", "", files)
suffix <- c("FC", "adj.pval")
new_col_names <- paste(rep(prefix, each = 2), suffix, sep = "_")
colnames(df) <- new_col_names
DF_siCtrl1 <- cbind(genes,df)

setwd("~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis")
write.table(DF_siCtrl1, "DF_siCtrl1.txt", quote = F, sep = "\t", row.names = T)

# cleaning
rm(list = ls())

# Normalisation vs. siCtrl2
setwd("~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis/siCtrl_2")

files <- list.files(path="~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis/siCtrl_2")
genes <- read.table(files[1], row.names = 1, header=TRUE, sep="\t")[,0] # gene names
df <- do.call(cbind,
              lapply(files, function(fn) 
                read.table(fn, header=TRUE, sep="\t")[,c(2,6)]))
prefix <- gsub(".txt", "", files)
suffix <- c("FC", "adj.pval")
new_col_names <- paste(rep(prefix, each = 2), suffix, sep = "_")
colnames(df) <- new_col_names
DF_siCtrl2 <- cbind(genes,df)

setwd("~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis")
write.table(DF_siCtrl2, "DF_siCtrl2.txt", quote = F, sep = "\t", row.names = T)


# cleaning
rm(list = ls())

# manually edit datasheets by shifting 1st row to the right and adding to the A1 position RefSeq
# save files with the xlsx extension

setwd("~/Desktop/RNAseq Vps4/2018.11.07 R 3.4.4 Complete reanalysis")

library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
# install.packages("imputeTS")
library(imputeTS)

# merging data

DF_NT <- read_excel("DF_NT.xlsx")
DF_siCtrl1 <- read_excel("DF_siCtrl1.xlsx")

DF <- merge(DF_NT, DF_siCtrl1, by = "RefSeq", all = TRUE)

DF_siCtrl2 <- read_excel("DF_siCtrl2.xlsx")

DF <- merge(DF, DF_siCtrl2, by = "RefSeq", all = TRUE)


# add gene symbols to RefSeq

DF_Symbol <- bitr(DF$RefSeq, fromType = "REFSEQ", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

new_col_names2 <- c("RefSeq", "Symbol")

colnames(DF_Symbol) <- new_col_names2

DF <- merge(DF, DF_Symbol, by = "RefSeq", all = TRUE)
DF <- DF[, c(1,50,2:49)]

# change class of columns

DF <- mutate_at(DF, vars(-c(RefSeq, Symbol)), funs(as.numeric(.)))

# change log2FC to FC

DF2 <- mutate_at(DF, vars(-c(RefSeq, Symbol,
                             siCtrl_1.vs.NT_adj.pval,
                             siCtrl_2.vs.NT_adj.pval,
                             siVps4A_66.vs.NT_adj.pval,
                             siVps4A_68.vs.NT_adj.pval,
                             siVps4AB_66_73.vs.NT_adj.pval,
                             siVps4AB_68_72.vs.NT_adj.pval,
                             siVps4B_72.vs.NT_adj.pval,
                             siVps4B_73.vs.NT_adj.pval,
                             NT.vs.siCtrl_1_adj.pval,
                             siCtrl_2.vs.siCtrl_1_adj.pval,
                             siVps4A_66.vs.siCtrl_1_adj.pval,
                             siVps4A_68.vs.siCtrl_1_adj.pval,
                             siVps4AB_66_73.vs.siCtrl_1_adj.pval,
                             siVps4AB_68_72.vs.siCtrl_1_adj.pval,
                             siVps4B_72.vs.siCtrl_1_adj.pval,
                             siVps4B_73.vs.siCtrl_1_adj.pval,
                             NT.vs.siCtrl_2_adj.pval,
                             siCtrl_1.vs.siCtrl_2_adj.pval,
                             siVps4A_66.vs.siCtrl_2_adj.pval,
                             siVps4A_68.vs.siCtrl_2_adj.pval,
                             siVps4AB_66_73.vs.siCtrl_2_adj.pval,
                             siVps4AB_68_72.vs.siCtrl_2_adj.pval,
                             siVps4B_72.vs.siCtrl_2_adj.pval,
                             siVps4B_73.vs.siCtrl_2_adj.pval
)), funs(2^(.)))

# round up numbers

DF2 <- mutate_at(DF2, vars(-c(RefSeq, Symbol,
                              siCtrl_1.vs.NT_adj.pval,
                              siCtrl_2.vs.NT_adj.pval,
                              siVps4A_66.vs.NT_adj.pval,
                              siVps4A_68.vs.NT_adj.pval,
                              siVps4AB_66_73.vs.NT_adj.pval,
                              siVps4AB_68_72.vs.NT_adj.pval,
                              siVps4B_72.vs.NT_adj.pval,
                              siVps4B_73.vs.NT_adj.pval,
                              NT.vs.siCtrl_1_adj.pval,
                              siCtrl_2.vs.siCtrl_1_adj.pval,
                              siVps4A_66.vs.siCtrl_1_adj.pval,
                              siVps4A_68.vs.siCtrl_1_adj.pval,
                              siVps4AB_66_73.vs.siCtrl_1_adj.pval,
                              siVps4AB_68_72.vs.siCtrl_1_adj.pval,
                              siVps4B_72.vs.siCtrl_1_adj.pval,
                              siVps4B_73.vs.siCtrl_1_adj.pval,
                              NT.vs.siCtrl_2_adj.pval,
                              siCtrl_1.vs.siCtrl_2_adj.pval,
                              siVps4A_66.vs.siCtrl_2_adj.pval,
                              siVps4A_68.vs.siCtrl_2_adj.pval,
                              siVps4AB_66_73.vs.siCtrl_2_adj.pval,
                              siVps4AB_68_72.vs.siCtrl_2_adj.pval,
                              siVps4B_72.vs.siCtrl_2_adj.pval,
                              siVps4B_73.vs.siCtrl_2_adj.pval)), funs(round(., 3)))

# add avaeraged number of reads per gene

reads <- read.csv("cts2.csv", header = TRUE)
reads <- dplyr::rename(reads, RefSeq = gene_id)
reads <- mutate(reads, 
                NT_reads = (reads$X2_NT + reads$X3_NT + reads$X4_NT)/3,
                siCtrl1_reads = (reads$X2_siCtrl_1 + reads$X3_siCtrl_1 + reads$X4_siCtrl_1)/3,
                siCtrl2_reads = (reads$X1_siCtrl_2 + reads$X2_siCtrl_2 + reads$X3_siCtrl_2 + reads$X4_siCtrl_2)/4,
                siVps4A_66_reads = (reads$X1_siVps4A_66 + reads$X3_siVps4A_66 + reads$X4_siVps4A_66 + reads$X5_siVps4A_66)/4,
                siVps4A_68_reads = (reads$X2_siVps4A_68 + reads$X3_siVps4A_68 + reads$X4_siVps4A_68)/3,
                siVps4B_72_reads = (reads$X1_siVps4B_72 + reads$X2_siVps4B_72 + reads$X3_siVps4B_72 + reads$X4_siVps4B_72)/4,
                siVps4B_73_reads = (reads$X1_siVps4B_73 + reads$X2_siVps4B_73 + reads$X3_siVps4B_73 + reads$X4_siVps4B_73)/4,
                siVps4AB_68_72_reads = (reads$X1_siVps4AB_68_72 + reads$X2_siVps4AB_68_72 + reads$X3_siVps4AB_68_72 + reads$X4_siVps4AB_68_72)/4,
                siVps4AB_66_73_reads = (reads$X1_siVps4AB_66_73 + reads$X2_siVps4AB_66_73 + reads$X3_siVps4AB_66_73 + reads$X4_siVps4AB_66_73)/4)
reads <- reads[, c(1,36:44)]
reads <- mutate_at(reads, vars(-c(RefSeq)), funs(round(., digits = 3)))

DF2 <- left_join(DF2, reads, by = "RefSeq")

write.table(DF2, "ES_Vps4_project_full.txt", quote = F, sep = "\t", row.names = F)

DF2 <- dplyr::mutate(DF2, readsSum = rowSums(DF2[, c(51:59)]))

DF_filtered <- filter(DF2, readsSum >= 10)

WriteXLS(DF_filtered, "ES_Vps4_project_full.xlsx")

DF_filtered_over100 <- filter(DF2, readsSum >= 100)

WriteXLS(DF_filtered_over100, "ES_Vps4_project_full_over100.xlsx")
