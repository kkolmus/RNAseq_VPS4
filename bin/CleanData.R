rm(list = ls())

############
### LIBS ###

# ipak function: install and load multiple R packages.
# credits: https://gist.github.com/stevenworthington
# Check to see if packages are installed. 
# Install them if they are not, then load them into the R session.

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("readxl", "WriteXLS", "dplyr", "futile.logger")
ipak(packages)

############
### MAIN ###

flog.threshold(DEBUG)

flog.debug("Install and load the DESeq2 package")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
suppressPackageStartupMessages(library(DESeq2))


flog.debug("Setting up working directory for the analysis")
data.dir = "~/Desktop/HT projects/RNAseq Vps4"
proj.dir = "~/Desktop/HT projects/RNAseq_VPS4"
setwd(proj.dir)

flog.debug("Prepare count matrix")
DF <- read_excel(file.path(data.dir, "rawReads.Vps4.18_03_17_cleaned.xlsx"))
# remove library ctrls, statistics and the gene symbol column
DF <- DF[-c(1:10,20813:20817), -2]
# fill up missing values
DF[is.na(DF)] <- 0
# remove duplicates
DF <- DF[!duplicated(DF$gene_id),] # no duplicates

# assign column names
new_colname <- c("gene_id",
                 "1_NT", "1_siCtrl_1", "1_siCtrl_2",
                 "1_siVps4A_66", "1_siVps4A_68", "1_siVps4B_72", "1_siVps4B_73",
                 "1_siVps4AB_68_72", "1_siVps4AB_66_73",
                 "2_NT", "2_siCtrl_1", "2_siCtrl_2",
                 "2_siVps4A_66", "2_siVps4A_68", "2_siVps4B_72", "2_siVps4B_73",
                 "2_siVps4AB_68_72", "2_siVps4AB_66_73",
                 "3_NT", "3_siCtrl_1", "3_siCtrl_2",
                 "3_siVps4A_66", "3_siVps4A_68", "3_siVps4B_72", "3_siVps4B_73",
                 "3_siVps4AB_68_72", "3_siVps4AB_66_73",
                 "4_NT", "4_siCtrl_1", "5_siCtrl_1", "4_siCtrl_2",
                 "4_siVps4A_66", "5_siVps4A_66", "4_siVps4A_68", "4_siVps4B_72", "4_siVps4B_73",
                 "4_siVps4AB_68_72", "4_siVps4AB_66_73")

colnames(DF) <- new_colname

# exporting data frame
write.csv(DF, "cts.csv", row.names = FALSE)

# loading data on experimental design
coldata <- read_excel("coldata.xlsx")

# exporting data frame
write.csv(coldata, "coldata.csv", row.names = FALSE)

# loading data frames for DESeq2 analysis
cts <- as.matrix(read.csv("cts.csv", row.names = 1))

coldata <- read.csv("coldata.csv", row.names = 1)
coldata$batch <- as.factor(coldata$batch)
coldata <- coldata[, c("batch", "condition")]

# confirming coherence

colnames(cts) <- sub("X", "", colnames(cts))
all(rownames(coldata) %in% colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# removed columns for reanalysis

DF2 <- DF[, -c(2,3,6,14)]

write.csv(DF2, "cts2.csv", row.names = FALSE)

cts2 <- as.matrix(read.csv("cts2.csv", row.names = 1))

coldata2 <- coldata[-c(1,2,5,13),]

write.csv(coldata2, "coldata2.csv", row.names = TRUE)


# confirming coherence

colnames(cts2) <- sub("X", "", colnames(cts2))
all(rownames(coldata2) %in% colnames(cts2))

cts2 <- cts2[, rownames(coldata2)]
all(rownames(coldata2) == colnames(cts2))