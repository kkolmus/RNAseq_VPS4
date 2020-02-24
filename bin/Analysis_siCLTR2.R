# cleaning
rm(list = ls())

############
### LIB ###

library(futile.logger)
library(readxl)
library(WriteXLS)
library(dplyr)
library(DESeq2)


############
### MAIN ###

flog.threshold(DEBUG)

flog.debug("Set project directory")
proj.dir = "~/Desktop/HT projects/RNAseq_VPS4/data"

cts <- as.matrix(read.csv("cts2.csv", row.names = 1))

coldata <- read.csv("coldata2.csv", row.names = 1)
coldata$batch <- as.factor(coldata$batch)
coldata <- coldata[, c("batch", "condition")]

colnames(cts) <- sub("X", "", colnames(cts))
all(rownames(coldata) %in% colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# analysis with DESeq2 taking into account batch and condition
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ batch + condition)
dds

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Note on factor levels
dds$condition <- relevel(dds$condition, ref = "siCtrl_2")

# Differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

# Extracting transformed values

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# Effects of transformations on the variance

# this gives log2(n + 1)
ntd <- normTransform(dds)

par(mfrow = c(3,1))

meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))


flog.debug("Heatmap of the count matrix")

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("batch", "condition")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


flog.debug("Heatmap of the sample-to-sample distances")

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$batch, vsd$condition, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$batch, vsd$condition, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col = colors)


flog.debug("Principal component plot of the samples")

plotPCA(vsd, intgroup=c("batch", "condition"))

pcaData <- plotPCA(vsd, intgroup=c("batch", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  ggtitle("PCA Vps4 project data") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


flog.debug("Contrasts for differential expression analysis")

NT.vs.siCtrl_2 <- results(dds, 
                          contrast = c("condition", "NT", "siCtrl_2"),
                          pAdjustMethod = "BH",
                          test = "Wald", alpha = 0.05)

write.table(as.data.frame(NT.vs.siCtrl_2),
            file = "NT.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

siCtrl_2.vs.siCtrl_2 <- results(dds, 
                                contrast = c("condition", "siCtrl_2", "siCtrl_2"),
                                pAdjustMethod = "BH",
                                test = "Wald", alpha = 0.05)

write.table(as.data.frame(siCtrl_2.vs.siCtrl_2),
            file = "siCtrl_2.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

siVps4A_66.vs.siCtrl_2 <- results(dds, 
                                  contrast = c("condition", "siVps4A_66", "siCtrl_2"),
                                  pAdjustMethod = "BH",
                                  test = "Wald", alpha = 0.05)

write.table(as.data.frame(siVps4A_66.vs.siCtrl_2),
            file = "siVps4A_66.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

siVps4A_68.vs.siCtrl_2 <- results(dds, 
                                  contrast = c("condition", "siVps4A_68", "siCtrl_2"),
                                  pAdjustMethod = "BH",
                                  test = "Wald", alpha = 0.05)

write.table(as.data.frame(siVps4A_68.vs.siCtrl_2),
            file = "siVps4A_68.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

siVps4B_72.vs.siCtrl_2 <- results(dds, 
                                  contrast = c("condition", "siVps4B_72", "siCtrl_2"),
                                  pAdjustMethod = "BH",
                                  test = "Wald", alpha = 0.05)

write.table(as.data.frame(siVps4B_72.vs.siCtrl_2),
            file = "siVps4B_72.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

siVps4B_73.vs.siCtrl_2 <- results(dds, 
                                  contrast = c("condition", "siVps4B_73", "siCtrl_2"),
                                  pAdjustMethod = "BH",
                                  test = "Wald", alpha = 0.05)

write.table(as.data.frame(siVps4B_73.vs.siCtrl_2),
            file = "siVps4B_73.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)


siVps4AB_68_72.vs.siCtrl_2 <- results(dds, 
                                      contrast = c("condition", "siVps4AB_68_72", "siCtrl_2"),
                                      pAdjustMethod = "BH",
                                      test = "Wald", alpha = 0.05)

write.table(as.data.frame(siVps4AB_68_72.vs.siCtrl_2),
            file = "siVps4AB_68_72.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)

siVps4AB_66_73.vs.siCtrl_2 <- results(dds, 
                                      contrast = c("condition", "siVps4AB_66_73", "siCtrl_2"),
                                      pAdjustMethod = "BH",
                                      test = "Wald", alpha = 0.05)

write.table(as.data.frame(siVps4AB_66_73.vs.siCtrl_2),
            file = "siVps4AB_66_73.vs.siCtrl_2.txt",
            quote = FALSE,
            sep = "\t",
            row.names = TRUE)


###################
### SessionInfo ###

sessionInfo()