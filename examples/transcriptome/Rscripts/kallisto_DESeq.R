library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(fdrtool)
library(geneplotter)
#library(EDASeq)
#library(apeglm)
library(rhdf5)

samples <- read.csv("samples.csv",header=TRUE)
exprnames <- do.call(paste,c(samples[c("Condition","Replicate")],sep="."))
exprnames <- sub(".([123])$",".r\\1",exprnames,perl=TRUE)
print(exprnames)

files <- file.path("results","kallisto",exprnames,"abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- exprnames
colnames(txi.kallisto$abundance) <- exprnames
write.csv(txi.kallisto$abundance,"reports/kallisto.TPM.csv")
write.csv(txi.kallisto$counts,"reports/kallisto.counts.csv")

# DEseq2 analyses

cond = factor( samples$Condition)
rep = factor( samples$Replicate)
#treatment = factor (samples$Condition)

#condition = treatment,
sampleTable <- data.frame(condition = cond, replicate = rep)
rownames(sampleTable) = exprnames
sampleTable



dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable, ~ condition)
dds$condition <- relevel(dds$condition, ref = "3hr")
#dds$condition <- factor(dds$condition, levels = c("24hr")

dds$condition
#nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
#nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/DESeq.pdf")

plotDispEsts(dds)

multidensity( counts(dds, normalized = T),
              xlab="mean counts", xlim=c(0, 1000))
multiecdf( counts(dds, normalized = T),
           xlab="mean counts", xlim=c(0, 1000))

#MA.idx = t(combn(1:4, 2))
#for( i in  seq_along( MA.idx[,1])){ 
#  MDPlot(counts(dds, normalized = T), 
#         c(MA.idx[i,1],MA.idx[i,2]), 
#         main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
#                       colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
#}

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]


df2 <- as.data.frame(colData(dds)[,c("condition")])
rownames(df2) = exprnames
colnames(df2) = c("Condition")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="VSD")

topVar <- head(order(rowVars(assay(vsd)),
          decreasing=TRUE),60)
mat  <- assay(vsd)[ topVar, ]
mat  <- mat - rowMeans(mat)

mat2  <- assay(vsd)[ topVar, ]

pheatmap(mat, show_rownames=TRUE,
        fontsize_row = 7,fontsize_col = 7,
        cluster_cols=FALSE, annotation_col=df2,main="VSD Most different")

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD")

topVar <- head(order(rowVars(assay(rld)),
    decreasing=TRUE),60)
mat  <- assay(rld)[ topVar, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, show_rownames=TRUE,
         fontsize_row = 7,fontsize_col = 7,
         cluster_cols=FALSE, annotation_col=df2,main="RLD Most different")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()

dds <- DESeq(dds)
res <- results(dds)
res
resOrdered <- res[order(res$pvalue,decreasing=TRUE),]

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


resTop <- subset(res,abs(res$log2FoldChange) > 4 & res$padj < 0.01)
topVar <- head(resTop[order(resTop$log2FoldChange,decreasing=FALSE),],60)
mat  <- assay(rld)[ rownames(topVar), ]

#reorder the columns
pheatmap(mat,method="complete",main = "4 fold up/down, p < 0.01 using RLD ordered by fold change", show_rownames = T, show_colnames=T,
         annotation_legend = FALSE, legend=T, cluster_rows=FALSE, cexRow=0.3,
         fontsize_row = 7,fontsize_col = 7)
pheatmap(mat,method="complete",main = "4 fold up/down, p < 0.01 using RLD, clustered by expression", 
         show_rownames = T, show_colnames=T,
         annotation_legend = FALSE, legend=T, cluster_rows=TRUE, cexRow=0.3,
         fontsize_row = 7,fontsize_col = 7)

topChange <- head(resOrdered[order(resOrdered$log2FoldChange,decreasing=TRUE),],100)
topChangePlot =  data.frame(log2FoldChange = topChange$log2FoldChange)
rownames(topChangePlot) = rownames(topChange)
pheatmap(topChangePlot,method="complete",main = "Fold Change Plot Up", show_rownames = T,show_colnames=F,
         cluster_rows=FALSE,cluster_cols=FALSE,cexRow=0.3,legend=T, cexRow=0.3,
         fontsize_row = 6)

topChange

topChange <- head(resOrdered[order(resOrdered$log2FoldChange,decreasing=FALSE),],100)
topChangePlot =  data.frame(log2FoldChange = topChange$log2FoldChange)
rownames(topChangePlot) = rownames(topChange)
pheatmap(topChangePlot,method="complete",main = "Fold Change Plot Down", show_rownames = T,show_colnames=F,
         cluster_rows=FALSE,cluster_cols=FALSE,cexRow=0.3,legend=T, cexRow=0.3,
         fontsize_row = 6)


dev.off()
