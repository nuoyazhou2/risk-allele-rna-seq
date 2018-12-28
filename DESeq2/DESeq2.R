#!/usr/usc/R/3.3.2/bin/R

library('DESeq2')
library("ggplot2")
library("RColorBrewer")
library("pheatmap")
library("stringr")

dir = "/staging/wl/caim/analysis/22Rv1_allele_RNA_Seq"
sampleDirs = list.files(dir)
sampleFiles = sapply(sampleDirs, function(x) paste0(dir,"/",x,"/htseq_",x,".txt"))
sampleCondition = str_match(sampleDirs, ".*-[0-9]+([C|T])_.*")[,2]
sampleName = sampleDirs
sampleTable <- data.frame(sampleName = sampleName, 
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "/",
                                       design= ~ condition)
#ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ] #pre-filtering to remove rows with only 0 or 1 read
# note that the "results" function below automatically performs independent filtering, so no need to do the filtering as above.
#ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref="NEGATIVE")
analysis = DESeq(ddsHTSeq,betaPrior=FALSE)
resultsNames(analysis)

write.table(as.data.frame(counts(analysis,normalized=TRUE)),file="normalized_counts.txt",quote=FALSE,row.names=TRUE,sep="\t")
write.table(as.data.frame(counts(analysis)),file="raw_counts.txt",quote=FALSE,row.names=TRUE,sep="\t")
write.table(as.data.frame(results(analysis, contrast=c("condition","T","C"))),file="DESeq2_T_vs_C_.txt",quote=FALSE,row.names=TRUE,sep="\t")

rld <- rlog(ddsHTSeq, blind=FALSE)
write.table(assay(rld),file="rld_transformed_normalized_counts.txt",quote=FALSE,row.names=TRUE,sep="\t")

save.image(file="DESeq2.RData")


# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)), method="euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- colnames(rld)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("sample-to-sample distances (euclidean_complete).pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, clustering_method = "complete")
dev.off()

pdf("sample-to-sample distances (euclidean average).pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, clustering_method = "average")
dev.off()


# Principal component plot of the samples
pdf("PCA.pdf",height=3,width=6)
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
p = ggplot(data, aes(PC1, PC2, color=condition)) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(axis.title=element_text(size=14)) + 
  theme_bw() + 
  scale_color_manual( values=c("#E69F00", "#56B4E9")) +
  labs(color='allele')
p
dev.off()

