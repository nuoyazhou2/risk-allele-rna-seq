#!/usr/bin/env Rscript

#'
#' Compare the effect of *risk allele/major allele* versus *KO/WT* in prostate cancer
#' 
#' Author: Mingyang Cai
#' 

library(ggplot2)

# read in DESeq2 output file
read_df = function(file, suffix, pval_cutoff=0.1, fc_cutoff=1.5) {
  df = read.delim(file, as.is=T)
  up_df = subset(df, (padj<pval_cutoff) & (log2FoldChange>log2(fc_cutoff)))
  names(up_df) = paste0(names(up_df), suffix)
  down_df = subset(df, (padj<pval_cutoff) & (log2FoldChange<log2(1/fc_cutoff)))
  names(down_df) = paste0(names(down_df), suffix)
  list(up_df, down_df)
}

# main
risk_file = "./DESeq2_T_vs_C_.txt"
ko_file = "./DESeq2_KO_vs_WT_.txt"
risk_df = read_df(risk_file, "_risk")
ko_df = read_df(ko_file, "_ko")
risk_up_df = risk_df[[1]]
risk_down_df = risk_df[[2]]
ko_up_df = ko_df[[1]]
ko_down_df = ko_df[[2]]

# print gene names
# cat(paste(row.names(risk_up_df), collapse="\n"))
# cat(paste(row.names(risk_down_df), collapse="\n"))

dat = merge(risk_up_df, ko_down_df, by=0)
names(dat)[which(names(dat)=="Row.names")] = "gene"
candidate_genes = c("CDH23","SIPA1","CNTN1","KRT8","FAIM2","ITGA5","KRT18")
p = ggplot(dat, aes(x=log2FoldChange_risk, y=log2FoldChange_ko)) + 
  geom_point(alpha=0.5, size=3, color="red") +
  geom_text(aes(label=ifelse(gene %in% candidate_genes, gene, "")), hjust=1, vjust=0, fontface="bold") +
  xlab("log2 fold change risk allele(T)/major allele(C)") +
  ylab("log2 fold change KO/WT") +
  theme_bw() +
  theme(
    plot.margin=unit(c(0.5,0.5,0.25,0.25),"in"),
    axis.text=element_text(size=15),
    axis.title=element_text(size=18, face="bold")
  )
pdf("risk_and_ko_comparison.pdf",8,8)
print(p)
dev.off()
