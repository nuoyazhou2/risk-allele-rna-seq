#!/usr/bin/env Rscript

#'
#' Combine RNA-Seq and capture-C results
#' 
#' Author: Mingyang Cai
#' 

library(ggplot2)
library(ggrepel)
library(dplyr)

set.seed(1234)
file = "/Users/Mingyang/Google Drive/Lu_lab/PCa/plot_clean/22Rv1_only/cis_trans/prep_rs55958994.txt"
candidate_genes = c("CDH23","SIPA1","CNTN1","KRT8","FAIM2","ITGA5","KRT18")
cutoff = log2(1.5)

# read in file
read_dat = function(file) {
  df = read.table(file, as.is=T, header=F)
  names(df) = c("gene", "count_captureC", "log2FoldChange_RNAseq", "adjusted-p")
  df = unique(df)
  
  dat = df %>%
    group_by(gene, log2FoldChange_RNAseq, `adjusted-p`) %>%
    summarize(capturec_value=max(count_captureC, na.rm=T))
  
  dat = dat %>%
    mutate(significant = capturec_value > cutoff & log2FoldChange_RNAseq > cutoff | capturec_value > cutoff & log2FoldChange_RNAseq < -cutoff)
  dat
}

dat = read_dat(file)

p_inset = ggplot(dat, aes(capturec_value, log2FoldChange_RNAseq)) + 
  geom_point(aes(size=-log10(`adjusted-p`), color=significant), alpha=0.5) + 
  geom_label_repel(aes(label=ifelse(gene %in% candidate_genes, gene, "")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_bw() + 
  xlab("") + 
  ylab("") + 
  xlim(0,5) + ylim(-7,-1) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = cutoff, linetype="dashed") + 
  geom_hline(yintercept = c(-cutoff, cutoff), linetype="dashed") + 
  theme(axis.title=element_text(size=18, face="bold"),
        axis.text=element_text(size=15),
        legend.position="none")

g = ggplotGrob(p_inset)

p = ggplot(dat, aes(capturec_value, log2FoldChange_RNAseq)) + 
  geom_point(aes(size=-log10(`adjusted-p`),color=significant), alpha=0.5) + 
  geom_label_repel(aes(label=ifelse(gene %in% candidate_genes, gene, "")),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_bw() + 
  xlab("capture-C counts") + 
  ylab("log2 fold change KO/WT") + 
  xlim(0,70) + ylim(-8,8) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = cutoff, linetype="dashed") + 
  geom_hline(yintercept = c(-cutoff, cutoff), linetype="dashed") + 
  theme(plot.margin=unit(c(0.5,0.5,0.25,0.25),"in"),
    axis.title=element_text(size=18, face="bold"),
    axis.text=element_text(size=15), 
    legend.title=element_text(size=12), 
    legend.text=element_text(size=10)) +
  annotation_custom(
    grob = g,
    xmin = 20,
    xmax = Inf,
    ymin = 1,
    ymax = Inf
  )

pdf("capture-C_versus_RNA-Seq.pdf", width = 10, height = 8)
p
dev.off()