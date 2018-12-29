#!/usr/bin/env Rscript

#'
#' Make circos plot and spider plot
#' 
#' Author: Mingyang Cai
#' 

source("spider_funcs.R")
library(RCircos)

# make circos plot
make_circos = function(dat, outfile) {
  data(UCSC.HG19.Human.CytoBandIdeogram)
  cyto.info = UCSC.HG19.Human.CytoBandIdeogram
  
  circos = dat[1:3]
  circos$V4 = "chr12"
  circos$V5 = 53301238
  circos$V6 = 53301985
  
  RCircos.Set.Core.Components(cyto.info,chr.exclude=NULL, tracks.inside=5, tracks.outside=0)
  pdf(file=outfile,height=8,width=8)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  RCircos.Link.Plot(circos,track.num=2,by.chromosome=TRUE)
  dev.off()
}

# prepare gene data
read_genes = function(file) {
  dat = read.table(file, header=F)
  dat[c(2:3,6)]
}

# make spider plot
make_spider = function(dat, gene_fle, outfile) {
  sub_dat <- subset(dat, dat$V1=="chr12")
  dom = sub_dat[2:3]
  genes<- read_genes(gene_file)
  
  pdf(outfile, 26, 3)
  makeSpiderGramSingle(dom=dom, vp.loc=53301238, color='orange', gene=genes, chrom.len=133851895)
  dev.off()
}


####################
# main 
####################
cov_file = "/Users/Mingyang/Google Drive/Lu_lab/PCa/captureC_clean/rs55958994_multiinter_cov_genes.txt"
dat = read.table(cov_file)
sub_dat = subset(dat, (dat$V6>0 & dat$V7>0 | dat$V6>0 & dat$V8>0 | dat$V7>0 & dat$V8>0) & dat$V1 != "chrM")
cat(sprintf("Among those doubly reproducible sites, %s are located in cis, perc: %s", sum(sub_dat$V1=="chr12"), sum(sub_dat$V1=="chr12")/dim(sub_dat)[1]))
sub_sub_dat = subset(dat, dat$V6>0 & dat$V7>0 & dat$V8>0 & dat$V1 != "chrM")
cat(sprintf("Among those triply reproducible sites, %s are located in cis, perc: %s", sum(sub_sub_dat$V1=="chr12"), sum(sub_sub_dat$V1=="chr12")/dim(sub_sub_dat)[1]))

make_circos(sub_dat, "circos.pdf")
gene_file = "/Users/Mingyang/Google Drive/Lu_lab/Shi/RNA-Seq-CaptureC-2018-May/chr12_genes"
make_spider(sub_dat, gene_file, "spider.pdf")
