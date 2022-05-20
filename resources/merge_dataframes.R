#!/usr/bin/env Rscript --vanilla

# suppressPackageStartupMessages(library("argparse"))
# 
# # create parser object
# parser <- ArgumentParser()
# 
# # specify our desired options 
# # by default ArgumentParser will add an help option 
# parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
#                     help="Print extra output [default]")
# parser$add_argument("--df1", 
#                     dest="df1", help="dataframe1")
# parser$add_argument("--df1_colname", 
#                     help="dataframe1 columnname to merge by")
# parser$add_argument("--df2", 
#                     dest="df2", help="dataframe2")
# parser$add_argument("--df2_colname", 
#                     help="dataframe2 columnname to merge by")
# parser$add_argument("--out", 
#                     dest="out", help="out filename")
# 
# # get command line options, if help option encountered print help and exit,
# # otherwise if options not found on command line then set defaults, 
# args <- parser$parse_args()

setwd("~/Ziegelbauer_lab/circRNADetection/scripts/circRNA/resources")
mm39=read.csv("mm39.circBase.collapsed.bed",
                      header=FALSE,
                      sep="\t",
                      comment.char = "#",
                      check.names = FALSE)
colnames(mm39)=c("chrom","start","end","name","score","strand")
# mm39$mm39ID=paste0(mm39$chrom,"##",mm39$start,"##",mm39$end,"##",mm39$strand)
mid2gene=read.csv("mID_to_gene.tsv",
                  header=FALSE,
                  sep="\t",
                  comment.char = "#",
                  check.names = FALSE)
colnames(mid2gene)=c("mID","NMid","gene_name")
df=merge(mm39,mid2gene,by.x="name",by.y="mID",all.x = TRUE)
mmu_mm9_circRNA_txt=read.csv("mmu_mm9_circRNA.txt",
                         header=TRUE,
                         sep="\t",
                         check.names = FALSE)
mmu_mm9_circRNA_txt=mmu_mm9_circRNA_txt[,-1]
mmu_mm9_circRNA_txt=mmu_mm9_circRNA_txt[,-1]
mmu_mm9_circRNA_txt=mmu_mm9_circRNA_txt[,-1]
mmu_mm9_circRNA_txt=mmu_mm9_circRNA_txt[,-1]
df=merge(df,mmu_mm9_circRNA_txt,by.x="name",by.y="circRNA ID",all.x = TRUE)
mm9_circID_to_name=read.csv("mm9_circID_to_name.txt",
                            header=TRUE,
                            sep="\t",
                            check.names = FALSE)
colnames(mm9_circID_to_name)=c("circID","mm9_name")
df=merge(df,mm9_circID_to_name,by.x="name",by.y="circID",all.x=TRUE)
