rm(list=ls())
library(ggVennDiagram)
library(ggplot2)
library(argparse)


# setwd("/Users/kopardevn/Desktop/Temp/Venn2")

parser <- ArgumentParser()

parser$add_argument("-i", "--ciriout", type="character", required=TRUE,
                    help="ciri2 ciri.out outputfile")

parser$add_argument("-e", "--circExplorerout", type="character", required=TRUE,
                    help="circExporer2 circularRNA_known.txt outputfile")

parser$add_argument("-p", "--plot", type="character", required=TRUE,
                    help="output Venn PNG file")

parser$add_argument("-l", "--cirionly", type="character", required=TRUE,
                    help="ciri only list of circRNAs")

parser$add_argument("-r", "--circExploreronly", type="character", required=TRUE,
                    help="circExplorer only list of circRNAs")

parser$add_argument("-c", "--common", type="character", required=TRUE,
                    help="common list of circRNAs")

args <- parser$parse_args()

# args=list()
# args$ciriout="Rep2_KO_72h.ciri.out"
# args$circExplorerout="Rep2_KO_72h.circularRNA_known.txt"

ciritable=read.csv(args$ciriout,header = TRUE,sep="\t")
circExplorertable=read.csv(args$circExplorerout,header=FALSE,"\t")
circExplorertable$circRNA_ID=paste0(circExplorertable$V1,":",circExplorertable$V2+1,"|",circExplorertable$V3)


xx=list(CIRI=ciritable$circRNA_ID,circExplorer=circExplorertable$circRNA_ID)

png(args$plot)
ggVennDiagram(xx)+ scale_fill_gradient(low="blue",high = "red")
dev.off()

regions=ggVennDiagram::get_region_items(xx)

write(regions$A,args$cirionly)
write(regions$B,args$circExploreronly)
write(regions$AB,args$common)
