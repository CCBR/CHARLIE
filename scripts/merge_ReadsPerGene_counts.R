rm(list=ls())

files=list.files(pattern="*ReadsPerGene*",recursive=TRUE)

read_counts<-function(fname,sname,id){
  x=read.csv(fname,header=FALSE,sep="\t")[,c(1,id)]
  colnames(x)=c("Gene",sname)
  return(x)
}

datasets_unstranded=list()
datasets_stranded=list()
datasets_revstranded=list()

for (i in 1:length(files)){
  sname=unlist(strsplit(basename(files[i]),"_p2"))[1]
  datasets_unstranded[[sname]]=read_counts(files[i],sname,2)
  datasets_stranded[[sname]]=read_counts(files[i],sname,3)
  datasets_revstranded[[sname]]=read_counts(files[i],sname,4)	
}


x=Reduce(function(d1, d2) merge(d1, d2, by = "Gene", all.x = TRUE, all.y = FALSE), 
       datasets_unstranded)
y=Reduce(function(d1, d2) merge(d1, d2, by = "Gene", all.x = TRUE, all.y = FALSE), 
         datasets_stranded)
z=Reduce(function(d1, d2) merge(d1, d2, by = "Gene", all.x = TRUE, all.y = FALSE), 
         datasets_revstranded)

write.table(x,file="unstranded_STAR_GeneCounts.tsv",quote = FALSE,row.names = FALSE,sep="\t")
write.table(y,file="stranded_STAR_GeneCounts.tsv",quote = FALSE,row.names = FALSE,sep="\t")
write.table(z,file="revstranded_STAR_GeneCounts.tsv",quote = FALSE,row.names = FALSE,sep="\t")
