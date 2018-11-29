#Author: Muhammad Zohaib Anwar
#License: GPL v3.0\n\n

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("File & Evalue must be provided", call.=FALSE)
} 

e=as.numeric(args[4])
setwd(args[1])
Sword=read.table(args[2],sep="\t",header=T)
Sword5 = Sword[Sword$e.value<=e,]
Swordsorted = Sword5[order(Sword5[,'e.value'],-Sword5[,'score']),]
Swordbest = Swordsorted[!duplicated(substr(Swordsorted[,1],1,nchar(as.character(Swordsorted[,1]))-2)),]
write.table(Swordbest,args[3],sep="\t",quote=F,row.names=F)
