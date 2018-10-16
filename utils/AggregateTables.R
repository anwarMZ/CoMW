#Author: Muhammad Zohaib Anwar
#License: GPL v3.0\n\n

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("File & Evalue must be provided", call.=FALSE)
} 

library(plyr)
library(reshape2)
setwd(args[1])
DBID=read.table(args[2],sep="\t",header=T)
DBID_agg <- ddply(DBID, "ContigID", numcolwise(sum))
write.table(DBID_agg,args[2],sep="\t",quote=F,row.names=F)

