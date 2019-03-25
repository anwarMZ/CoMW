#Author: Anders Lanzen
#License: GPL v3.0\n\n

args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Not enough Parameters", call.=FALSE)
} 

setwd(args[1])
cag = read.table(args[2],sep="\t",header=T,row.names=1)

svg('./TempFiles/MappedReadsperdataset.svg')
hist(colSums(cag),col="grey",breaks=10,main="mapped reads per dataset",xlab="reads")
dev.off()
svg('./TempFiles/MappedReadsperContig.svg')
hist(log10(rowSums(cag)),col="grey",main="mapped reads per contig",xlab="log10(reads)")
dev.off()

minReads=min(colSums(cag))


require(vegan)
cagt = t(cag)
cagr = t(decostand(cagt,method="total"))

svg('./TempFiles/MeanRealtiveContigExpression.svg')
hist(log10(rowMeans(cagr)),col="grey",main="Mean relative contig expression")
dev.off()

exp=as.double(args[3])
included_contigs = row.names(cagr)[(rowMeans(cagr)>=exp/minReads)]
write.table(included_contigs,sep="\t",file=paste("TempFiles/",args[4],"_IncludedContigs.txt",sep=""),quote = FALSE)
cagKeep = cag[(rowMeans(cagr)>=exp/minReads),]
cagKeep = cbind(ContigID=rownames(cagKeep), cagKeep)
write.table(cagKeep,sep="\t",file=paste(args[4],"_AbundanceFiltered.tsv",sep=""),quote = FALSE, row.names=FALSE)
cagrKeep = cagr[(rowMeans(cagr)>=exp/minReads),]
write.table(cagrKeep,sep="\t",file=paste("TempFiles/",args[4],"_RelativeAbundances.tsv",sep=""),quote = FALSE)
