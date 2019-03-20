"""

Author: Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This script aligns quality filtered mRNA (merged or paired-end reads) against the assembled
contigs from RNA-Seq de novo transcriptome assemblers (e.g. Trinity). Given a directory 
with FASTQ files merged or paired-end and a FASTA file consisting the assembled contigs,
the script aligns using BWA mapper and produces an abundance table. This script can be 
parallelized using the threads option -t. 

Dependencies:
1. $CoMW/utils/MapReads_to_contigs.sh

Example:
python map_reads_to_contigs.py -f $Contigs.fasta -i $Fastq_dir -o $Output_dir -t 12 -m n 
aligns paired-end fastq reads present in $Fastq_dir against contigs.fasta using 12 threads
and producing the abundance table in $Output_dir

python map_reads_to_contigs.py -f $Contigs.fasta -i $Fastq_dir -o $Output_dir -t 16 -m y 
aligns merged fastq reads present in $Fastq_dir against contigs.fasta using 16 threads and
producing the abundance table in $Output_dir

"""

import subprocess
import sys
import argparse
import os
import csv
import pandas
import os.path as path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-t", "--threads", help= "Number of Threads", type=int, default=1)
parser.add_argument("-m", "--merged", help='Single or Paired-end, default = paired]', default = 'paired')
requiredName = parser.add_argument_group('required arguments')
requiredName.add_argument("-f", "--fastafile", help= "Fasta file of contigs")
requiredName.add_argument("-i", "--readsdir", help= "Fastq file directory")
requiredName.add_argument("-o", "--outputfile", help= "Output file")


args = parser.parse_args()
if not args.fastafile: print("No contigs file provided")
if not args.readsdir: print("No reads directory provided")
if not args.outputfile: print("No output file provided")
if not (args.fastafile or args.outputfile or args.readsdir): sys.exit(1)



def mapping(filename, directory, cpus):
	command = str(utildir+"/MapReads_to_contigs.sh "+str(filename)+" " +str(directory) +" "+ str(cpus))
	os.system(command)


def merge(filename, fileout):
	for i in range(0,len(names)):
		names[i]=names[i][:names[i][:names[i].index("contig_abundances")].rfind("_")-3]

	snames=set(names)
	csv=pandas.read_csv(filename,sep="\t")
	for n in snames:
		csv[n]=csv[n+"_R1"]+csv[n+"_R2"]
	for n in snames:
		csv=csv.drop([n+"_R1"], axis=1)
		csv=csv.drop([n+"_R2"], axis=1)
	csv.to_csv(fileout, sep='\t',index_label=False, index=False)

CoMWdir = os.path.realpath(__file__)
dbdir =  path.abspath(path.join(__file__ ,"../../databases"))
utildir = path.abspath(path.join(__file__ ,"../../utils"))
contigsdir, contigsfile = os.path.split(args.fastafile)
outputdir, outputfile = os.path.split(args.outputfile)
readsdir = str(args.readsdir)+"/"
seq_type = str(args.merged)
threads = str(args.threads)

if __name__ == "__main__":
	
	mapping(filename = contigsdir+"/"+contigsfile , directory = readsdir, cpus = threads)
	dicts=[]
	names=[]
	for i in os.listdir(readsdir):
		if "contig_abundances.txt" in i:
			f=open(readsdir+i,'r')
			names.append(str(i))
			r = csv.reader(f, delimiter='\t')
			dict1 = {row[0]: row[2] for row in r if "*" not in row}
			dicts.append(dict1)

	keys = set().union(*(d.keys() for d in dicts))
	with open("temp_combined.tsv", 'w') as f:
		w = csv.writer(f, delimiter='\t')
		w.writerow(["ContigID"]+[n[:n.index("contig_abundances")-1]for n in names])
		for key in keys:
			w.writerow([key] + [d.get(key, "0") for d in dicts])

	if seq_type == 'paired':
		merge(filename = "temp_combined.tsv", fileout = outputdir+"/"+outputfile)	
		os.remove("temp_combined.tsv")
	else:
		csv=pandas.read_csv("temp_combined.tsv",sep="\t")
		csv.to_csv(outputdir+"/"+outputfile, sep='\t',index_label=False, index=False)
