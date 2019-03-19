"""

Authors: Anders Lanzen & Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This is a script that uses Infernal (a secondary-structure-aware aligner). cmsearch module
of the infernal predicts the secondary structure of RNA sequences and similarities based 
on the consensus structure models of RFam. This script uses the $utils/parsecm.py then to
parse the oputput and filter the non-coding RNA contigs based on the confidence threshold
of alignment.

Dependencies:
1. $CoMW/utils/parsecm.py

Example:
python filter_ncRNA.py -f $Contigs.fasta -e 3 -t 16 -o $Contigs_ncrna_filtered.fasta -r n
Given an input fasta file $Contigs.fasta the script uses infernal to align the RNA contigs
using 16 threads in parallel against RFam and filter the non-coding RNAs with a confidence
threshold of 1E-1 and write to $Contigs_ncRNA_filtered.fasta

"""

import subprocess
import sys
import argparse
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os.path as path
import shutil
from Bio.SeqIO import FastaIO



parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-e", "--evalue", help= "E-value in integar", type=int, default=1)
parser.add_argument("-t", "--threads", help= "Number of threads", type=int, default=1)
parser.add_argument("-r", "--remove", help= "Delete temporary files created [y/n], default y", default = 'n')
requiredName = parser.add_argument_group('required arguments')
requiredName.add_argument("-f", "--fastafile", help= "Input fasta file")
requiredName.add_argument("-o", "--outputfile", help= "Output fasta file")



args = parser.parse_args()
if not args.fastafile: print("No fasta file provided")
if not args.outputfile: print("No output prefix provided")
if not (args.outputfile or args.fastafile): sys.exit(1)
if args.remove not in ['y','n']:sys.exit(1)


def multi2linefasta(fastafile):
	mfasta = fastafile.replace(tail,"_formatted"+tail)
	ifile = open(fastafile,'rU')
	with open(mfasta, 'w') as ofile:
		for record in SeqIO.parse(ifile, "fasta"):
			sequence = str(record.seq)
			ofile.write('>'+record.id+'\n'+sequence+'\n')


def filter_fasta(fastafile_formatted, idfile):
	f=open(idfile,'r')
	lines=f.readlines()
	f.close()
	ids=[]
	for i in range(1,len(lines)):
		lines[i]=lines[i].strip()
		ids.append(lines[i])
	f=open(outputdir+"/"+outputfile,'w')
	for record in SeqIO.parse(fastafile_formatted,'fasta'):
		if record.id not in ids:
			f.write(">"+str(record.id)+"\n")
			f.write(str(record.seq)+"\n")
	f.close()


CoMWdir = os.path.realpath(__file__)
dbdir =  path.abspath(path.join(__file__ ,"../../databases"))
utildir = path.abspath(path.join(__file__ ,"../../utils"))
fastadir, fastaf = os.path.split(args.fastafile)
outputdir, outputfile = os.path.split(args.outputfile)
cpus = args.threads
e = str(args.evalue)
head, tail = os.path.splitext(fastaf)


if __name__ == "__main__":

	if not os.path.exists(outputdir+"/TempFiles"):
		os.makedirs(outputdir+"/TempFiles")
	print("Running Infernal to detect ncRNAs\n")
	print(outputdir+"/"+fastaf.replace(".fasta","_cmsearch.out"))
	command=["cmsearch", "--cpu", str(cpus), "-o", outputdir+"/"+fastaf.replace(".fasta","_cmsearch.out"), dbdir+"/"+"Rfam.cm" , fastadir+"/"+fastaf]
	#subprocess.call(command)
	print(fastadir+"/"+fastaf)
	print("Parsing Output\n")
	command = ["python", utildir+"/parsecm.py", outputdir+"/"+fastaf.replace(".fasta","_cmsearch.out"), "1E-"+e]
	subprocess.call(command)
	print("ncRNAs predicted\n")
	print("Filtering FASTA file provided\n")
	print(fastadir+"/"+fastaf)
	multi2linefasta(fastafile = fastadir+"/"+fastaf)
	print(fastadir+"/"+fastaf.replace(tail,"_formatted"+tail))
	print( outputdir +"/"+fastaf.replace(".fasta","_cmsearchncRNA.txt"))
	filter_fasta(fastafile_formatted = fastadir+"/"+fastaf.replace(tail,"_formatted"+tail), idfile = outputdir +"/"+fastaf.replace(".fasta","_cmsearchncRNA.txt"))
	print("FASTA file filtered\n")
	agree= args.remove
	if agree is 'y':
		shutil.rmtree(fastadir+"/TempFiles/")