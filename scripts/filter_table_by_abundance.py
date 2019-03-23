"""

Authors: Anders Lanzen & Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This is an optional script filters the contigs less than a given threshold of relative 
expression. eg if e=1 \n Only Contigs with sum > 1/sum(Minimum Reads) are selected. Filters
out contigs from both count table [output from map_reads_to_contigs.py] and fasta file of 
contigs assembled.

Dependencies:
1. $CoMW/utils/Filteration.R

Example:
python filter_table_by_abundance.py -i $Abundance_table.tsv -f $Contigs.fasta -e 1 -o out_prefix -r y 
Given an abundance table the script filters $Abundance_table.tsv and $Contigs.fasta using
expression 1% and producing the new abundance table and contigs file with output prefix in
same directory

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


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-e", "--expression", help= "Relative expression in integars", default=1)
parser.add_argument("-o", "--outputprefix", help= "Output prefix for filtered table and fasta file")
parser.add_argument("-r", "--remove", help= "Delete temporary files created [y/n], default y", default = 'n' )
requiredName = parser.add_argument_group('required arguments')
requiredName.add_argument("-i", "--inputfile", help= "Table file from BWA mapper output")
requiredName.add_argument("-f", "--fastafile", help= "Fasta file")


args = parser.parse_args()
if not args.inputfile: print("No input file provided")
if not args.fastafile: print("No fasta file provided")
if not args.outputprefix: print("No output prefix provided")
if not (args.inputfile or args.outputprefix or args.fastafile): sys.exit(1)
if args.remove not in ['y','n']:sys.exit(1)



def normalize_tab(path,tabfile,exp,prefix):
	subprocess.call(['Rscript', utildir+'/Filteration.R',path,tabfile,exp,prefix])	

def filter_fasta(fastafile, prefix):
	f=open(inputdir+"/TempFiles/"+prefix+"_IncludedContigs.txt",'r')
	lines=f.readlines()
	f.close()

	ids=[]
	for i in range(1,len(lines)):
		lines[i]=lines[i].strip()
		ids.append(lines[i].split("\t")[1])

	f=open(fastadir+"/"+prefix+"_AbundanceFiltered_"+fastaf,'w')
	for record in SeqIO.parse(fastafile,'fasta'):
		if record.id in ids:
			f.write(">"+str(record.id)+"\n")
			f.write(str(record.seq)+"\n")

	f.close()


CoMWdir = os.path.realpath(__file__)
dbdir =  path.abspath(path.join(__file__ ,"../../databases"))
utildir = path.abspath(path.join(__file__ ,"../../utils"))
fastadir, fastaf = os.path.split(args.fastafile)
outprefix= str(args.outputprefix)
inputdir, inputfile = os.path.split(args.inputfile)
e=str(args.expression)


if __name__ == "__main__":

	if not os.path.exists(inputdir+"/TempFiles"):
		os.makedirs(inputdir+"/TempFiles")
	print("Normalaizing and filtering abudnance table based on the relative expression\n")
	normalize_tab(path = inputdir, tabfile = inputfile, exp = e, prefix = outprefix)
	print("Abundance table normalized\n")
	print("Filtering FASTA file provided\n")
	filter_fasta(fastadir+"/"+fastaf, prefix = outprefix)
	print("FASTA file filtered")
	agree= args.remove
	if agree is 'y':
		shutil.rmtree(inputdir+"/TempFiles/")