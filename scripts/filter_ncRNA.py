"""

Authors: Anders Lanzen & Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This is an optional script filters the contigs less than a given threshold of relative expression. eg if e=1 \n Only Contigs with sum > 1/sum(Minimum Reads) are selected.
Filters out contigs from both count table [output from map_reads_to_contigs.py] and fasta file of contigs assembled.


Example
Given an input count table and FATSA file generates a new count table and FASTA file that includes only contigs that have a relative expression of higher than the threshold specified by the user. 


Dependencies:
1. $CoMW/utils/Filteration.R
2. Bio.Seq http://biopython.org/DIST/docs/api/Bio.Seq-module.html from biopython http://biopython.org  


Example:
python filter_table_by_abundance.py -i abundance_table.tsv -f contigs.fasta -e 1 -o out_prefix -r y 
filters abundance_table.tsv and contigs.fasta using expression 1% and producing the new abundance table and contigs file with output prefix in same directory

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
parser.add_argument("-f", "--fastafile", help= "Fasta file")
parser.add_argument("-e", "--evalue", help= "Evalue in integar", type=int, default=1)
parser.add_argument("-t", "--threads", help= "Number of Threads", type=int, default=1)
parser.add_argument("-o", "--outputfile", help= "Output file for fasta file")
parser.add_argument("-r", "--remove", help= "Delete temporary files created [y/n], default y", default = 'y' )

args = parser.parse_args()
if not args.fastafile: print "No fasta file provided"
if not args.outputfile: print "No output prefix provided"
if not (args.outputfile and args.fastafile): sys.exit(1)
if args.remove not in ['y','n']:sys.exit(1)


def filter_fasta(fastafile, prefix):
	f=open(fastadir+"/TempFiles/"+prefix+"_IncludedContigs.txt",'r')
	lines=f.readlines()
	f.close()

	ids=[]
	for i in range(1,len(lines)):
		lines[i]=lines[i].strip()
		ids.append(lines[i].split("\t")[1])

	f=open(fastadir+"/"+prefix+"_ncRNA_Removed_"+fastaf,'w')
	for record in SeqIO.parse(fastafile,'fasta'):
		if record.id in ids:
			f.write(">"+str(record.id)+"\n")
			f.write(str(record.seq)+"\n")

	f.close()


CoMWdir = os.path.realpath(__file__)
dbdir =  path.abspath(path.join(__file__ ,"../../databases"))
utildir = path.abspath(path.join(__file__ ,"../../utils"))
fastadir, fastaf = os.path.split(args.fastafile)
outputdir, outputfile = os.path.split(args.outputfile)
cpus = args.threads
e=str(args.evalue)


if __name__ == "__main__":

	if not os.path.exists(fastadir+"/TempFiles"):
		os.makedirs(fastadir+"/TempFiles")
	print "Running Infernal to detect ncRNAs\n"
	command=["cmsearch", "--cpu "+str(cpus), "-o "+fastadir+"/cmsearch_"+fastaf.replace(".fasta",".out"), dbdir+"/"+"Rfam.cm" , fastadir+"/"+fastaf]
	subprocess.check_output(command)
	print "Parsing Output\n"
	command =["python",utildir+"/cmsearch/parsecm.py", fastadir+"/cmsearch_"+fastaf.replace(".fasta",".out"), " 1E-" +e]
	print "ncRNAs predicted\n"
	print "Filtering FASTA file provided\n"
	#filter_fasta(fastadir+"/"+fastaf, prefix = outprefix)
	print "FASTA file filtered\n"
	agree= args.remove
	if agree is 'y':
		shutil.rmtree(fastadir+"/TempFiles/")
	


