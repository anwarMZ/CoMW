"""

Authors: Anders Lanzen & Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This is an optional script that filters the non-coding RNAs from the assembled contigs based on the alignment against RFam database.


Example
Given an input FATSA file generates a new FASTA file that includes only contigs that have filtered non-coding RNAs from the contigs based on the confidence threshold provided by the user. 


Dependencies:
1. $CoMW/utils/parsecm.py
2. Bio.Seq http://biopython.org/DIST/docs/api/Bio.Seq-module.html from biopython http://biopython.org  
3. Infernal in path http://eddylab.org/infernal/


Example:
python filter_ncRNA.py -f contigs.fasta -e 3 -t 16 -o contigs_ncrna_filtered.fasta -r y
 -e 1 -o out_prefix -r y 
filters contigs.fasta using ncRNA database and produces filtered contigs file

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
e=str(args.evalue)
head, tail = os.path.splitext(fastaf)


if __name__ == "__main__":

	if not os.path.exists(outputdir+"/TempFiles"):
		os.makedirs(outputdir+"/TempFiles")
	print "Running Infernal to detect ncRNAs\n"
	command=["cmsearch", "--cpu "+str(cpus), "-o "+outputdir+"/"+fastaf.replace(".fasta","_cmsearch.out"), dbdir+"/"+"Rfam.cm" , fastadir+"/"+fastaf]
	subprocess.check_output(command)
	print "Parsing Output\n"
	command = ["python", utildir+"/parsecm.py", outputdir+"/"+fastaf.replace(".fasta","_cmsearch.out"), "1E-"+e]
	subprocess.check_output(command)
	print "ncRNAs predicted\n"
	print "Filtering FASTA file provided\n"
	multi2linefasta(fastafile = fastadir+"/"+fastaf)
	filter_fasta(fastafile_formatted = fastadir+"/"+fastaf.replace(tail,"_formatted"+tail), idfile = outputdir +"/"+fastaf.replace(".fasta","_cmsearchncRNA.txt"))
	print "FASTA file filtered\n"
	agree= args.remove
	if agree is 'y':
		shutil.rmtree(fastadir+"/TempFiles/")
	
