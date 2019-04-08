"""

Author: Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This script will use SWORD to align the assembled contigs from previous step against
database of choice from following options
1. Md5nr https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-141 and
eggNOG annotation http://eggnogdb.embl.de/#/app/home
2. CAZy https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2686590/
3. NCyc https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty741/5085377
to provide alignment results in BM9 format using multiple threads. 

Dependencies:
1. Databases in $CoMW/databases

Example:
python align_contigs_to_database.py -f $Contigs.fasta -s 12 -n 6 -o $SWORD_result.tsv -t 12 -d 1 -r y
Given an input FASTA file $Contigs.fasta  is aligned against Md5nr using 12 threads and 6 
possible ORFs generated an alignment file $[*].tsv. The input file is splitted into 12 parts
after translation in order to save running memory 

python align_contigs_to_database.py -f $Contigs.fasta -s 12 -n 1 -o $SWORD_result.tsv -t 12 -d 2 -r y
Given an input FASTA file $Contigs.fasta  is aligned against CAZy using 12 threads and 1 
possible ORFs generated an alignment file $[*].tsv. The input file is splitted into 12 parts
after translation in order to save running memory 

python align_contigs_to_database.py -f $Contigs.fasta -s 12 -n 3 -o $SWORD_result.tsv -t 12 -d 3 -r y
Given an input FASTA file $Contigs.fasta  is aligned against NCyc using 12 threads and 3
possible ORFs generated an alignment file $[*].tsv. The input file is splitted into 12 
parts after translation in order to save running memory 


"""

import subprocess
import sys
import argparse
import os
from subprocess import Popen, PIPE
import shlex
from pyfasta import Fasta
import pyfasta
import os.path as path
import shutil


##Reading parapeters
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-s", "--splitsize", help= "Number of parts Fasta file to be splitted in", required=True, type=int, default = 1)
parser.add_argument("-n", "--ORFs", help= "Number of ORFs (1-6) to be calculated for alignment",required=True, type=int, default=1)
parser.add_argument("-t", "--threads", help= "Number of threads to be run",required=True, type=int, default = 1)
parser.add_argument("-d", "--database", help= "Alignment database of choice 1: Md5nr, 2: CAZy, 3: NCyc", required=True, type=int, default = 1 )
parser.add_argument("-r", "--remove", help= "Remove temporary files [y/n]", default = 'n' )
requiredName = parser.add_argument_group('required arguments')
requiredName.add_argument("-f", "--inputfastafile", help= "Fasta file of assembled contigs, output from Trinity", required=True )
requiredName.add_argument("-o", "--outputfile", help= "Output file .tsv format",required=True)

args = parser.parse_args()
if not args.inputfastafile: print("No input file provided")
if not args.outputfile: print("No output file provided")
if not (args.inputfastafile or args.outputfile): sys.exit(1)
if args.database not in [1,2,3]:sys.exit(1)
if args.remove not in ['y','n']:sys.exit(1)


def translation(filename, orfs):
	transeqcommand = ["transeq","-sformat","pearson", "-clean", "-frame", str(orfs), "-sequence", str(filename), "-outseq", str(outputdir+"/TempFiles/"+transout)]
	subprocess.call(transeqcommand)


def split(filename, size):
	splitcommand = ["pyfasta", "split", "-n", str(size), str(outputdir+"/TempFiles/"+transout)]
	subprocess.call(splitcommand)
	

def batchsword(fastalist, Subjectdatabase, threads):	
	for i in fastalist:
		command = ["sword", "-i", str(outputdir+"/TempFiles/"+i), "-t", str(threads), "-o", str(outputdir+"/TempFiles/"+i[:i.index(".fasta")]+".result.tsv"), "-f", "bm9",  "-j", str(Subjectdatabase), "-c", str(30000)]
		subprocess.call(command)

CoMWdir = os.path.realpath(__file__)
dbdir =  path.abspath(path.join(__file__ ,"../../databases"))
utildir = path.abspath(path.join(__file__ ,"../../utils"))
outputdir, outputfile = os.path.split(args.outputfile)
inputdir, inputfile = os.path.split(args.inputfastafile)


if __name__ == "__main__":
	
	if not os.path.exists(outputdir+"/TempFiles"):
		os.makedirs(outputdir+"/TempFiles")
	n = args.ORFs
	transout="Translated_"+inputfile
	print("\n\nTranslating the query sequence to proteins in " +str(args.ORFs) + " ORF(s)")
	translation(filename = inputdir+"/"+inputfile, orfs = n)
	print("\n\nTranslation done")

	n=args.splitsize
	print("\n\nSplitting translated fasta file into "+str(n)+ " to use less memory")
	split(filename = outputdir+"/"+transout, size=n)
	print("\n\nSplitting done")

	files=[]
	x=os.listdir(outputdir+"/TempFiles/")
	for i in x:
		if ".fasta" in i:
			files.append(i)
	files.remove(transout)
	files.remove(transout+".gdx")
	files.remove(transout+".flat")
	t = int(args.threads)
	db=int(args.database)
	if db is 1:
		database=dbdir+"/M5NR_protien.fasta"
		batchsword(fastalist = files, Subjectdatabase =  database, threads = t)	
	elif db is 2:
		database=dbdir+"/CAZyDB.07202017.fa"
		batchsword(fastalist = files, Subjectdatabase =  database, threads = t)
	elif db is 3:
		database=dbdir+"/NCyc_100.faa"
		batchsword(fastalist = files, Subjectdatabase =  database, threads = t)
	else:
		print("Wrong input for database: Only options 1,2 & 3, see help with -h")

	for f in files:
		os.system("cat "+outputdir+"/TempFiles/"+f.replace(".fasta",".result.tsv")+" >> " + outputdir+"/"+outputfile)
	agree= args.remove
	if agree is 'y':
		shutil.rmtree(outputdir+"/TempFiles/")
