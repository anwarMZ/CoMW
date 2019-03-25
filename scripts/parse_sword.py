"""

Author: Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This script is used for parsing BM9 output file from SWORD alignement to using a specific threshold e.g. 1E-5 against a database of choice from following
1. Md5nr https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-141 and eggNOG annotation http://eggnogdb.embl.de/#/app/home
2. CAZy https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2686590/
3. NCyc https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty741/5085377
to produce a map of identified genes or orthologs against each contig.

Dependencies:
1. Databases and annotations in $CoMW/databases

Example:
python parse_sword.py -i Sword_result.BM9 -e 3 -o parsed_SWORD_result.tsv -d 2 
Given an input SWord_output in BM9 this script parse BM9 file to produce a readable format parsed_SWord_result.tsv and a map file against the CAZy database

"""

import subprocess
import sys
import argparse
import os
import shlex
import os.path as path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-e", "--Evalue", help= "Evalue for threshold eg: 5,6", type=int, default=5)
parser.add_argument("-d", "--database", help= "1: Md5nr, 2: CAZy, 3: NCyc",  type=int, default = 1)
requiredName = parser.add_argument_group('required arguments')
requiredName.add_argument("-i", "--inputfile", help= "SWORD output in bm9 format", required=True,)
requiredName.add_argument("-o", "--outputfile", help= "Parsed Result file in .tsv format", required=True,)

args = parser.parse_args()
if not args.inputfile: print("No input file provided")
if not args.outputfile: print("No output file provided")
if not args.Evalue: print("No evalue file provided")
if not (args.inputfile or args.outputfile): sys.exit(1)
db=int(args.database)
if args.database not in [1,2,3]:sys.exit(1)




def Md5nr_map(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()

	md5sumDatabase=dbdir+"/eggNOG.md52id2ont" 
	f=open(md5sumDatabase,"r")
	db=f.readlines()
	f.close()

	dbdic={}
	md5dic={}


	for i in db:
		i=i.strip()
		key=i.split("\t")[0]
		s=i.split("\t")[1]
		element=s
		dbdic[key] = s

	for line in lines:
		line=line.strip()
		key1=line.split("\t")[0]
		element1=line.split("\t")[1]
		md5dic[key1]=element1

	f=open(filename.replace(".tsv","eggNOG.map"),"w")
	for key in md5dic:
		if md5dic[key] in dbdic:
			id_desc = dbdic[md5dic[key]]+"\t"
			md5 = md5dic[key].strip()
			s = key.strip() + "\t"  + id_desc.split("\t")[0].strip()
			f.write(s+"\n")
	f.close()	


def CAZy_map(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()

	f=open(filename.replace(".tsv","CAZy.map"),"w")
	for i in range(1,len(lines)):
		lines[i]=lines[i].strip()
		s = lines[i].split("\t")[0].strip()+"\t"+lines[i].split("\t")[1].strip()
		f.write(s+"\n")
	f.close()

		

def NCyc_map(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()

	md5sumDatabase=dbdir+"/id2gene.map.txt" 
	f=open(md5sumDatabase,"r")
	db=f.readlines()
	f.close()

	dbdic={}
	md5dic={}


	for i in db:
		i=i.strip()
		key=i.split("\t")[0]
		s=i.split("\t")[1]
		element=s
		dbdic[key] = s
	

	for line in lines:
		line=line.strip()
		key1=line.split("\t")[0]
		element1=line.split("\t")[1]
		md5dic[key1]=element1
	
	f=open(filename.replace(".tsv","NCyc.map"),"w")
	for key in md5dic:
		if md5dic[key] in dbdic:
			md5 = md5dic[key].strip()
			s = key.strip() + "\t" +md5
			f.write(s+"\n")

	f.close()	



CoMWdir = os.path.realpath(__file__)
dbdir =  path.abspath(path.join(__file__ ,"../../databases"))
utildir = path.abspath(path.join(__file__ ,"../../utils"))
outputdir, outputfile = os.path.split(args.outputfile)
inputdir, inputfile = os.path.split(args.inputfile)
e="1E-"+str(args.Evalue)
e=format(float(e),"."+str(args.Evalue)+"f")
e=str(e)


if __name__ == "__main__":
	
	SWORD= open(inputdir+"/"+inputfile,'r')
	parser=SWORD.readlines()
	SWORD.close()
	header = "Query id\tSubject id\t% identity\talignment length\tmismatches\tgap openings\tq. start\tq. end\ts. start\ts. end\te-value\tscore\n"
	SWORDparsed=open(outputdir+"/Tempfile.tsv",'w')
	SWORDparsed.write(header)
	for line in parser:
		if "Fields" not in line and "Query id" not in line:
			SWORDparsed.write(line)
	SWORDparsed.close()
	subprocess.call(['Rscript', utildir+'/ParsingSword.R',outputdir,"Tempfile.tsv",outputfile,e])
	
	if db is 1:
		Md5nr_map(filename = outputdir+"/"+outputfile)
	elif db is 2:
		CAZy_map(filename = outputdir+"/"+outputfile)
	elif db is 3:
		NCyc_map(filename = outputdir+"/"+outputfile)
	else:
		print("Wrong input for database: Only options 1,2 & 3, see help with -h")
