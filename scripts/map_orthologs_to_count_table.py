"""

Author: Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This script will map the aligned genes to the count table using the map generated in parse_sword.py 

Dependencies:
1. $CoMW/utils/AggregateTables.R

Example:
python map_orthologs_to_count_table.py -i $Abundance_table.tsv -m $SWORD_result_eggNOG.map -o $ggNOG_Counttable.tsv 
Given an input abundance table $abundance_table.tsv this script maps the identified genes using the map generated in parse_sword.py

"""

import subprocess
import sys
import argparse
import os
import os.path as path

 
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
requiredName = parser.add_argument_group('required arguments')
requiredName.add_argument("-i", "--inputfile", help= "Table file from BWA mapper output")
requiredName.add_argument("-m", "--mapfile", help= "Map file from SWORD parsed output")
requiredName.add_argument("-o", "--outputfile", help= "Output file in tsv file")

args = parser.parse_args()
if not args.inputfile: print("No input file provided")
if not args.outputfile: print("No output file provided")
if not args.mapfile: print("No map file provided")
if not (args.inputfile or args.outputfile or args.mapfile): sys.exit(1)

CoMWdir = os.path.realpath(__file__)
utildir = path.abspath(path.join(__file__ ,"../../utils"))
outputdir, outputfile = os.path.split(args.outputfile)
inputdir, inputfile = os.path.split(args.inputfile)
mapdir, mapfile = os.path.split(args.mapfile)



def aggregate(filedir, filename):
	subprocess.call(['Rscript', utildir+'/AggregateTables.R',filedir, filename])

if __name__ == "__main__":
	
	table = open(inputdir+"/"+inputfile,'r')
	text = table.readlines()
	table.close()

	contigs = open(mapdir+"/"+mapfile, 'r')
	lines=contigs.readlines()
	contigs.close()

	counttable=open(outputdir+"/"+outputfile,'w')
	counttable.write(text[0])
	for i in range(0,len(lines)):
		contigID=lines[i].split("\t")[0].strip()
		contigIDx=contigID[:len(contigID)-2]
		db_id=lines[i].split("\t")[1].strip()
		for items in range(1,len(text)):
			if contigIDx in text[items]:
				count = text[items][text[items].index("\t"):].strip()
				counttable.write(db_id+"\t"+count+"\n")
	counttable.close()
	aggregate(filedir = outputdir, filename = outputfile)
