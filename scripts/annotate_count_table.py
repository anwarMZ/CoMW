"""
Author: Muhammad Zohaib Anwar
License: GPL v3.0\n\n


Description:
This script will annotate a given countatble against the database of choice from the following
1. Md5nr https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-141 and eggNOG annotation http://eggnogdb.embl.de/#/app/home
2. CAZy https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2686590/
3. NCyc https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty741/5085377


Dependencies:
1. Databases and annotations in $CoMW/databases


Example:
python annotate_count_table.py -i counttable.tsv -o counttable_annotated.tsv -d 1 
Given an input count table counttable.tsv is annotated using eggNOG hierarchial annotation

python annotate_count_table.py -i counttable.tsv -o counttable_annotated.tsv -d 2 
Given an input count table counttable.tsv is annotated using CAZy hierarchial annotation

python annotate_count_table.py -i counttable.tsv -o counttable_annotated.tsv -d 3 
Given an input count table counttable.tsv is annotated using NCyc hierarchial annotation


"""

import subprocess
import sys
import argparse
import os
import os.path as path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i", "--inputfile", help= "Table file from mapping output", required=True )
parser.add_argument("-o", "--outputfile", help= "Output file .tsv format",required=True)
parser.add_argument("-d", "--database", help= "1: Md5nr, 2: CAZy, 3: NCyc",  type=int, default = 1 )
args = parser.parse_args()
if not args.inputfile: print "No input file provided"
if not args.outputfile: print "No output file provided"
if not (args.inputfile and args.outputfile): sys.exit(1)
db=int(args.database)
if args.database not in [1,2,3]:sys.exit(1)



def annotate_eggNOG(tabfile, outfile):

	f=open(tabfile,"r")
	infile=f.readlines()
	f.close()

	##Reading Databases
	f=open(dbdir+"/all.funcat.eggnogv3.txt",'r') #Reading Sword COG_cat dictionary
	lines=f.readlines()
	f.close()

	fun_cat={}
	for items in lines:
		items=items.strip()
		fun_cat[items.split("\t")[0]]=items.split("\t")[1]

	##1st Level of Function
	levelI={}
	levelI["J"]="INFORMATION STORAGE AND PROCESSING"
	levelI["A"]="INFORMATION STORAGE AND PROCESSING"
	levelI["K"]="INFORMATION STORAGE AND PROCESSING"
	levelI["L"]="INFORMATION STORAGE AND PROCESSING"
	levelI["B"]="INFORMATION STORAGE AND PROCESSING"

	levelI["D"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["Y"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["V"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["T"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["M"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["N"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["Z"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["W"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["U"]="CELLULAR PROCESSES AND SIGNALING"
	levelI["O"]="CELLULAR PROCESSES AND SIGNALING"

	levelI["C"]="METABOLISM"
	levelI["G"]="METABOLISM"
	levelI["E"]="METABOLISM"
	levelI["F"]="METABOLISM"
	levelI["H"]="METABOLISM"
	levelI["I"]="METABOLISM"
	levelI["P"]="METABOLISM"
	levelI["Q"]="METABOLISM"

	levelI["R"]="POORLY CHARACTERIZED"
	levelI["S"]="POORLY CHARACTERIZED"

	levelI["X"]="MOBILOME"
	
	##2nd Level of Function
	f=open(dbdir+"/fun2003-2014.tab",'r') #Reading Sword COG_cat dictionary
	lines=f.readlines()[1:]
	f.close()

	levelII={}
	for items in lines:
		items=items.strip()
		levelII[items.split("\t")[0]]=items.split("\t")[1]
	

	##3rd Level of Function
	f=open(dbdir+"/eggNOG.md52id2ont",'r') #Reading Sword COG_cat dictionary
	lines=f.readlines()
	f.close()
	
	levelIII={}
	for items in lines:
		items=items.strip()
		levelIII[items.split("\t")[1]]=items.split("\t")[2]


	infile[0]=infile[0].strip()
	infile[0]="#"+infile[0]+"\ttaxonomy"

	for i in range(1,len(infile)):
		infile[i]=infile[i].strip()
		if '"' in infile[i]:
			infile[i]=infile[i].replace('"','')
		First_level=""
		Second_level=""
		Third_level=""
		if infile[i].split("\t")[0] in fun_cat:
			cat=fun_cat[infile[i].split("\t")[0]]
		else:
			cat="S"
		for j in cat:
			First_level=First_level+str(levelI[j])+","
		First_level=First_level[:len(First_level)-1]
		first=First_level.split(",")
		First_level=""
		for k in set(first):
			First_level=First_level+k+".."
		First_level=First_level[:len(First_level)-2]
	
		for j in cat:
			Second_level=Second_level+str(levelII[j])+".."
		Second_level=Second_level[:len(Second_level)-2]
		Second=Second_level.split("..")
		Second_level=""
		for k in set(Second):
			Second_level=Second_level+k+".."
		Second_level=Second_level[:len(Second_level)-2]
	
		if infile[i].split("\t")[0] in levelIII:
			Third_level=levelIII[infile[i].split("\t")[0]]
		else:
			Third_level="Function unknown"
		infile[i]=infile[i]+ "\tk__COG; p__"+str(First_level)+"; c__"+str(Second_level)+"; o__"+str(Third_level)+"; f__; g__; s__"

	out=open(outfile,'w')
	out.write(infile[0]+"\n")
	for i in range(1,len(infile)):
		if ".." in infile[i]:
			if ".." in infile[i][infile[i].index("c__") : infile[i].index("; o__")] and ".." not in infile[i][infile[i].index("p__") : infile[i].index("; c__")]:
				taxonomy = infile[i][infile[i].index("c__")+3 : infile[i].index("; o__")]
				inheritence=taxonomy.split("..")
				for j in range(0,len(inheritence)):
					ID=infile[i][:infile[i].index("\t")]+"_"+str(j+1)
					x=infile[i][infile[i].index("\t"):infile[i].index("c__")]+"c__"+inheritence[j]+infile[i][infile[i].index("; o__"):]
					out.write(ID+x+"\n")
			if ".." in infile[i][infile[i].index("c__") : infile[i].index("; o__")] and ".." in infile[i][infile[i].index("p__") : infile[i].index("; c__")]:
				taxonomyII = infile[i][infile[i].index("c__") +3 : infile[i].index("; o__")]
				inheritenceII=taxonomyII.split("..")
				for l in range(0,len(inheritenceII)):
					list_values = [key for key,val in levelII.items() if val==inheritenceII[l]]
					ID=infile[i][:infile[i].index("\t")]+"_"+str(l+1)
					x=infile[i][infile[i].index("\t"):infile[i].index("p__")]+"p__"+levelI[list_values[0]]+infile[i][infile[i].index("; c__"):]
					x=x.replace(x[x.index("c__"):x.index("o__")],"c__"+inheritenceII[l]+"; ")
					out.write(ID+x+"\n")
		else:	
			out.write(infile[i]+"\n")
	out.close()


def annotate_CAZy(tabfile,outfile):
	
	f=open(dbdir+"/CAZY_hierarchy.txt",'r')
	lines1=f.readlines()
	f.close()

	values=['Archaea', 'Bacteria', 'Eukaryota', 'unclassified']

	CAZy_id=[]
	CAZy_org=[]
	key=""
	value=""
	for i in range(0,len(lines1)):
		lines1[i]=lines1[i].strip()
		tokens=lines1[i].split("\t")
		if tokens[0] in values:
			value=tokens[0]
		if tokens[0] not in values and "." in tokens[0]:
			key=tokens[0]
			CAZy_id.append(key)
			CAZy_org.append(value)
	
	f=open(tabfile,'r')
	lines=f.readlines()
	f.close()
	
	Level2=""
	Level3=""
	for i in range(1,len(lines)):
		lines[i]=lines[i].strip()
		tokens=lines[i].split("\t")
		if "GH" in tokens[0].split("|")[1]:
			EnzymeClass="Glycoside Hydrolases"
		if "GT" in tokens[0].split("|")[1]:
			EnzymeClass="GlycosylTransferases"
		if "PL" in tokens[0].split("|")[1]:
			EnzymeClass="Polysaccharide Lyases"
		if "CE" in tokens[0].split("|")[1]:
			EnzymeClass="Carbohydrate Esterases"
		if "AA" in tokens[0].split("|")[1]:
			EnzymeClass="Auxiliary Activities"
		if "CBM" in tokens[0].split("|")[1]:
			EnzymeClass="Carbohydrate-Binding Modules"
		#print	tokens[0].split("|")[1]
		Level2=tokens[0].split("|")[1]
		lines[i]=lines[i]+"\tk__CAZy; p__"+EnzymeClass+"; c__"+Level2+"; o__; f__; g__; s__\n"
		if len(tokens[0].split("|")) > 2:
			Level3=tokens[0].split("|")[2]
			lines[i]=lines[i].replace("o__;","o__"+Level3+";")

	f=open(outfile,'w')
	f.write(lines[0].strip()+"\ttaxonomy\n")
	for i in range(1,len(lines)):
		f.write(lines[i])
	f.close()
	

def annotate_NCyc(tabfile, outfile):

	f=open(tabfile,"r")
	infile=f.readlines()
	f.close()
	

	f=open(dbdir+"/NCyc_Cat.txt",'r') 
	lines=f.readlines()
	f.close()

	
	f=open(dbdir+"/id2gene.map.txt",'r')
	lines1=f.readlines()
	f.close()


	fun_cat={}
	for items in lines:
		items=items.strip()
		fun_cat[items.split("\t")[0]]=items.split("\t")[1]
	
	fun_desc={}
	for items in lines:
		items=items.strip()
		fun_desc[items.split("\t")[0]]=items.split("\t")[2]
	
	fun_ID={}
	for items in lines1:
		items=items.strip()
		fun_ID[items.split("\t")[0]]=items.split("\t")[1]

	f=open(dbdir+"/NCyc_100.faa",'r')
	dbfile=f.readlines()
	f.close()	
	
	source = {}
	for dbheader in dbfile:
		dbheader=dbheader.strip()
		if ">" in dbheader:
			if "[" in dbheader:
				key = dbheader[1:dbheader.index("description")-2]
				desc = dbheader[dbheader.index("description")+12:]
				item=desc.split(" ")[0]
				source[key]=item
	
	infile[0]=infile[0].strip()
	infile[0]="#"+infile[0]+"\ttaxonomy"
	out=open(outfile,'w')
	out.write(infile[0].strip()+"\n")

	for i in range(1,len(infile)):
		infile[i]=infile[i].strip()
		if '"' in infile[i]:
			infile[i]=infile[i].replace('"','')
		First_level_1=""
		First_level_2="NULL"
		Second_level=""
		if infile[i].split("\t")[0] in fun_ID:
			Second_level = fun_ID[infile[i].split("\t")[0]]
		else:		
			Second_level = source[infile[i].split("\t")[0]]


		
		if Second_level in fun_cat:
			#print Second_level
			First_level_1=fun_cat[Second_level]
			#print First_level_1
			if "," in First_level_1:
				First_level_2=First_level_1.split(",")[1]
				First_level_1=First_level_1.split(",")[0]

		else:
			First_level_1="Unknown"
		
		
		if not "NULL" in First_level_2:
			infile[i]=infile[i]+ "\tk__NCyc; p__"+str(First_level_1)+"; c__"+str(Second_level)+"; o__; f__; g__; s__"
			out.write(infile[i].strip()+"\n")
			infile[i]="1_"+infile[i].split("k__")[0]+ "k__NCyc; p__"+str(First_level_2)+"; c__"+str(Second_level)+"; o__; f__; g__; s__"
			out.write(infile[i].strip()+"\n")
		else:
			infile[i]=infile[i]+ "\tk__NCyc; p__"+str(First_level_1)+"; c__"+str(Second_level)+"; o__; f__; g__; s__"
			out.write(infile[i].strip()+"\n")
	out.close()



CoMWdir = os.path.realpath(__file__)
dbdir =  path.abspath(path.join(__file__ ,"../../databases"))
utildir = path.abspath(path.join(__file__ ,"../../utils"))
outputdir, outputfile = os.path.split(args.outputfile)
inputdir, inputfile = os.path.split(args.inputfile)

if __name__ == "__main__":

	if db is 1:
		annotate_eggNOG(tabfile = inputdir+"/"+inputfile, outfile = outputdir+"/"+outputfile)
	elif db is 2:
		annotate_CAZy(tabfile = inputdir+"/"+inputfile, outfile = outputdir+"/"+outputfile)
	elif db is 3:
		annotate_NCyc(tabfile = inputdir+"/"+inputfile, outfile = outputdir+"/"+outputfile)
	else:
		print "Wrong input for database: Only options 1,2 & 3, see help with -h"
