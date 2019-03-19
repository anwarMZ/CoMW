'''
Created on Jan 5, 2011
Updated 29 April 2016 to deal with infernal 1.1.1
Updated 10 October 2016 to fix bug failing to identify abundant ncRNAs (>1000)

@author: Anders Lanzen

Used to parse output from cmsearch to give output identifiers with each Model
If a read looks like more than one Model - the lowest bitscore wins it
'''

import sys
import os

class Model:
	def __init__(self,name):
		self.name = name
		self.reads = []
    
	def printReadNumber(self,ampnoise_us_annotation=False):
		print("%s\t%s" % (self.name,self.getTotalReads(ampnoise_us_annotation)))
    
	def getTotalReads(self,ampnoise_us_annotation=False):
		total=0
		for read in self.reads:
			if ampnoise_us_annotation and ("_" in read): 
				total+=int(read[read.find("_")+1:])
			else: total+=1
		return total
        
	def addRead(self, readname):
		self.reads.append(readname)

class CmsearchOut:
	def __init__(self, cmOutFile):
		self.cmOut = cmOutFile
		self.readModel = {}
		self.readEValue = {}
		self.models = {}
        
	def parseOld(self,evalue):
		tRead = open(self.cmOut,"r")
		model = None
		read = None
		for line in tRead:
			if line.startswith("CM:"):
				model = line[4:-1]
			elif line[0]=='>':
				read = line[1:-1]
			elif line.startswith(" Score = "):
				try: e=float(line[line.find("E = ")+4:line.find(", P = ")])
				except: e=1.0
				if e <= evalue:
					if read in self.readModel.keys():
						if self.readEValue[read] > e:
							self.readEValue[read] = e
							self.readModel[read] = model
                            #print "%s reassigned to %s" % (read,model)
					else:
						self.readModel[read] = model
						self.readEValue[read] = e
                        #print "%s assigned to %s" %(read, model)
		for read in self.readModel.keys():
			modelname = self.readModel[read]
			if modelname in self.models.keys():
				self.models[modelname].addRead(read)
			else:
				m = Model(modelname)
				self.models[modelname]=m
				m.addRead(read)
                
	def parse(self,evalue):
		"""Updated method dealing with 1.1+ output"""
        
		inScoreTable = False
		tRead = open(self.cmOut,"r")
		model = None
		read = None
		for line in tRead:
			items = line.split()
			if line.startswith("Query:"):
				model = items[1]
                #print("DEBUG: New model: %s" %model)
    
			elif len(items)>4 and items[0].startswith("-") and items[1].startswith("-"):
				inScoreTable=True
                #print("DEBUG: Parsing hits to %s" %model)
			elif  "inclusion threshold" in line or (inScoreTable and len(items)<=4):
				inScoreTable=False
                #print("DEBUG: Finished parsing hits to %s" %model)
			elif inScoreTable and len(items)>4:
				read= items[5]
				e = float(items[2])
                
				if e <= evalue:
					if read in self.readModel.keys():
						if self.readEValue[read] > e:
							self.readEValue[read] = e
							self.readModel[read] = model
                            #print "DEBUG: %s reassigned to %s" % (read,model)
					else:
						self.readModel[read] = model
						self.readEValue[read] = e
                        #print "DEBUG: %s assigned to %s" %(read, model)
        
		for read in self.readModel.keys():
			modelname = self.readModel[read]
			if modelname in self.models.keys():
				self.models[modelname].addRead(read)
			else:
				m = Model(modelname)
				self.models[modelname]=m
				m.addRead(read)
        
    
if __name__=="__main__":
	if len(sys.argv) < 3: 
		print("Use parsecm.py cmsearch-output e-value-cutoff [options] \n \
                Options:\n \ -underscore : Fasta file has underscore AmpliconNoise style annotation for abundance\n \
                -old : support for older infernal (before v1.1)\n")
	else:
		inputdir, inputfile = os.path.split(sys.argv[1])	
		mt = CmsearchOut(inputdir+"/"+inputfile)
	        
		if (len(sys.argv)>3 and sys.argv[3]=="-underscore") or (len(sys.argv)>4 and sys.argv[4]=="-underscore"):
			ampnoise_us_annotation = True
		else:
			ampnoise_us_annotation = False
		if (len(sys.argv)>3 and sys.argv[3]=="-old") or (len(sys.argv)>4 and sys.argv[4]=="-old"):
			mt.parseOld(float(sys.argv[2]))
		else:
			mt.parse(float(sys.argv[2]))
        #print("\n----------------------\n")
    
	ncRNA = open(inputdir+"/"+inputfile.replace(".out","ncRNA.txt"),"w")
	for Model in mt.models.values():
		Model.printReadNumber(ampnoise_us_annotation)
		fof = open(inputdir+"/TempFiles/"+inputfile.replace(".out","_")+Model.name+".fof","w")
		for read in Model.reads:
			fof.write(read+"\n")
		ncRNA.write(read + "\n")
		fof.close()
	ncRNA.close()
        #print("\n----------------------\n")
        #for rm in mt.readModel.keys():
        #    print("%s\t%s" % (rm, mt.readModel[rm]))