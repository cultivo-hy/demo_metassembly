#!/usr/bin/env python

import sys
import re
import os
import string

def PrintTable(outfile,data):
	#data=dictionary
	outfile.write("File")
	for value in sys.argv[2:]:
		outfile.write("\t"+value)
	outfile.write("\n")
	for key in data.keys():
		outfile.write(key)
		for value in data[key]:
			outfile.write("\t"+value)
		outfile.write("\n")
	outfile.write("\n")


USAGE="USAGE: "+sys.argv[0]+" outfile inFile1 inFile2 ..."

if(len(sys.argv) < 2):
	sys.exit(USAGE)

outf=sys.argv[1]
if os.access(outf,os.F_OK):
	iout=open(outf,'r')
	line=iout.readline()
	while line == '\n':
		line=iout.readline()
	rei=re.search("WARNING executed HMMER version has not been tested !!!",line)
	if rei:
		sys.exit("ERROR: Firs argument is a CEGMA completeness report file\nyou forgot to specify the out file")
	rei=re.search("Statistics of the completeness of the genome based on 248 CEGs",line)
	if rei:
		sys.exit("ERROR: Firs argument is a CEGMA completeness report file\nyou forgot to specify the out file")
	iout.close()


metrics={'Complete_Prots':[],
	 'Complete_Percent_completeness':[],
	 'Complete_Total_CEGS':[],
	 'Partial_Prots':[],
	 'Partial_Percent_completeness':[],
	 'Partial_Total_CEGS':[]}

for file in sys.argv[2:]:
	f=open(file,'r')
	while 1:		
		line=f.readline()
		if not line : break
		                                 #Prots  %Complete       Total    Average     Ortho
		recomp=re.search("^\s+Complete\s+(\d+)\s+([\d\.]+)\s+-\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)\s*",line)
		if recomp:
			metrics['Complete_Prots'].append(recomp.group(1))
			metrics['Complete_Percent_completeness'].append(recomp.group(2))
			metrics['Complete_Total_CEGS'].append(recomp.group(3))

		repart=re.search("^\s+Partial\s+(\d+)\s+([\d\.]+)\s+-\s+(\d+)\s+([\d\.]+)\s+([\d\.]+)\s*",line)
		if repart:
			metrics['Partial_Prots'].append(repart.group(1))
                        metrics['Partial_Percent_completeness'].append(repart.group(2))
                        metrics['Partial_Total_CEGS'].append(repart.group(3))

out=open(outf,'w')
PrintTable(out,metrics)


