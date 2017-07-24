#!/usr/bin/env python

import sys
import re
import os
import string

def PrintTable(outfile,data):
	int_metrics=['gaps','gaps_total_length','gaps_gapBrk','percentage_gaps_broken','errors_within_contigs',
			'error_free_bases','percent_error_free','longest_scf','longest_scf_gapBrk',
			'N50','N50_gapBrk','N50_gapBrk/N50','mean_scf','mean_scf_gapBrk']
	#data=dictionary:
	outfile.write('filename')
	for value in data['filename']:
		value.rstrip()
		outfile.write("\t"+value)
	outfile.write("\n")
	for key in int_metrics:
		if not data.has_key(key):
			sys.exit("ERROR: metric %s not found" %(key))
		if key == 'filename':
			continue
		outfile.write(key)
		for value in data[key]:
			outfile.write("\t"+value)
		outfile.write("\n")
	outfile.write("\n")


def ParseColNames(instr):
	instr=re.sub("#filename","filename",instr)
	instr=re.sub("_br","_gapBrk",instr)
	instr=re.sub("longest","longest_scf",instr)
	instr=re.sub("mean_length","mean_scf",instr)
	instr=re.sub("bases","total_length",instr)
	instr=re.sub("error_free","error_free_bases",instr)
	instr=re.sub("FCD","errors_within_contigs",instr)
	return(instr)
	

USAGE="USAGE: "+sys.argv[0]+" outfile inFile1 inFile2 ..."

if(len(sys.argv) < 2):
	sys.exit(USAGE)

outf=sys.argv[1]
if os.access(outf,os.F_OK):
	iout=open(outf,'r')
	line=iout.readline()
	while line == '\n':
		line=iout.readline()
	rei=re.search("^#filename\tbases\tsequences\tmean_length\tlongest\tN50\tN50_n",line)
	if rei:
		sys.exit("ERROR: Firs argument is a REAPR 05.summary.report.tsv file\nyou forgot to specify the out file")
	rei=re.search("^Stats for original assembly \'00\.assembly\.fa\':",line)
	if rei:
		sys.exit("ERROR: Firs argument is a REAPR 05.summary.report.txt file\nyou forgot to specify the out file")
	iout.close()


metrics={}

for file in sys.argv[2:]:
	f=open(file,'r')
	lstr=f.readline()
	lstr.rstrip()
	lstr=ParseColNames(lstr)
	lstr=re.split("\s+",lstr)
	for metric in lstr:
		if not metrics.has_key(metric): 
			metrics[metric]=[]
	
	valstr=f.readline()
	valstr.rstrip()
	valstr=re.split("\s+",valstr)
	for idx in range(0,len(valstr)):
		metrics[lstr[idx]].append(valstr[idx])

	if not metrics.has_key('percentage_gaps_broken'):
		metrics['percentage_gaps_broken']=[]
	gapsBrk=metrics['gaps_gapBrk'][len(metrics['gaps_gapBrk'])-1]
	Igaps=metrics['gaps'][len(metrics['gaps'])-1]
	brokenGaps=float(Igaps)-float(gapsBrk)
	metrics['percentage_gaps_broken'].append(str(float(brokenGaps)/float(Igaps)))

	metric='percent_error_free'
	if not metrics.has_key(metric):
		metrics[metric]=[]
	lastv=len(metrics['error_free_bases'])-1
	metrics[metric].append(str(float(metrics['error_free_bases'][lastv])/float(metrics['total_length'][lastv])))
	
	metric='N50_gapBrk/N50'
	if not metrics.has_key(metric):
		metrics[metric]=[]
	lastv=len(metrics['N50_gapBrk'])-1
	metrics[metric].append(str(float(metrics['N50_gapBrk'][lastv])/float(metrics['N50'][lastv])))

	

out=open(outf,'w')
PrintTable(out,metrics)


