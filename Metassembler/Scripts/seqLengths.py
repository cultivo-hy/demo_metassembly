#!/usr/bin/env python2.7
import sys
import re
import os

USAGE="USAGE: "+os.path.basename(sys.argv[0])+" in_fasta"

if(len(sys.argv) < 2):
	sys.exit(USAGE)

fasta=sys.argv[1]

f=open(fasta,'r')

chrom=""
count=0

while 1:
	line=f.readline()
	if not line: break

	line=line.rstrip('\n')	

	header=re.search("^>(.*)",line)
	if header:
		if chrom:
			print "%s\t%lu" %(chrom,count)
		chrom=header.group(1)
		count=0
	else:
		count+=len(line)

if chrom:
	print "%s\t%lu" %(chrom,count)
else :
	print "Did not find any chromosome header"
