#!/usr/bin/env python2.7
import sys
import re
import os
import string

def PrintRow(outfile,data):
	outfile.write(data[0])
	for x in data[1:]:
		outfile.write("\t"+x)
	outfile.write("\n")

USAGE="USAGE: "+sys.argv[0]+" outFile inFile1 inFile2 ...";

if(len(sys.argv) < 2):
	sys.exit(USAGE)

outf=sys.argv[1]
if os.access(outf,os.F_OK):
	iout=open(outf,'r')
	line=iout.readline()
	rei=re.search("n=(\d+)\s\[(\d+),\s(\d+)\]\s([\d\.]+)\s\+/-\s([\d\.]+)\ssum=(\d+)\sn50=(\d+)\sn50cnt=(\d+)\s",line)
	if rei:
		sys.exit("ERROR:First argument is a stats file,\n you forgot to specify the out file")
	iout.close()

out=open(outf,'w')

names=["Genome/rank"]
ns=["N"]
sums=["Sum"]	
mins=["Min"]
maxs=["Max"]
means=["Mean"]
sds=["Sd"]
n50s=["N50"]
n50cnts=["N50_count"]

for file in sys.argv[2:]:
	names.append(os.path.basename(file))
	f=open(file,'r');
	line=f.readline();
	if not line:
		sys.exit("Couldn't open %s" %(file));
	print line
	reg=re.search("n=(\d+)\s\[(\d+),\s(\d+)\]\s([\d\.]+)\s\+/-\s([\d\.]+)\ssum=(\d+)\sn50=(\d+)\sn50cnt=(\d+)\s",line)
	n=reg.group(1)
	min=reg.group(2)
	max=reg.group(3)
	mean=reg.group(4)
	sd=reg.group(5)
	sum=reg.group(6)
	n50=reg.group(7)
	n50cnt=reg.group(8)

	ns.append(n)
	mins.append(min)
	maxs.append(max)
	means.append(mean)
	sds.append(sd)
	sums.append(sum)
	n50s.append(n50)
	n50cnts.append(n50cnt)
	f.close()

PrintRow(out,names)
PrintRow(out,ns)
PrintRow(out,mins)
PrintRow(out,maxs)
PrintRow(out,means)
PrintRow(out,sds)
PrintRow(out,sums)
PrintRow(out,n50s)
PrintRow(out,n50cnts)
out.close()

