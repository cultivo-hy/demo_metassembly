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

USAGE="USAGE: "+sys.argv[0]+" outfile inFile1 inFile2 ..."

if(len(sys.argv) < 2):
	sys.exit(USAGE)

outf=sys.argv[1]
if os.access(outf,os.F_OK):
	iout=open(outf,'r')
	line=iout.readline()
	while line == '\n':
		line=iout.readline()
	rei=re.search("-----Metassembly .+ statistics-----",line)
	if rei:
		sys.exit("ERROR: Firs argument is a assem report file\nyou forgot to specify the out file")
	iout.close()

out=open(outf,'w')

names=["Genome/stat"]
ScfCnts=["ScfCnt"]
ScfMeans=["ScfMean"]
ScfSds=["ScfSd"]
ScfMins=["ScfMin"]
ScfMaxs=["ScfMax"]
ScfMeds=["ScfMed"]
ScfN50s=["ScfN50"]
ScfTots=["ScfTOT"]
NEvents=["NEvents"]
N1s=["N1"]
NAs=["NA"]
NBs=["NB"]
NCs=["NC"]
BKPs=["BKP"]
INVs=["INV"]
QIs=["INS"]
LNs=["LNKs"]

for file in sys.argv[2:]:
	names.append(os.path.basename(file))
	f=open(file,'r')
	line=f.readline()
	while line == '\n':
		line=f.readline()
	hre=re.search("-----Metassembly .+ statistics-----",line)
	if hre:
		line=f.readline()
	else:
		sys.exit("ERROR: %s is not a asseMerge report, no --Metassembly header" %(file))

	hre=re.search("^-+$",line)
	if hre:
		line=f.readline()
		scre=re.search("Scaffolds:\s+(\d+)",line)
		if scre:
			ScfCnts.append(scre.group(1))
		else:
			sys.exit("No Scaffolds in %s" %(file))

		line=f.readline()
		smre=re.search("length mean:\s+(\d+)",line)
		if smre:
			ScfMeans.append(smre.group(1))
		else:
			sys.exit("No length mean in "+file)

		line=f.readline()
		ssd=re.search("length sd:\s+(\d+)",line)
		if ssd:
			ScfSds.append(ssd.group(1))
		else:
			sys.exit("No length sd in "+file)

		line=f.readline()
		smin=re.search("length min:\s+(\d+)",line)
		if smin:
			ScfMins.append(smin.group(1))
		else:
			sys.exit("No length min in"+file)
		
		line=f.readline()
		smax=re.search("length max:\s+(\d+)",line)
		if smax:
			ScfMaxs.append(smax.group(1))
		else:
			sys.exit("No length max in "+file)

		line=f.readline()
		smed=re.search("length med:\s+(\d+)",line)
		if smed:
			ScfMeds.append(smed.group(1))
		else:
			sys.exit("No length med in "+file)

		line=f.readline()
		s50=re.search("N50:\s+(\d+)",line)
		if s50:
			ScfN50s.append(s50.group(1))
		else:
			sys.exit("No N50 in "+file)
		
		line=f.readline()
		stot=re.search("total length:\s+(\d+)",line)
		if stot:
			ScfTots.append(stot.group(1))
		else:
			sys.exit("No total length in "+file)

		line=f.readline()
		while line == '\n':
			line=f.readline()
		if not line  == "-----Events---------------\n":
			sys.exit("ERROR parsing %s no ----Events--- line" %(file))
		line=f.readline()
		h2re=re.search("^-+$",line)
		if not h2re:
			sys.exit("ERROR parsing %s, no \'---\' line"%(file))
		line=f.readline()
		nev=re.search("----Ns:Total number of events:\s+(\d+)",line)
		if nev:
			NEvents.append(nev.group(1))
		else:
			sys.exit("No ---Ns:total... in "+file)
		
		line=f.readline()
		n1=re.search("Case 1:\s+(\d+)",line)
		if n1:
			N1s.append(n1.group(1))
		else:
			sys.exit("No Case 1 in "+file)

		line=f.readline()
		nA=re.search("Case 2\.A:\s+(\d+)",line)
		if nA:
			NAs.append(nA.group(1))
		else:
			sys.exit("No Case A in "+file)

		line=f.readline()
		nB=re.search("Case 2\.B:\s+(\d+)",line)
		if nB:
			NBs.append(nB.group(1))
		else:
			sys.exit("No Case B in "+file)

		line=f.readline()
		nC=re.search("Case 2\.C:\s+(\d+)",line)
		if nC:
			NCs.append(nC.group(1))
		else:
			sys.exit("No Case C in "+file)

		line=f.readline()
		rel=re.search('^-+$',line)
		if not rel:
			sys.exit("ERROR parshin %s no ---- line"%(file))

		line=f.readline()
		bkp=re.search("BKP:\s+(\d+)",line)
		if bkp:
			BKPs.append(bkp.group(1))
		else:
			sys.exit("No BKP in "+file)

		line=f.readline()
		inv=re.search("INV:\s+(\d+)",line)
		if inv:
			INVs.append(inv.group(1))
		else:
			sys.exit("No INV in "+file)

		line=f.readline()
		qi=re.search("Insertions:\s+(\d+)",line)
		if qi:
			QIs.append(qi.group(1))
		else:
			sys.exit("No Insertions in "+file)

		line=f.readline()
		ln=re.search("Scf Joined-Links:\s+\d+\s+-\s+(\d+)",line)
		if ln:
			LNs.append(ln.group(1))
		else:
			sys.exit("No Scf Joined-Links in "+file)


	else:
		sys.exti("ERROR %s is not a asseMerge report, no \'----\' line" %(file))

PrintRow(out,names)		
PrintRow(out,ScfCnts)
PrintRow(out,ScfMeans)
PrintRow(out,ScfSds)
PrintRow(out,ScfMins)
PrintRow(out,ScfMaxs)
PrintRow(out,ScfMeds)
PrintRow(out,ScfN50s)
PrintRow(out,ScfTots)
PrintRow(out,NEvents)
PrintRow(out,N1s)
PrintRow(out,NAs)
PrintRow(out,NBs)
PrintRow(out,NCs)
PrintRow(out,BKPs)
PrintRow(out,INVs)
PrintRow(out,QIs)
PrintRow(out,LNs)










