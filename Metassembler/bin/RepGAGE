#!/usr/bin/env python2.7

import os
import sys
import re

USAGE="USAGE: "+sys.argv[0]+" outName inGage1 inGage2 ..."

if len(sys.argv) < 3:
	sys.exit(USAGE)

def PrintRow(outfile,list):
	outfile.write(list[0])
	for x in list[1:]:
		outfile.write("\t"+x)
	outfile.write("\n")

out=sys.argv[1]
if os.access(out,os.F_OK):
	outi=open(out,'r')
	line=outi.readline()
	if line == "Contig Stats\n":
		sys.exit("ERROR: %s is a gage report file,\n you forgot to specify the out file" %(out))
	outi.close()

outfile=open(out,'w')

names=["File/Stat"]
gss=["Genome Size"]
Ass=["Assembly Size"]
MRBs=["Missing Ref Bases"]
MRBps=["Missing Ref Bases %"]
MABs=["Missing Assem Bases"]
MABps=["Missing Assem Bases %"]
MACs=["Missing Assem Ctgs"]
MACps=["Missing Assem Ctgs %"]
DRBs=["Dup Ref Bases"]
CRBs=["Comp Ref Bases"]
SNPs=["SNPs"]
IND5s=["Indels <5"]
IND6s=["Indels >=5"]
INVs=["Inversions"]
RELs=["Relocations"]
TRAs=["Translocations"]

CN50s=["Contig N50"]
CN50counts=["Contig N50 count"]
CMins=["Contig Min"]
CMaxs=["Contig Max"]
CEsizes=["Contig Mean"]

CCSN50s=["Corr Contig N50"]
CCSN50counts=["Corr Contig N50 count"]
CCSMins=["Corr Contig Min"]
CCSMaxs=["Corr Contig Max"]
CCEsizes=["Corr Contig Mean"]

SSN50s=["Scf N50"]
SSN50counts=["Scf N50 count"]
SSMins=["Scf Min"]
SSMaxs=["Scf Max"]
SSEsizes=["Scf Mean"]


CSN50s=["Corr Scf N50"]
CSN50counts=["Corr Scf N50 count"]
CSMaxs=["Corr Scf Max"]
CSEsizes=["Corr Scf Mean"]

for filename in sys.argv[2:]:

        names.append(os.path.basename(filename))
	#Rearrangement stats
	gs=None
	As=None
	MRB=None
	MRBp=None
	MAB=None
	MABp=None
	MAC=None
	MACp=None
	DRB=None
	CRB=None
	SNP=None
	IND5=None
	IND6=None
	INV=None
	REL=None
	TRA=None

	#Contig stats
	CN50=None
	CN50count=None
	CMin=None
	CMax=None
	CEsize=None
	
	#Corrected Contig Stats
	CCSN50=None
	CCSN50count=None
	CCSMin=None
	CCSMax=None
	CCEsize=None

	#Scaffold stats
	SSN50=None
	SSN50count=None
	SSMin=None
	SSMax=None
	SSEsize=None

	#Corrected Scaffold Stats
	CSN50=None
	CSN50count=None
	CSMax=None
	CSEsize=None

	contigFlag=False
	scfFlag=False

	file=open(filename,'r')

	while 1:
		line=file.readline()
		if not line: break
		
		if not gs: 
			gs=re.search("Genome Size:\s(.*)",line)
			if gs:	
				gs=gs.group(1)
				continue
		
		if not As:
			As=re.search("Assembly Size:\s(.*)",line)
			if As: 
				As=As.group(1)
				continue

		if not MRB:
			MRB=re.search("Missing Reference Bases:\s(\d+)\(([\d\.]+)%\)",line)
			if MRB:
				MRBp=MRB.group(2)
				MRB=MRB.group(1)
				continue

		if not MAB:
			MAB=re.search("Missing Assembly Bases:\s(\d+)\(([\d\.]+)%\)",line)
			if MAB:
				MABp=MAB.group(2)
				MAB=MAB.group(1)
				continue
		
		if not MAC:
			MAC=re.search("Missing Assembly Contigs:\s(\d+)\(([\d\.]+)%\)",line)
			if MAC:
				MACp=MAC.group(2)
				MAC=MAC.group(1)
				continue

		if not DRB:
			DRB=re.search("Duplicated Reference Bases:\s(.*)",line)
			if DRB:
				DRB=DRB.group(1)
				continue

		if not CRB:
			CRB=re.search("Compressed Reference Bases:\s(.*)",line)
			if CRB:
				CRB=CRB.group(1)
				continue

		if not SNP:
			SNP=re.search("SNPs:\s(.*)",line)
			if SNP:
				SNP=SNP.group(1)
				continue

		if not IND5:
			IND5=re.search("Indels < 5bp:\s(.*)",line)
			if IND5:
				IND5=IND5.group(1)
				continue

		if not IND6:
			IND6=re.search("Indels >= 5:\s(.*)",line)
			if IND6:
				IND6=IND6.group(1)
				continue

		if not INV:
			INV=re.search("Inversions:\s(.*)",line)
			if INV:
				INV=INV.group(1)
				continue

		if not REL:
			REL=re.search("Relocation:\s(.*)",line)
			if REL:
				REL=REL.group(1)
				continue

		if not TRA:
			TRA=re.search("Translocation:\s(.*)",line)
			if TRA:
				TRA=TRA.group(1)
				continue

		x=re.search("Contig Stats",line)
		if x:
			set="Ctg"

		x=re.search("Corrected Contig Stats",line)
		if x: 
			set="Corr Ctg"
	
		x=re.search("Scaffold Stats",line)
		if x:
			set="Scf"

		x=re.search("Corrected Scaffold Stats",line)
		if x:
			set="Corr Scf"

		if not CN50 and not CN50count and set == "Ctg":
			reCN50=re.search("N50:\s(.*)\sCOUNT:\s(\d*)",line)
			if reCN50:
				CN50=reCN50.group(1)
				CN50count=reCN50.group(2)
				continue


		if not CMin and set == "Ctg":
			CMin=re.search("Min:\s(.*)",line)
			if CMin:
				CMin=CMin.group(1)
				continue

		if not CMax and set == "Ctg":
			CMax=re.search("Max:\s(.*)",line)
			if CMax:
				CMax=CMax.group(1)
				continue

		if not CEsize and set == "Ctg":
			CEsize=re.search("E-size:([\d\.]+)",line)
			if CEsize:
				CEsize=CEsize.group(1)
				continue

		if not CCSN50 and not CCSN50count and set == "Corr Ctg":
			reCCSN50=re.search("N50:\s(.*)\sCOUNT:\s(\d*)",line)
			if reCCSN50:
				CCSN50=reCCSN50.group(1)
				CCSN50count=reCCSN50.group(2)
				continue

		if not CCSMin and set == "Corr Ctg":
			CCSMin=re.search("Min:\s(.*)",line)
			if CCSMin:
				CCSMin=CCSMin.group(1)
				continue

		if not CCSMax and set == "Corr Ctg":
			CCSMax=re.search("Max:\s(.*)",line)
			if CCSMax:
				CCSMax=CCSMax.group(1)
				continue

		if not CCEsize and set == "Corr Ctg":
			CCEsize=re.search("E-size:([\d\.]+)",line)
			if CCEsize:
				CCEsize=CCEsize.group(1)
				continue

		if not SSN50 and not SSN50count and set == "Scf":
			reSSN50=re.search("N50:\s(.*)\sCOUNT:\s(\d*)",line)
			if reSSN50:
				SSN50=reSSN50.group(1)
				SSN50count=reSSN50.group(2)
				continue

		if not SSMin and set == "Scf":
			SSMin=re.search("Min:\s(.*)",line)
			if SSMin:
				SSMin=SSMin.group(1)
				continue

		if not SSMax and set == "Scf":
			SSMax=re.search("Max:\s(.*)",line)
			if SSMax:
				SSMax=SSMax.group(1)
				continue

		if not SSEsize and set == "Scf":
			SSEsize=re.search("E-size:([\d\.]+)",line)
			if SSEsize:
				SSEsize=SSEsize.group(1)
				continue

		if not CSN50 and not CSN50count and not CSEsize and set == "Corr Scf":
			corrScfStats=re.search("([\d\.]+,){6}(\d+),(\d+),(\d+)",line)
			if corrScfStats:
				#CSMax=corrScfStats.group(3)
				CSN50=corrScfStats.group(3)
				CSN50count=corrScfStats.group(2)
				CSEsize=corrScfStats.group(4)


	if gs:	gss.append(gs)
	else: sys.exit("Couldn't find \"%s\" in %s" %(gss[0],file))

	if As:	Ass.append(As)
	else: sys.exit("Couldn't find \"%s\" in %s" %(Ass[0],file))

	if MRB: MRBs.append(MRB)
	else: sys.exit("Couldn't find \"%s\" in %s" %(MRBs[0],file))
	
	if MRBp: MRBps.append(MRBp)
	else: sys.exit("Couldn't find \"%s\" in %s" %(MRBps[0]),file)

	if MAB: MABs.append(MAB)
	else: sys.exit("Couldn't find \"%s\" in %s" %(MABs[0],file))

	if MABp: MABps.append(MABp)
	else: sys.exit("Couldn't find \"%s\" in %s" %(MABps[0],file))
	
	if MAC: MACs.append(MAC)
	else: sys.exit("Couldn't find \"%s\" in %s" %(MACs[0],file))

 	if MACp: MACps.append(MACp)
	else: sys.exit("Couldn't find \"%s\" in %s" %(MACps[0],file))

	if DRB: DRBs.append(DRB)
        else: sys.exit("Couldn't find \"%s\" in %s" %(DRBs[0],file))

	if CRB: CRBs.append(CRB)
        else: sys.exit("Couldn't find \"%s\" in %s" %(CRBs[0],file))

	if SNP: SNPs.append(SNP)
        else: sys.exit("Couldn't find \"%s\" in %s" %(SNPs[0],file))

	if IND5: IND5s.append(IND5)
        else: sys.exit("Couldn't find \"%s\" in %s" %(IND5s[0],file))
	
	if IND6: IND6s.append(IND6)
        else: sys.exit("Couldn't find \"%s\" in %s" %(IND6s[0],file))

	if INV: INVs.append(INV)
        else: sys.exit("Couldn't find \"%s\" in %s" %(INVs[0],file))

	if REL: RELs.append(REL)
        else: sys.exit("Couldn't find \"%s\" in %s" %(RELs[0],file))

	if TRA: TRAs.append(TRA)
        else: sys.exit("Couldn't find \"%s\" in %s" %(TRAs[0],file))

	if CN50: CN50s.append(CN50)
        else: sys.exit("Couldn't find \"%s\" in %s" %(CN50s[0],file))

	if CN50count: CN50counts.append(CN50count)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CN50counts[0],file))
	
	if CMin: CMins.append(CMin)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CMins[0],file))

	if CMax: CMaxs.append(CMax)
        else: sys.exit("Couldn't find \"%s\" in %s" %(CMaxs[0],file))

	if CEsize: CEsizes.append(CEsize)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CEsize[0],file))

	if CCSN50: CCSN50s.append(CCSN50)
        else: sys.exit("Couldn't find \"%s\" in %s" %(CCSN50s[0],file))

	if CCSN50count: CCSN50counts.append(CCSN50count)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CCSN50counts[0],file))
	
	if CCSMin: CCSMins.append(CCSMin)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CCSMins[0],file))

	if CCSMax: CCSMaxs.append(CCSMax)
        else: sys.exit("Couldn't find \"%s\" in %s" %(CCSMaxs[0],file))
	
	if CCEsize: CCEsizes.append(CCEsize)
	else: sys.exit("Couldn't fine \"%s\" in %s" %(CCEsizes[0],file))

	if SSN50: SSN50s.append(SSN50)
        else: sys.exit("Couldn't find \"%s\" in %s" %(SSN50s[0],file))

	if SSN50count: SSN50counts.append(SSN50count)
	else: sys.exit("Couldn't find \"%s\" in %s" %(SSN50counts[0],file))

	if SSMin: SSMins.append(SSMin)
        else: sys.exit("Couldn't find \"%s\" in %s" %(SSMins[0],file))

	if SSMax: SSMaxs.append(SSMax)
        else: sys.exit("Couldn't find \"%s\" in %s" %(SSMaxs[0],file))

	if SSEsize: SSEsizes.append(SSEsize)
	else: sys.exit("Couldn't find \"%s\" in %s" %(SSEsizes[0]),file)

#	if CSMax: CSMaxs.append(CSMax)
#	else: sys.exit("Couldn't find \"%s\" in %s" %(CSMaxs[0],file))

	if CSN50: CSN50s.append(CSN50)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CSN50s[0],file))
	
	if CSN50count: CSN50counts.append(CSN50count)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CSN50counts[0],file))

	if CSEsize: CSEsizes.append(CSEsize)
	else: sys.exit("Couldn't find \"%s\" in %s" %(CSEsizes[0],file))

	file.close()

PrintRow(outfile,names)
PrintRow(outfile,gss)
PrintRow(outfile,Ass)
PrintRow(outfile,MRBs)
PrintRow(outfile,MRBps)
PrintRow(outfile,MABs)
PrintRow(outfile,MABps)
PrintRow(outfile,MACs)
PrintRow(outfile,MACps)
PrintRow(outfile,DRBs)
PrintRow(outfile,CRBs)
PrintRow(outfile,SNPs)
PrintRow(outfile,IND5s)
PrintRow(outfile,IND6s)
PrintRow(outfile,INVs)
PrintRow(outfile,RELs)
PrintRow(outfile,TRAs)
PrintRow(outfile,CN50s)
PrintRow(outfile,CN50counts)
PrintRow(outfile,CMins)
PrintRow(outfile,CMaxs)
PrintRow(outfile,CEsizes)
PrintRow(outfile,CCSN50s)
PrintRow(outfile,CCSN50counts)
PrintRow(outfile,CCSMins)
PrintRow(outfile,CCSMaxs)
PrintRow(outfile,CCEsizes)
PrintRow(outfile,SSN50s)
PrintRow(outfile,SSN50counts)
PrintRow(outfile,SSMins)
PrintRow(outfile,SSMaxs)
PrintRow(outfile,SSEsizes)
PrintRow(outfile,CSN50s)
PrintRow(outfile,CSN50counts)
#PrintRow(outfile,CSMaxs)
PrintRow(outfile,CSEsizes)

outfile.close()

