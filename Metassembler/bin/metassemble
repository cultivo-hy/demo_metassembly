#!/usr/bin/env python2.7
import subprocess
import sys
import re
import os
import string
import argparse

#------------------------------------ Bin Paths --------------------------------#
Ncoords="__Ncoords_PATH"
nucmer="__nucmer_PATH"
deltafilter="__delta-filter_PATH"
showcoords="__show-coords_PATH"
asseMerge="__asseMerge_PATH"
mateAn="__mateAn_PATH"
bwtbuild="__bwtbuild_PATH"
bwtaln="__bwtaln_PATH"
gage="__gage_PATH"
splitscaf="__splitscaf_PATH"
stats="__stats_PATH"
meta2fasta="__meta2fasta_PATH"
seqLengths="__seqLengths_PATH"
reportCE="__reportCE_PATH"
samtools="__samtools_PATH"
RepStats="__repstats_PATH"
RepMetassem="__repmerge_PATH"

#------------------------------------ Functions --------------------------------#
#------------------- runCmd ----------------------------------------------------#
def runCmd(desc,cmd,outf=None,errf=None,stdinf=None):
	
	print "\n---------- Run bash command ----------\n"
	print "%s:\n%s...\n" %(desc,string.join(cmd+(" ",)))
	try:
		rc=subprocess.call(cmd,stdout=outf,stderr=errf,stdin=stdinf)
	except OSError as e:
		if outf: outf.close()
		if errf: errf.close()
		if stdinf: stdinf.close()
		sys.exit("Error running %s\nExit status: %s\nError message: %s" %(string.join(cmd+(' ',)),e[0],e[1]))
	except:
		if outf: outf.close()
		if errf: errf.close()
		if stdinf: stdinf.close()
		sys.exit("Error running %s" %(string.join(cmd+(' ',))))
	print "----------------------------------------\n"

#------------------- parseConfig -----------------------------------------------#
def parseConfig(cfg):
        global Ncoords
        global nucmer
        global deltafilter
        global showcoords
        global asseMerge
        global mateAn
        global bwtbuild
        global bwtaln
        global gage
        global splitscaf
        global stats
        global meta2fasta
        global seqLengths
        global reportCE
	global samtools
	global RepStats
	global RepMetassem

        cfgp={}
	cfgp['global']={'bowtie2_threads':'1',
                        'bowtie2_read1':None,
                        'bowtie2_read2':None,
                        'bowtie2_maxins':None,
                        'bowtie2_minins':None,
                        'genomeLength':None,
                        'mateAn_file':'','mateAn_sam':'',
			'mateAn_A':'','mateAn_B':'',
			'mateAn_s':'','mateAn_m':'',
			'mateAn_b':'','mateAn_l':'',
			'mateAn_e':'0','mateAn_z':'3','mateAn_q':'20',
			'mateAn_N':'1','mateAn_f':'6','mateAn_c':'0.05',
			'mateAn_o':True,'mateAn_n':'30',
			'nucmer_l':'50','nucmer_c':'300',
			'asseMerge_z':'3','asseMerge_d':'2','asseMerge_i':'13',
			'asseMerge_s':'0','asseMerge_c':'15','asseMerge_e':'200',
			'asseMerge_l':'60','asseMerge_a':'60','asseMerge_L':'20',
			'asseMerge_t':'50','asseMerge_p':'0.65',
			'fasta':'','ID':'','delta':'','1delta':'','1coords':'',
			'meta2fasta_keepUnaligned':'0','meta2fasta_sizeUnaligned':['0','0'],
			'meta2fasta_keepDF':'0', 'meta2fasta_sizeDF':['0','0'],
			'done':False
			}
	cfgp['N_asms']=0

	ipars=''
	aparms=cfgp['global'].keys()
	bparms=['Ncoords','nucmer','delta-filter','show-coords','asseMerge','mateAn','bowtie2-build','bowtie2','gage','splitscaf','stats','meta2fasta','seqLenghts','reportCE','samtools','RepStats','RepMetassem']
	while 1:
		line=cfg.readline()
		if not line:
			break
		line=line.strip()
		line=line.rstrip('\n')
		if not line:
			continue
		if line[0] == '#':
			continue
		iasm=re.search('^\[(global|\d+)\]$',line)
		if iasm:
			ipars=iasm.group(1)
			if ipars != 'global':
				cfgp[ipars]=cfgp['global'].copy()
				cfgp['N_asms']+=1
			continue
		posp=re.search('^([^=]+)=.*',line) # Captures possible parametrs, may contain unrecognized parameters
		inp=re.search('^('+string.join(aparms+bparms,'|')+')=(.*)$',line) #Captures true parameters
		if inp:
			parm=inp.group(1)
			if parm == 'Ncoords':
				Ncoords=inp.group(2)
			elif parm == 'nucmer':
				nucmer=inp.group(2)
			elif parm == 'delta-filter':
				deltafilter=inp.group(2)
			elif parm == 'show-coords':
				showcoords=inp.group(2)
			elif parm == 'asseMerge':
				asseMerge=inp.group(2)
			elif parm == 'mateAn':
				mateAn=inp.group(2)
			elif parm == 'bowtie2-build':
				bwtbuild=inp.group(2)	
			elif parm == 'bowtie2':
				bwtaln=inp.group(2)
			elif parm == 'gage':
				gage=inp.group(2)
			elif parm == 'splitscaf':
				splitscaf=inp.group(2)
			elif parm == 'stats':
				stats=inp.group(2)
			elif parm == 'meta2fasta':
				meta2fasta=inp.group(2)
			elif parm == 'seqLengths':
				seqLengths=inp.group(2)
			elif parm == 'reportCE':
				reportCE=inp.group(2)
			elif parm == 'samtools':
				samtools=inp.group(2)
			elif parm == 'RepStats':
				RepStats=inp.group(2)
			elif parm == 'RepMetassem':
				RepMetassem=inp.group(2)
			elif aparms.count(parm) and ipars:
				if parm == 'mateAn_o':
					if inp.group(2) == 'False':
						cfgp[ipars]['mateAn_o']=False
					else:
						cfgp[ipars]['mateAn_o']=True
				elif parm == 'done':
					if inp.group(2) == 'False':
						cfgp[ipars]['done']=False
					elif inp.group(2) == 'True':
						cfgp[ipars]['done']=True
					else: 
						sys.exit('ERROR: %s is not a valid value for \'done\', it can only be set to \'True\' or \'False\'\n' %(inp.group(2)))
				elif parm == 'meta2fasta_sizeUnaligned' or parm=='meta2fasta_sizeDF' :
					cfgp[ipars][parm]=inp.group(2).strip().split()
				else:
					cfgp[ipars][parm]=inp.group(2)
					if ((parm == 'meta2fasta_keepUnaligned' or parm == 'meta2fasta_keepDF') and
                                            cfgp[ipars][parm] != '0' and cfgp[ipars][parm] != '1' and cfgp[ipars][parm] != '2'
                                            and cfgp[ipars][parm] != '3'):
						sys.exit('ERROR: %s %s is not a valid value for %s, possible values are:0,1,2, and 3\n' %(cfgp[ipars][parm][0],cfgp[ipars][parm][1],parm))
			elif aparms.count(parm) and not ipars:
				print "WARNING: %s was provided but could not be assigned to any header of [global], [1], ..." %(parm)

				
		elif posp:
			sys.exit("ERROR: Unrecognized parameter: \"%s\" in line:\n\"%s\"" %(posp.group(1),line))
		else:
			sys.exit("ERROR: Please check your configuration file. Remember that parameters are specified using the sintax:\nparm_name=value\nThe following line raised an error:\n\"%s\"\n" %(line))
			

	IDs={}
	for set in cfgp.keys():
		if (set == 'global' or
                    set == 'N_asms'):
			continue
		if(not cfgp[set]['bowtie2_read1'] or not cfgp[set]['bowtie2_read2']):
			sys.exit('ERROR: In assembly [%s]. You must specify the paths to the mate-pair libraries:\nbowtie2_read1=<path1>[...,<pathN>]\nbowtie2_read2=<path1>[...,<pathN>]\nin your configuration file.\nThese two parameters will be fed to bowtie2, pleas check the -1 and -2 options in their manual\n' %(set))
		if(not cfgp[set]['bowtie2_maxins'] or not cfgp[set]['bowtie2_minins']):
			sys.exit('ERROR: In assembly [%s]. You must specify the max and min insert length for the mate-pair libraries:\nmaxins=<int>\nminins=<int>\nin your configuration file.\nThese two parameters will be fed to bowtie2, please check the options --maxins and --minins in their manual\n' %(set))
		if(not cfgp[set]['delta'] and cfgp[set]['1delta']):
			print "WARNING: In assembly [%s]. A 1delta file was provided but the delta\nfile was not; the merging process needs a delta file so the\n WGA step will still have to be run.\nIf you want to skip the WGA you must provide a delta file.\n" %(set)		
		if(not cfgp[set]['delta'] and cfgp[set]['1coords']):
			print "WARNING: In assembly [%s]. A 1coords file was provided but the delta\nfile was not; the merging process needs a delta so the\n WGA step will still have to be run.\nIf you want to skip the WGA you must provide a delta file.\n" %(set)

		if not ( 
		    cfgp[set]['mateAn_file'] or 
		    (cfgp[set]['mateAn_A'] and cfgp[set]['mateAn_B']) or 
		    (cfgp[set]['mateAn_s'] and cfgp[set]['mateAn_m'])  
		   ):
			sys.exit('ERROR: There is not enough information to compute the CE-statistic for  [%s],\nyou must specify one of the following:\nmateAn_file=<path> or\nmateAn_A=<int> and mateAn_B=<int> or\nmateAn_m=<float> and mateAn_s=<float>\nin your configuration file\n' %(set))

		if not cfgp[set]['fasta']:
			sys.exit('ERROR: You must specify the path to assembly [%s]:\nfasta=<path>\nin your configuration file\n' %(set))
		if not cfgp[set]['ID']:
			cfgp[set]['ID']=os.path.basename(cfgp[set]['fasta'])
		if IDs.has_key(cfgp[set]['ID']):
			sys.exit('ERROR: duplicated ID.\nBy default the ID for assembly \'i\' is set to basename(<fastai>),\nhowever, you can specify any ID using:\nID=<string>\nin your configuration file\n')

	if not cfgp.has_key('1') or not cfgp.has_key('2'):
		sys.exit('ERROR: At least two assemblies must be merged, [1] and [2]\n')

	for x in range(1,cfgp['N_asms']+1):
		if not cfgp.has_key(str(x)):
			sys.exit('ERROR: The labels [1],[2],...,[N] should be consecutive, you are missing %s\n' %(str(x)))
	return(cfgp)
	
#------------------- CEstatistic -----------------------------------------------#
def CEstatistic(cfg_asm,outd='.',threads='1',read1='',read2='',maxins='',minins=''):
	outd.rstrip()
#	asmdir=outd+'/'+cfg_asm['ID']
	zdir=outd+'/CEstat'
	bwtdir=zdir+'/BWTaln'
	cedir=zdir+'/'+cfg_asm['ID']+'.ce'
	ceop=cedir+'/'+cfg_asm['ID']
	if cfg_asm['mateAn_file']:
		return(cfg_asm['mateAn_file'])
	if not os.path.isdir(asmdir):
		os.mkdir(asmdir)
	if not os.path.isdir(zdir):
		os.mkdir(zdir)
	if not os.path.isdir(cedir):
		os.mkdir(cedir)	
	cmdm=(mateAn,)
	if(cfg_asm['mateAn_A'] and cfg_asm['mateAn_B']):
		cmdm+=('-A',cfg_asm['mateAn_A'],'-B',cfg_asm['mateAn_B'])
	else:
		cmdm+=('-m',cfg_asm['mateAn_m'],'-s',cfg_asm['mateAn_s'])
	if not cfg_asm['mateAn_o']:
		cmdm+=('-o')
	cmdm+=(  '-p',ceop,'-e',cfg_asm['mateAn_e'],'-z',cfg_asm['mateAn_z'],
		'-q',cfg_asm['mateAn_q'],'-N',cfg_asm['mateAn_N'],'-f',cfg_asm['mateAn_f'],
		'-c',cfg_asm['mateAn_c'],'-n',cfg_asm['mateAn_n']
             )

	ceout=open(ceop+'.out','w')
	ceerr=open(ceop+'.err','w')

	if cfg_asm['mateAn_b'] and cfg_asm['mateAn_l']:
		cmdm+=('-b',cfg_asm['mateAn_b'],'-l',cfg_asm['mateAn_l'])
		runCmd(desc='mateAn -b -l',outf=out,errf=err,cmd=cmdm)
		ceout.close()
		ceerr.close()
		return(ceop+'.mateAn')

	if cfg_asm['mateAn_sam']:
		cmdm+=('-a',cfg_asm['mateAn_sam'])
		runCmd(desc='mateAn -a sam',outf=out,errf=err,cmd=cmdm)
		ceout.close()
		ceerr.close()

		#Gzip *.sort.bedpe file
		gzerr=open(ceop+".sort.bedpe.gz.err",'w')
		runCmd(desc="gzip %s.sort.bedpe file" %(ceop+".sort.bedpe"),cmd=("gzip","-q","-f",ceop+".sort.bedpe"))
		gzerr.close()

		#Convert sam to bam
		bamout=open(bwtp+'.mtp.bam','w')
		bamerr=open(bwtp+'.mtp.bam.err','w')
		runCmd(desc='convert sam to bam',outf=bamout,errf=bamerr,cmd=(samtools,'view','-Sb',bwtp+'.mtp.sam'))
		bamout.close()
		bamerr.close()
		os.remove(bwtp+'.mtp.sam')
		
		return(ceop+".mateAn")

	if not os.path.isdir(bwtdir):
		os.mkdir(bwtdir)

	if not (read1 and read2 and maxins and minins) or (maxins == '0' or minins == '0'):
		sys.exit('ERROR: Not enough information to compute CE-stat for %s.\nread1=%s\nread2=%s\nmaxins=%s\nminins=%s\n' %(read1,read2,maxins,minins))
	#Construct bowtie2 index
	bwtp=bwtdir+"/"+cfg_asm['ID']
	bldout=open(bwtp+'.bld.out','w')
	blderr=open(bwtp+'.bld.err','w')
	runCmd(desc='bowtie2-build',outf=bldout,errf=blderr,cmd=(bwtbuild,cfg_asm['fasta'],bwtp+'.bld'))
	bldout.close()
	blderr.close()
	#Align reads
	balout=open(bwtp+'.mtp.sam','w')
	balerr=open(bwtp+'.mtp.err','w')
	runCmd(desc='bowtie2',outf=balout,errf=balerr,cmd=(bwtaln,'-x',bwtp+'.bld','-1',read1,'-2',read2,
							   '--maxins',maxins,'--minins',minins,'--rf','--threads',threads))
	balout.close()
	balerr.close()
	#Compute CE-statistic
	cmdm+=('-a',bwtp+'.mtp.sam')
	runCmd(desc='mateAn -a sam',outf=ceout,errf=ceerr,cmd=cmdm)
	ceout.close()
	ceerr.close()

	#Gzip *.sort.bedpe file
	gzerr=open(ceop+".sort.bedpe.gz.err",'w')
	runCmd(desc="gzip %s.sort.bedpe file" %(ceop+".sort.bedpe"),cmd=("gzip","-q","-f",ceop+".sort.bedpe"))
	gzerr.close()

	#Convert sam to bam
	bamout=open(bwtp+'.mtp.bam','w')
	bamerr=open(bwtp+'.mtp.bam.err','w')
	runCmd(desc='convert sam to bam',outf=bamout,errf=bamerr,cmd=(samtools,'view','-Sb',bwtp+'.mtp.sam'))
	bamout.close()
	bamerr.close()
	os.remove(bwtp+'.mtp.sam')

	return(ceop+'.mateAn')

#------------------- isReadable ------------------------------------------------#
def isReadable(file):
	try:
		f=open(file,'r')
	except:
		return(False)
	else:
		f.close()
		return(True)

#------------------- wFastaStats ------------------------------------------------#
def wFastaStats(fasta=None,genomeLength=None):
	if not isReadable(fasta):
		sys.exit("wFastaStats: Cannot read %s" %(fasta))		
	mlengths=open(fasta+".lengths",'w')
	runCmd(desc="seqLengths",outf=mlengths,cmd=(seqLengths,fasta))
	mlengths.close()

	mlens2=open(fasta+".lengths.2",'w')
	runCmd(desc="cut seqLengths",outf=mlens2,cmd=("cut", "-f", "2",fasta+".lengths"))
	mlens2.close()

	mstats=open(fasta+".lengths.stats",'w')
	if genomeLength:
		runCmd(desc="stats",outf=mstats,cmd=(stats,"-n50",genomeLength,fasta+".lengths.2"))
	else:
		runCmd(desc="stats",outf=mstats,cmd=(stats,fasta+".lengths.2"))
	mstats.close()
	os.remove(fasta+".lengths.2")

#------------------- getScriptDir ------------------------------------------------#
def getScriptDir():
	path = os.path.realpath(sys.argv[0])
	if os.path.isdir(path):
		return path
	else:
		return os.path.dirname(path)

#------------------------------------ Main -------------------------------------#
sys.stdout=os.fdopen(sys.stdout.fileno(),'w',0)
parser=argparse.ArgumentParser(description='Metassemble a set of input assemblies.')
reqs=parser.add_argument_group('required arguments')
reqs.add_argument('--conf',type=file,metavar='<configuration file>',help="Path to configuration file. (Required)",required=True)
parser.add_argument('--outd',type=str,metavar='<out directory>',default='.')
parser.add_argument('--mgage',type=str,metavar='<mail>',help='For each intermediate metassembly send an email to <mail> when the merging step is done and the GAGE evaluation tool can be run by the user. This will only work if the Unix mail command is available and functioning', default='')
parser.add_argument('--dogage',type=str,metavar='<ref path>',help='For each intermediate metassembly run the GAGE evaluation tool using <ref path> as the reference.',default='')
args=parser.parse_args()
args.outd.rstrip()
if os.path.isfile(args.outd):
        sys.exit(args.outd+" already exists and is not a directory!")
if not os.path.isdir(args.outd):
	os.mkdir(args.outd)
GAGEan=args.mgage or args.dogage
cfg=parseConfig(args.conf)

initial_scfStats=[]
initial_ctgStats=[]

#Check binary paths
metassembleDir=getScriptDir()+"/"
if Ncoords == "__Ncoords_PATH": Ncoords=metassembleDir+"Ncoords"
if nucmer == "__nucmer_PATH": nucmer="nucmer"
if deltafilter == "__delta-filter_PATH": deltafilter="delta-filter"
if showcoords == "__show-coords_PATH": showcoords="show-coords"
if asseMerge == "__asseMerge_PATH": asseMerge=metassembleDir+"asseMerge"
if mateAn == "__mateAn_PATH": mateAn=metassembleDir+"mateAn"
if bwtbuild == "__bwtbuild_PATH": bwtbuild="bowtie2-build"
if bwtaln == "__bwtaln_PATH": bwtaln="bowtie2"
if gage == "__gage_PATH": gage=metassembleDir+"gage"
if splitscaf == "__splitscaf_PATH": splitscaf=metassembleDir+"splitscafffa"
if stats == "__stats_PATH": stats=metassembleDir+"stats"
if meta2fasta == "__meta2fasta_PATH": meta2fasta=metassembleDir+"meta2fasta"
if seqLengths == "__seqLengths_PATH": seqLengths=metassembleDir+"seqLengths"
if reportCE == "__reportCE_PATH": reportCE=metassembleDir+"reportCE"
if samtools == "__samtools_PATH": samtools="samtools"
if RepStats == "__repstats_PATH": RepStats=metassembleDir+"RepStats"
if RepMetassem == "__repmerge_PATH": RepMetassem=metassembleDir+"RepMetassem"

## Get permutation string, this will be used for some file names
permS=cfg['1']['ID']
for x in range(2,cfg['N_asms']+1):
	i=str(x)
	permS=cfg[i]['ID']+"."+permS

## Compute input assemblies CE-stat and contiguity stats
print("\n\nNumber of input assemblies: %d" %(cfg['N_asms']))
for x in range(1,cfg['N_asms']+1):
	i=str(x)
	asmdir=args.outd+'/'+cfg[i]['ID']
	asmp=asmdir+'/'+cfg[i]['ID']
	asmfa=asmp+'.fa'
	asmNs=asmfa+'.Ns'
	if not os.path.isdir(asmdir): os.mkdir(asmdir)
	if not isReadable(cfg[i]['fasta']): sys.exit('Cannot read %s' %(cfg[i]['fasta']))
	cfg[i]['fasta_nucmer']=cfg[i]['fasta']
	if isReadable(cfg[i]['fasta']+".Ns"):
		subprocess.call(('ln','-s',cfg[i]['fasta']+'.Ns',asmNs))
	subprocess.call(('ln','-s',cfg[i]['fasta'],asmfa))
	if not isReadable(cfg[i]['fasta']): sys.exit('Broken link for %s' %(cfg[i]['fasta']))
	if not isReadable(asmNs):
		runCmd('Ncoords',cmd=(Ncoords,'-N',asmfa,'-D',asmdir+'/'))
	cfg[i]['fasta']=asmfa
	cfg[i]['Ns']=asmNs
	### Get/Compute Ce-stat data
	zdir=asmdir+'/CEstat'
	cedir=zdir+'/'+cfg[i]['ID']+'.ce'
	if not os.path.isdir(zdir): os.mkdir(zdir)
	if not os.path.isdir(cedir): os.mkdir(cedir)
	if cfg[i]['mateAn_file']:
		cefile=cedir+'/'+cfg[i]['ID']+'.mateAn'
		if not isReadable(cfg[i]['mateAn_file']): sys.exit('ERROR: cannot read %s' %(cfg[i]['mateAn_file']))
		subprocess.call(('ln','-s',cfg[i]['mateAn_file'],cefile))
		cfg[i]['mateAn_file']=cefile
	elif cfg[i]['done']:
		cfg[i]['mateAn_file']=cedir+'/'+cfg[i]['ID']+'.mateAn'
	else:
		print('\n\n\n---- Computing CEstatistic for %s ----' %(cfg[i]['ID']))
		cfg[i]['mateAn_file']=CEstatistic(outd=asmdir,threads=cfg[i]['bowtie2_threads'],
						  read1=cfg[i]['bowtie2_read1'],read2=cfg[i]['bowtie2_read2'],
						  maxins=cfg[i]['bowtie2_maxins'],minins=cfg[i]['bowtie2_minins'],
						  cfg_asm=cfg[i])
	## Compute scaffold and coting stats
	if not (cfg[i]['done'] and isReadable(cfg[i]['fasta']+".lengths.stats")):
		wFastaStats(fasta=cfg[i]['fasta'],genomeLength=cfg[i]['genomeLength'])
	initial_scfStats.append(cfg[i]['fasta']+".lengths.stats")

	mctgs=open(asmp+".ctgs.fasta",'w')
	runCmd(desc="Get contigs",outf=mctgs,cmd=(splitscaf,cfg[i]['fasta']))
	mctgs.close()
	if not (cfg[i]['done'] and isReadable(asmp+".ctgs.fasta.lengths.stats")):
		wFastaStats(fasta=asmp+".ctgs.fasta",genomeLength=cfg[i]['genomeLength'])
	initial_ctgStats.append(asmp+".ctgs.fasta.lengths.stats")

## Report Initial Assemblies Stats
# Report Scf Stats
scf_error=open(args.outd+'/'+permS+'_InitialScfStats.err','w')
runCmd(desc="scaffold stats",errf=scf_error,cmd=(RepStats,args.outd+'/'+permS+'_InitialScfStats') + tuple(initial_scfStats))
scf_error.close()
# Report Ctg Stats
ctg_error=open(args.outd+'/'+permS+'_InitialCtgStats.err','w')
runCmd(desc='contig stats',errf=ctg_error,cmd=(RepStats,args.outd+'/'+permS+'_InitialCtgStats') + tuple(initial_ctgStats) )
ctg_error.close()

##Compute metassembly
#Initialize metassembly stats
merge_scfStats=[cfg['1']['fasta']+".lengths.stats"]
merge_ctgStats=[args.outd+"/"+cfg['1']['ID']+"/"+cfg['1']['ID']+".ctgs.fasta.lengths.stats"]
merge_mergeStats=[]

outmerge=args.outd+"/Metassembly"
if not os.path.isdir(outmerge):os.mkdir(outmerge)

#Compute pairwise merges
Rpath=cfg['1']['fasta']
RZstat=cfg['1']['mateAn_file']
if not cfg['2']['done'] and not isReadable(RZstat):
	sys.exit("No initial primary CE-stat")
	
rA=cfg['1']['ID']

for x in range(2,cfg['N_asms']+1):
	i=str(x)
	qA=cfg[i]['ID']
	mergeDir=outmerge+"/Q"+qA+"."+rA
	mergeBase="Q"+qA+"."+rA
	print("\n\n\n---- Merging %s and %s ==> %s" %(qA,rA,mergeBase))
	Qpath=cfg[i]['fasta']
	QZstat=cfg[i]['mateAn_file']
	
	if not isReadable(Qpath+".Ns"):
		runCmd(desc="Get Ncoords for %s"%(Qpath),
		       cmd=(Ncoords,"-N",Qpath,"-D",os.path.dirname(Qpath)+"/"))

	if not os.path.isdir(mergeDir):
		runCmd(desc="Create "+mergeDir,cmd=("mkdir",mergeDir))

	nucpre=mergeDir+"/"+mergeBase

	if not cfg[i]['delta'] and not cfg[i]['done']:
		runCmd(desc="nucmer "+Rpath+" "+Qpath,cmd=(nucmer,"--maxmatch","-l",cfg[i]['nucmer_l'],"-c",cfg[i]['nucmer_c'],
							    "-p",nucpre,
							    Rpath,Qpath))
		cfg[i]['delta']=nucpre+'.delta'
	
	if not (cfg[i]['1delta']) and not cfg[i]['done']:
		filt=open(nucpre+".1delta",'w')
		runCmd(desc="delta-filter %s.delta"%(nucpre),outf=filt,
			cmd=(deltafilter,"-1",cfg[i]['delta']))
		filt.close()
		cfg[i]['1delta']=nucpre+'.1delta'
	
	if not cfg[i]['1coords'] and not cfg[i]['done']:
		coords=nucpre+'.1coords'
		sco=open(coords,'w')
		runCmd(desc="show-coords %s.1delta"%(nucpre), outf=sco,
			cmd=(showcoords,"-cHlTr",cfg[i]['1delta']))
		sco.close()
		cfg[i]['1coords']=nucpre+'.1coords'

	if not os.path.isdir(mergeDir+"/M1"): os.mkdir(mergeDir+"/M1")

        OutMD=mergeDir+"/M1"
	OutM=OutMD+'/'+mergeBase

	if not cfg[i]['done']:

		for file in (RZstat, Rpath+'.Ns', QZstat, Qpath+'.Ns', cfg[i]['delta'], cfg[i]['1coords']):
			if not isReadable(RZstat): print "RZstat is None"
			if not isReadable(Rpath+'.Ns') : print "Rpath.Ns is None"
			if not isReadable(cfg[i]['delta']): print "cfg[%s]['delta'] is none" %(i)
			if not isReadable(cfg[i]['1coords']): print "cfg[%s]['1coords'] is None" %(i)
			if not isReadable(file):
				sys.exit("ERROR: cannot open %s.\n" %(file))
	        mout=open(OutM+".mout",'w')
	        merr=open(OutM+".merr",'w')
	        runCmd(desc="asseMerge ",outf=mout,errf=merr,cmd=(asseMerge, "-R",RZstat,"-r",Rpath+".Ns","-Q",QZstat,"-q",Qpath+".Ns","-D",cfg[i]['delta'],"-M",cfg[i]['1coords'],"-O",OutM,"-x",mergeBase,"-c",cfg[i]['asseMerge_c'],"-L",cfg[i]['asseMerge_L'], "-p",cfg[i]['asseMerge_p'],"-z",cfg[i]['asseMerge_z'],"-d",cfg[i]['asseMerge_d'],"-i",cfg[i]['asseMerge_i'],"-s",cfg[i]['asseMerge_s'],"-e",cfg[i]['asseMerge_e'],"-l",cfg[i]['asseMerge_l'],"-a",cfg[i]['asseMerge_a'],"-t",cfg[i]['asseMerge_t']))
		mout.close()
		merr.close()
	
		mfasta=open(OutM+".m2f.out",'w')
		m2fcmd=[]
		m2fcmd.append('--delta')
		m2fcmd.append(cfg[i]['delta']) 
		m2fcmd.append('--1delta')
                m2fcmd.append(cfg[i]['1delta'])
		m2fcmd.append('--keepUnaligned')
		m2fcmd.append(cfg[i]['meta2fasta_keepUnaligned'])
		m2fcmd.append('--sizeUnaligned')
		m2fcmd=m2fcmd + cfg[i]['meta2fasta_sizeUnaligned']
		m2fcmd.append('--keepDF')
		m2fcmd.append(cfg[i]['meta2fasta_keepDF'])
		m2fcmd.append('--sizeDF')
		m2fcmd=m2fcmd + cfg[i]['meta2fasta_sizeDF']
		m2fcmd=tuple(m2fcmd)
	        runCmd(desc="meta2fasta",outf=mfasta,cmd=(meta2fasta,'--out',OutM,
								     '--meta',OutM+".metassem",
								     '--fastas',Rpath,Qpath)+m2fcmd)
	        mfasta.close()
	
		runCmd(desc="Ncoords",cmd=(Ncoords,'-N',OutM+'.fasta','-D',OutMD+'/'))
	
		if GAGEan:
			gageDir=OutMD+'/GAGEan'
			if not os.path.isdir(gageDir): os.mkdir(gageDir)
			if args.mgage:
				mailout=open(gageDir+'/MAILGAGE','w')
				mailout.write("Run GAGEan for %s\n" %(mergeBase))
				mailout.write("cd %s; getCorrectnessStats.sh %s %s %s\n" %(gageDir,'<ref fasta>',OutM+'.ctgs.fasta',OutM+'.fasta'))
				mailout.close()
				mailin=open(gageDir+'/MAILGAGE','r')
				runCmd(desc='Send GAGE mail',stdinf=mailin,cmd=('mail','-s','Run GAGE for %s' %(mergeBase),args.mgage))
				mailin.close()
			if args.dogage:
				print("-- Running GAGE for %s" %(mergeBase))
				curDir=os.getcwd()
				os.chdir(gageDir)
				gagep=gageDir+'/'+mergeBase+'.gage'
				gageout=open(gagep+'.report','w')
				gageerr=open(gagep+'.err','w')
				runCmd(desc='GAGEan for %s' %(mergeBase),outf=gageout,errf=gageerr,cmd=(gage,args.dogage,OutM+'.ctgs.fasta',OutM+'.fasta'))
				gageout.close()
				gageerr.close()
				os.chdir(curDir)
	
	if not(cfg[i]['done'] and isReadable(OutM+".fasta.lengths.stats")):
		wFastaStats(fasta=OutM+".fasta",genomeLength=cfg[i]['genomeLength'])
	
	if not(cfg[i]['done'] and isReadable(OutM+".ctgs.fasta")):
		mctgs=open(OutM+".ctgs.fasta",'w')
		runCmd(desc="Get contigs",outf=mctgs,cmd=(splitscaf,OutM+".fasta"))
		mctgs.close()
	if not(cfg[i]['done'] and isReadable(OutM+".ctgs.fasta.lengths.stats")):
		wFastaStats(fasta=OutM+".ctgs.fasta",genomeLength=cfg[i]['genomeLength'])

	cfg[mergeBase]=cfg['1']
	cfg[mergeBase]['mateAn_file']=''
	cfg[mergeBase]['mateAn_b']=''
	cfg[mergeBase]['mateAn_l']=''
	cfg[mergeBase]['mateAn_a']=''
	cfg[mergeBase]['mateAn_sam']=''
	cfg[mergeBase]['delta']=''
	cfg[mergeBase]['1delta']=''
	cfg[mergeBase]['1coords']=''
	cfg[mergeBase]['ID']=mergeBase
	cfg[mergeBase]['fasta']=OutM+'.fasta'

	merge_scfStats.append(OutM+".fasta.lengths.stats")
	merge_ctgStats.append(OutM+".ctgs.fasta.lengths.stats")
	merge_mergeStats.append(OutM+".report")

	if (int(i) == cfg['N_asms'] and not isReadable(OutMD+'/CEstat/'+mergeBase+'.ce/'+mergeBase+'.mateAn') ) or \
	    (int(i) < cfg['N_asms'] and not cfg[str(int(i)+1)]['done'] and \
             not isReadable(OutMD+'/CEstat/'+mergeBase+'.ce/'+mergeBase+'.mateAn')):
		print("-- Computing CE-stat for %s" %(mergeBase))
	
		cfg[mergeBase]['mateAn_file']=CEstatistic( outd=OutMD,threads=cfg[mergeBase]['bowtie2_threads'],
						read1=cfg[mergeBase]['bowtie2_read1'],read2=cfg[mergeBase]['bowtie2_read2'],
						maxins=cfg[mergeBase]['bowtie2_maxins'],minins=cfg[mergeBase]['bowtie2_minins'],
						cfg_asm=cfg[mergeBase])

	elif (int(i) == cfg['N_asms'] and isReadable(OutMD+'/CEstat/'+mergeBase+'.ce/'+mergeBase+'.mateAn') ) or \
              (not cfg[str(int(i)+1)]['done'] and isReadable(OutMD+'/CEstat/'+mergeBase+'.ce/'+mergeBase+'.mateAn')):
		print "WARNING: Skipping CEstat computation for %s and using already existant %s file." %(mergeBase, \
								OutMD+'/CEstat/'+mergeBase+'.ce/'+mergeBase+'.mateAn')
		cfg[mergeBase]['mateAn_file']=OutMD+'/CEstat/'+mergeBase+'.ce/'+mergeBase+'.mateAn'
		
	Rpath=cfg[mergeBase]['fasta']
	RZstat=cfg[mergeBase]['mateAn_file']

	if int(i) == cfg['N_asms'] or not cfg[i]['done']:
		runCmd(desc="reportCE",cmd=(reportCE,   '-P',OutM+'.BKPs',\
						'-P',OutM+'.LNs', \
						'-P',OutM+'.N1s', \
						'-P',OutM+'.NAs', \
						'-P',OutM+'.NBs', \
						'-P',OutM+'.NCs', \
						'-P',OutM+'.QIs', \
						'-C',cfg[mergeBase]['mateAn_file']))

	rA=qA+'.'+rA

## Report Metassemblies stats
#Report Ctg Stats
m_ctgErr=open(args.outd+'/Metassembly/'+permS+"_metaCtgStats.err",'w')
runCmd(desc="meta ctg stats",errf=m_ctgErr,cmd=(RepStats,args.outd+'/Metassembly/'+permS+'_metaCtgStats')+tuple(merge_ctgStats))
m_ctgErr.close()
#Report Scf Stats
m_scfErr=open(args.outd+'/Metassembly/'+permS+'_metaScfStats.err','w')
runCmd(desc="meta scf stats",errf=m_scfErr,cmd=(RepStats,args.outd+'/Metassembly/'+permS+'_metaScfStats')+tuple(merge_scfStats))
m_scfErr.close()
#Report Merge Stats
m_mergeErr=open(args.outd+'/Metassembly/'+permS+'_metaMergeStats.err','w')
runCmd(desc="meta merge stats",errf=m_mergeErr,cmd=(RepMetassem,args.outd+'/Metassembly/'+permS+'_metaMergeStats')+tuple(merge_mergeStats))
m_mergeErr.close()
