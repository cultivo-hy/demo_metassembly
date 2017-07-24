#!/usr/bin/env python2.7
import subprocess
import sys
import re
import os
import string
import argparse

#------------------------------------ Bin Paths --------------------------------#
mateAn="__mateAn_PATH"
bwtbuild="__bwtbuild_PATH"
bwtaln="__bwtaln_PATH"
samtools="__samtools_PATH"

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

#------------------- getScriptDir ------------------------------------------------#
def getScriptDir():
	path = os.path.realpath(sys.argv[0])
	if os.path.isdir(path):
		return path
	else:
		return os.path.dirname(path)

#------------------------------------ Main -------------------------------------#
parser=argparse.ArgumentParser(description='Run pipeline for CE-stat computation')
bwt=parser.add_argument_group('Bowtie2 arguments.\n Please refer to bowtie2 -h for further help')
bwt.add_argument('--fas',type=str,metavar='<str>',required=True,help="fasta file")
bwt.add_argument('--max',type=str,metavar='<int>',required=True,help="bowtie2 --maxins")
bwt.add_argument('--min',type=str,metavar='<int>',required=True,help="bowtie2 --minins")
bwt.add_argument('--r1',type=str,metavar='<str>',required=True,help="bowtie2 -1")
bwt.add_argument('--r2',type=str,metavar='<str>',required=True,help="bowtie2 -2")
bwt.add_argument('--thr',type=str,metavar='<int>',required=False,default='1',help="bowtie2 --threads")
mateAnArgs=parser.add_argument_group('mateAn arguments.\n Either --A and --B, or --m and --s must be specified.\n Please refer to mateAn -h for further help')
mateAnArgs.add_argument('--A',type=str,metavar='<mateAn A>',required=False)
mateAnArgs.add_argument('--B',type=str,metavar='<mateAn B>',required=False)
mateAnArgs.add_argument('--m',type=str,metavar='<mateAn m>',required=False)
mateAnArgs.add_argument('--s',type=str,metavar='<mateAn s>',required=False)
mateAnArgs.add_argument('--e',type=str,metavar='<mateAn e>',required=False,default='0', help='only print the CE-statistic for chromosome positions that are at least --e bases\napart from the edges')
mateAnArgs.add_argument('--z',type=str,metavar='<mateAn z>',required=False,default='3', help='zstat threshold for considering deviations from the mean as significant')
mateAnArgs.add_argument('--q',type=str,metavar='<mateAn q>',required=False,default='20', help='minimum mapping quality')
mateAnArgs.add_argument('--N',type=str,metavar='<mateAn N>',required=False,default='1', help='only print the CE-statistic for positions with coverage >= N')
mateAnArgs.add_argument('--f',type=str,metavar='<mateAn f>',required=False,default='6', help='filter out inserts with at least --f deviations from the mean')
mateAnArgs.add_argument('--o',action='store_true',required=False,help='If specified, do not adjust mu0 and sd0')
mateAnArgs.add_argument('--c',type=str,metavar='<mateAn c>',required=False,default='0.05', help='adjust mu0 and sd0 until |mean CE-stat| <= c')
mateAnArgs.add_argument('--n',type=str,metavar='<mateAn n>',required=False,default='30', help='when computing global CE-statistic mean and sd, only consider positions with coverage >= n')
out=parser.add_argument_group('Output directory and prefix')
out.add_argument('--outp',type=str,metavar='<str>',help='--outp will be used to name some output files',required=True)
out.add_argument('--outd',type=str,metavar='<str>',help='Output will be written under outd/CEstat',default='.')
args=parser.parse_args()
outd=args.outd
outd.rstrip()

if not ( (args.A and args.B) or (args.m and args.s) ):
	sys.exit('Either --A and --B, or --m and --s must be specified')

if (outd != '.' or outd != './') and not os.path.isdir(outd):
	os.mkdir(outd)

#Check binary paths
metassembleDir=getScriptDir()+"/"
if mateAn == "__mateAn_PATH": mateAn=metassembleDir+"mateAn"
if bwtbuild == "__bwtbuild_PATH": bwtbuild="bowtie2-build"
if bwtaln == "__bwtaln_PATH": bwtaln="bowtie2"
if samtools == "__samtools_PATH": samtools="samtools"

## CEstatistic 
if outd == '.' or outd == './': outd=os.getcwd()
outd.rstrip()
zdir=outd+'/CEstat'
bwtdir=zdir+'/BWTaln'
cedir=zdir+'/'+args.outp+'.ce'
ceop=cedir+'/'+args.outp
if not os.path.isdir(zdir):
	os.mkdir(zdir)
if not os.path.isdir(cedir):
	os.mkdir(cedir)	
cmdm=(mateAn,)
if(args.A and args.B):
	cmdm+=('-A',args.A,'-B',args.B)
else:
	cmdm+=('-m',args.m,'-s',args.s)
if args.o:
	cmdm+=('-o')
cmdm+=(  '-p',ceop,'-e',args.e,'-z',args.z,
	'-q',args.q,'-N',args.N,'-f',args.f,
	'-c',args.c,'-n',args.n
            )
ceout=open(ceop+'.out','w')
ceerr=open(ceop+'.err','w')
if not os.path.isdir(bwtdir):
	os.mkdir(bwtdir)
if not (args.r1 and args.r2 and args.max and args.min) or (args.max == '0' or args.min == '0'):
	sys.exit('ERROR: Not enough information to compute CE-stat for %s.\nr1=%s\nr2=%s\nmax=%s\nmin=%s\n' %(args.r1,args.r2,args.max,args.min))
#Construct bowtie2 index
bwtp=bwtdir+"/"+args.outp

bldcmd=open(bwtp+'.bld.cmd','w')
bldcmd.write(string.join((bwtbuild,args.fas,bwtp+'.bld',' ',)))
bldcmd.close()

bldout=open(bwtp+'.bld.out','w')
blderr=open(bwtp+'.bld.err','w')
runCmd(desc='bowtie2-build',outf=bldout,errf=blderr,cmd=(bwtbuild,args.fas,bwtp+'.bld'))
bldout.close()
blderr.close()
#Align reads
balcmd=open(bwtp+'.mtp.cmd','w')
balcmd.write(string.join((bwtaln,'-x',bwtp+'.bld','-1',args.r1,'-2',args.r2,
                          '--maxins',args.max,'--minins',args.min,'--rf','--threads',args.thr,' ',)))
balcmd.close()
balout=open(bwtp+'.mtp.sam','w')
balerr=open(bwtp+'.mtp.err','w')
runCmd(desc='bowtie2',outf=balout,errf=balerr,cmd=(bwtaln,'-x',bwtp+'.bld','-1',args.r1,'-2',args.r2,
						   '--maxins',args.max,'--minins',args.min,'--rf','--threads',args.thr))
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

