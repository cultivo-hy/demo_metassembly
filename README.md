demo_metassembly
=====

demo_metassembly is a de novo sequence metassembler intended for short paired-end reads and large genomes. you can install this tool simply by using ./install.sh (only provided now) it consist of abyss1 , abyss2 , SOAPdenovo2 , Metassembly and another related tool(nvBowtie,kind of machine learning -not packaged yet)

Contents
========

* [Quick Start](#quick-start)
	* [Install demo_metassembly on Debian or Ubuntu](#install-demo_metassembly-on-debian-or-ubuntu) : currently not provided
	* [Install demo_metassembly on Mac OS X](#install-demo_metassembly-on-mac-os-x) : currently not provided
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Usage](#usage)

Quick Start
===========

## Install demo_metassembly on Debian or Ubuntu

Run the command

	sudo apt-get install demo_metassembly  : currently not provided

or download and install the
[Debian package](http://www.bcgsc.ca/platform/bioinfo/software/demo_metassembly).


Dependencies
============

demo_metassembly requires the following libraries:

* [Boost](http://www.boost.org/)
* [Open MPI](http://www.open-mpi.org)
* [sparsehash](https://code.google.com/p/sparsehash/)
* [SOAPdenovo2](http://soap.genomics.org.cn/soapdenovo.html)
* [Metassember](https://sourceforge.net/projects/metassembler/)
* [MUMmer](http://mummer.sourceforge.net/)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtool](http://www.htslib.org/)

Option
* [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)
* [Quast](http://quast.bioinf.spbau.ru/)

Installation
============
This pipeline is executed in Linux system. so i recommend to use in linux and use AWS instance becuast it has many depedency between program, so it makes many error. you just copy install.sh and make file to run. it will automatically install pipeline needed. or see install.sh and follow the script.

Usage
=====
If you have a data, skip this part. using sratoolkit, download genome data. 
install sratoolkit in ~/tool:
	
	tar -xzf sratoolkit.current-centos_linux64.tar.gz
	cd ~/tool/sratoolkit.2.8.2-1-centos_linux64/bin
	./fastq-dump -I --split-files SRR800798		## it make paired end sequencing data
	
and it need to change format of data to use. you can change format using changePairEndFormat.py in Format folder.

* [Abyss](https://github.com/bcgsc/abyss)

If pipeline is installed well, you can use abyss 1 & abyss 2. 
basic commands to run abyss 1 is :

	abyss-pe name=**** k=** in="/home/ubuntu/data/genome_1 /home/ubuntu/data/genome_2"
	
and if you want to run it parallel, add 'np' command 'np=16' :

	abyss-pe name=**** k=** np=16 in="/home/ubuntu/data/genome_1 /home/ubuntu/data/genome_2"
	
It will be parallelized automatically using mpi and multithreading. 
also you can run abyss2 by adding commands "B=100M H=3 kc=3":

	abyss-pe name=**** k=** np=16 in="/home/ubuntu/data/genome_1 /home/ubuntu/data/genome_2" B=100M H=3 kc=3 v=-v
	
then you can get result of scaffolds. check details of parameter in Abyss github

* [SOAPdenovo2](https://github.com/aquaskyline/SOAPdenovo2)

If pipeline is installed well, you can use SOAPdenovo63mer & SOAPdenovo127mer. 
127mer need more memory than 63mer because of detailed execution. 
before you run SOAPdenovo2, You make your configure file like this:

	max_rd_len=110 
	[LIB]
	avg_ins=200
	asm_flags=3
	reverse_seq=0
	rank=1
	q1=/home/ubuntu/data/SRR800798_1.fastq
	q2=/home/ubuntu/data/SRR800798_2.fastq
	
and then, run using SOAPdenovo***mer :
	
	./SOAPdenovo***mer all -s configureName.config -o outputName 1>ass.log 2>ass.err
	
if you want to run it parllel, add command '-p' :

	./SOAPdenovo***mer all -s configureName.config -p 16 -o outputName 1>ass.log 2>ass.err
	
then you can get result of scaffolds. check details of parameter in SOAPdenovo2 github

* [Metassembler](https://sourceforge.net/p/metassembler)

if pipeline is installed well, you can use Metassembler. it can be executed by shell script file. 
See file 'executeMetaseembly.sh'. you just change some value related with your data repository and tool repository.   
If you want to run it parallel on bottlenecks, open file and add this command "-p 16" in bowtie2 :

	bowtie2 -x MergePipeline/A.CEstat/BWTaln/A -1 jump2k.1.fq -2 jump2k.2.fq -p 16 --minins 200 ~~.....
	
then you can get result of Metassembler. check details of code in Metasembler Sourceforge
