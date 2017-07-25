#!/bin/bash
## It is an install that assumes a new instance of AWS.
## If you want to install what you need, just select and install.
#==============================base==============================#
sudo apt-get install python
sudo apt-get install make
sudo apt-get update
sudo apt-get install g++
sudo apt-get install c++

#============================== Abyss 2.0 ==============================#
wget http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/2.0.2/abyss-2.0.2.tar.gz
tar -xvf abyss-2.0.2.tar.gz
rm abyss-2.0.2.tar.gz

#============================== sparsehash ==============================#
mkdir tool
cd tool
git clone https://github.com/sparsehash/sparsehash.git
cd ./sparsehash
./configure
make 
sudo make install  

#============================== openMPI ==============================#
cd ~/tool
wget https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz
tar -xvf openmpi-2.1.1.tar.gz
rm openmpi-2.1.1.tar.gz

cd openmpi-2.1.1
./configure 
sudo make install
sudo ldconfig

#============================== Boost ==============================#
cd ~/abyss-2.0.2
wget wget http://downloads.sourceforge.net/project/boost/boost/1.56.0/boost_1_56_0.tar.bz2
tar -xvf boost_1_56_0.tar.bz2 
rm boost_1_56_0.tar.bz2 
sudo apt-get update
#============================== Abyss 2.0 ==============================#
./configure
make
sudo make install
#============================== Abyss 1 ==============================#
sudo apt-get install abyss
sudo apt-get update

cd ~/tool
# ==============================option_SRAtoolKit ==============================#
wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
tar -xzf sratoolkit.current-centos_linux64.tar.gz
rm sratoolkit.current-centos_linux64.tar.gz
#cd ./sratoolkit.2.8.2-1-centos_linux64/bin
#./fastq-dump -I --split-files SRR800798
#./changePairEndFormat.py
mkdir ~/data
mv SRR800798_* ~/data
# you need to change format to use pair-end format.
# check github : https://github.com/cultivo-hy/demo_metassembly/dataFormat/changePairEndFormat.py
# change mode and run ./changePairEndFormat.py. will get pair-end data

cd ~
#============================== SOAPdenovo2 ==============================#
git clone https://github.com/aquaskyline/SOAPdenovo2.git 
cd SOAPdenovo2
make
#to run, you make config file " ****.config " see the sample of website or github
# configure file format	
# max_rd_len=110
# [LIB]
# avg_ins=200
# asm_flags=3
# reverse_seq=0
# rank=1
# q1=/home/ubuntu/data/SRR800798_1.fastq
# q2=/home/ubuntu/data/SRR800798_2.fastq

cd ~
#============================== Metassembly ==============================#
wget https://downloads.sourceforge.net/project/metassembler/v1.5/Metassembler.1.5.tar.gz
tar -xvf Metassembler.1.5.tar.gz
rm Metassembler.1.5.tar.gz
#============================== MUMmer ==============================#
cd ~/tool
wget https://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
tar -xvf MUMmer3.23.tar.gz
rm MUMmer3.23.tar.gz
cd ./MUMmer3.23
make check
make all
cd scripts
make
cd ~/tool/MUMmer3.23
sudo make install

#============================== bowtie2 ==============================#
sudo apt-get install bowtie2

#============================== samtool ==============================#
sudo apt-get install samtools


cd ~/Metassembler
make
sudo make install 

# now run using script. see the sample of website or github
