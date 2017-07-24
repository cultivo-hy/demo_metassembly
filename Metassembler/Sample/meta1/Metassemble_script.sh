#!/bin/bash
## Script that uses the metassemble.py wrapper to metassemble assembly A and B, using
## assembly A as the primary assembly, and a 2kb mate-pair simulated library.

### Generate configuration file that will be input for metassemble.py
echo -e "\
############################################\n\
###   Metassemble A and B configuration file\n\
############################################\n\
[global]\n\
\n\
bowtie2_threads=1\n\
bowtie2_read1=$(pwd)/jump2k.1.fq\n\
bowtie2_read2=$(pwd)/jump2k.2.fq\n\
bowtie2_maxins=3000\n\
bowtie2_minins=1000\n\
\n\

\n\
mateAn_A=1300\n\
mateAn_B=2300\n\
\n\
[1]\n\
\n\
fasta=$(pwd)/A.fa\n\
ID=A\n\
\n\
[2]\n\
\n\
fasta=$(pwd)/B.fa\n\
ID=B\n\
" > B.A.metassemble.config

### Run metassemble
../../bin/metassemble --conf B.A.metassemble.config --outd ./MergeMetassemble  

mv B.A.metassemble.config ./MergeMetassemble
