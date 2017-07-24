#!/bin/bash
## Step by step script to metassemble assembly A and assembly B, using 
## assembly A as the primary assembly, and a 2kb mate-pair simulated library.

####### Create output directory "MergePipeline"

mkdir MergePipeline

####### Compute CE statistic for assembly A.fa and get positions where ambiguous bases are found (i.e. rows of Ns)

mkdir MergePipeline/A.CEstat
mkdir MergePipeline/A.CEstat/BWTaln
mkdir MergePipeline/A.CEstat/A.ce

bowtie2-build A.fa MergePipeline/A.CEstat/BWTaln/A
bowtie2 -x MergePipeline/A.CEstat/BWTaln/A -1 jump2k.1.fq -2 jump2k.2.fq --minins 1000 --maxins 3000 --rf 1> MergePipeline/A.CEstat/BWTaln/A.mtp.2k.sam 2> MergePipeline/A.CEstat/BWTaln/A.mtp.2k.err

../../bin/mateAn -A 1300 -B 2300 -a MergePipeline/A.CEstat/BWTaln/A.mtp.2k.sam -p MergePipeline/A.CEstat/A.ce/A

../../bin/Ncoords -N A.fa -D MergePipeline
####### Compute CE statistic for assembly B.fa and get positions where ambiguous bases are found (i.e. rows of Ns)

mkdir MergePipeline/B.CEstat
mkdir MergePipeline/B.CEstat/BWTaln
mkdir MergePipeline/B.CEstat/B.ce

bowtie2-build B.fa MergePipeline/B.CEstat/BWTaln/B
bowtie2 -x MergePipeline/B.CEstat/BWTaln/B -1 jump2k.1.fq -2 jump2k.2.fq --minins 1000 --maxins 3000 --rf 1> MergePipeline/B.CEstat/BWTaln/B.mtp.2k.sam 2> MergePipeline/B.CEstat/BWTaln/B.mtp.2k.err

../../bin/mateAn -A 1300 -B 2300 -a MergePipeline/B.CEstat/BWTaln/B.mtp.2k.sam -p MergePipeline/B.CEstat/B.ce/B

../../bin/Ncoords -N B.fa -D MergePipeline
####### Perform pairwise whole genome alignment and filtering

nucmer --maxmatch -c 20 -l 50 -p MergePipeline/B.A A.fa B.fa
delta-filter -1 MergePipeline/B.A.delta > MergePipeline/B.A.1delta
show-coords -cHlTr MergePipeline/B.A.1delta > MergePipeline/B.A.1coords

####### Metassemble

mkdir MergePipeline/M1
../../bin/asseMerge -R MergePipeline/A.CEstat/A.ce/A.mateAn -r MergePipeline/A.Ns -Q MergePipeline/B.CEstat/B.ce/B.mateAn -q MergePipeline/B.Ns -D MergePipeline/B.A.1delta -M MergePipeline/B.A.1coords -O MergePipeline/M1/B.A -x B.A -c 5

####### Construct Metassembly sequence

../../bin/meta2fasta --out MergePipeline/M1/B.A --meta MergePipeline/M1/B.A.metassem --fastas A.fa B.fa

