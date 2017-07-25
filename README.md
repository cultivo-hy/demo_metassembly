demo_metassembly
=====

demo_metassembly is a *de novo* sequence metassembler intended for short paired-end
reads and large genomes. you can install this tool simply by using ./install.sh (only provided now)
it can install abyss1 , abyss2 , SOAPdenovo2 , Metassembly and another related tool.

Contents
========

* [Quick Start](#quick-start)
	* [Install demo_metassembly on Debian or Ubuntu](#install-demo_metassembly-on-debian-or-ubuntu)
	* [Install demo_metassembly on Mac OS X](#install-demo_metassembly-on-mac-os-x) : not currently provided
* [Dependencies](#dependencies)
* [Usage](#usage)
Quick Start
===========

## Install demo_metassembly on Debian or Ubuntu

Run the command

	sudo apt-get install demo_metassembly  : not currently provided

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

Usage
=====
I will update each detailed method to execute.

Abyss 
See here https://github.com/bcgsc/abyss
SOAPdenovo2 
See here https://github.com/aquaskyline/SOAPdenovo2
Metassembler
see here https://sourceforge.net/p/metassembler/code/ci/master/tree/Sample/ and follow.
