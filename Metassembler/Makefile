#--------------------------------------- 
# Top level Makefile for metassembly.
#
# CURDIR should be Metassembly/
#
#---------------------------------------


SHELL = /bin/sh

BIN_DIR := $(CURDIR)/bin
SRC_DIR := $(CURDIR)/src
SCRP_DIR := $(CURDIR)/Scripts
GAGE_DIR := $(CURDIR)/GAGE

BASH := $(filter /%,$(shell /bin/sh -c 'type bash'))
GC := $(filter /%,$(shell /bin/sh -c 'type g++'))
PERL := $(filter /%,$(shell /bin/sh -c 'type perl'))
SED := $(filter /%,$(shell /bin/sh -c 'type sed'))
NUC := $(filter /%,$(shell /bin/sh -c 'type nucmer'))
DELTA_FI := $(filter /%,$(shell /bin/sh -c 'type delta-filter'))
SHOW_CO := $(filter /%,$(shell /bin/sh -c 'type show-coords'))
DNADIF := $(filter /%,$(shell /bin/sh -c 'type dnadiff'))
MUMPLOT := $(filter /%,$(shell /bin/sh -c 'type mummerplot'))
#GAGE := $(filter /%,$(shell /bin/sh -c 'type getCorrectnessStats.sh'))
BWT2 := $(filter /%,$(shell /bin/sh -c 'type bowtie2'))
BWTBLD := $(filter /%,$(shell /bin/sh -c 'type bowtie2-build'))
#SPLITSCF := $(filter /%,$(shell /bin/sh -c 'type splitscafffa'))
PYTHON := $(filter /%,$(shell /bin/sh -c 'type python'))
SAMTOOLS := $(filter /%,$(shell /bin/sh -c 'type samtools'))


#-- Export these variable to other Makefiles
export BIN_DIR SHELL GC PERL SED NUC DELTA_FI SHOW_CO \
       DNADIF MUMPLOT BWT2 BWTBLD PYTHON BASH \
       GAGE_DIR SAMTOOLS #SPLITSCF GAGE

#-- PHONY rules --#
.PHONY: all check clean uninstall

all: Src scripts GAGEan | $(BIN_DIR)

check:
ifndef BIN_DIR
	@echo "ERROR: could not find bin directory"
endif
ifndef GC
	@echo "ERROR: 'g++' compiler not found"
endif
ifndef NUC
	@echo "ERROR: 'nucmer' not found"
endif
ifndef DELTA_FI
	@echo "ERROR: 'delta-filter' not found"
endif
ifndef SHOW_CO
	@echo "ERROR: 'show-coords' not found"
endif
#ifndef GAGE
#	@echo "ERROR: 'getCorrectnessStats.sh' not found"
#endif
ifndef PERL
	@echo "ERROR: 'perl' not found"
endif
ifndef SED
	@echo "ERROR: 'sed' not found"
endif
ifndef BWT2
	@echo "ERROR: 'bowtie2' not found"
endif
ifndef BWTBLD
	@echo "ERROR: 'bowtie2-build' not found"
endif
#ifndef SPLITSCF
#	@echo "ERROR: 'splitscafffa' not found"
#endif
ifndef SAMTOOLS
	@echo "ERROR: 'samtools' not found"
endif
	@echo "check complete"

clean:| $(BIN_DIR)
	-rm -rf $(BIN_DIR)

install: all

Src: | $(BIN_DIR)
	cd $(SRC_DIR); $(MAKE) all

scripts: | $(BIN_DIR)
	cd $(SCRP_DIR); $(MAKE) all

GAGEan: | $(BIN_DIR)
	cd $(GAGE_DIR); $(MAKE) all

uninstall: clean

$(BIN_DIR):
	-mkdir $(BIN_DIR)

#-- End of Makefile --#

