#--------------------------------------------
# Makefile for Scripts
#
# Imported variables from top level Makefile:
# BIN_DIR ANMUM_DIR PIPE_DIR ZSTAT_DIR GC

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif 

ALL := RepGAGE RepMetassem RepStats RepREAPR RepCEGMA splitscafffa stats seqLengths

#-- PHONY rules --#
.PHONY: all clean

all: $(ALL) | $(BIN_DIR)

clean: | $(BIN_DIR)
	-rm -f $(ALL)
	cd $(BIN_DIR); rm -f $(ALL)

RepGAGE: RepGAGE.py | $(BIN_DIR)
	cp RepGAGE.py $(BIN_DIR)/RepGAGE
#	$(SED) -e 's?/usr/bin/python?$(PYTHON)?g' \
#		RepGAGE.py > $(BIN_DIR)/RepGAGE
	chmod 755 $(BIN_DIR)/RepGAGE

RepMetassem: RepMetassem.py | $(BIN_DIR)
	cp RepMetassem.py $(BIN_DIR)/RepMetassem
#	$(SED) -e 's?/usr/bin/python?$(PYTHON)?g' \
#		RepMetassem.py > $(BIN_DIR)/RepMetassem
	chmod 755 $(BIN_DIR)/RepMetassem

RepStats: RepStats.py | $(BIN_DIR)
	cp RepStats.py $(BIN_DIR)/RepStats
#	$(SED) -e 's?/usr/bin/python?$(PYTHON)?g' \
#		RepStats.py > $(BIN_DIR)/RepStats
	chmod 755 $(BIN_DIR)/RepStats

RepREAPR: RepREAPR.py | $(BIN_DIR)
	cp RepREAPR.py $(BIN_DIR)/RepREAPR
	chmod 755 $(BIN_DIR)/RepREAPR

RepCEGMA: RepCEGMA.py | $(BIN_DIR)
	cp RepCEGMA.py $(BIN_DIR)/RepCEGMA
	chmod 755 $(BIN_DIR)/RepCEGMA


stats: stats.pl | $(BIN_DIR)
	$(SED) -e 's?/usr/bin/perl -w?$(PERL) -w?g' \
	       stats.pl >$(BIN_DIR)/stats
#	cp stats.pl $(BIN_DIR)/stats
	chmod 755 $(BIN_DIR)/stats

splitscafffa: splitscafffa.pl | $(BIN_DIR)
#	cp splitscafffa.pl $(BIN_DIR)/splitscafffa
	$(SED) -e 's?/usr/bin/perl -w?$(PERL) -w?g' \
		splitscafffa.pl > $(BIN_DIR)/splitscafffa
	chmod 755 $(BIN_DIR)/splitscafffa

seqLengths: seqLengths.py | $(BIN_DIR)
	cp seqLengths.py $(BIN_DIR)/seqLengths
#	$(SED) -e 's?/usr/bin/python?$(PYTHON)?g' \
#		seqLengths.py > $(BIN_DIR)/seqLengths
	chmod 755 $(BIN_DIR)/seqLengths

$(BIN_DIR):
	-mkdir $(BIN_DIR)
