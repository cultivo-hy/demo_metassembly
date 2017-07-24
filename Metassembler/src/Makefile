#------------------------------------------------
# Makefile for src/
#

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif

BIN_RULE = $(GC) $^ -o $(BIN_DIR)/$@; chmod 755 $(BIN_DIR)/$@
OBJ_RULE = $(GC) $< -c -o $@

ALL := mateAn asseMerge reportCE metassemble meta2fasta Ncoords computeCE

### PHONY rules :) ###
.PHONY: all clean

all: $(ALL) | $(BIN_DIR)

clean: | $(BIN_DIR)
	-rm -f *.o
	cd $(BIN_DIR) 
	-rm -f $(ALL)

### not PHONY :( ###

Ncoords: Ncoords.cpp | $(BIN_DIR)
	$(BIN_RULE)

CEstat.o: CEstat.cpp CEstat.hh
	$(OBJ_RULE)

merge.o: merge.cpp merge.hh
	$(OBJ_RULE)

mateAn: mateAn.cpp CEstat.o
	$(BIN_RULE)

delta.o: delta.cc delta.hh
	$(OBJ_RULE)

asseMerge: asseMerge.cpp CEstat.o merge.o delta.o
	$(BIN_RULE)

reportCE: reportCE.cpp CEstat.o
	$(BIN_RULE)

metassemble: metassemble.py
	cp metassemble.py $(BIN_DIR)/metassemble
	chmod 755 $(BIN_DIR)/metassemble

computeCE: computeCE.py
	cp computeCE.py $(BIN_DIR)/computeCE
	chmod 755 $(BIN_DIR)/computeCE

meta2fasta: meta2fasta.pl | $(BIN_DIR)
	$(SED) -e 's?/usr/bin/perl?$(PERL)?g' \
               meta2fasta.pl > $(BIN_DIR)/meta2fasta
	chmod 755 $(BIN_DIR)/meta2fasta

$(BIN_DIR):
	-mkdir $(BIN_DIR)
