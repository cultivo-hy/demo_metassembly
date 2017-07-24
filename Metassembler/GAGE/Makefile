#--------------------------------------------
# Makefile for Scripts
#
# Imported variables from top level Makefile:
# BIN_DIR ANMUM_DIR PIPE_DIR ZSTAT_DIR GC

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif 

ALL := gage

#-- PHONY rules --#
.PHONY: all clean

all: $(ALL) | $(BIN_DIR)

clean: | $(BIN_DIR)
	-rm -f $(ALL)
	cd $(BIN_DIR); rm -f $(ALL)

gage: gage.sh | $(BIN_DIR)
	$(SED)  -e 's?/bin/bash?$(BASH)?g' \
		-e 's?__GAGE_DIR?$(GAGE_DIR)?g' \
		gage.sh > $(BIN_DIR)/gage
	chmod 755 $(BIN_DIR)/gage

$(BIN_DIR):
	-mkdir $(BIN_DIR)
