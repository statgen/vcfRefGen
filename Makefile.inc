VERSION=0.1.3

########################
# Library Path:
#
# Default path to the statgen library.  You can either update your copy of
# this Makefile to be where you want to get the Makefile from or you can 
# overwrite LIB_PATH_GENERAL or LIB_PATH_VCF_REF_GEN.
CURRENT_PATH := $(dir $(lastword $(MAKEFILE_LIST)))
LIB_PATH_GENERAL ?=../libStatGen
LIB_PATH_VCF_REF_GEN ?= $(LIB_PATH_GENERAL)

# add any additional ../ as necessary if it is a relative path
ACTUAL_PATH =  $(patsubst ../%, $(CURRENT_PATH)../%, $(LIB_PATH_VCF_REF_GEN))
INCLUDE_MAKE = $(ACTUAL_PATH)/Makefiles/$(PARENT_MAKE)

########################
# Include the base Makefile
include $(INCLUDE_MAKE)

.phony: libStatGen

########################
# This is to handle the case where the statgen library is not already there.
$(INCLUDE_MAKE):
	@echo Unable to locate: $(INCLUDE_MAKE)
	@echo To change the location, set LIB_PATH_GENERAL or LIB_PATH_VCF_REF_GEN to the appropriate path to libStatGen.  Or specify \"make LIB_PATH_GENERAL=yourPath\" or \"make LIB_PATH_VCF_REF_GEN=yourPath\"
	@echo Use make libStatGen if you have git and want to download and build the current libStatGen at that location


libStatGen:
	@if test -d $(LIB_PATH_VCF_REF_GEN); \
	then echo $(LIB_PATH_VCF_REF_GEN) already exists; \
	else git clone git://github.com/statgen/libStatGen.git $(LIB_PATH_VCF_REF_GEN); fi
	$(MAKE) -C $(LIB_PATH_VCF_REF_GEN) --no-print-directory all; \
