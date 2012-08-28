EXE=vcfRefGen
########################
TOOLBASE = 
SRCONLY = Main.cpp
HDRONLY = 

COMPILE_ANY_CHANGE = Main
########################
DATE=$(shell date)
USER=$(shell whoami)
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DVERSION="\"${VERSION}\"" -DUSER="\"${USER}\""

########################
PARENT_MAKE = Makefile.src
include ../Makefile.inc

