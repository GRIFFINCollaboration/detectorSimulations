# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

# Extensions des fichiers
ObjSuf        = o
SrcSuf        = cpp
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

# Compilateur et options
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC 
LD            = g++
SOFLAGS       = -shared 
LDFLAGS       = -O 

# Options de ROOT
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

#CXXFLAGS     += $(ROOTCFLAGS) $(OPTFLAGS)
CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS) 

# Specify Shared libraries
#-------------------------
#SHARELIB      = libDetector.so libNewData.so libNew2Data.so

SHARELIB      = libSpiceData.so  libS3Data.so libGriffinData.so

all:            $(SHARELIB)

# Create dictionaries
#-------------------------
libSpiceData.so:	TSpiceData.o   TSpiceDataDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TSpiceDataDict.cpp:	TSpiceData.h
			rootcint -f $@ -c $^

libS3Data.so:	TS3Data.o   TS3DataDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TS3DataDict.cpp:	TS3Data.h
			rootcint -f $@ -c $^
			
libGriffinData.so:	TGriffinData.o   TGriffinDataDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TGriffinDataDict.cpp:	TGriffinData.h
			rootcint -f $@ -c $^			

#libNewData.so:	TNewData.o   TNewDataDict.o
#		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
#
#TNewDataDict.cpp:	TNewData.h
#			rootcint -f $@ -c $^

#------------------------------------------------------------------------------
clean:
	@rm -f core *~ *.o *Dict*

distclean:
	  make clean; rm -f *.so

.SUFFIXES: .$(SrcSuf)

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<


# Manage dependances
#-------------------------
TSpiceData.o:	TSpiceData.cpp	TSpiceData.h
TS3Data.o:	TS3Data.cpp	TS3Data.h
TGriffinData.o:	TGriffinData.cpp	TGriffinData.h
#TNewData.o:	TNewData.cpp	TNewData.h