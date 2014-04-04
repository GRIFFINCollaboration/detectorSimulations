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

SHARELIB      = libSpiceData.so  libS3Data.so libGriffinData.so\
				libFragment.so  libTigFragment.so

all:            $(SHARELIB)

# Create dictionaries
#-------------------------

# S P I C E
libSpiceData.so:	TSpiceData.o   TSpiceDataDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TSpiceDataDict.cpp:	TSpiceData.h
			rootcint -f $@ -c $^

# S 3
libS3Data.so:	TS3Data.o   TS3DataDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TS3DataDict.cpp:	TS3Data.h
			rootcint -f $@ -c $^

# G R I F F I N 
libGriffinData.so:	TGriffinData.o   TGriffinDataDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TGriffinDataDict.cpp:	TGriffinData.h
			rootcint -f $@ -c $^			

# F R A G M E N T 	(G R S I S P O O N) 
libFragment.so:	TFragment.o   TFragmentDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TFragmentDict.cpp:	TFragment.h
			rootcint -f $@ -c $^

# T I G  F R A G M E N T 	(G R S I S P O O N) 
libTigFragment.so:	TTigFragment.o   TTigFragmentDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TTigFragmentDict.cpp:	TTigFragment.h
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


# Manage dependancies
#-------------------------
TFragment.o:		TFragment.cpp		TFragment.h
TTigFragment.o:		TTigFragment.cpp	TTigFragment.h
TSpiceData.o:		TSpiceData.cpp		TSpiceData.h
TS3Data.o:			TS3Data.cpp			TS3Data.h
TGriffinData.o:		TGriffinData.cpp	TGriffinData.h
