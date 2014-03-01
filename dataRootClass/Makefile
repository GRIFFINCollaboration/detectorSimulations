#  Mhd le 26 Feb 2014 : compile TSpiceData.

# Makefile for the ROOT test programs.
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

#------------------------------------------------------------------------------
#SHARELIB      = libDetector.so \
#		libmyDetector1Data.so libmyDetector1.so \
#		libmyDetector2Data.so libmyDetector2.so \

SHARELIB      = libSpiceData.so 


all:            $(SHARELIB)

#------------------------------------------------------------------------------

libSpiceData.so:	TSpiceData.o   TSpiceDataDict.o
		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@

TSpiceDataDict.cpp:	TSpiceData.h
			rootcint -f $@ -c $^

#libmyDetector1Data.so:	TmyDetector1Data.o   TmyDetector1DataDict.o
#		$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
#
#TmyDetector1DataDict.cpp:	TmyDetector1Data.h
#			rootcint -f $@ -c $^
			
#------------------------------------------------------------------------------
clean:
	@rm -f core *~ *.o *Dict*

distclean:
	  make clean; rm -f *.so

.SUFFIXES: .$(SrcSuf)

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<


# dependances
TSpiceData.o:	TSpiceData.cpp	TSpiceData.h
#TmyDetector1Data.o:	TmyDetector1Data.cpp	TmyDetector1Data.h
#TmyDetector2Data.o:	TmyDetector2Data.cpp	TmyDetector2Data.h


