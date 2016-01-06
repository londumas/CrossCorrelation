#===================================================================================
#
#         FILE: Makefile
#
#        USAGE: ---
#
#  DESCRIPTION: ---
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Hélion du Mas des Bourboux, Helion.du-Mas-des-Bourboux@cea.fr
#      COMPANY: CEA (France)
#      VERSION: ---
#      CREATED: ---
#     REVISION: ---
#===================================================================================

#WORKDIR = /home/helion/Documents/Thèse
WORKDIR = /home/usr201/mnt/hdumasde

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -O3 -Wall -fPIC 
#CXXFLAGS      = -g -Wall -fPIC 
LD            = g++
#LDFLAGS       = -g
LDFLAGS       = -O3
SOFLAGS       = -shared

#CFITSIO (to be defined in .cshrc)
CFITSIO = /home/nfs/manip/mnt/bao/EXTLibs/cfitsio

CXXFLAGS     += $(ROOTCFLAGS) -I $(CFITSIO)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS) -lMinuit  -lTreePlayer -lMLP  -L$(CFITSIO) -lcfitsio

# ---------------------------------------------------------------------

BIN = ../bin
OBJ = ../obj
ROOT = $(WORKDIR)/Code/Root/Library
MATH = $(WORKDIR)/Code/Cpp/Library

INC = Constants.h GetDelta.h $(ROOT)/RootHistoFunctions.h $(MATH)/mathFunctions.h 
#TTreeSpectra.h
MYLIBS = $(LIBS) -lc $(GLIBS)

# ---------------------------------------------------------------------
# lc : library of c libc.so
# lm : library math

all: main

echo:
	@echo GLIBS = $(GLIBS)
	@echo LIBS = $(LIBS)
	@echo CXXFLAGS = $(CXXFLAGS)

main: $(BIN)/main.exe
	@echo "building main.exe"

$(BIN)/main.exe: $(OBJ)/main.o $(OBJ)/GetDelta.o $(OBJ)/RootHistoFunctions.o $(OBJ)/mathFunctions.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(MYLIBS)

$(OBJ)/main.o: main.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ main.cpp

$(OBJ)/GetDelta.o: GetDelta.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ GetDelta.cpp
	
$(OBJ)/RootHistoFunctions.o: $(ROOT)/RootHistoFunctions.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ $(ROOT)/RootHistoFunctions.cpp
	
$(OBJ)/mathFunctions.o: $(MATH)/mathFunctions.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ $(MATH)/mathFunctions.cpp

clean:
	rm -f $(OBJ)/*.o $(BIN)/main.exe
	
	
