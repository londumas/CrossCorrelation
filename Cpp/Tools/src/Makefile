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
ROOT = ../../../Root/Library
MATH = ../../../Cpp/Library
CONSTANTS = ../../../Constants

INC = tools.h $(ROOT)/RootHistoFunctions.h $(MATH)/mathFunctions.h $(CONSTANTS)/globalValues.h ../../../chain_annalys_delta/Correlation/src/Cosmology.h
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

$(BIN)/main.exe: $(OBJ)/main.o $(OBJ)/tools.o $(OBJ)/RootHistoFunctions.o $(OBJ)/mathFunctions.o $(OBJ)/Cosmology.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(MYLIBS)

$(OBJ)/main.o: main.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ main.cpp

$(OBJ)/tools.o: tools.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ tools.cpp

$(OBJ)/Cosmology.o: ../../../chain_annalys_delta/Correlation/src/Cosmology.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ ../../../chain_annalys_delta/Correlation/src/Cosmology.cpp
	
$(OBJ)/RootHistoFunctions.o: $(ROOT)/RootHistoFunctions.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ $(ROOT)/RootHistoFunctions.cpp
	
$(OBJ)/mathFunctions.o: $(MATH)/mathFunctions.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c -o $@ $(MATH)/mathFunctions.cpp

clean:
	rm -f $(OBJ)/*.o $(BIN)/main.exe
	
	
