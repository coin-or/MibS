#===========================================================================#
# This file is part of a Mixed Integer Bilevel Solver                       #
# developed using the BiCePS Linear Integer Solver (BLIS).                  #
#                                                                           #
# Authors: Scott DeNegre, Lehigh University                                 #
#          Ted Ralphs, Lehigh University                                    #
#                                                                           #
# Copyright (C) 2007-2015 Lehigh University, Scott DeNegre, and Ted Ralphs. #
# All Rights Reserved.                                                      #
#                                                                           #
# This software is licensed under the Eclipse Public License. Please see    #
# accompanying file for terms.                                              #
#===========================================================================#

# CHANGEME: This should be the name of your executable
EXE = mibs

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement

OBJS =  MibSMain.o \
	MibSCutGenerator.o \
	MibSModel.o \
	MibSBilevel.o \
	MibSHeuristic.o \
	MibSSolution.o \
	MibSParams.o \
	MibSTreeNode.o \
	MibSBranchStrategyMaxInf.o \
	MibSBranchStrategyPseudo.o \
	MibSBranchStrategyStrong.o \
	MibSObjectInt.o \
	MibSBranchObjectInt.o

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I../include

# CHANGEME: Directory to the sources for the (example) problem definition
# files
SRCDIR = .

##########################################################################
#  CPLEX directories                                                     #
##########################################################################

USE_CPLEX = TRUE
CPLEXDIR      = /usr/local/cplex/
CPLEXINCDIR   = $(CPLEXDIR)/include
CPLEXLIBDIR   = $(CPLEXDIR)/lib/x86-64_linux/static_pic
CPLNFLAGS  = -lcplex -lm -lpthread -lreadline -lncurses

##########################################################################
#  Compiler and options                                                  #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -g -pipe -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion
ifeq ($(USE_CPLEX),TRUE)
CXXFLAGS += -DUSE_CPLEX -I$(CPLEXINCDIR)
endif

# additional C++ Compiler options for linking
CXXLINKFLAGS =

# Directory with COIN libraries
COININSTDIR = /home/ted/COIN/trunk/build-polyps-debug

###########################################################################
# Libraries necessary to link with COIN
###########################################################################

COIN_HAS_PKGCONFIG = TRUE

ifeq ($(COIN_HAS_PKGCONFIG),TRUE)
LIBS = `PKG_CONFIG_PATH=$(COININSTDIR)/lib/pkgconfig:$(COININSTDIR)/share/pkgconfig: pkg-config --libs blis osi-cbc osi-sym`
ifeq ($(USE_CPLEX),TRUE)
LIBS += `PKG_CONFIG_PATH=$(COININSTDIR)/lib/pkgconfig:$(COININSTDIR)/share/pkgconfig: pkg-config --libs osi-cplex`
endif
else
LIBS =  `cat $(COINDOCDIR)/Blis/blis_addlibs.txt` \
        `cat $(COINDOCDIR)/Cbc/cbc_addlibs.txt` \
        `cat $(COINDOCDIR)/SYMPHONY/sym_addlibs.txt`
ifeq ($(USE_CPLEX),TRUE)
LIBS += -L$(CPLEXLIBDIR) $(CPLNFLAGS)
endif
endif

###########################################################################
# Include directorie
###########################################################################

ifeq ($(COIN_HAS_PKGCONFIG),TRUE)
CXXFLAGS += `PKG_CONFIG_PATH=$(COININSTDIR)/lib/pkgconfig:$(COININSTDIR)/share/pkgconfig: pkg-config --cflags blis osi-cbc osi-sym` $(ADDINCFLAGS)
else
CXXFLAGS += $(ADDINCFLAGS)
endif

###########################################################################
# The following is necessary under cygwin, if native compilers are used
###########################################################################

CYGPATH_W = echo

###########################################################################
# Targets
###########################################################################

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
