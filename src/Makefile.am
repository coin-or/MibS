#===========================================================================#
# This file is part of a Mixed Integer Bilevel Solver                       #
#                                                                           #
# Authors: Scott DeNegre, Lehigh University                                 #
#          Ted Ralphs, Lehigh University                                    #
#          Sahar Tahernajad, Lehigh University                              #
# Copyright (C) 2007-2019 Lehigh University, Scott DeNegre, Ted Ralphs,     #
#                         and Sahar Tahernajad                              #
# All Rights Reserved.                                                      #
#                                                                           #
# This software is licensed under the Eclipse Public License. Please see    #
# accompanying file for terms.                                              #
#===========================================================================#

AUTOMAKE_OPTIONS = foreign

########################################################################
#                           Common Flags                               #
########################################################################

LDADD = $(MIBSLIB_LIBS)
DEPENDENCIES = $(MIBSLIB_DEPENDENCIES)
AM_CPPFLAGS = $(MIBSLIB_CFLAGS)

########################################################################
#                               libMibS                                #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libMibs.la

# List all source files for this library, including headers
libMibs_la_SOURCES = \
	MibSBranchStrategyStrong.cpp \
	MibSObjectInt.cpp \
	MibSBilevel.cpp \
	MibSBranchStrategyStrong.hpp \
	MibSObjectInt.hpp \
	MibSBilevel.hpp \
	MibSConfig.hpp \
	MibSParams.cpp \
	MibSBranchObjectInt.cpp \
	MibSConstants.hpp \
	MibSParams.hpp \
	MibSBranchObjectInt.hpp \
	MibSCutGenerator.cpp \
	MibSSolTypes.hpp \
	MibSBranchStrategyMaxInf.cpp \
	MibSCutGenerator.hpp \
	MibSSolution.cpp \
	MibSBranchStrategyMaxInf.hpp \
	MibSHeuristic.cpp \
	MibSSolution.hpp \
	MibSBranchStrategyPseudo.cpp \
	MibSHeuristic.hpp \
	MibSTreeNode.cpp \
	MibSBranchStrategyPseudo.hpp \
	MibSTreeNode.hpp \
	MibSBranchStrategyRel.cpp \
	MibSModel.cpp \
	MibSBranchStrategyRel.hpp \
	MibSModel.hpp \
	MibSHelper.hpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libMibs_la_LIBADD = $(MIBSLIB_LIBS)
endif

# This is for libtool (on Windows)
libMibs_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                        MibS binary                                   #
########################################################################

bin_PROGRAMS = mibs
mibs_SOURCES = MibSMain.cpp
mibs_LDADD = libMibs.la $(MIBSLIB_LIBS)
mibs_DEPENDENCIES = libMibs.la $(MIBSLIB_DEPENDENCIES)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# AM_CPPFLAGS = $(MIBSLIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	MibSBilevel.hpp \
	MibSConfig.hpp \
	MibSParams.hpp \
	MibSBranchObjectInt.hpp \
	MibSConstants.hpp \
	MibSSolTypes.hpp \
	MibSBranchStrategyMaxInf.hpp \
	MibSCutGenerator.hpp \
	MibSSolution.hpp \
	MibSBranchStrategyPseudo.hpp \
	MibSHeuristic.hpp \
	MibSTreeNode.hpp \
	MibSBranchStrategyRel.hpp \
	MibSModel.hpp \
	MibSBranchStrategyStrong.hpp \
	MibSObjectInt.hpp \
	MibSHelper.hpp

########################################################################
#                     Installing manifest (MSVC++)                     #
########################################################################

if COIN_CXX_IS_CL
install-exec-hook:
	test -s mibs.exe.manifest && \
	mt -manifest mibs.exe.manifest -outputresource:mibs.exe;\
	cp mibs.exe $(bindir)/mibs.exe
endif

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

install-exec-local:
	$(install_sh_DATA) config_mibs.h $(DESTDIR)$(includecoindir)/MibsConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/MibSConfig.h
