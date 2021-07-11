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

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

pkgconfiglib_DATA = mibs.pc 

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = 

########################################################################
#                           Extra Targets                              #
########################################################################

#test: all
#	cd test; $(MAKE) test

#unitTest: test

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc uninstall-doxygen-docs

.PHONY: test unitTest doxydoc

