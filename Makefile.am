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
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

#if ALWAYS_FALSE
#  SUBDIRS += test
#endif

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

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )


clean-local: clean-doxydoc
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest doxydoc

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = mibs.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/MibS

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@:$(DESTDIR)$(pkgconfiglibdir) \
	$(PKG_CONFIG) --libs mibs > $(addlibsdir)/mibs_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libMibs.lib @MIBSLIB_LIBS_INSTALLED@" > $(addlibsdir)/mibs_addlibs.txt
else
	echo -L@abs_lib_dir@ -lMibs @MIBSLIB_LIBS_INSTALLED@ > $(addlibsdir)/mibs_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/mibs_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
