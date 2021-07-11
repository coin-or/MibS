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

#############################################################################
#                       Names and other basic things                        #
#############################################################################

m4_define([gitversion], m4_esyscmd([sh -c "git describe --always | cut -d '=' -f 2 | head --bytes -1"]))
AC_INIT([MibS],[devel],[https://github.com/coin-or/MibS/issues/new],[coin-or-mibss],[https://github.com/coin-or/MibS])

AC_DEFINE([MIBS_GIT_VERSION], [gitversion], [MibS git version])

AC_COPYRIGHT([
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

])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/MibSBilevel.hpp)

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.
AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Check for C++11 support
AX_CXX_COMPILE_STDCXX(11,, [optional])

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                              COIN components                              #
#############################################################################

m4_foreach_w([myvar],[Alps Bcps Blis CoinUtils Osi Cgl Clp SYMPHONY],
  [ AC_COIN_CHK_PKG(myvar,[MibSLib])
    if test $coin_has_[]m4_tolower(myvar) != yes ; then
      AC_MSG_ERROR([Required package myvar is not available.])
    fi
  ])

AC_COIN_CHK_PKG(Sample,,[coindatasample],,dataonly)

AC_COIN_CHK_PKG(Cbc, [MibSLib])

AC_COIN_CHK_PKG(OsiClp,[MibSLib],[osi-clp])
AC_COIN_CHK_PKG(OsiCbc,[MibSLib],[osi-cbc])
AC_COIN_CHK_PKG(OsiCpx,[MibSLib],[osi-cplex])
AC_COIN_CHK_PKG(OsiSym,[MibSLib],[osi-sym])

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi Clp Cgl Alps Bcps Blis OsiCplex)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/mibs.par
                 mibs.pc])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_mibs.h])

AC_COIN_FINALIZE_FLAGS([MibSLib])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
