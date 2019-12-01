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

AC_PREREQ(2.59)

m4_define([gitversion], m4_esyscmd([sh -c "git describe --always | cut -d '=' -f 2 | head --bytes -1"]))
AC_INIT([MibS],[trunk],[https://github.com/coin-or/MibS/issues/new])

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

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(MibS)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(MibS)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# We filter out the -pedantic-errors flag for C and C++ compiler in
# case the user wants to use MPI, since otherwise the compilation
# fails
if test "${with_mpi_lib+set}" = set; then
  CFLAGS=`echo $CFLAGS | sed -e s/-pedantic-errors//g`
  CXXFLAGS=`echo $CXXFLAGS | sed -e s/-pedantic-errors//g`
fi

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_CHECK_PACKAGE(CoinDepend, [alps bcps blis coinutils osi cgl clp osi-clp], [MibSLib])
if test $coin_has_coindepend != yes ; then
  AC_MSG_ERROR([Required packages not available.])
fi

AC_COIN_CHECK_PACKAGE(SYMPHONY, [symphony osi-sym], [MibSLib])
if test $coin_has_symphony != yes ; then
  AC_MSG_ERROR([Required package SYMPHONY not available.])
fi

AC_COIN_CHECK_PACKAGE(Cbc, [cbc osi-cbc], [MibSLib])

#############################################################################
#                              check MPI libraries                          #
#############################################################################

AC_COIN_CHECK_USER_LIBRARY([Mpi], [MPI], [mpi.h], [MPI_Irecv], [], [AlpsLib])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available
AC_COIN_CHECK_PACKAGE(Cplex,  [osi-cplex],  [MibSLib])

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

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
