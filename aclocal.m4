# generated automatically by aclocal 1.16.3 -*- Autoconf -*-

# Copyright (C) 1996-2020 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
m4_if(m4_defn([AC_AUTOCONF_VERSION]), [2.71],,
[m4_warning([this file was generated for autoconf 2.71.
You have another version of autoconf.  It may work, but is not guaranteed to.
If you have problems, you may need to regenerate the build system entirely.
To do so, use the procedure documented by the package, typically 'autoreconf'.])])

# Copyright (C) 2002-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_AUTOMAKE_VERSION(VERSION)
# ----------------------------
# Automake X.Y traces this macro to ensure aclocal.m4 has been
# generated from the m4 files accompanying Automake X.Y.
# (This private macro should not be called outside this file.)
AC_DEFUN([AM_AUTOMAKE_VERSION],
[am__api_version='1.16'
dnl Some users find AM_AUTOMAKE_VERSION and mistake it for a way to
dnl require some minimum version.  Point them to the right macro.
m4_if([$1], [1.16.3], [],
      [AC_FATAL([Do not call $0, use AM_INIT_AUTOMAKE([$1]).])])dnl
])

# _AM_AUTOCONF_VERSION(VERSION)
# -----------------------------
# aclocal traces this macro to find the Autoconf version.
# This is a private macro too.  Using m4_define simplifies
# the logic in aclocal, which can simply ignore this definition.
m4_define([_AM_AUTOCONF_VERSION], [])

# AM_SET_CURRENT_AUTOMAKE_VERSION
# -------------------------------
# Call AM_AUTOMAKE_VERSION and AM_AUTOMAKE_VERSION so they can be traced.
# This function is AC_REQUIREd by AM_INIT_AUTOMAKE.
AC_DEFUN([AM_SET_CURRENT_AUTOMAKE_VERSION],
[AM_AUTOMAKE_VERSION([1.16.3])dnl
m4_ifndef([AC_AUTOCONF_VERSION],
  [m4_copy([m4_PACKAGE_VERSION], [AC_AUTOCONF_VERSION])])dnl
_AM_AUTOCONF_VERSION(m4_defn([AC_AUTOCONF_VERSION]))])

# Copyright (C) 2011-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_AR([ACT-IF-FAIL])
# -------------------------
# Try to determine the archiver interface, and trigger the ar-lib wrapper
# if it is needed.  If the detection of archiver interface fails, run
# ACT-IF-FAIL (default is to abort configure with a proper error message).
AC_DEFUN([AM_PROG_AR],
[AC_BEFORE([$0], [LT_INIT])dnl
AC_BEFORE([$0], [AC_PROG_LIBTOOL])dnl
AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([ar-lib])dnl
AC_CHECK_TOOLS([AR], [ar lib "link -lib"], [false])
: ${AR=ar}

AC_CACHE_CHECK([the archiver ($AR) interface], [am_cv_ar_interface],
  [AC_LANG_PUSH([C])
   am_cv_ar_interface=ar
   AC_COMPILE_IFELSE([AC_LANG_SOURCE([[int some_variable = 0;]])],
     [am_ar_try='$AR cru libconftest.a conftest.$ac_objext >&AS_MESSAGE_LOG_FD'
      AC_TRY_EVAL([am_ar_try])
      if test "$ac_status" -eq 0; then
        am_cv_ar_interface=ar
      else
        am_ar_try='$AR -NOLOGO -OUT:conftest.lib conftest.$ac_objext >&AS_MESSAGE_LOG_FD'
        AC_TRY_EVAL([am_ar_try])
        if test "$ac_status" -eq 0; then
          am_cv_ar_interface=lib
        else
          am_cv_ar_interface=unknown
        fi
      fi
      rm -f conftest.lib libconftest.a
     ])
   AC_LANG_POP([C])])

case $am_cv_ar_interface in
ar)
  ;;
lib)
  # Microsoft lib, so override with the ar-lib wrapper script.
  # FIXME: It is wrong to rewrite AR.
  # But if we don't then we get into trouble of one sort or another.
  # A longer-term fix would be to have automake use am__AR in this case,
  # and then we could set am__AR="$am_aux_dir/ar-lib \$(AR)" or something
  # similar.
  AR="$am_aux_dir/ar-lib $AR"
  ;;
unknown)
  m4_default([$1],
             [AC_MSG_ERROR([could not determine $AR interface])])
  ;;
esac
AC_SUBST([AR])dnl
])

# AM_AUX_DIR_EXPAND                                         -*- Autoconf -*-

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# For projects using AC_CONFIG_AUX_DIR([foo]), Autoconf sets
# $ac_aux_dir to '$srcdir/foo'.  In other projects, it is set to
# '$srcdir', '$srcdir/..', or '$srcdir/../..'.
#
# Of course, Automake must honor this variable whenever it calls a
# tool from the auxiliary directory.  The problem is that $srcdir (and
# therefore $ac_aux_dir as well) can be either absolute or relative,
# depending on how configure is run.  This is pretty annoying, since
# it makes $ac_aux_dir quite unusable in subdirectories: in the top
# source directory, any form will work fine, but in subdirectories a
# relative path needs to be adjusted first.
#
# $ac_aux_dir/missing
#    fails when called from a subdirectory if $ac_aux_dir is relative
# $top_srcdir/$ac_aux_dir/missing
#    fails if $ac_aux_dir is absolute,
#    fails when called from a subdirectory in a VPATH build with
#          a relative $ac_aux_dir
#
# The reason of the latter failure is that $top_srcdir and $ac_aux_dir
# are both prefixed by $srcdir.  In an in-source build this is usually
# harmless because $srcdir is '.', but things will broke when you
# start a VPATH build or use an absolute $srcdir.
#
# So we could use something similar to $top_srcdir/$ac_aux_dir/missing,
# iff we strip the leading $srcdir from $ac_aux_dir.  That would be:
#   am_aux_dir='\$(top_srcdir)/'`expr "$ac_aux_dir" : "$srcdir//*\(.*\)"`
# and then we would define $MISSING as
#   MISSING="\${SHELL} $am_aux_dir/missing"
# This will work as long as MISSING is not called from configure, because
# unfortunately $(top_srcdir) has no meaning in configure.
# However there are other variables, like CC, which are often used in
# configure, and could therefore not use this "fixed" $ac_aux_dir.
#
# Another solution, used here, is to always expand $ac_aux_dir to an
# absolute PATH.  The drawback is that using absolute paths prevent a
# configured tree to be moved without reconfiguration.

AC_DEFUN([AM_AUX_DIR_EXPAND],
[AC_REQUIRE([AC_CONFIG_AUX_DIR_DEFAULT])dnl
# Expand $ac_aux_dir to an absolute path.
am_aux_dir=`cd "$ac_aux_dir" && pwd`
])

# AM_CONDITIONAL                                            -*- Autoconf -*-

# Copyright (C) 1997-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_CONDITIONAL(NAME, SHELL-CONDITION)
# -------------------------------------
# Define a conditional.
AC_DEFUN([AM_CONDITIONAL],
[AC_PREREQ([2.52])dnl
 m4_if([$1], [TRUE],  [AC_FATAL([$0: invalid condition: $1])],
       [$1], [FALSE], [AC_FATAL([$0: invalid condition: $1])])dnl
AC_SUBST([$1_TRUE])dnl
AC_SUBST([$1_FALSE])dnl
_AM_SUBST_NOTMAKE([$1_TRUE])dnl
_AM_SUBST_NOTMAKE([$1_FALSE])dnl
m4_define([_AM_COND_VALUE_$1], [$2])dnl
if $2; then
  $1_TRUE=
  $1_FALSE='#'
else
  $1_TRUE='#'
  $1_FALSE=
fi
AC_CONFIG_COMMANDS_PRE(
[if test -z "${$1_TRUE}" && test -z "${$1_FALSE}"; then
  AC_MSG_ERROR([[conditional "$1" was never defined.
Usually this means the macro was only invoked conditionally.]])
fi])])

# Copyright (C) 1999-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# There are a few dirty hacks below to avoid letting 'AC_PROG_CC' be
# written in clear, in which case automake, when reading aclocal.m4,
# will think it sees a *use*, and therefore will trigger all it's
# C support machinery.  Also note that it means that autoscan, seeing
# CC etc. in the Makefile, will ask for an AC_PROG_CC use...


# _AM_DEPENDENCIES(NAME)
# ----------------------
# See how the compiler implements dependency checking.
# NAME is "CC", "CXX", "OBJC", "OBJCXX", "UPC", or "GJC".
# We try a few techniques and use that to set a single cache variable.
#
# We don't AC_REQUIRE the corresponding AC_PROG_CC since the latter was
# modified to invoke _AM_DEPENDENCIES(CC); we would have a circular
# dependency, and given that the user is not expected to run this macro,
# just rely on AC_PROG_CC.
AC_DEFUN([_AM_DEPENDENCIES],
[AC_REQUIRE([AM_SET_DEPDIR])dnl
AC_REQUIRE([AM_OUTPUT_DEPENDENCY_COMMANDS])dnl
AC_REQUIRE([AM_MAKE_INCLUDE])dnl
AC_REQUIRE([AM_DEP_TRACK])dnl

m4_if([$1], [CC],   [depcc="$CC"   am_compiler_list=],
      [$1], [CXX],  [depcc="$CXX"  am_compiler_list=],
      [$1], [OBJC], [depcc="$OBJC" am_compiler_list='gcc3 gcc'],
      [$1], [OBJCXX], [depcc="$OBJCXX" am_compiler_list='gcc3 gcc'],
      [$1], [UPC],  [depcc="$UPC"  am_compiler_list=],
      [$1], [GCJ],  [depcc="$GCJ"  am_compiler_list='gcc3 gcc'],
                    [depcc="$$1"   am_compiler_list=])

AC_CACHE_CHECK([dependency style of $depcc],
               [am_cv_$1_dependencies_compiler_type],
[if test -z "$AMDEP_TRUE" && test -f "$am_depcomp"; then
  # We make a subdir and do the tests there.  Otherwise we can end up
  # making bogus files that we don't know about and never remove.  For
  # instance it was reported that on HP-UX the gcc test will end up
  # making a dummy file named 'D' -- because '-MD' means "put the output
  # in D".
  rm -rf conftest.dir
  mkdir conftest.dir
  # Copy depcomp to subdir because otherwise we won't find it if we're
  # using a relative directory.
  cp "$am_depcomp" conftest.dir
  cd conftest.dir
  # We will build objects and dependencies in a subdirectory because
  # it helps to detect inapplicable dependency modes.  For instance
  # both Tru64's cc and ICC support -MD to output dependencies as a
  # side effect of compilation, but ICC will put the dependencies in
  # the current directory while Tru64 will put them in the object
  # directory.
  mkdir sub

  am_cv_$1_dependencies_compiler_type=none
  if test "$am_compiler_list" = ""; then
     am_compiler_list=`sed -n ['s/^#*\([a-zA-Z0-9]*\))$/\1/p'] < ./depcomp`
  fi
  am__universal=false
  m4_case([$1], [CC],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac],
    [CXX],
    [case " $depcc " in #(
     *\ -arch\ *\ -arch\ *) am__universal=true ;;
     esac])

  for depmode in $am_compiler_list; do
    # Setup a source with many dependencies, because some compilers
    # like to wrap large dependency lists on column 80 (with \), and
    # we should not choose a depcomp mode which is confused by this.
    #
    # We need to recreate these files for each test, as the compiler may
    # overwrite some of them when testing with obscure command lines.
    # This happens at least with the AIX C compiler.
    : > sub/conftest.c
    for i in 1 2 3 4 5 6; do
      echo '#include "conftst'$i'.h"' >> sub/conftest.c
      # Using ": > sub/conftst$i.h" creates only sub/conftst1.h with
      # Solaris 10 /bin/sh.
      echo '/* dummy */' > sub/conftst$i.h
    done
    echo "${am__include} ${am__quote}sub/conftest.Po${am__quote}" > confmf

    # We check with '-c' and '-o' for the sake of the "dashmstdout"
    # mode.  It turns out that the SunPro C++ compiler does not properly
    # handle '-M -o', and we need to detect this.  Also, some Intel
    # versions had trouble with output in subdirs.
    am__obj=sub/conftest.${OBJEXT-o}
    am__minus_obj="-o $am__obj"
    case $depmode in
    gcc)
      # This depmode causes a compiler race in universal mode.
      test "$am__universal" = false || continue
      ;;
    nosideeffect)
      # After this tag, mechanisms are not by side-effect, so they'll
      # only be used when explicitly requested.
      if test "x$enable_dependency_tracking" = xyes; then
	continue
      else
	break
      fi
      ;;
    msvc7 | msvc7msys | msvisualcpp | msvcmsys)
      # This compiler won't grok '-c -o', but also, the minuso test has
      # not run yet.  These depmodes are late enough in the game, and
      # so weak that their functioning should not be impacted.
      am__obj=conftest.${OBJEXT-o}
      am__minus_obj=
      ;;
    none) break ;;
    esac
    if depmode=$depmode \
       source=sub/conftest.c object=$am__obj \
       depfile=sub/conftest.Po tmpdepfile=sub/conftest.TPo \
       $SHELL ./depcomp $depcc -c $am__minus_obj sub/conftest.c \
         >/dev/null 2>conftest.err &&
       grep sub/conftst1.h sub/conftest.Po > /dev/null 2>&1 &&
       grep sub/conftst6.h sub/conftest.Po > /dev/null 2>&1 &&
       grep $am__obj sub/conftest.Po > /dev/null 2>&1 &&
       ${MAKE-make} -s -f confmf > /dev/null 2>&1; then
      # icc doesn't choke on unknown options, it will just issue warnings
      # or remarks (even with -Werror).  So we grep stderr for any message
      # that says an option was ignored or not supported.
      # When given -MP, icc 7.0 and 7.1 complain thusly:
      #   icc: Command line warning: ignoring option '-M'; no argument required
      # The diagnosis changed in icc 8.0:
      #   icc: Command line remark: option '-MP' not supported
      if (grep 'ignoring option' conftest.err ||
          grep 'not supported' conftest.err) >/dev/null 2>&1; then :; else
        am_cv_$1_dependencies_compiler_type=$depmode
        break
      fi
    fi
  done

  cd ..
  rm -rf conftest.dir
else
  am_cv_$1_dependencies_compiler_type=none
fi
])
AC_SUBST([$1DEPMODE], [depmode=$am_cv_$1_dependencies_compiler_type])
AM_CONDITIONAL([am__fastdep$1], [
  test "x$enable_dependency_tracking" != xno \
  && test "$am_cv_$1_dependencies_compiler_type" = gcc3])
])


# AM_SET_DEPDIR
# -------------
# Choose a directory name for dependency files.
# This macro is AC_REQUIREd in _AM_DEPENDENCIES.
AC_DEFUN([AM_SET_DEPDIR],
[AC_REQUIRE([AM_SET_LEADING_DOT])dnl
AC_SUBST([DEPDIR], ["${am__leading_dot}deps"])dnl
])


# AM_DEP_TRACK
# ------------
AC_DEFUN([AM_DEP_TRACK],
[AC_ARG_ENABLE([dependency-tracking], [dnl
AS_HELP_STRING(
  [--enable-dependency-tracking],
  [do not reject slow dependency extractors])
AS_HELP_STRING(
  [--disable-dependency-tracking],
  [speeds up one-time build])])
if test "x$enable_dependency_tracking" != xno; then
  am_depcomp="$ac_aux_dir/depcomp"
  AMDEPBACKSLASH='\'
  am__nodep='_no'
fi
AM_CONDITIONAL([AMDEP], [test "x$enable_dependency_tracking" != xno])
AC_SUBST([AMDEPBACKSLASH])dnl
_AM_SUBST_NOTMAKE([AMDEPBACKSLASH])dnl
AC_SUBST([am__nodep])dnl
_AM_SUBST_NOTMAKE([am__nodep])dnl
])

# Generate code to set up dependency tracking.              -*- Autoconf -*-

# Copyright (C) 1999-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_OUTPUT_DEPENDENCY_COMMANDS
# ------------------------------
AC_DEFUN([_AM_OUTPUT_DEPENDENCY_COMMANDS],
[{
  # Older Autoconf quotes --file arguments for eval, but not when files
  # are listed without --file.  Let's play safe and only enable the eval
  # if we detect the quoting.
  # TODO: see whether this extra hack can be removed once we start
  # requiring Autoconf 2.70 or later.
  AS_CASE([$CONFIG_FILES],
          [*\'*], [eval set x "$CONFIG_FILES"],
          [*], [set x $CONFIG_FILES])
  shift
  # Used to flag and report bootstrapping failures.
  am_rc=0
  for am_mf
  do
    # Strip MF so we end up with the name of the file.
    am_mf=`AS_ECHO(["$am_mf"]) | sed -e 's/:.*$//'`
    # Check whether this is an Automake generated Makefile which includes
    # dependency-tracking related rules and includes.
    # Grep'ing the whole file directly is not great: AIX grep has a line
    # limit of 2048, but all sed's we know have understand at least 4000.
    sed -n 's,^am--depfiles:.*,X,p' "$am_mf" | grep X >/dev/null 2>&1 \
      || continue
    am_dirpart=`AS_DIRNAME(["$am_mf"])`
    am_filepart=`AS_BASENAME(["$am_mf"])`
    AM_RUN_LOG([cd "$am_dirpart" \
      && sed -e '/# am--include-marker/d' "$am_filepart" \
        | $MAKE -f - am--depfiles]) || am_rc=$?
  done
  if test $am_rc -ne 0; then
    AC_MSG_FAILURE([Something went wrong bootstrapping makefile fragments
    for automatic dependency tracking.  If GNU make was not used, consider
    re-running the configure script with MAKE="gmake" (or whatever is
    necessary).  You can also try re-running configure with the
    '--disable-dependency-tracking' option to at least be able to build
    the package (albeit without support for automatic dependency tracking).])
  fi
  AS_UNSET([am_dirpart])
  AS_UNSET([am_filepart])
  AS_UNSET([am_mf])
  AS_UNSET([am_rc])
  rm -f conftest-deps.mk
}
])# _AM_OUTPUT_DEPENDENCY_COMMANDS


# AM_OUTPUT_DEPENDENCY_COMMANDS
# -----------------------------
# This macro should only be invoked once -- use via AC_REQUIRE.
#
# This code is only required when automatic dependency tracking is enabled.
# This creates each '.Po' and '.Plo' makefile fragment that we'll need in
# order to bootstrap the dependency handling code.
AC_DEFUN([AM_OUTPUT_DEPENDENCY_COMMANDS],
[AC_CONFIG_COMMANDS([depfiles],
     [test x"$AMDEP_TRUE" != x"" || _AM_OUTPUT_DEPENDENCY_COMMANDS],
     [AMDEP_TRUE="$AMDEP_TRUE" MAKE="${MAKE-make}"])])

# Do all the work for Automake.                             -*- Autoconf -*-

# Copyright (C) 1996-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This macro actually does too much.  Some checks are only needed if
# your package does certain things.  But this isn't really a big deal.

dnl Redefine AC_PROG_CC to automatically invoke _AM_PROG_CC_C_O.
m4_define([AC_PROG_CC],
m4_defn([AC_PROG_CC])
[_AM_PROG_CC_C_O
])

# AM_INIT_AUTOMAKE(PACKAGE, VERSION, [NO-DEFINE])
# AM_INIT_AUTOMAKE([OPTIONS])
# -----------------------------------------------
# The call with PACKAGE and VERSION arguments is the old style
# call (pre autoconf-2.50), which is being phased out.  PACKAGE
# and VERSION should now be passed to AC_INIT and removed from
# the call to AM_INIT_AUTOMAKE.
# We support both call styles for the transition.  After
# the next Automake release, Autoconf can make the AC_INIT
# arguments mandatory, and then we can depend on a new Autoconf
# release and drop the old call support.
AC_DEFUN([AM_INIT_AUTOMAKE],
[AC_PREREQ([2.65])dnl
dnl Autoconf wants to disallow AM_ names.  We explicitly allow
dnl the ones we care about.
m4_pattern_allow([^AM_[A-Z]+FLAGS$])dnl
AC_REQUIRE([AM_SET_CURRENT_AUTOMAKE_VERSION])dnl
AC_REQUIRE([AC_PROG_INSTALL])dnl
if test "`cd $srcdir && pwd`" != "`pwd`"; then
  # Use -I$(srcdir) only when $(srcdir) != ., so that make's output
  # is not polluted with repeated "-I."
  AC_SUBST([am__isrc], [' -I$(srcdir)'])_AM_SUBST_NOTMAKE([am__isrc])dnl
  # test to see if srcdir already configured
  if test -f $srcdir/config.status; then
    AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
  fi
fi

# test whether we have cygpath
if test -z "$CYGPATH_W"; then
  if (cygpath --version) >/dev/null 2>/dev/null; then
    CYGPATH_W='cygpath -w'
  else
    CYGPATH_W=echo
  fi
fi
AC_SUBST([CYGPATH_W])

# Define the identity of the package.
dnl Distinguish between old-style and new-style calls.
m4_ifval([$2],
[AC_DIAGNOSE([obsolete],
             [$0: two- and three-arguments forms are deprecated.])
m4_ifval([$3], [_AM_SET_OPTION([no-define])])dnl
 AC_SUBST([PACKAGE], [$1])dnl
 AC_SUBST([VERSION], [$2])],
[_AM_SET_OPTIONS([$1])dnl
dnl Diagnose old-style AC_INIT with new-style AM_AUTOMAKE_INIT.
m4_if(
  m4_ifdef([AC_PACKAGE_NAME], [ok]):m4_ifdef([AC_PACKAGE_VERSION], [ok]),
  [ok:ok],,
  [m4_fatal([AC_INIT should be called with package and version arguments])])dnl
 AC_SUBST([PACKAGE], ['AC_PACKAGE_TARNAME'])dnl
 AC_SUBST([VERSION], ['AC_PACKAGE_VERSION'])])dnl

_AM_IF_OPTION([no-define],,
[AC_DEFINE_UNQUOTED([PACKAGE], ["$PACKAGE"], [Name of package])
 AC_DEFINE_UNQUOTED([VERSION], ["$VERSION"], [Version number of package])])dnl

# Some tools Automake needs.
AC_REQUIRE([AM_SANITY_CHECK])dnl
AC_REQUIRE([AC_ARG_PROGRAM])dnl
AM_MISSING_PROG([ACLOCAL], [aclocal-${am__api_version}])
AM_MISSING_PROG([AUTOCONF], [autoconf])
AM_MISSING_PROG([AUTOMAKE], [automake-${am__api_version}])
AM_MISSING_PROG([AUTOHEADER], [autoheader])
AM_MISSING_PROG([MAKEINFO], [makeinfo])
AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
AC_REQUIRE([AM_PROG_INSTALL_STRIP])dnl
AC_REQUIRE([AC_PROG_MKDIR_P])dnl
# For better backward compatibility.  To be removed once Automake 1.9.x
# dies out for good.  For more background, see:
# <https://lists.gnu.org/archive/html/automake/2012-07/msg00001.html>
# <https://lists.gnu.org/archive/html/automake/2012-07/msg00014.html>
AC_SUBST([mkdir_p], ['$(MKDIR_P)'])
# We need awk for the "check" target (and possibly the TAP driver).  The
# system "awk" is bad on some platforms.
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([AC_PROG_MAKE_SET])dnl
AC_REQUIRE([AM_SET_LEADING_DOT])dnl
_AM_IF_OPTION([tar-ustar], [_AM_PROG_TAR([ustar])],
	      [_AM_IF_OPTION([tar-pax], [_AM_PROG_TAR([pax])],
			     [_AM_PROG_TAR([v7])])])
_AM_IF_OPTION([no-dependencies],,
[AC_PROVIDE_IFELSE([AC_PROG_CC],
		  [_AM_DEPENDENCIES([CC])],
		  [m4_define([AC_PROG_CC],
			     m4_defn([AC_PROG_CC])[_AM_DEPENDENCIES([CC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
		  [_AM_DEPENDENCIES([CXX])],
		  [m4_define([AC_PROG_CXX],
			     m4_defn([AC_PROG_CXX])[_AM_DEPENDENCIES([CXX])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJC],
		  [_AM_DEPENDENCIES([OBJC])],
		  [m4_define([AC_PROG_OBJC],
			     m4_defn([AC_PROG_OBJC])[_AM_DEPENDENCIES([OBJC])])])dnl
AC_PROVIDE_IFELSE([AC_PROG_OBJCXX],
		  [_AM_DEPENDENCIES([OBJCXX])],
		  [m4_define([AC_PROG_OBJCXX],
			     m4_defn([AC_PROG_OBJCXX])[_AM_DEPENDENCIES([OBJCXX])])])dnl
])
AC_REQUIRE([AM_SILENT_RULES])dnl
dnl The testsuite driver may need to know about EXEEXT, so add the
dnl 'am__EXEEXT' conditional if _AM_COMPILER_EXEEXT was seen.  This
dnl macro is hooked onto _AC_COMPILER_EXEEXT early, see below.
AC_CONFIG_COMMANDS_PRE(dnl
[m4_provide_if([_AM_COMPILER_EXEEXT],
  [AM_CONDITIONAL([am__EXEEXT], [test -n "$EXEEXT"])])])dnl

# POSIX will say in a future version that running "rm -f" with no argument
# is OK; and we want to be able to make that assumption in our Makefile
# recipes.  So use an aggressive probe to check that the usage we want is
# actually supported "in the wild" to an acceptable degree.
# See automake bug#10828.
# To make any issue more visible, cause the running configure to be aborted
# by default if the 'rm' program in use doesn't match our expectations; the
# user can still override this though.
if rm -f && rm -fr && rm -rf; then : OK; else
  cat >&2 <<'END'
Oops!

Your 'rm' program seems unable to run without file operands specified
on the command line, even when the '-f' option is present.  This is contrary
to the behaviour of most rm programs out there, and not conforming with
the upcoming POSIX standard: <http://austingroupbugs.net/view.php?id=542>

Please tell bug-automake@gnu.org about your system, including the value
of your $PATH and any error possibly output before this message.  This
can help us improve future automake versions.

END
  if test x"$ACCEPT_INFERIOR_RM_PROGRAM" = x"yes"; then
    echo 'Configuration will proceed anyway, since you have set the' >&2
    echo 'ACCEPT_INFERIOR_RM_PROGRAM variable to "yes"' >&2
    echo >&2
  else
    cat >&2 <<'END'
Aborting the configuration process, to ensure you take notice of the issue.

You can download and install GNU coreutils to get an 'rm' implementation
that behaves properly: <https://www.gnu.org/software/coreutils/>.

If you want to complete the configuration process using your problematic
'rm' anyway, export the environment variable ACCEPT_INFERIOR_RM_PROGRAM
to "yes", and re-run configure.

END
    AC_MSG_ERROR([Your 'rm' program is bad, sorry.])
  fi
fi
dnl The trailing newline in this macro's definition is deliberate, for
dnl backward compatibility and to allow trailing 'dnl'-style comments
dnl after the AM_INIT_AUTOMAKE invocation. See automake bug#16841.
])

dnl Hook into '_AC_COMPILER_EXEEXT' early to learn its expansion.  Do not
dnl add the conditional right here, as _AC_COMPILER_EXEEXT may be further
dnl mangled by Autoconf and run in a shell conditional statement.
m4_define([_AC_COMPILER_EXEEXT],
m4_defn([_AC_COMPILER_EXEEXT])[m4_provide([_AM_COMPILER_EXEEXT])])

# When config.status generates a header, we must update the stamp-h file.
# This file resides in the same directory as the config header
# that is generated.  The stamp files are numbered to have different names.

# Autoconf calls _AC_AM_CONFIG_HEADER_HOOK (when defined) in the
# loop where config.status creates the headers, so we can generate
# our stamp files there.
AC_DEFUN([_AC_AM_CONFIG_HEADER_HOOK],
[# Compute $1's index in $config_headers.
_am_arg=$1
_am_stamp_count=1
for _am_header in $config_headers :; do
  case $_am_header in
    $_am_arg | $_am_arg:* )
      break ;;
    * )
      _am_stamp_count=`expr $_am_stamp_count + 1` ;;
  esac
done
echo "timestamp for $_am_arg" >`AS_DIRNAME(["$_am_arg"])`/stamp-h[]$_am_stamp_count])

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_SH
# ------------------
# Define $install_sh.
AC_DEFUN([AM_PROG_INSTALL_SH],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
if test x"${install_sh+set}" != xset; then
  case $am_aux_dir in
  *\ * | *\	*)
    install_sh="\${SHELL} '$am_aux_dir/install-sh'" ;;
  *)
    install_sh="\${SHELL} $am_aux_dir/install-sh"
  esac
fi
AC_SUBST([install_sh])])

# Copyright (C) 2003-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# Check whether the underlying file-system supports filenames
# with a leading dot.  For instance MS-DOS doesn't.
AC_DEFUN([AM_SET_LEADING_DOT],
[rm -rf .tst 2>/dev/null
mkdir .tst 2>/dev/null
if test -d .tst; then
  am__leading_dot=.
else
  am__leading_dot=_
fi
rmdir .tst 2>/dev/null
AC_SUBST([am__leading_dot])])

# Add --enable-maintainer-mode option to configure.         -*- Autoconf -*-
# From Jim Meyering

# Copyright (C) 1996-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MAINTAINER_MODE([DEFAULT-MODE])
# ----------------------------------
# Control maintainer-specific portions of Makefiles.
# Default is to disable them, unless 'enable' is passed literally.
# For symmetry, 'disable' may be passed as well.  Anyway, the user
# can override the default with the --enable/--disable switch.
AC_DEFUN([AM_MAINTAINER_MODE],
[m4_case(m4_default([$1], [disable]),
       [enable], [m4_define([am_maintainer_other], [disable])],
       [disable], [m4_define([am_maintainer_other], [enable])],
       [m4_define([am_maintainer_other], [enable])
        m4_warn([syntax], [unexpected argument to AM@&t@_MAINTAINER_MODE: $1])])
AC_MSG_CHECKING([whether to enable maintainer-specific portions of Makefiles])
  dnl maintainer-mode's default is 'disable' unless 'enable' is passed
  AC_ARG_ENABLE([maintainer-mode],
    [AS_HELP_STRING([--]am_maintainer_other[-maintainer-mode],
      am_maintainer_other[ make rules and dependencies not useful
      (and sometimes confusing) to the casual installer])],
    [USE_MAINTAINER_MODE=$enableval],
    [USE_MAINTAINER_MODE=]m4_if(am_maintainer_other, [enable], [no], [yes]))
  AC_MSG_RESULT([$USE_MAINTAINER_MODE])
  AM_CONDITIONAL([MAINTAINER_MODE], [test $USE_MAINTAINER_MODE = yes])
  MAINT=$MAINTAINER_MODE_TRUE
  AC_SUBST([MAINT])dnl
]
)

# Check to see how 'make' treats includes.	            -*- Autoconf -*-

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MAKE_INCLUDE()
# -----------------
# Check whether make has an 'include' directive that can support all
# the idioms we need for our automatic dependency tracking code.
AC_DEFUN([AM_MAKE_INCLUDE],
[AC_MSG_CHECKING([whether ${MAKE-make} supports the include directive])
cat > confinc.mk << 'END'
am__doit:
	@echo this is the am__doit target >confinc.out
.PHONY: am__doit
END
am__include="#"
am__quote=
# BSD make does it like this.
echo '.include "confinc.mk" # ignored' > confmf.BSD
# Other make implementations (GNU, Solaris 10, AIX) do it like this.
echo 'include confinc.mk # ignored' > confmf.GNU
_am_result=no
for s in GNU BSD; do
  AM_RUN_LOG([${MAKE-make} -f confmf.$s && cat confinc.out])
  AS_CASE([$?:`cat confinc.out 2>/dev/null`],
      ['0:this is the am__doit target'],
      [AS_CASE([$s],
          [BSD], [am__include='.include' am__quote='"'],
          [am__include='include' am__quote=''])])
  if test "$am__include" != "#"; then
    _am_result="yes ($s style)"
    break
  fi
done
rm -f confinc.* confmf.*
AC_MSG_RESULT([${_am_result}])
AC_SUBST([am__include])])
AC_SUBST([am__quote])])

# Fake the existence of programs that GNU maintainers use.  -*- Autoconf -*-

# Copyright (C) 1997-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_MISSING_PROG(NAME, PROGRAM)
# ------------------------------
AC_DEFUN([AM_MISSING_PROG],
[AC_REQUIRE([AM_MISSING_HAS_RUN])
$1=${$1-"${am_missing_run}$2"}
AC_SUBST($1)])

# AM_MISSING_HAS_RUN
# ------------------
# Define MISSING if not defined so far and test if it is modern enough.
# If it is, set am_missing_run to use it, otherwise, to nothing.
AC_DEFUN([AM_MISSING_HAS_RUN],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([missing])dnl
if test x"${MISSING+set}" != xset; then
  MISSING="\${SHELL} '$am_aux_dir/missing'"
fi
# Use eval to expand $SHELL
if eval "$MISSING --is-lightweight"; then
  am_missing_run="$MISSING "
else
  am_missing_run=
  AC_MSG_WARN(['missing' script is too old or missing])
fi
])

# Helper functions for option handling.                     -*- Autoconf -*-

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_MANGLE_OPTION(NAME)
# -----------------------
AC_DEFUN([_AM_MANGLE_OPTION],
[[_AM_OPTION_]m4_bpatsubst($1, [[^a-zA-Z0-9_]], [_])])

# _AM_SET_OPTION(NAME)
# --------------------
# Set option NAME.  Presently that only means defining a flag for this option.
AC_DEFUN([_AM_SET_OPTION],
[m4_define(_AM_MANGLE_OPTION([$1]), [1])])

# _AM_SET_OPTIONS(OPTIONS)
# ------------------------
# OPTIONS is a space-separated list of Automake options.
AC_DEFUN([_AM_SET_OPTIONS],
[m4_foreach_w([_AM_Option], [$1], [_AM_SET_OPTION(_AM_Option)])])

# _AM_IF_OPTION(OPTION, IF-SET, [IF-NOT-SET])
# -------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
AC_DEFUN([_AM_IF_OPTION],
[m4_ifset(_AM_MANGLE_OPTION([$1]), [$2], [$3])])

# Copyright (C) 1999-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_CC_C_O
# ---------------
# Like AC_PROG_CC_C_O, but changed for automake.  We rewrite AC_PROG_CC
# to automatically call this.
AC_DEFUN([_AM_PROG_CC_C_O],
[AC_REQUIRE([AM_AUX_DIR_EXPAND])dnl
AC_REQUIRE_AUX_FILE([compile])dnl
AC_LANG_PUSH([C])dnl
AC_CACHE_CHECK(
  [whether $CC understands -c and -o together],
  [am_cv_prog_cc_c_o],
  [AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])
  # Make sure it works both with $CC and with simple cc.
  # Following AC_PROG_CC_C_O, we do the test twice because some
  # compilers refuse to overwrite an existing .o file with -o,
  # though they will create one.
  am_cv_prog_cc_c_o=yes
  for am_i in 1 2; do
    if AM_RUN_LOG([$CC -c conftest.$ac_ext -o conftest2.$ac_objext]) \
         && test -f conftest2.$ac_objext; then
      : OK
    else
      am_cv_prog_cc_c_o=no
      break
    fi
  done
  rm -f core conftest*
  unset am_i])
if test "$am_cv_prog_cc_c_o" != yes; then
   # Losing compiler, so override with the script.
   # FIXME: It is wrong to rewrite CC.
   # But if we don't then we get into trouble of one sort or another.
   # A longer-term fix would be to have automake use am__CC in this case,
   # and then we could set am__CC="\$(top_srcdir)/compile \$(CC)"
   CC="$am_aux_dir/compile $CC"
fi
AC_LANG_POP([C])])

# For backward compatibility.
AC_DEFUN_ONCE([AM_PROG_CC_C_O], [AC_REQUIRE([AC_PROG_CC])])

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_RUN_LOG(COMMAND)
# -------------------
# Run COMMAND, save the exit status in ac_status, and log it.
# (This has been adapted from Autoconf's _AC_RUN_LOG macro.)
AC_DEFUN([AM_RUN_LOG],
[{ echo "$as_me:$LINENO: $1" >&AS_MESSAGE_LOG_FD
   ($1) >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD
   ac_status=$?
   echo "$as_me:$LINENO: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   (exit $ac_status); }])

# Check to make sure that the build environment is sane.    -*- Autoconf -*-

# Copyright (C) 1996-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SANITY_CHECK
# ---------------
AC_DEFUN([AM_SANITY_CHECK],
[AC_MSG_CHECKING([whether build environment is sane])
# Reject unsafe characters in $srcdir or the absolute working directory
# name.  Accept space and tab only in the latter.
am_lf='
'
case `pwd` in
  *[[\\\"\#\$\&\'\`$am_lf]]*)
    AC_MSG_ERROR([unsafe absolute working directory name]);;
esac
case $srcdir in
  *[[\\\"\#\$\&\'\`$am_lf\ \	]]*)
    AC_MSG_ERROR([unsafe srcdir value: '$srcdir']);;
esac

# Do 'set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   am_has_slept=no
   for am_try in 1 2; do
     echo "timestamp, slept: $am_has_slept" > conftest.file
     set X `ls -Lt "$srcdir/configure" conftest.file 2> /dev/null`
     if test "$[*]" = "X"; then
	# -L didn't work.
	set X `ls -t "$srcdir/configure" conftest.file`
     fi
     if test "$[*]" != "X $srcdir/configure conftest.file" \
	&& test "$[*]" != "X conftest.file $srcdir/configure"; then

	# If neither matched, then we have a broken ls.  This can happen
	# if, for instance, CONFIG_SHELL is bash and it inherits a
	# broken ls alias from the environment.  This has actually
	# happened.  Such a system could not be considered "sane".
	AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
  alias in your environment])
     fi
     if test "$[2]" = conftest.file || test $am_try -eq 2; then
       break
     fi
     # Just in case.
     sleep 1
     am_has_slept=yes
   done
   test "$[2]" = conftest.file
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
AC_MSG_RESULT([yes])
# If we didn't sleep, we still need to ensure time stamps of config.status and
# generated files are strictly newer.
am_sleep_pid=
if grep 'slept: no' conftest.file >/dev/null 2>&1; then
  ( sleep 1 ) &
  am_sleep_pid=$!
fi
AC_CONFIG_COMMANDS_PRE(
  [AC_MSG_CHECKING([that generated files are newer than configure])
   if test -n "$am_sleep_pid"; then
     # Hide warnings about reused PIDs.
     wait $am_sleep_pid 2>/dev/null
   fi
   AC_MSG_RESULT([done])])
rm -f conftest.file
])

# Copyright (C) 2009-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_SILENT_RULES([DEFAULT])
# --------------------------
# Enable less verbose build rules; with the default set to DEFAULT
# ("yes" being less verbose, "no" or empty being verbose).
AC_DEFUN([AM_SILENT_RULES],
[AC_ARG_ENABLE([silent-rules], [dnl
AS_HELP_STRING(
  [--enable-silent-rules],
  [less verbose build output (undo: "make V=1")])
AS_HELP_STRING(
  [--disable-silent-rules],
  [verbose build output (undo: "make V=0")])dnl
])
case $enable_silent_rules in @%:@ (((
  yes) AM_DEFAULT_VERBOSITY=0;;
   no) AM_DEFAULT_VERBOSITY=1;;
    *) AM_DEFAULT_VERBOSITY=m4_if([$1], [yes], [0], [1]);;
esac
dnl
dnl A few 'make' implementations (e.g., NonStop OS and NextStep)
dnl do not support nested variable expansions.
dnl See automake bug#9928 and bug#10237.
am_make=${MAKE-make}
AC_CACHE_CHECK([whether $am_make supports nested variables],
   [am_cv_make_support_nested_variables],
   [if AS_ECHO([['TRUE=$(BAR$(V))
BAR0=false
BAR1=true
V=1
am__doit:
	@$(TRUE)
.PHONY: am__doit']]) | $am_make -f - >/dev/null 2>&1; then
  am_cv_make_support_nested_variables=yes
else
  am_cv_make_support_nested_variables=no
fi])
if test $am_cv_make_support_nested_variables = yes; then
  dnl Using '$V' instead of '$(V)' breaks IRIX make.
  AM_V='$(V)'
  AM_DEFAULT_V='$(AM_DEFAULT_VERBOSITY)'
else
  AM_V=$AM_DEFAULT_VERBOSITY
  AM_DEFAULT_V=$AM_DEFAULT_VERBOSITY
fi
AC_SUBST([AM_V])dnl
AM_SUBST_NOTMAKE([AM_V])dnl
AC_SUBST([AM_DEFAULT_V])dnl
AM_SUBST_NOTMAKE([AM_DEFAULT_V])dnl
AC_SUBST([AM_DEFAULT_VERBOSITY])dnl
AM_BACKSLASH='\'
AC_SUBST([AM_BACKSLASH])dnl
_AM_SUBST_NOTMAKE([AM_BACKSLASH])dnl
])

# Copyright (C) 2001-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# AM_PROG_INSTALL_STRIP
# ---------------------
# One issue with vendor 'install' (even GNU) is that you can't
# specify the program used to strip binaries.  This is especially
# annoying in cross-compiling environments, where the build's strip
# is unlikely to handle the host's binaries.
# Fortunately install-sh will honor a STRIPPROG variable, so we
# always use install-sh in "make install-strip", and initialize
# STRIPPROG with the value of the STRIP variable (set by the user).
AC_DEFUN([AM_PROG_INSTALL_STRIP],
[AC_REQUIRE([AM_PROG_INSTALL_SH])dnl
# Installed binaries are usually stripped using 'strip' when the user
# run "make install-strip".  However 'strip' might not be the right
# tool to use in cross-compilation environments, therefore Automake
# will honor the 'STRIP' environment variable to overrule this program.
dnl Don't test for $cross_compiling = yes, because it might be 'maybe'.
if test "$cross_compiling" != no; then
  AC_CHECK_TOOL([STRIP], [strip], :)
fi
INSTALL_STRIP_PROGRAM="\$(install_sh) -c -s"
AC_SUBST([INSTALL_STRIP_PROGRAM])])

# Copyright (C) 2006-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_SUBST_NOTMAKE(VARIABLE)
# ---------------------------
# Prevent Automake from outputting VARIABLE = @VARIABLE@ in Makefile.in.
# This macro is traced by Automake.
AC_DEFUN([_AM_SUBST_NOTMAKE])

# AM_SUBST_NOTMAKE(VARIABLE)
# --------------------------
# Public sister of _AM_SUBST_NOTMAKE.
AC_DEFUN([AM_SUBST_NOTMAKE], [_AM_SUBST_NOTMAKE($@)])

# Check how to create a tarball.                            -*- Autoconf -*-

# Copyright (C) 2004-2020 Free Software Foundation, Inc.
#
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# _AM_PROG_TAR(FORMAT)
# --------------------
# Check how to create a tarball in format FORMAT.
# FORMAT should be one of 'v7', 'ustar', or 'pax'.
#
# Substitute a variable $(am__tar) that is a command
# writing to stdout a FORMAT-tarball containing the directory
# $tardir.
#     tardir=directory && $(am__tar) > result.tar
#
# Substitute a variable $(am__untar) that extract such
# a tarball read from stdin.
#     $(am__untar) < result.tar
#
AC_DEFUN([_AM_PROG_TAR],
[# Always define AMTAR for backward compatibility.  Yes, it's still used
# in the wild :-(  We should find a proper way to deprecate it ...
AC_SUBST([AMTAR], ['$${TAR-tar}'])

# We'll loop over all known methods to create a tar archive until one works.
_am_tools='gnutar m4_if([$1], [ustar], [plaintar]) pax cpio none'

m4_if([$1], [v7],
  [am__tar='$${TAR-tar} chof - "$$tardir"' am__untar='$${TAR-tar} xf -'],

  [m4_case([$1],
    [ustar],
     [# The POSIX 1988 'ustar' format is defined with fixed-size fields.
      # There is notably a 21 bits limit for the UID and the GID.  In fact,
      # the 'pax' utility can hang on bigger UID/GID (see automake bug#8343
      # and bug#13588).
      am_max_uid=2097151 # 2^21 - 1
      am_max_gid=$am_max_uid
      # The $UID and $GID variables are not portable, so we need to resort
      # to the POSIX-mandated id(1) utility.  Errors in the 'id' calls
      # below are definitely unexpected, so allow the users to see them
      # (that is, avoid stderr redirection).
      am_uid=`id -u || echo unknown`
      am_gid=`id -g || echo unknown`
      AC_MSG_CHECKING([whether UID '$am_uid' is supported by ustar format])
      if test $am_uid -le $am_max_uid; then
         AC_MSG_RESULT([yes])
      else
         AC_MSG_RESULT([no])
         _am_tools=none
      fi
      AC_MSG_CHECKING([whether GID '$am_gid' is supported by ustar format])
      if test $am_gid -le $am_max_gid; then
         AC_MSG_RESULT([yes])
      else
        AC_MSG_RESULT([no])
        _am_tools=none
      fi],

  [pax],
    [],

  [m4_fatal([Unknown tar format])])

  AC_MSG_CHECKING([how to create a $1 tar archive])

  # Go ahead even if we have the value already cached.  We do so because we
  # need to set the values for the 'am__tar' and 'am__untar' variables.
  _am_tools=${am_cv_prog_tar_$1-$_am_tools}

  for _am_tool in $_am_tools; do
    case $_am_tool in
    gnutar)
      for _am_tar in tar gnutar gtar; do
        AM_RUN_LOG([$_am_tar --version]) && break
      done
      am__tar="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$$tardir"'
      am__tar_="$_am_tar --format=m4_if([$1], [pax], [posix], [$1]) -chf - "'"$tardir"'
      am__untar="$_am_tar -xf -"
      ;;
    plaintar)
      # Must skip GNU tar: if it does not support --format= it doesn't create
      # ustar tarball either.
      (tar --version) >/dev/null 2>&1 && continue
      am__tar='tar chf - "$$tardir"'
      am__tar_='tar chf - "$tardir"'
      am__untar='tar xf -'
      ;;
    pax)
      am__tar='pax -L -x $1 -w "$$tardir"'
      am__tar_='pax -L -x $1 -w "$tardir"'
      am__untar='pax -r'
      ;;
    cpio)
      am__tar='find "$$tardir" -print | cpio -o -H $1 -L'
      am__tar_='find "$tardir" -print | cpio -o -H $1 -L'
      am__untar='cpio -i -H $1 -d'
      ;;
    none)
      am__tar=false
      am__tar_=false
      am__untar=false
      ;;
    esac

    # If the value was cached, stop now.  We just wanted to have am__tar
    # and am__untar set.
    test -n "${am_cv_prog_tar_$1}" && break

    # tar/untar a dummy directory, and stop if the command works.
    rm -rf conftest.dir
    mkdir conftest.dir
    echo GrepMe > conftest.dir/file
    AM_RUN_LOG([tardir=conftest.dir && eval $am__tar_ >conftest.tar])
    rm -rf conftest.dir
    if test -s conftest.tar; then
      AM_RUN_LOG([$am__untar <conftest.tar])
      AM_RUN_LOG([cat conftest.dir/file])
      grep GrepMe conftest.dir/file >/dev/null 2>&1 && break
    fi
  done
  rm -rf conftest.dir

  AC_CACHE_VAL([am_cv_prog_tar_$1], [am_cv_prog_tar_$1=$_am_tool])
  AC_MSG_RESULT([$am_cv_prog_tar_$1])])

AC_SUBST([am__tar])
AC_SUBST([am__untar])
]) # _AM_PROG_TAR

# Copyright (C) 2020 COIN-OR Foundation
# All Rights Reserved.
# This file is distributed under the Eclipse Public License 2.0.
# See LICENSE for details.
#
# This file defines the common autoconf macros for COIN-OR

# Check requirements
AC_PREREQ([2.71])

#                           COIN_CHECK_VPATH                              #

# This macro sets the variable coin_vpath_config to true if this is a VPATH
# configuration, otherwise it sets it to false. `VPATH' just means we're
# following best practices and not building in the source directory.

AC_DEFUN([AC_COIN_CHECK_VPATH],
[
  AC_MSG_CHECKING(whether this is a VPATH configuration)
  if test `cd $srcdir ; pwd` != `pwd`; then
    coin_vpath_config=yes
  else
    coin_vpath_config=no
  fi
  AC_MSG_RESULT($coin_vpath_config)
])


#                            COIN_VPATH_LINK                              #

# This macro ensures that the given files are available in a VPATH
# configuration, using the same name and relative path as in the source
# tree. It expects a white-space separated list of files.
# This macro is a small wrapper around AC_CONFIG_LINKS.

AC_DEFUN([AC_COIN_VPATH_LINK],
[
  m4_foreach_w(linkvar,[$1],[AC_CONFIG_LINKS(linkvar:linkvar)])
])


#                          COIN_PROJECTVERSION                            #

# This macro is used by COIN_INITIALIZE and sets up variables related to
# versioning numbers of the project.
#   COIN_PROJECTVERSION([libtool_version_string])
#
# If libtool_version_string is given, then defines AC_COIN_LIBTOOL_VERSIONINFO,
# which will be picked up by AC_COIN_PROG_LIBTOOL to set libtools -version-info flag.
#
# Further, breaks up AC_PACKAGE_VERSION into AC_PACKAGE_VERSION_MAJOR,
# AC_PACKAGE_VERSION_MINOR, AC_PACKAGE_VERSION_RELEASE, assuming it has
# the form major.minor.release, but use 9999 for missing values.
# AC_PACKAGE_VERSION is defined by AC_INIT

AC_DEFUN([AC_COIN_PROJECTVERSION],
[
  m4_define(AC_PACKAGE_VERSION_MAJOR, m4_bregexp(AC_PACKAGE_VERSION, [^\([0-9]*\).*], [\1]))
  m4_define(AC_PACKAGE_VERSION_MINOR, m4_bregexp(AC_PACKAGE_VERSION, [^[0-9]*\.\([0-9]*\).*], [\1]))
  m4_define(AC_PACKAGE_VERSION_RELEASE, m4_bregexp(AC_PACKAGE_VERSION, [^[0-9]*\.[0-9]*\.\([0-9]*\).*], [\1]))

  AC_DEFINE_UNQUOTED(m4_toupper(AC_PACKAGE_NAME)_VERSION,
    "AC_PACKAGE_VERSION",
    [Version number of project])
  AC_DEFINE_UNQUOTED(m4_toupper(AC_PACKAGE_NAME)_VERSION_MAJOR,
    m4_ifnblank(AC_PACKAGE_VERSION_MAJOR,AC_PACKAGE_VERSION_MAJOR,9999),
    [Major version number of project.])
  AC_DEFINE_UNQUOTED(m4_toupper(AC_PACKAGE_NAME)_VERSION_MINOR,
    m4_ifnblank(AC_PACKAGE_VERSION_MINOR,AC_PACKAGE_VERSION_MINOR,9999),
    [Minor version number of project.])
  AC_DEFINE_UNQUOTED(m4_toupper(AC_PACKAGE_NAME)_VERSION_RELEASE,
    m4_ifnblank(AC_PACKAGE_VERSION_RELEASE,AC_PACKAGE_VERSION_RELEASE,9999),
    [Release version number of project.])

  m4_define(AC_COIN_LIBTOOL_VERSIONINFO,$1)
])


#                          COIN_ENABLE_MSVC                               #

# This macro is invoked by any PROG_compiler macro to establish the
# --enable-msvc option.
# If unset but we have a Windows environment, look for some known C compilers
# and set enable_msvc if (i)cl is picked up (or has been set via CC by user)

AC_DEFUN([AC_COIN_ENABLE_MSVC],
[
  AC_CANONICAL_BUILD

  AC_ARG_ENABLE([msvc],
    [AS_HELP_STRING([--enable-msvc],[look for and allow only Intel/Microsoft compilers on MinGW/MSys/Cygwin])],
    [enable_msvc=$enableval],
    [enable_msvc=no
     case $build in
       *-mingw* | *-cygwin* | *-msys* )
         AC_CHECK_PROGS(CC, [gcc clang icl cl])
         case "$CC" in *cl ) enable_msvc=yes ;; esac
       ;;
     esac])
])


#                        COIN_COMPFLAGS_DEFAULTS                          #

# Overwrite default compiler flags, since the autoconf defaults don't work
# well for cl/icl and we want to take --enable-debug and --enable-msvc into
# account:
# - debugging enabled: enable debug symbols (-g/-Z7)
# - debugging disabled: disable debug code (-DNDEBUG); enable (more) optimization (-O2)
# - enable exceptions for (i)cl

AC_DEFUN([AC_COIN_COMPFLAGS_DEFAULTS],
[
dnl We want --enable-msvc setup and checked
  AC_REQUIRE([AC_COIN_ENABLE_MSVC])

  AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],[build debugging symbols and turn off compiler optimization])],
    [enable_debug=$enableval],
    [enable_debug=no])

dnl This macro should run before the compiler checks (does not seem to be sufficient for CFLAGS)
  AC_BEFORE([$0],[AC_PROG_CXX])
  AC_BEFORE([$0],[AC_PROG_CC])
  AC_BEFORE([$0],[AC_PROG_F77])
  AC_BEFORE([$0],[AC_PROG_FC])

  if test "$enable_debug" = yes ; then
    if test "$enable_msvc" = yes ; then
      : ${FFLAGS:="-nologo -fpp -Z7 -MDd $ADD_FFLAGS"}
      : ${FCFLAGS:="-nologo -fpp -Z7 -MDd $ADD_FCFLAGS"}
      : ${CFLAGS:="-nologo -Z7 -MDd $ADD_CFLAGS"}
      : ${CXXFLAGS:="-nologo -EHs -Z7 -MDd $ADD_CXXFLAGS"}
    else
      : ${FFLAGS:="-g $ADD_FFLAGS"}
      : ${FCFLAGS:="-g $ADD_FCFLAGS"}
      : ${CFLAGS:="-g $ADD_CFLAGS"}
      : ${CXXFLAGS:="-g $ADD_CXXFLAGS"}
    fi
  else
    if test "$enable_msvc" = yes ; then
      : ${FFLAGS:="-nologo -fpp -O2 -MD $ADD_FFLAGS"}
      : ${FCFLAGS:="-nologo -fpp -O2 -MD $ADD_FCFLAGS"}
      : ${CFLAGS:="-nologo -DNDEBUG -O2 -MD $ADD_CFLAGS"}
      : ${CXXFLAGS:="-nologo -EHs -DNDEBUG -O2 -MD $ADD_CXXFLAGS"}
    else
      : ${FFLAGS:="-O2 $ADD_FFLAGS"}
      : ${FCFLAGS:="-O2 $ADD_FCFLAGS"}
      : ${CFLAGS:="-O2 -DNDEBUG $ADD_CFLAGS"}
      : ${CXXFLAGS:="-O2 -DNDEBUG $ADD_CXXFLAGS"}
    fi
  fi
])


#                            COIN_DEBUGLEVEL                              #

# This macro makes the switches --with-prjct-verbosity and
# --with-prjct-checklevel available, which define the preprocessor macros
# COIN_PRJCT_VERBOSITY and COIN_PRJCT_CHECKLEVEL to the specified value
# (default is 0).

AC_DEFUN([AC_COIN_DEBUGLEVEL],
[
  AC_ARG_WITH(m4_tolower(AC_PACKAGE_NAME)-verbosity,
    AS_HELP_STRING([--with-m4_tolower(AC_PACKAGE_NAME)-verbosity],[specify the debug verbosity level]),
    [if test "$withval" = yes; then withval=1 ; fi
     coin_verbosity=$withval],
    [coin_verbosity=0])
  AC_DEFINE_UNQUOTED(m4_toupper(AC_PACKAGE_NAME)_VERBOSITY,
                     $coin_verbosity,
                     [Define to the debug verbosity level (0 is no output)])

  AC_ARG_WITH(m4_tolower(AC_PACKAGE_NAME)-checklevel,
    AS_HELP_STRING([--with-m4_tolower(AC_PACKAGE_NAME)-checklevel],[specify the sanity check level]),
    [if test "$withval" = yes; then withval=1 ; fi
     coin_checklevel=$withval],
    [coin_checklevel=0])
  AC_DEFINE_UNQUOTED(m4_toupper(AC_PACKAGE_NAME)_CHECKLEVEL,
                     $coin_checklevel,
                     [Define to the debug sanity check level (0 is no test)])
])


#                            COIN_INITIALIZE                               #

# AC_COIN_INITIALIZE(version)

# This macro does everything that is required in the early part of the
# configure script, such as defining a few variables.
#   version (optional): the libtool library version
#                       project version number is used if not available
# - enforces the required autoconf version
# - sets the project's version numbers
# - changes the default compiler flags
# - get and patch the build and host types
# - Make silent build rules the default (https://www.gnu.org/software/automake/
#   manual/html_node/Automake-Silent-Rules.html)
# - disable automake maintainer mode
# - Initialize automake
#   - do not be as strict as for GNU projects
#   - don't AC_DEFINE PACKAGE or VERSION (but there're still defined as shell
#     variables in configure, and as make variables).
#   - disable dist target
#   - place objects from sources in subdirs into corresponding subdirs
#   - enable all automake warnings
# - Figure out the path where libraries are installed.
# - Add a configure flag to indicate whether .pc files should be made relocatable.
#   This is off by default, as it creates libtool warnings and doesn't account for
#   cases where --libdir is used.

AC_DEFUN([AC_COIN_INITIALIZE],
[
  AC_PREREQ([2.71])

  AC_COIN_PROJECTVERSION($1)

  AC_REQUIRE([AC_COIN_COMPFLAGS_DEFAULTS])

dnl this needs to used before AC_CANONICAL_BUILD
  AC_CANONICAL_BUILD
  AC_CANONICAL_HOST

  # libtool has some magic for host_os and build_os being mingw, but doesn't know about msys
  if test $host_os = msys ; then
    host_os=mingw
    host=`echo $host | sed -e 's/msys/mingw/'`
  fi

  if test $build_os = msys ; then
    build_os=mingw
    build=`echo $build | sed -e 's/msys/mingw/'`
  fi

dnl this needs to be used before AM_INIT_AUTOMAKE
  AM_SILENT_RULES([yes])

dnl without this, automake adds targets that may try to rerun the autotools
dnl on the users side, which asks for trouble
  AM_MAINTAINER_MODE([disable])

  AM_INIT_AUTOMAKE([foreign no-define no-dist subdir-objects -Wall])

  # Figure out the path where libraries are installed.
  # Unless the user specifies --prefix, prefix is set to NONE until the
  # end of configuration, at which point it will be set to $ac_default_prefix.
  # Unless the user specifies --exec-prefix, exec-prefix is set to NONE until
  # the end of configuration, at which point it's set to '${prefix}'.
  # Sheesh.  So do the expansion, then back it out.
  save_prefix=$prefix
  save_exec_prefix=$exec_prefix
  if test "x$prefix" = xNONE ; then
    prefix=$ac_default_prefix
  fi
  if test "x$exec_prefix" = xNONE ; then
    exec_prefix=$prefix
  fi
  expanded_libdir=$libdir
  while expr "$expanded_libdir" : '.*$.*' >/dev/null 2>&1 ; do
    eval expanded_libdir=$expanded_libdir
  done
  prefix=$save_prefix
  exec_prefix=$save_exec_prefix

  AC_ARG_ENABLE([relocatable],
    [AS_HELP_STRING([--enable-relocatable],[whether prefix in installed .pc files should be setup relative to pcfiledir])],
    [coin_enable_relocatable=$enableval],
    [coin_enable_relocatable=no])
  AM_CONDITIONAL([COIN_RELOCATABLE], [test $coin_enable_relocatable = yes])
])


#                           COIN_PROG_LIBTOOL                             #

# Set up libtool parameters and create and patchup libtool
# (https://www.gnu.org/software/libtool/manual/html_node/LT_005fINIT.html)
# Pass no-win32-dll to omit passing win32-dll to LT_INIT.
# Checks whether ar-lib wrapper needs to be used.
#
# Set up LT_LDFLAGS variable, which user can initialize and configure augments:
# - Add version info flag using the libtool library info, if defined,
#   otherwise use the project version if a full major.minor.release number is available.
# - Add -no-undefined as shared libraries should have no undefined symbols.
#   For Windows DLLs, it is mandatory to add this.
#
# Defines automake conditional COIN_STATIC_BUILD on whether we build
# shared or static. Used in .pc files.
#
# Patch libtool to circumvent some issues when using MSVC and MS lib.
# This needs to be run after config.status has created libtool.
# 1. Relax check which libraries can be used when linking a DLL.
#    libtool's func_mode_link() would reject linking a .lib file when building a DLL,
#    even though this .lib file may just be the one that eventually loads a depending DLL,
#    e.g., mkl_core_dll.lib. Setting deplibs_check_method=pass_all will still print a
#    warning, but the .lib is still passed to the linker.
# 2. Ensure always_export_symbols=no if win32-dll. Even when we pass win32-dll,
#    libtool forces always_export_symbols=yes for --tag=CXX if using MS compiler.
#    This leads to a nm call that collects ALL C-functions from a library
#    and explicitly dll-exporting them, leading to warnings about duplicates
#    regarding those that are properly marked for dll-export in the source.
# 3. Do not add mkl_*.lib to old_deplibs, which can result in trying to unpack and repack
#    the MKL libraries (which are pretty big). Instead, treat them like other -l<...> libs.
# 4. Add MKL libraries to dependency_libs in .la file, which I guess should be
#    the case due to point 5.
#
# Patch libtool also to circumvent some issues when using MinGW (Msys+GCC).
# 1. Relax check which libraries can be used when linking a DLL.
#    libtool's func_mode_link() would reject linking MinGW system libraries,
#    e.g., -lmingw32, when building a DLL, because it does not find this
#    library in the installation path, and then falls back to build only
#    static libraries. Setting deplibs_check_method=pass_all will avoid
#    this faulty check.

AC_DEFUN([AC_COIN_PROG_LIBTOOL],
[
dnl checkout AR and decide whether to use ar-lib wrapper
  AM_PROG_AR

dnl create libtool
  LT_INIT(disable-static pic-only m4_bmatch($1,no-win32-dll,,win32-dll))

  case "$am_cv_ar_interface" in
    lib )
      AC_CONFIG_COMMANDS([libtoolclpatch],
        [sed -e '/^deplibs_check_method/s/.*/deplibs_check_method="pass_all"/g' \
             m4_bmatch($1,no-win32-dll,,[-e 's|always_export_symbols=yes|always_export_symbols=no|g']) \
             -e '/func_append old_deplibs/s/\(.*\)/case $arg in *mkl_*.lib) ;; *) \1 ;; esac/g' \
             -e '/static library .deplib is not portable/a case $deplib in *mkl_*.lib) newdependency_libs="$deplib $newdependency_libs" ;; esac' \
         libtool > libtool.tmp
         mv libtool.tmp libtool
         chmod 755 libtool])
      ;;
    * )
      case $build in
        *-mingw* )
          AC_CONFIG_COMMANDS([libtoolmingwpatch],
            [sed -e '/^deplibs_check_method/s/.*/deplibs_check_method="pass_all"/g' libtool > libtool.tmp
             mv libtool.tmp libtool
             chmod 755 libtool])
        ;;
      esac
      ;;
  esac

  AC_ARG_VAR(LT_LDFLAGS,[Flags passed to libtool when building libraries or executables that are installed])

  m4_ifnblank(AC_COIN_LIBTOOL_VERSIONINFO,
    LT_LDFLAGS="$LT_LDFLAGS -version-info AC_COIN_LIBTOOL_VERSIONINFO"
    AC_MSG_NOTICE(libtool version info: -version-info AC_COIN_LIBTOOL_VERSIONINFO),
    m4_ifnblank(AC_PACKAGE_VERSION_MAJOR,
     m4_ifnblank(AC_PACKAGE_VERSION_MINOR,
      m4_ifnblank(AC_PACKAGE_VERSION_RELEASE,
       LT_LDFLAGS="$LT_LDFLAGS -version-number AC_PACKAGE_VERSION_MAJOR:AC_PACKAGE_VERSION_MINOR:AC_PACKAGE_VERSION_RELEASE"
       AC_MSG_NOTICE(libtool version info: -version-number AC_PACKAGE_VERSION_MAJOR:AC_PACKAGE_VERSION_MINOR:AC_PACKAGE_VERSION_RELEASE)
    ))))

  LT_LDFLAGS="$LT_LDFLAGS -no-undefined"

  AM_CONDITIONAL([COIN_STATIC_BUILD],[test "$enable_shared" = no])
])


#                    COIN_PROG_CC   COIN_PROG_CXX                         #

# Macros to find C and C++ compilers according to specified lists of compiler
# names. For Fortran compilers, see coin_fortran.m4.
#
# If enable-msvc, then test for Intel (on Windows) and MS C/C++ compiler
# explicitly and add the compile wrapper before AC_PROG_CC/CXX. The compile
# wrapper works around issues related to finding MS link.exe. (Unix link.exe
# occurs first in PATH, which causes compile and link checks to fail.) For
# the same reason, set LD to use the compile wrapper. If CC/CXX remains unset
# (neither icl or cl was found, and CC was not set by the user), stop with
# an error.
#
# Declares ADD_C(XX)FLAGS variables for additional compiler flags.
#
# Note that automake redefines AC_PROG_CC to invoke _AM_PROG_CC_C_O (an
# equivalent to AC_PROG_CC_C_O, plus wrap CC in the compile script if needed)
# and _AM_DEPENDENCIES (`dependency style'). Look in aclocal.m4 to see this.

AC_DEFUN_ONCE([AC_COIN_PROG_CC],
[
  AC_REQUIRE([AC_COIN_ENABLE_MSVC])

dnl Setting up libtool with LT_INIT will AC_REQUIRE AC_PROG_CC, but we want
dnl to invoke it from this macro first so that we can supply a parameter.
  AC_BEFORE([$0],[LT_INIT])

  if test $enable_msvc = yes ; then
    AC_CHECK_PROGS(CC, [icl cl])
    if test -n "$CC" ; then
      CC="$am_aux_dir/compile $CC"
      ac_cv_prog_CC="$CC"
      LD="$CC"
      : ${AR:=lib}
    else
      AC_MSG_ERROR([Neither MS nor Intel C compiler found in PATH and CC is unset.])
    fi
  fi

dnl Look for some C compiler and check that it works. If the user has set CC
dnl or we found icl/cl above, this should not overwrite CC. Unlike the macros
dnl that establish C++ or Fortran compilers, PROG_CC also takes care of adding
dnl the compile wrapper if needed.
  AC_PROG_CC([gcc clang cc icc icl cl cc xlc xlc_r pgcc])

  AC_ARG_VAR(ADD_CFLAGS,[Additional C compiler options (if not overwriting CFLAGS)])
])

# Note that automake redefines AC_PROG_CXX to invoke _AM_DEPENDENCIES
# (`dependency style') but does not add an equivalent for AC_PROG_CXX_C_O,
# hence we need an explicit call.

AC_DEFUN_ONCE([AC_COIN_PROG_CXX],
[
  AC_REQUIRE([AC_COIN_ENABLE_MSVC])
  AC_REQUIRE([AC_COIN_PROG_CC])

  if test $enable_msvc = yes ; then
    AC_CHECK_PROGS(CXX, [icl cl])
    if test -n "$CXX" ; then
      CXX="$am_aux_dir/compile $CXX"
      ac_cv_prog_CXX="$CXX"
      LD="$CXX"
      : ${AR:=lib}
    else
      AC_MSG_ERROR([Neither MS nor Intel C++ compiler found in PATH and CXX is unset.])
    fi
  fi

dnl Look for some C++ compiler and check that it works. If the user has set
dnl CXX or we found icl/cl above, this should not overwrite CXX.
  AC_PROG_CXX([g++ clang++ c++ pgCC icpc gpp cxx cc++ icl cl FCC KCC RCC xlC_r aCC CC])

dnl If the compiler does not support -c -o, wrap it with the compile script.
  AC_PROG_CXX_C_O
  if test $ac_cv_prog_cxx_c_o = no ; then
    CXX="$am_aux_dir/compile $CXX"
  fi

  AC_ARG_VAR(ADD_CXXFLAGS,[Additional C++ compiler options (if not overwriting CXXFLAGS)])
])

#                             COIN_CXXLIBS                                #

# Determine the C++ runtime libraries required for linking a C++ library
# with a Fortran or C compiler.  The result is available in CXXLIBS.

AC_DEFUN([AC_COIN_CXXLIBS],
[AC_REQUIRE([AC_COIN_PROG_CXX])dnl
AC_LANG_PUSH(C++)
AC_ARG_VAR(CXXLIBS,[Libraries necessary for linking C++ code with non-C++ compiler])
if test -z "$CXXLIBS"; then
  if test "$GXX" = "yes"; then
    case "$CXX" in
      icpc* | */icpc*)
        CXXLIBS="-lstdc++"
        ;;
      *)
        # clang uses libc++ as the default standard C++ library, not libstdc++
        # this test is supposed to recognize whether the compiler is clang
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <ciso646>]], [[
#ifndef _LIBCPP_VERSION
       choke me
#endif
          ]])],
          [CXXLIBS="-lc++"],
          [CXXLIBS="-lstdc++ -lm"])  dnl The -lm works around an issue with libtool removing -lm from the dependency_libs in an .la file (look for postdeps= in libtool)
        ;;
    esac
  else
    case $build in
     *-mingw* | *-cygwin* | *-msys* )
       if test "$enable_msvc" = yes ; then
         CXXLIBS=nothing
       fi;;
     *-linux-*)
      case "$CXX" in
      icpc* | */icpc*)
        CXXLIBS="-lstdc++"
             ;;
      pgCC* | */pgCC*)
        CXXLIBS="-lstd -lC -lc"
             ;;
      esac;;
    *-ibm-*)
      CXXLIBS="-lC -lc"
      ;;
    *-hp-*)
      CXXLIBS="-L/opt/aCC/lib -l++ -lstd_v2 -lCsup_v2 -lm -lcl -lc"
      ;;
    *-*-solaris*)
      CXXLIBS="-lCstd -lCrun"
    esac
  fi
fi
if test -z "$CXXLIBS"; then
  AC_MSG_WARN([Could not automatically determine CXXLIBS (C++ link libraries; necessary if main program is in Fortran or C).])
else
  AC_MSG_NOTICE([Assuming that CXXLIBS is $CXXLIBS.])
fi
if test x"$CXXLIBS" = xnothing; then
  CXXLIBS=
fi
AC_LANG_POP(C++)
]) # AC_COIN_CXXLIBS


#                            COIN_RPATH_FLAGS                             #

# This macro, in case shared objects are used, defines a variable
# RPATH_FLAGS that can be used by the linker to hardwire the library
# search path for the given directories.  This is useful for example
# Makefiles

AC_DEFUN([AC_COIN_RPATH_FLAGS],
[RPATH_FLAGS=

if test $enable_shared = yes; then
  case $build in
    *-linux-*)
      if test "$GCC" = "yes"; then
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -Wl,--rpath -Wl,$dir"
        done
      fi ;;
    *-darwin*)
        RPATH_FLAGS=nothing ;;
    *-ibm-*)
      case "$CC" in
      xlc* | */xlc* | mpxlc* | */mpxlc*)
        RPATH_FLAGS=nothing ;;
      esac ;;
    *-hp-*)
        RPATH_FLAGS=nothing ;;
    *-mingw* | *-msys* )
        RPATH_FLAGS=nothing ;;
    *-*-solaris*)
        RPATH_FLAGS=
        for dir in $1; do
          RPATH_FLAGS="$RPATH_FLAGS -R$dir"
        done
  esac

  if test "$RPATH_FLAGS" = ""; then
    AC_MSG_WARN([Could not automatically determine how to tell the linker about automatic inclusion of the path for shared libraries.  The test examples might not work if you link against shared objects.  You will need to set the LD_LIBRARY_PATH, DYLP_LIBRARY_PATH, or LIBDIR variable manually.])
  fi
  if test "$RPATH_FLAGS" = "nothing"; then
    RPATH_FLAGS=
  fi
fi

AC_SUBST(RPATH_FLAGS)
]) # AC_COIN_RPATH_FLAGS


#                   COIN_DEFINENAMEMANGLING                               #

# COIN_DEFINENAMEMANGLING (name,namemangling)
# -------------------------------------------------------------------------
# Determine C/C++ name mangling to allow linking with header-less libraries.
#  name ($1) an identifier to prefix C macro names
#  namemangling ($2) the name mangling scheme as determined by COIN_NAMEMANGLING
#                    or COIN_TRY_LINK
#
# Defines the C macros $1_FUNC and $1_FUNC_ (in uppercase) to be used for
# declaring functions from library $1.

# -------------------------------------------------------------------------

AC_DEFUN([AC_COIN_DEFINENAMEMANGLING],
[
  AH_TEMPLATE($1_FUNC, [Define to a macro mangling the given C identifier (in lower and upper case).])
  AH_TEMPLATE($1_FUNC_, [As $1_FUNC, but for C identifiers containing underscores.])
  case "$2" in
   "lower case, no underscore, no extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [name])
      AC_DEFINE($1_FUNC_[(name,NAME)], [name]) ;;
   "lower case, no underscore, extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [name])
      AC_DEFINE($1_FUNC_[(name,NAME)], [name [##] _]) ;;
   "lower case, underscore, no extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [name [##] _])
      AC_DEFINE($1_FUNC_[(name,NAME)], [name [##] _]) ;;
   "lower case, underscore, extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [name [##] _])
      AC_DEFINE($1_FUNC_[(name,NAME)], [name [##] __]) ;;
   "upper case, no underscore, no extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [NAME])
      AC_DEFINE($1_FUNC_[(name,NAME)], [NAME]) ;;
   "upper case, no underscore, extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [NAME])
      AC_DEFINE($1_FUNC_[(name,NAME)], [NAME [##] _]) ;;
   "upper case, underscore, no extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [NAME [##] _])
      AC_DEFINE($1_FUNC_[(name,NAME)], [NAME [##] _]) ;;
   "upper case, underscore, extra underscore")
      AC_DEFINE($1_FUNC[(name,NAME)],  [NAME [##] _])
      AC_DEFINE($1_FUNC_[(name,NAME)], [NAME [##] __]) ;;
   *)
      AC_MSG_WARN([Unsupported or unknown name-mangling scheme: $2])
      ;;
  esac
])


#                            COIN_TRY_LINK                                #

# This is a helper macro for checking if a library can be linked based on
# a function name only.
#   COIN_TRY_LINK(func,lflags,pcfiles,action-if-success,action-if-failed,msg)
#  func ($1) the name of the function to try to link
#  lflags ($2) linker flags to use
#  pcfiles ($3) pc files to query for additional linker flags
#  action-if-success ($4) commands to execute if any linking was successful
#  action-if-failed ($5) commands to execute if no linking was successful
#  msg ($6) 'no' to suppress 'checking for' messages
#
# The macro tries different name mangling schemes and expands into
# action-if-success for the first one that succeeds.
# It sets variable func_namemangling according to the found name mangling
# scheme, which can be used as input for COIN_DEFINENAMEMANGLING.

AC_DEFUN([AC_COIN_TRY_LINK],
[
dnl setup LIBS by adding $2 and those from $3
  ac_save_LIBS="$LIBS"
  m4_ifnblank([$2], [LIBS="$2 $LIBS"])
  m4_ifnblank([$3],
  [ AC_REQUIRE([AC_COIN_HAS_PKGCONFIG])
    if test -n "$3" ; then
      temp_LFLAGS=`PKG_CONFIG_PATH="$COIN_PKG_CONFIG_PATH" $PKG_CONFIG --libs $pkg_static $3`
      LIBS="$temp_LFLAGS $LIBS"
    fi
  ])

  $1_namemangling=unknown

  for ac_extra in "no extra underscore" "extra underscore" ; do
    for ac_case in "lower case" "upper case" ; do
      for ac_trail in "underscore" "no underscore" ; do
        case $ac_case in
          "lower case")
            ac_name=m4_tolower($1)
            ;;
          "upper case")
            ac_name=m4_toupper($1)
            ;;
        esac
        if test "$ac_trail" = "underscore" ; then
          ac_name=${ac_name}_
        fi
        if test "$ac_extra" = "extra underscore" ; then
          ac_name=${ac_name}_
        fi
        m4_if([$6],[no],[],[AC_MSG_CHECKING([for function $ac_name in $LIBS])])
        AC_TRY_LINK_FUNC([$ac_name],
          [$1_namemangling="${ac_case}, ${ac_trail}, ${ac_extra}"
           ac_success=yes],
          [ac_success=no])
        m4_if([$6],[no],[],[AC_MSG_RESULT([$ac_success])])
        if test $ac_success = yes ; then
          break 3
        fi
      done
    done
  done
  LIBS=$ac_save_LIBS

  if test $ac_success = yes ; then
    m4_ifblank([$4],[:],[$4])
    m4_ifnblank([$5],[else $5])
  fi
])


#                           COIN_HAS_PKGCONFIG                            #

# Check whether a suitable pkg-config tool is available.  pkgconf is the
# up-and-coming thing, replacing pkg-config, so it is preferred.  The default
# minimal version number is 0.16.0 because COIN-OR .pc files include a
# URL field which breaks pkg-config version <= 0.15. A more recent minimum
# version can be specified as a parameter.  Portions of the macro body are
# derived from macros in pkg.m4.

# If a pkg-config tool is found, then the variable PKGCONFIG is set
# to its path, otherwise it is set to "". Further, the AM_CONDITIONAL
# COIN_HAS_PKGCONFIG is set and PKGCONFIG is AC_SUBST'ed.

# Finally, the search path for .pc files is assembled from the value of
# $prefix and $PKG_CONFIG_PATH in a variable COIN_PKG_CONFIG_PATH, which is
# also AC_SUBST'ed.  The search path will include the installation directory
# for .pc files for COIN-OR packages. COIN-OR .pc files are installed
# in ${libdir}/pkgconfig and COIN_INITIALIZE takes care of setting up
# $expanded_libdir based on $libdir.

# Of course, this whole house of cards balances on the shaky assumption that
# the user is sane and has installed all packages in the same place and does
# not change that place when make executes. If not, well, it's their
# responsibility to augment PKG_CONFIG_PATH in the environment.

AC_DEFUN([AC_COIN_HAS_PKGCONFIG],
[
  AC_ARG_VAR([PKG_CONFIG],[path to pkg-config utility])

dnl This is a modified version of PKG_PROG_PKG_CONFIG from pkg.m4.
  if test -z "$PKG_CONFIG" ; then
   AC_CHECK_TOOLS([PKG_CONFIG],[pkgconf pkg-config])
  fi
  if test -n "$PKG_CONFIG" ; then
    pkg_min_version=m4_default([$1],[0.16.0])
    AC_MSG_CHECKING([$PKG_CONFIG is at least version $pkg_min_version])
    if $PKG_CONFIG --atleast-pkgconfig-version $pkg_min_version ; then
      pkg_version=`$PKG_CONFIG --version`
      AC_MSG_RESULT([yes: $pkg_version])
    else
     AC_MSG_RESULT([no])
     PKG_CONFIG=""
   fi
  fi

dnl Check if PKG_CONFIG supports the short-errors flag.
dnl This is a modified version of _PKG_SHORT_ERRORS_SUPPORTED from pkg.m4.
  if test -n "$PKG_CONFIG" &&
     $PKG_CONFIG --atleast-pkgconfig-version 0.20 ; then
    pkg_short_errors=" --short-errors "
  else
    pkg_short_errors=""
  fi

dnl Check whether -static option of pkg-config should be used when requesting
dnl libs
  pkg_static=
  if test -n "$PKG_CONFIG" ; then
    case "$LDFLAGS" in "-static" | "* -static*" ) pkg_static=--static ;; esac
  fi

dnl Create a automake conditional and PKG_CONFIG variable
  AM_CONDITIONAL([COIN_HAS_PKGCONFIG], [test -n "$PKG_CONFIG"])
  AC_SUBST(PKG_CONFIG)

  COIN_PKG_CONFIG_PATH="${PKG_CONFIG_PATH}"
  AC_SUBST(COIN_PKG_CONFIG_PATH)

  COIN_PKG_CONFIG_PATH="${expanded_libdir}/pkgconfig:${COIN_PKG_CONFIG_PATH}"
  if test -n "$PKG_CONFIG"; then
    AC_MSG_NOTICE([$PKG_CONFIG path is "$COIN_PKG_CONFIG_PATH"])
  fi
])  # COIN_HAS_PKGCONFIG


#                      COIN_CHK_MOD_EXISTS                                #

# COIN_CHK_MOD_EXISTS(MODULE, PACKAGES,
#                     [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])

# Check to see whether a particular module (set of packages) exists. Derived
# from PKG_CHECK_MODULES() from pkg.m4, but set only the variables $1_VERSIONS
# and $1_PKG_ERRORS. PACKAGES ($2) is a space-separated list of .pc file names
# (without the .pc suffix). Note that variable values will have one line per
# package (embedded newlines) if more than one package is given in PACKAGES.

AC_DEFUN([AC_COIN_CHK_MOD_EXISTS],
[
  AC_REQUIRE([AC_COIN_HAS_PKGCONFIG])

  if test -n "$PKG_CONFIG" ; then
    if PKG_CONFIG_PATH="$COIN_PKG_CONFIG_PATH" $PKG_CONFIG --exists "$2" ; then
      m4_toupper($1)_VERSIONS=`PKG_CONFIG_PATH="$COIN_PKG_CONFIG_PATH" $PKG_CONFIG --modversion "$2" 2>/dev/null | tr '\n' ' '`
      $3
    else
      m4_toupper($1)_PKG_ERRORS=`PKG_CONFIG_PATH="$COIN_PKG_CONFIG_PATH" $PKG_CONFIG $pkg_short_errors --errors-to-stdout --print-errors "$2"`
      $4
    fi
  else
    AC_MSG_ERROR("Cannot check for existence of module $1 without pkgconf")
  fi
])


#                          COIN_CHK_HERE                                  #

# COIN_CHK_HERE([prim],[client packages],[pcfile])

# Augment the _LFLAGS, _CFLAGS, and _PCFILES variables of the client
# packages with the values from PRIM_LFLAGS_PUB, PRIM_CFLAGS_PUB, and
# PRIM_PCFILES_PUB. This macro is intended for the case where a single project
# builds several objects and one object includes another. For example,
# the various OsiXxxLib solvers, which depend on OsiLib. We can't consult
# osi.pc (it's not installed yet) but the relevant variables are ready at
# hand. The name of prim is often different from the name of the .pc file
# ($3), hence the separate parameter. If $3 is not given, it defaults to
# tolower($1).

# This macro should be called after FINALIZE_FLAGS is invoked for the
# client packages, for two reasons: First, COIN-OR packages tend to use
# .pc files, so we're probably adding a package to _PCFILES that isn't yet
# installed. Also, FINALIZE_FLAGS will have accessed the .pc files already
# in _PCFILES and expanded them into _LIBS and _CFLAGS, saving the original
# _LIBS and _CFLAGS in _LIBS_NOPC and _CFLAGS_NOPC.

AC_DEFUN([AC_COIN_CHK_HERE],
[
dnl Make sure the necessary variables exist for each client package.
  m4_foreach_w([myvar],[$2],
    [AC_SUBST(m4_toupper(myvar)_LFLAGS)
     AC_SUBST(m4_toupper(myvar)_CFLAGS)
     AC_SUBST(m4_toupper(myvar)_PCFILES)
    ])

dnl Add the .pc file and augment LFLAGS and CFLAGS.
dnl From the CFLAGS of $1, remove the -D$1_BUILD, though.
    m4_foreach_w([myvar],[$2],
      [m4_toupper(myvar)_PCFILES="$m4_toupper(myvar)_PCFILES m4_default($3,m4_tolower($1))"
       m4_toupper(myvar)_LFLAGS="$m4_toupper(myvar)_LFLAGS $m4_toupper($1)_LFLAGS"
       m4_toupper(myvar)_CFLAGS="$m4_toupper(myvar)_CFLAGS `echo $m4_toupper($1)_CFLAGS | sed -e s/-D[]m4_toupper($1)_BUILD//`"

       # Define BUILDTOOLS_DEBUG to enable debugging output
       if test "$BUILDTOOLS_DEBUG" = 1 ; then
         AC_MSG_NOTICE([CHK_HERE adding $1 to myvar:])
         AC_MSG_NOTICE([m4_toupper(myvar)_PCFILES: "${m4_toupper(myvar)_PCFILES}"])
         AC_MSG_NOTICE([m4_toupper(myvar)_LFLAGS: "${m4_toupper(myvar)_LFLAGS}"])
         AC_MSG_NOTICE([m4_toupper(myvar)_CFLAGS: "${m4_toupper(myvar)_CFLAGS}"])
       fi
      ])
])   # COIN_CHK_HERE


#                      COIN_DEF_PRIM_ARGS                                 #

# COIN_DEF_PRIM_ARGS([prim],[base],[lflags],[cflags],[dflags],[build])

# This is a utility macro to handle the standard arguments that COIN-OR
# configuration files supply for a component (package or library):
#   --with-prim: use primitive (yes / no / special)
#   --with-prim-lflags: linker flags for the primitive
#   --with-prim-cflags: preprocessor & compiler flags for the primitive
#   --with-prim-data: data directory for the primitive
# These arguments allow the user to override default macro behaviour from the
# configure command line.
# The string prim, lower-cased, is used in the flag name.
# The parameters base, lflags, cflags, and dflags have the value yes or no and
# determine whether a --with option will be defined for prim, lflags, cflags,
# and data, respectively. They must be literals, as the macro uses them to cut
# out unused options. To use the results, construct the name of the shell
# variable as specified in the autoconf doc'n for ARG_WITH.

# Setting the final parameter to 'default_build' will cause the phrase
# "'build' will look for a COIN ThirdParty package" to be inserted in the
# documentation for --with-prim.

AC_DEFUN([AC_COIN_DEF_PRIM_ARGS],
[
  m4_define([extraHelp],[
    m4_normalize(Use $1. [If an argument is given,]
      ['yes' is equivalent to --with-m4_tolower($1),]
      m4_case($6,[default_build],
      ['no' is equivalent to --without-m4_tolower($1)[,]
       'build' will look for a COIN-OR ThirdParty package.],
      ['no' is equivalent to --without-m4_tolower($1).])
      m4_case($3$4$5,nonono,,
        nonoyes,
        Any other argument is applied as for --with-m4_tolower($1)-data,
        noyesno,
        Any other argument is applied as for --with-m4_tolower($1)-cflags,
        noyesyes,
        Any other argument is applied as for --with-m4_tolower($1)-cflags,
        Any other argument is applied as for --with-m4_tolower($1)-lflags))])

  m4_if($2,yes,
    [AC_ARG_WITH([m4_tolower($1)],
       AS_HELP_STRING([--with-m4_tolower($1)],extraHelp))])

  m4_if($3,yes,
    [AC_ARG_WITH([m4_tolower($1)-lflags],
       AS_HELP_STRING([--with-m4_tolower($1)-lflags],
         [Linker flags for $1 appropriate for your environment.
          (Most often, -l specs for libraries.)]))])

  m4_if($4,yes,
    [AC_ARG_WITH([m4_tolower($1)-cflags],
       AS_HELP_STRING([--with-m4_tolower($1)-cflags],
         [Compiler flags for $1 appropriate for your environment.
          (Most often, -I specs for header file directories.)]))])

  m4_if($5,yes,
    [AC_ARG_WITH([m4_tolower($1)-data],
       AS_HELP_STRING([--with-m4_tolower($1)-data],
         [A data directory specification for $1 appropriate for your
          environment.]))])
])   # COIN_DEF_PRIM_ARGS

#                     FIND_PRIM_LIB (obsolete)                            #

# COIN_FIND_PRIM_LIB([prim],[lflgs],[cflgs],[dflgs],
#                    [func],[other libraries],
#                    [dfltaction],[cmdlineopts])

# Determine whether we can use primary library prim ($1) and assemble
# information on the required linker flags (prim_lflags), compiler flags
# (prim_cflags), and data directories (prim_data) as specified by cmdlineopts.
# Run a link check if the user provides [func]. Linker flags for the link are
# the concatenation of [lflgs] and [other libraries].

# cmdlineopts ($8) specifies the set of configure command line options
# defined and processed: 'nodata' produces --with-prim, --with-prim-lflags,
# and --with-prim-cflags; 'dataonly' produces only --with-prim and
# --with-prim-data; anything else ('all' works well) produces all four
# command line options. Shell code produced by the macro is tailored based
# on cmdlineopts. `nodata' is the default.

# --with-prim is interpreted as follows:
#   * --with-prim=no is equivalent to --without-prim.
#   * --with-prim or --with-prim=yes is equivalent to
#       --with-prim-lflags=-lprim
#       --with-prim-data=/usr/local/share
#   * --with-prim=build attempts to invent something that will find a COIN
#     ThirdParty library or data
#       --with-prim-lflags="-L\$(libdir) -lcoinprim"
#       --with-prim-cflgs="-I\$(pkgincludedir)/ThirdParty"
#       --with-prim-data="\$(pkgdatadir)"
#   * Any other value is taken as equivalent to
#       --with-prim-data=value (dataonly) or
#       --with-prim-lflags=value (anything else)

# The algorithm first checks for a user-specified value of --with-prim;
# if this is no, prim is skipped. Next, it looks for user specified values
# given with command line parameters --with-prim-lflags, --with-prim-cflags,
# and --with-prim-data. If none of these are specified, it will use the values
# passed as parameters. It's all or none: any command line options disable all
# parameters.

# Default action ($7) (no, yes, build) is the default value of --with-prim
# if the user offers no guidance via command line parameters. The (hardwired)
# default is yes. `build' doesn't have a hope of working except for COIN
# ThirdParty packages, and even then it's pretty shaky. You should be
# using CHK_PKG, because COIN packaging for ThirdParty software creates a .pc
# file.

# The macro doesn't test that the specified values actually work unless
# [func] is given as a parameter. This is deliberate --- there's no guarantee
# that the specified library can be accessed just yet with the specified
# flags. Except for the link check, all we're doing here is filling in
# variables using a complicated algorithm.

AC_DEFUN([AC_COIN_FIND_PRIM_LIB_OBS],
[
  dflt_action=m4_default([$7],[yes])

  # Initialize variables for the primary library.
  m4_tolower(coin_has_$1)=noInfo
  m4_tolower($1_lflags)=
  m4_tolower($1_cflags)=
  m4_tolower($1_data)=

  # --with-prim is always present.
dnl If the client specified dataonly, its value is assigned to prim_data.
  withval="$m4_tolower(with_$1)"
  if test -n "$withval" ; then
    case "$withval" in
      no )
        m4_tolower(coin_has_$1)=skipping
        ;;
      yes )
        m4_tolower(coin_has_$1)=requested
        ;;
      build )
        m4_tolower(coin_has_$1)=build
        ;;
      * )
        m4_tolower(coin_has_$1)=yes
        m4_if(m4_default($8,nodata),dataonly,
          [m4_tolower($1_data)="$withval"],
          [m4_tolower($1_lflags)="$withval"])
        ;;
    esac
  fi

  m4_if(m4_default($8,nodata),dataonly,[],
    [# --with-prim-lflags and --with-prim-cflags are present
     # Specifying --with-prim=no overrides the individual options for lflags and cflags.
     if test "$m4_tolower(coin_has_$1)" != skipping ; then
       withval="$m4_tolower(with_$1_lflags)"
       if test -n "$withval" ; then
         case "$withval" in
           build | no | yes )
             AC_MSG_ERROR(["$withval" is not valid here; please specify linker flags appropriate for your environment.])
             ;;
           * )
             m4_tolower(coin_has_$1)=yes
             m4_tolower($1_lflags)="$withval"
             ;;
         esac
       fi

       withval="$m4_tolower(with_$1_cflags)"
       if test -n "$withval" ; then
         case "$withval" in
           build | no | yes )
             AC_MSG_ERROR(["$withval" is not valid here; please specify compiler flags appropriate for your environment.])
             ;;
           * )
             m4_tolower(coin_has_$1)=yes
             m4_tolower($1_cflags)="$withval"
             ;;
         esac
       fi
     fi])

  m4_if(m4_default($8,nodata),nodata,[],
    [# --with-prim-data will be present
     # Specifying --with-prim=no overrides the individual option for data.
     if test "$m4_tolower(coin_has_$1)" != skipping ; then
       withval="$m4_tolower(with_$1_data)"
       if test -n "$withval" ; then
         case "$withval" in
           build | no | yes )
             AC_MSG_ERROR(["$withval" is not valid here; please give a data directory specification appropriate for your environment.])
             ;;
           * )
             m4_tolower(coin_has_$1)=yes
             m4_tolower($1_data)="$withval"
             ;;
         esac
       fi
     fi])

  # At this point, coin_has_prim can be one of
  # - noInfo (no user options specified),
  # - skipping (user said no),
  # - requested,
  # - build (user said yes or build and gave no further guidance), or
  # - yes (user specified one or more --with-prim options).
  # If we're already at yes or skipping, we're done looking.

  # If there are no user options (noInfo) and the default is no, we're skipping.
  # Otherwise, the default must be yes or build; consider the package requested.
  # A default action we don't recognise defaults to yes.
  if test "$m4_tolower(coin_has_$1)" = noInfo ; then
    case $dflt_action in
      no )
        m4_tolower(coin_has_$1)=skipping
        ;;
      build )
        m4_tolower(coin_has_$1)=build
        ;;
      * )
        m4_tolower(coin_has_$1)=requested
        ;;
    esac
  fi

  # Now coin_has_prim can be one of skipping, yes, build, or requested.
  # For build or requested, use the parameter values or invent some.
  case $m4_tolower(coin_has_$1) in
    build | requested)
      m4_if(m4_default($8,nodata),dataonly,[],
        [m4_ifnblank([$2],
           [m4_tolower($1_lflags)="$2"],
           [if test "$m4_tolower(coin_has_$1)" = build ; then
              m4_tolower($1_lflags)="-L\$(libdir) -l[]m4_tolower(coin$1)"
            else
              m4_tolower($1_lflags)="-l[]m4_tolower($1)"
            fi])
         m4_ifnblank([$3],
           [m4_tolower($1_cflags)="$3"],
           [if test "$m4_tolower(coin_has_$1)" = build ; then
              m4_tolower($1_cflags)="-I\$(pkgincludedir)/ThirdParty"
            fi])])
      m4_if(m4_default($8,nodata),nodata,[],
        [m4_tolower($1_data)=m4_default([$3],
           [if test "$m4_tolower(coin_has_$1)" = build ; then
              m4_tolower($1_data)="\$(pkgdatadir)"
            else
              m4_tolower($1_data)="/usr/local/share"
            fi])])
dnl FIXME This catches user-supplied as well as made up. And not checking is
dnl a feature. We need to discuss  -lh, 210114-
dnl go to linkcheck if we have a symbol to look for ($5)
dnl otherwise skip project, since we cannot check whether made-up l/cflags work
      m4_ifnblank([$5],m4_tolower(coin_has_$1)=yes,m4_tolower(coin_has_$1)=skipping)
      ;;
    skipping | yes )
      ;;
    * )
      AC_MSG_WARN([Unexpected status "$m4_tolower(coin_has_$1)" in COIN_FIND_PRIM_LIB])
      ;;
  esac

  # At this point, coin_has_prim is yes or skipping.
  m4_ifnblank([$5],
    [if test $m4_tolower(coin_has_$1) != skipping ; then
       # Time to run a link check since we have a function name.
       # Use whatever we've collected for lflags, plus other libraries ($6).
       ac_save_LIBS=$LIBS
       LIBS="$m4_tolower($1_lflags) $6"
       AC_TRY_LINK_FUNC([$5],[],m4_tolower(coin_has_$1)=no)
       LIBS=$ac_save_LIBS
     fi])

  # The final value of coin_has_prim will be yes, no, or skipping. No means that
  # the link check failed. Yes means that we passed the link check, or no link
  # check was performed. Skipping means the user said `don't use.'

  # Define BUILDTOOLS_DEBUG to enable debugging output
  if test "$BUILDTOOLS_DEBUG" = 1 ; then
    AC_MSG_NOTICE([FIND_PRIM_LIB result for $1: "$m4_tolower(coin_has_$1)"])
    AC_MSG_NOTICE([Collected values for package '$1'])
    AC_MSG_NOTICE([m4_tolower($1_lflags) is "$m4_tolower($1_lflags)"])
    AC_MSG_NOTICE([m4_tolower($1_cflags) is "$m4_tolower($1_cflags)"])
    AC_MSG_NOTICE([m4_tolower($1_data) is "$m4_tolower($1_data)"])
    AC_MSG_NOTICE([m4_tolower($1_pcfiles) is "$m4_tolower($1_pcfiles)"])
  fi

])  # COIN_FIND_PRIM_LIB


#                          COIN_CHK_LIB (obsolete)                        #

# COIN_CHK_LIB([prim],[client packages],[lflgs],[cflgs],[dflgs],
#              [func],[other libraries],
#              [dfltaction],[cmdopts])

# Determine whether we can use primary library prim ($1) and assemble
# information on the required linker flags (prim_lflags), compiler flags
# (prim_cflags), and data directories (prim_data). A link check will be
# performed in COIN_FIND_PRIM_LIB if [func] is specified, using link flags
# formed by concatenating the values of [lflgs] and [other libraries].

# The configure command line options offered to the user are controlled
# by cmdopts ($9). 'nodata' offers --with-prim, --with-prim-lflags, and
# --with-prim-cflags; 'dataonly' offers --with-prim and --with-prim-data;
# 'all' offers all four. DEF_PRIM_ARGS and FIND_PRIM_LIB are tailored
# accordingly. The (hardwired) default is 'nodata'.

# Macro parameters lflgs ($3), cflgs ($4), and dflgs ($5) are used for
# --with-prim-lflags, --with-prim-cflags, and --with-prim-data if and only if
# there are no user-supplied values on the command line. It's all or nothing;
# any user-supplied value causes all macro parameters to be ignored.

# Default action ($8) (no, yes, build) is the default action if the user
# offers no guidance via command line parameters. Really, 'build' has no
# hope of working except for COIN ThirdParty packages. Don't use it for
# other COIN packages. You should be using CHK_PKG because COIN packaging
# for ThirdParty software creates a .pc file.

# Define an automake conditional COIN_HAS_PRIM to record the result. If we
# decide to use prim, also define a preprocessor symbol COIN_HAS_PRIM.

# Linker and compiler flag information will be propagated to the space-
# separated list of client packages ($2) using the _LFLAGS and _CFLAGS
# variables of client packages. These variables match Libs.private and
# Cflags.private, respectively, in a .pc file.

# Data directory information is used differently. Typically what's wanted is
# individual variables specifying the data directory for each primitive. Hence
# the macro defines PRIM_DATA for the primitive.

AC_DEFUN([AC_COIN_CHK_LIB_OBS],
[
  AC_MSG_CHECKING(m4_ifnblank($6,[for package $1 with function $6],[for package $1]))

dnl Make sure the necessary variables exist for each client package.
  m4_foreach_w([myvar],[$2],
    [AC_SUBST(m4_toupper(myvar)_LFLAGS)
     AC_SUBST(m4_toupper(myvar)_CFLAGS)
    ])

  # Check to see if the user has overridden configure parameters from the
  # environment.
  m4_tolower(coin_has_$1)=noInfo
  if test x"$COIN_SKIP_PROJECTS" != x ; then
    for pkg in $COIN_SKIP_PROJECTS ; do
      if test "$m4_tolower(pkg)" = "$m4_tolower($1)" ; then
        m4_tolower(coin_has_$1)=skipping
      fi
    done
  fi

dnl If we are not skipping this project, define and process the command line
dnl options according to the cmdopts parameter. Then invoke FIND_PRIM_PKG to do
dnl the heavy lifting.
  if test "$m4_tolower(coin_has_$1)" != skipping ; then
    m4_case(m4_default($9,nodata),
      nodata,[AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,no,$8)],
      dataonly,[AC_COIN_DEF_PRIM_ARGS([$1],yes,no,no,yes,$8)],
      [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,yes,$8)])
    AC_COIN_FIND_PRIM_LIB(m4_tolower($1),[$3],[$4],[$5],[$6],[$7],[$8],[$9])
    AC_MSG_RESULT([$m4_tolower(coin_has_$1)])
  else
    AC_MSG_RESULT([$m4_tolower(coin_has_$1) due to COIN_SKIP_PROJECTS])
  fi

  # Possibilities are `yes', `no', or `skipping'. Normalise to `yes' or `no'.
  if test "$m4_tolower(coin_has_$1)" != yes ; then
    m4_tolower(coin_has_$1)=no
  fi

dnl Create an automake conditional COIN_HAS_PRIM.
  AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),
                   [test $m4_tolower(coin_has_$1) = yes])

  # If we've located the package, define preprocessor symbol COIN_HAS_PRIM
  # and augment the necessary variables for the client packages.
  if test $m4_tolower(coin_has_$1) = yes ; then
    AC_DEFINE(m4_toupper(AC_PACKAGE_NAME)_HAS_[]m4_toupper($1),[1],
      [Define to 1 if the $1 package is available])
    m4_foreach_w([myvar],[$2],
      [m4_toupper(myvar)_LFLAGS="$m4_tolower($1_lflags) $m4_toupper(myvar)_LFLAGS"
       m4_toupper(myvar)_CFLAGS="$m4_tolower($1_cflags) $m4_toupper(myvar)_CFLAGS"
      ])

dnl Finally, set up PRIM_DATA, unless the user specified nodata.
    m4_if(m4_default([$9],nodata),nodata,[],
      [AC_SUBST(m4_toupper($1)_DATA)
       m4_toupper($1)_DATA=$m4_tolower($1_data)])
  fi
])   # COIN_CHK_LIB


#                           COIN_CHK_LIBM                           #

# COIN_CHK_LIBM([client1 client2 ...])
# Check if a library spec is needed for the math library. If something is
# needed, the macro adds the flags to CLIENT_LFLAGS for each client.

AC_DEFUN([AC_COIN_CHK_LIBM],
[
  AC_REQUIRE([AC_PROG_CC])

  m4_foreach_w([myvar],[$1],[AC_SUBST(m4_toupper(myvar)_LFLAGS)])

  coin_save_LIBS="$LIBS"
  LIBS=
  AC_SEARCH_LIBS([cos],[m],
    [if test "$ac_cv_search_cos" != 'none required' ; then
       m4_foreach_w([myvar],[$1],
         [m4_toupper(myvar)_LFLAGS="$ac_cv_search_cos $m4_toupper(myvar)_LFLAGS"
         ])
     fi])
  LIBS="$coin_save_LIBS"
]) # AC_COIN_CHK_LIBM


#                           COIN_CHK_ZLIB                               #

# COIN_CHK_ZLIB([client1 client2 ...])

# This macro checks for the libz library.  If found, it sets the automake
# conditional COIN_HAS_ZLIB and defines the C preprocessor variable
# COIN_HAS_ZLIB. The default is to use zlib, if it's present.  For each client,
# it adds the linker flags to the variable CLIENT_LFLAGS for each client.

AC_DEFUN([AC_COIN_CHK_ZLIB],
[
  AC_REQUIRE([AC_COIN_PROG_CC])

  m4_foreach_w([myvar],[$1],[AC_SUBST(m4_toupper(myvar)_LFLAGS)])

  coin_has_zlib=no

  AC_ARG_ENABLE([zlib],
    [AS_HELP_STRING([--disable-zlib],[do not compile with compression library zlib])],
    [coin_enable_zlib=$enableval],
    [coin_enable_zlib=yes])

  if test x$coin_enable_zlib = xyes ; then
    AC_CHECK_HEADER([zlib.h],[coin_has_zlib=yes])
    if test x$coin_has_zlib = xyes ; then
      AC_CHECK_LIB([z],[gzopen],[],[coin_has_zlib=no])
    fi
    if test x$coin_has_zlib = xyes ; then
      m4_foreach_w([myvar],[$1],
        [m4_toupper(myvar)_LFLAGS="-lz $m4_toupper(myvar)_LFLAGS"
        ])
      AC_DEFINE(m4_toupper(AC_PACKAGE_NAME)_HAS_ZLIB,[1],[Define to 1 if zlib is available])
    fi
  fi

  AM_CONDITIONAL([COIN_HAS_ZLIB],[test x$coin_has_zlib = xyes])
]) # AC_COIN_CHK_ZLIB


#                            COIN_CHK_BZLIB                             #

# COIN_CHK_BZLIB([client1 client2 ...])

# This macro checks for the libbz2 library.  If found, it defines the
# C preprocessor variable COIN_HAS_BZLIB and the automake conditional
# COIN_HAS_BZLIB.  Further, for a (space separated) list of clients, it adds
# the linker flag to the variable CLIENT_LFLAGS for each client.

AC_DEFUN([AC_COIN_CHK_BZLIB],
[
  AC_REQUIRE([AC_PROG_CC])

  m4_foreach_w([myvar],[$1],[AC_SUBST(m4_toupper(myvar)_LFLAGS)])

  coin_has_bzlib=no

  AC_ARG_ENABLE([bzlib],
    [AS_HELP_STRING([--disable-bzlib],[do not compile with compression library bzlib])],
    [coin_enable_bzlib=$enableval],
    [coin_enable_bzlib=yes])

  if test $coin_enable_bzlib = yes ; then
    AC_CHECK_HEADER([bzlib.h],[coin_has_bzlib=yes])
    if test $coin_has_bzlib = yes ; then
      AC_CHECK_LIB([bz2],[BZ2_bzReadOpen],[],[coin_has_bzlib=no])
    fi
    if test $coin_has_bzlib = yes ; then
      m4_foreach_w([myvar],[$1],
        [m4_toupper(myvar)_LFLAGS="-lbz2 $m4_toupper(myvar)_LFLAGS"
        ])
      AC_DEFINE(m4_toupper(AC_PACKAGE_NAME)_HAS_BZLIB,[1],[Define to 1 if bzlib is available])
    fi
  fi

  AM_CONDITIONAL([COIN_HAS_BZLIB],[test x$coin_has_bzlib = xyes])
]) # AC_COIN_CHK_BZLIB


#                              COIN_CHK_GMP                               #

# COIN_CHK_GMP([client1 client2 ...])

# This macro checks for the gmp library.  If found, it defines the C
# preprocessor variable COIN_HAS_GMP and the automake conditional COIN_HAS_GMP.
# Further, for a (space separated) list of clients, it adds the linker
# flag to the variable CLIENT_LFLAGS for each client.

AC_DEFUN([AC_COIN_CHK_GMP],
[
  AC_REQUIRE([AC_PROG_CC])

  m4_foreach_w([myvar],[$1],[AC_SUBST(m4_toupper(myvar)_LFLAGS)])

  coin_has_gmp=no

  AC_ARG_ENABLE([gmp],
    [AS_HELP_STRING([--disable-gmp],[do not compile with GNU multiple precision library])],
    [coin_enable_gmp=$enableval],
    [coin_enable_gmp=yes])

  if test $coin_enable_gmp = yes ; then
    AC_CHECK_HEADER([gmp.h],
      [AC_CHECK_LIB([gmp],[__gmpz_init],[coin_has_gmp=yes])])
    if test $coin_has_gmp = yes ; then
      m4_foreach_w([myvar],[$1],
        [m4_toupper(myvar)_LFLAGS="-lgmp $m4_toupper(myvar)_LFLAGS"
        ])
      AC_DEFINE(m4_toupper(AC_PACKAGE_NAME)_HAS_GMP,[1],[Define to 1 if GMP is available])
    fi
  fi

  AM_CONDITIONAL(COIN_HAS_GMP,[test x$coin_has_gmp = xyes])
]) # AC_COIN_CHK_GMP


#                         COIN_CHK_GNU_READLINE                         #

# COIN_CHK_GNU_READLINE([client1 client2 ...])

# This macro checks for GNU's readline.  It verifies that the header
# readline/readline.h is available, and that the -lreadline library contains
# "readline".  It is assumed that #include <stdio.h> is included in the source
# file before the #include<readline/readline.h> If found, it defines the C
# preprocessor variable COIN_HAS_READLINE.  Further, for a (space separated)
# list of clients, it adds the linker flag to the variable CLIENT_LFLAGS for
# each client.

AC_DEFUN([AC_COIN_CHK_GNU_READLINE],
[
  AC_REQUIRE([AC_PROG_CC])

  m4_foreach_w([myvar],[$1],[AC_SUBST(m4_toupper(myvar)_LFLAGS)])

  coin_has_readline=no

  AC_ARG_ENABLE([readline],
    [AS_HELP_STRING([--disable-readline],[do not compile with readline library])],
    [coin_enable_readline=$enableval],
    [coin_enable_readline=yes])

  if test $coin_enable_readline = yes ; then
    AC_CHECK_HEADER([readline/readline.h],
      [coin_has_readline=yes],[],[#include <stdio.h>])
    coin_save_LIBS="$LIBS"
    LIBS=
    if test $coin_has_readline = yes ; then
      AC_SEARCH_LIBS([tputs],
        [ncurses termcap curses],[],[coin_has_readline=no])
    fi
    if test $coin_has_readline = yes ; then
      AC_CHECK_LIB([readline],[readline],[],[coin_has_readline=no])
    fi
    if test $coin_has_readline = yes ; then
      m4_foreach_w([myvar],[$1],
        [m4_toupper(myvar)_LFLAGS="$LIBS $m4_toupper(myvar)_LFLAGS"
        ])
      AC_DEFINE(m4_toupper(AC_PACKAGE_NAME)_HAS_READLINE,[1],[Define to 1 if readline is available])
    fi
    LIBS="$coin_save_LIBS"
  fi

  AM_CONDITIONAL(COIN_HAS_READLINE,[test x$coin_has_readline = xyes])
]) # AC_COIN_CHK_GNU_READLINE



#                           COIN_DOXYGEN                                  #

# This macro determines the configuration information for doxygen, the tool
# used to generate online documentation of COIN code. It takes one parameter,
# a list of projects (mixed-case, to match the install directory names) that
# should be processed as external tag files. E.g., COIN_DOXYGEN([Clp Osi]).

# This macro will define the following variables:
#  coin_have_doxygen    Yes if doxygen is found, no otherwise
#  coin_doxy_usedot     Defaults to `yes'; --with-dot will still check to see
#                        if dot is available
#  coin_doxy_tagname    Name of doxygen tag file (placed in doxydoc directory)
#  coin_doxy_logname    Name of doxygen log file (placed in doxydoc directory)
#  coin_doxy_tagfiles   List of doxygen tag files used to reference other
#                       doxygen documentation

# It's not immediately obvious, but the code in this macro, configure-time
# substitions in doxygen.conf.in, and build-time edits of doxygen.conf in
# Makemain.inc combine to hardwire the assumptions that a tag file is named
# proj_doxy.tag, that PKG_TARNAME is coin-or-proj, and that the doxygen
# documentation is in the GNU default location $(docdir)/$PKG_TARNAME. Have
# a look over the complete machinery before you start changing things. The
# point of the build-time edits is to allow the user to redefine docdir at
# build time, as per GNU standards. Failure to use coin-or-proj as PKG_TARNAME
# will surely break linking of documentation with tag files.

AC_DEFUN([AC_COIN_DOXYGEN],
[
  AC_MSG_NOTICE([configuring doxygen documentation options])

  # Check to see if Doxygen and LaTeX are available.
  AC_CHECK_PROG([coin_have_doxygen],[doxygen],[yes],[no])
  AC_CHECK_PROG([coin_have_latex],[latex],[yes],[no])

  AC_ARG_WITH([dot],
    AS_HELP_STRING([--with-dot],[use dot (from graphviz) when creating documentation with doxygen if available; --without-dot to disable]),
    [],[withval=yes])

  # Look for the dot tool from the graphviz package, unless the user has disabled it.
  if test x"$withval" = xno ; then
    coin_doxy_usedot=NO
    AC_MSG_CHECKING([for dot ])
    AC_MSG_RESULT([disabled])
  else
    AC_CHECK_PROG([coin_doxy_usedot],[dot],[YES],[NO])
  fi

  # Generate a tag file name and a log file name.
  AC_SUBST([coin_doxy_tagname],[m4_tolower(AC_PACKAGE_NAME)_doxy.tag])
  AC_SUBST([coin_doxy_logname],[m4_tolower(AC_PACKAGE_NAME)_doxy.log])

  AM_CONDITIONAL(COIN_HAS_DOXYGEN, [test $coin_have_doxygen = yes])
  AM_CONDITIONAL(COIN_HAS_LATEX, [test $coin_have_latex = yes])

dnl Process the list of project names and massage each one into the name of
dnl a tag file. The value of coin_doxy_tagfiles is substituted for TAGFILES
dnl in doxygen.conf.in. Further substitution for @baredocdir_nosub@ will happen
dnl as an edit during make install. See comments in Makemain.inc.
  coin_doxy_tagfiles="m4_foreach_w([myvar],[$1],[@baredocdir_nosub@/coin-or-m4_tolower(myvar)/doxydoc/m4_tolower(myvar)_doxy.tag=@baredocdir_nosub@/coin-or-m4_tolower(myvar)/doxydoc/html ])"
  AC_SUBST([coin_doxy_tagfiles])
]) # AC_COIN_DOXYGEN


#                          COIN_EXAMPLE_FILES                             #

# This macro deals with `example' files. Its original use is for data files
# associated with the Data projects. It allows for uniform handling of files
# that are not processed but should be present in the build directory. Most of
# the work is done by make. This macro just sets up the necessary variables.
# The justification for this approach is to allow the use of file name wild
# cards in the argument. The underlying machinery of autoconf and automake
# kind of breaks down when the file names are not literal text.

# One tradeoff for allowing the use of wildcards is that it's not possible to
# give the created link a name that's different from the source file.

# The variables (file names accumulate over calls):
# EXAMPLE_DIST_FILES: The originals that should be distributed
# EXAMPLE_UNCOMP_FILES: like EXAMPLE_DIST_FILES, but names of compressed files
# are stripped of their .gz suffix.
# EXAMPLE_CLEAN_FILES: contains only the names of files that were stripped of
# .gz. Essential for distclean in a non-VPATH build.
# One of EXAMPLE_DIST_FILES or EXAMPLE_UNCOMP_FILES will be the list of file
# names for install / uninstall. See Makefile.am in a Data project for usage.

AC_DEFUN([AC_COIN_EXAMPLE_FILES],
[
  AC_REQUIRE([AC_COIN_CHECK_VPATH])[]dnl
  AC_REQUIRE([AC_PROG_LN_S])[]dnl
  AM_CONDITIONAL([COIN_VPATH_BUILD],[test x$coin_vpath_config = xyes])

  # Get the names of the files and sort them to the various variables. As
  # a convenient side-effect, the loop will remove newlines from the list
  # of files. Avoid adding duplicate file names. (Allowing wild cards makes
  # duplicates far too likely.)
  # The redirect of stderr on the ls is arguable. It avoids error messages
  # in the configure output if the package maintainer specifies files that
  # don't exist (for example, generic wild card specs that anticipate future
  # growth) but it will hide honest errors.
  files=`cd $srcdir; ls $1 2>/dev/null`
  for file in $files; do
    if expr " $EXAMPLE_DIST_FILES " : '.* '"$file"' .*' >/dev/null 2>&1 ; then
      continue ;
    fi
    EXAMPLE_DIST_FILES="$EXAMPLE_DIST_FILES $file"
    if expr "$file" : '.*\.gz$' >/dev/null 2>&1 ; then
      EXAMPLE_UNCOMP_FILES="$EXAMPLE_UNCOMP_FILES `expr "$file" : '\(.*\)\.gz$'`"
      EXAMPLE_CLEAN_FILES="$EXAMPLE_CLEAN_FILES `expr "$file" : '\(.*\)\.gz$'`"
    else
      EXAMPLE_UNCOMP_FILES="$EXAMPLE_UNCOMP_FILES $file"
    fi
  done

  AC_SUBST(EXAMPLE_DIST_FILES)
  AC_SUBST(EXAMPLE_UNCOMP_FILES)
  AC_SUBST(EXAMPLE_CLEAN_FILES)
]) # AC_COIN_EXAMPLE_FILES


#                           COIN_FINALIZE_FLAGS                           #

# COIN_FINALIZE_FLAGS([prim1 prim2 ...])

# For each PRIM, produce final versions of variables specifying linker and
# compiler flags.  PRIM_LFLAGS_NOPC and PRIM_CFLAGS_NOPC are appropriate for
# use in .pc files together with PRIM_PCFILES. PRIM_LFLAGS and PRIM_CFLAGS
# are augmented to contain flags obtained by invoking PKG_CONFIG on packages
# listed in PRIM_PCFILES and are appropriate for use in Makefile.am files.

AC_DEFUN([AC_COIN_FINALIZE_FLAGS],
[
  m4_foreach_w([myvar],[$1],
    [ if test "$BUILDTOOLS_DEBUG" = 1 ; then
        AC_MSG_NOTICE([FINALIZE_FLAGS for myvar:])
      fi
      AC_SUBST(m4_toupper(myvar)_LFLAGS_NOPC,[$m4_toupper(myvar)_LFLAGS])
      AC_SUBST(m4_toupper(myvar)_CFLAGS_NOPC,[$m4_toupper(myvar)_CFLAGS])
      if test -n "${m4_toupper(myvar)_PCFILES}" ; then
        temp_CFLAGS=`PKG_CONFIG_PATH="$COIN_PKG_CONFIG_PATH" $PKG_CONFIG --cflags ${m4_toupper(myvar)_PCFILES}`
        temp_LFLAGS=`PKG_CONFIG_PATH="$COIN_PKG_CONFIG_PATH" $PKG_CONFIG --libs $pkg_static ${m4_toupper(myvar)_PCFILES}`
        m4_toupper(myvar)_CFLAGS="$temp_CFLAGS ${m4_toupper(myvar)_CFLAGS}"
        m4_toupper(myvar)_LFLAGS="$temp_LFLAGS ${m4_toupper(myvar)_LFLAGS}"
      fi

      # setup XYZ_EXPORT symbol for library users
      libexport_attribute=
      if test "$enable_shared" = yes ; then
        case $build_os in
          cygwin* | mingw* | msys* | cegcc* )
            libexport_attribute="__declspec(dllimport)"
            if test "$enable_static" = yes ; then
              AC_MSG_ERROR([Cannot do DLL and static LIB builds simultaneously. Do not add --enable-static without --disable-shared.])
            fi
          ;;
        esac
      fi
      AC_DEFINE_UNQUOTED(m4_toupper(myvar)_EXPORT, [$libexport_attribute], [Library Visibility Attribute])

      # add -DXYZ_BUILD to XYZ_CFLAGS
      m4_toupper(myvar)_CFLAGS="${m4_toupper(myvar)_CFLAGS} -D[]m4_toupper(myvar)_BUILD"

      # Define BUILDTOOLS_DEBUG to enable debugging output
      if test "$BUILDTOOLS_DEBUG" = 1 ; then
        AC_MSG_NOTICE([m4_toupper(myvar)_LFLAGS_NOPC: "${m4_toupper(myvar)_LFLAGS_NOPC}"])
        AC_MSG_NOTICE([m4_toupper(myvar)_CFLAGS_NOPC: "${m4_toupper(myvar)_CFLAGS_NOPC}"])
        AC_MSG_NOTICE([adding "${m4_toupper(myvar)_PCFILES}"])
        AC_MSG_NOTICE([m4_toupper(myvar)_LFLAGS: "${m4_toupper(myvar)_LFLAGS}"])
        AC_MSG_NOTICE([m4_toupper(myvar)_CFLAGS: "${m4_toupper(myvar)_CFLAGS}"])
      fi
    ])
])


#                              COIN_FINALIZE                              #

# This macro should be called at the very end of the configure.ac file.
# It creates the output files (by using AC_OUTPUT), and might do some other
# things (such as generating links to data files in a VPATH configuration).
# It also prints the "success" message.
# Note: If the variable coin_skip_ac_output is set to yes, then no output
# files are written.

AC_DEFUN([AC_COIN_FINALIZE],
[
  AC_OUTPUT

dnl AC_MSG_NOTICE([In case of trouble, first consult the documentation at https://coin-or.github.io/user_introduction])
  AC_MSG_NOTICE([Configuration of $PACKAGE_NAME successful])
])


#                       DEPRECATE_BUILD_OPTION                            #

# COIN_DEPRECATE_BUILD_OPTION([macro name],[dfltaction])

# Macro to encapsulate the warning / fatal error for the deprecated 'build'
# option. From the LIBHDR chain this is fatal. From the PKG chain, it's a
# warning.

AC_DEFUN([AC_COIN_DEPRECATE_BUILD_OPTION],
[
  m4_bmatch([ $2 ],[.*build .*],
    [m4_bmatch([$1],
      [.*_LIB.*],
      [m4_fatal(m4_newline([m4_text_wrap(
         [$1:]
      [The 'default_build' or 'build' option is deprecated and not]
      [supported by COIN_CHK_LIB or COIN_CHK_LIBHDR.]
      [More generally, it is expected that COIN Third-Party packages]
      [will be phased out over time. Plan accordingly.]
      [If you really want to use 'default_build', use COIN_CHK_PKG.]
      [All COIN-OR ThirdParty packages produce .pc files.],[**],[**])]))],
      [.*_PKG.*],
      m4_errprintn(m4_text_wrap(
        [$1:]
        [The 'default_build' or 'build' option is deprecated.]
        [More generally, it is expected that COIN-OR ThirdParty packages will be]
        [phased out over time. Plan accordingly.],[**],[**])))])
])dnl         # COIN_DEPRECATE_BUILD_OPTION

#                         LIBHDR_DFLT_HACK                                #

# COIN_LIBHDR_DFLT_HACK([which],[dfltaction],[dfltuse],[dfltlink])

# A pure macro hack to get around POSIX m4's nine-parameter limit by
# overloading the [dfltaction] parameter to the PRIM_LIBHDR, CHK_LIBHDR
# and CHK_LIB macros. It's ugly but it's pure macro and will resolve at
# run_autotools.

# We're trying to control two distinct things:
# * should the library be used by default (default_use or default_skip)
#   or built as a third-party package (default_build)
# * should the link check use only the function (separate) or both the
#   includes and the function (together).

# This macro looks at the possibilities and returns the appropriate
# string. [which] ($1) controls whether it returns library usage ('usage')
# or the form of the link check ('link'). [dfltaction] ($2) is the parameter
# to be analysed. [dfltuse] ($3) specifies the default for usage, in case
# [dfltaction] doesn't specify it. [dfltlink] ($4) covers the default for
# the link check, in case [dfltaction] doesn't specify it.

# Specify [dfltaction] as a space-separated string of keywords, e.g.,
# [default_use separate]. m4_normalize deals with excess leading, interior,
# trailing spaces.

# The macro structure is expandable (new overloads can be added) and order
# insensitive provided all keywords are unique.

# It's critical that the macro produce a result with no leading or trailing
# spaces! Be careful if you modify it. In normal use, you want this macro to
# expand immediately, so don't quote it.

AC_DEFUN([AC_COIN_LIBHDR_DFLT_HACK],
[m4_case([$1],
   [usage],
     m4_bmatch([m4_normalize([ $2 ])],
                             [.* default_use .*],[default_use],
                             [.* default_skip .*],[default_skip],
                             [.* default_build .*],[default_build],
                             [$3]),
   [link],
     m4_bmatch([m4_normalize([ $2 ])],
                             [.* separate .*],[separate],
                             [.* together .*],[together],
                             [$4]))])dnl
# LIBHDR_DFLT_HACK


#                         FIND_PRIM_LIBHDR                                #

# COIN_FIND_PRIM_LIBHDR([prim],[lflgs],[cflgs],[dflgs],
#                       [function-body],[includes],
#                       [dfltaction],[cmdlineopts])

# Determine whether we can use primary library prim ($1) and assemble
# information on the required linker flags (prim_lflags), compiler flags
# (prim_cflags), and data directories (prim_data) as specified by cmdlineopts
# ($8).

# A compile check will be performed if [includes] ($6) is specified.  A link
# check will be performed if [function-body] ($5) is specified. The default is
# to use a program composed by concatenating [includes] and [function-body]
# with link flags specified by [lflgs] ($2), and compiler flags specified by
# [cflgs] ($3) (but note these can be overridden; see below). If [dfltaction]
# ($7) contains 'separate', the link check will be performed without the header
# (see LIBHDR_DFLT_HACK).

# There's really no reasonable assumptions that can be made about how a data
# directory is specified, so there's no attempt to test for existence.

# [cmdlineopts] ($8) specifies the set of configure command line options
# processed: 'nodata' processes --with-prim, --with-prim-lflags,
# and --with-prim-cflags; 'dataonly' processes only --with-prim and
# --with-prim-data; anything else ('all' works well) processes all four
# command line options. Shell code produced by the macro is tailored based
# on [cmdlineopts]. `nodata' is the default.

# --with-prim is interpreted as follows:
#   * --with-prim=no or --without-prim
#     prim status is set to skipping
#   * --with-prim or --with-prim=yes
#     prim status is set to requested
#   * Any other value is taken as equivalent to
#       --with-prim-data=value (dataonly) or
#       --with-prim-lflags=value (anything else)
#     prim status is set to requested

# If dfltaction contains default_skip, an explicit --with-prim or non-null
# --with-prim-{lflags,cflags} is required to override the default.

# The algorithm first checks for a user-specified value of --with-prim; values
# are interpreted as above.  Next, it looks for user specified values given
# with command line parameters --with-prim-lflags, --with-prim-cflags, and
# --with-prim-data. If any of these are specified, the value overrides the
# value passed to the macro as a parameter.

# The usage portion of [dfltaction] ($7) (default_skip, default_use) is
# used as the default value of --with-prim if the user offers no guidance
# via command line parameters. The (hardwired) default is default_use.
# See LIBHDR_DFLT_HACK for the full story.

# If you really wanted to use the obsolete 'default_build', use COIN_CHK_PKG
# instead. All COIN ThirdParty packages produce .pc files. This macro doesn't
# deal with it.

# This macro doesn't test that the specified values actually work unless
# [function-body] and/or [includes] are given as parameters. This
# is deliberate.  There's no guarantee that the specified library can be
# accessed just yet with the specified flags. Put another way, unless the
# user requests a compile and/or link check, all we're doing here is filling
# in variables using a complicated algorithm.

AC_DEFUN([AC_COIN_FIND_PRIM_LIBHDR],
[
dnl Trap default_build right up front and fail in run_autotools.
  AC_COIN_DEPRECATE_BUILD_OPTION([COIN_FIND_PRIM_LIBHDR],[$8])

dnl Set default values for flags, action, and status. These are taken from
dnl the macro parameters.

  m4_tolower($1_lflags)="$2"
  m4_tolower($1_cflags)="$3"
  m4_tolower($1_data)="$4"
  m4_if(AC_COIN_LIBHDR_DFLT_HACK([usage],[$7],[default_use]),default_use,
    [m4_tolower(coin_has_$1)=requested],
    [m4_tolower(coin_has_$1)=skipping
     m4_tolower($1_failmode)='default'])
  m4_tolower($1_userflags)='no'

dnl See if the user specified --with-prim.  If the value is something
dnl other than 'yes' or 'no' and the client specified dataonly, the value is
dnl assigned to prim_data, otherwise to prim_lflags. But in the absence of an
dnl explicit 'no', allow non-null lflags, cflags, or data to be equivalent to
dnl 'yes'.

  withval="$m4_tolower(with_$1)"
  if test -n "$withval" ; then
    case "$withval" in
      no )
        m4_tolower(coin_has_$1)=skipping
        m4_tolower($1_failmode)='command line'
        ;;
      yes )
        m4_tolower(coin_has_$1)=requested
        m4_tolower($1_failmode)=''
        ;;
      * )
        m4_tolower(coin_has_$1)=requested
        m4_tolower($1_failmode)=''
        m4_if(m4_default($8,nodata),dataonly,
          [m4_tolower($1_data)="$withval"],
          [m4_tolower($1_lflags)="$withval"])
        m4_tolower($1_userflags)='yes'
        ;;
    esac
  else
    if test -n "$m4_tolower(with_$1_lflags)" ||
       test -n "$m4_tolower(with_$1_cflags)" ||
       test -n "$m4_tolower(with_$1_data)" ; then
      m4_tolower(coin_has_$1)=requested
      m4_tolower($1_failmode)=''
      m4_tolower($1_userflags)='yes'
    fi
  fi

dnl As long as we are not dataonly and we are not skipping prim, check for
dnl --with-prim-lflags and --with-prim-cflags. Values will override parameter
dnl values.

  m4_if(m4_default($8,nodata),dataonly,[],
    [if test "$m4_tolower(coin_has_$1)" != skipping ; then
       withval="$m4_tolower(with_$1_lflags)"
       if test -n "$withval" ; then
         m4_tolower(coin_has_$1)=requested
         m4_tolower($1_lflags)="$withval"
         m4_tolower($1_userflags)='yes'
       fi

       withval="$m4_tolower(with_$1_cflags)"
       if test -n "$withval" ; then
         m4_tolower(coin_has_$1)=requested
         m4_tolower($1_cflags)="$withval"
         m4_tolower($1_userflags)='yes'
       fi
     fi])

dnl If we are not nodata and we are not skipping prim, check for
dnl --with-prim-data. A value will override the parameter value.

  m4_if(m4_default($8,nodata),nodata,[],
    [if test "$m4_tolower(coin_has_$1)" != skipping ; then
       withval="$m4_tolower(with_$1_data)"
       if test -n "$withval" ; then
         m4_tolower(coin_has_$1)=requested
         m4_tolower($1_data)="$withval"
         m4_tolower($1_userflags)='yes'
       fi
     fi])

dnl At this point, coin_has_prim can be one of skipping (user said no,
dnl or default was skip without override), or requested (user said yes,
dnl or default was use without override).

dnl If we have [includes], try to compile them.

  m4_ifnblank([$6],
    [if test $m4_tolower(coin_has_$1) != skipping ; then
       ac_save_CPPFLAGS=$CPPFLAGS
       CPPFLAGS="$m4_tolower($1_cflags)"
       AC_COMPILE_IFELSE([AC_LANG_PROGRAM([$6],[])],
         [],
         [m4_tolower(coin_has_$1)='no'
          m4_tolower($1_failmode)="header compile"])
       CPPFLAGS=$ac_save_CPPFLAGS
     fi])

dnl If we have a function-body, try to compile and link. Use both the
dnl [includes] and [function-body] unless the user has requested
dnl otherwise by '*separate' as [dfltaction].

  m4_ifnblank([$5],
    [if test "$m4_tolower(coin_has_$1)" != skipping ; then
       ac_save_LIBS=$LIBS
       ac_save_CPPFLAGS=$CPPFLAGS
       LIBS="$m4_tolower($1_lflags)"
       CPPFLAGS="$m4_tolower($1_cflags)"
       m4_if(
         AC_COIN_LIBHDR_DFLT_HACK([link],[$7],[unused],[together]),[separate],
         [AC_LINK_IFELSE([AC_LANG_SOURCE([$5])],
           [],
           [m4_tolower(coin_has_$1)='no'
            if test -n "$m4_tolower($1_failmode)" ; then
              m4_tolower($1_failmode)="$m4_tolower($1_failmode), bare link"
            else
              m4_tolower($1_failmode)="bare link"
            fi])],
         [AC_LINK_IFELSE([AC_LANG_PROGRAM([$6],[$5])],
           [],
           [m4_tolower(coin_has_$1)='no'
            if test -n "$m4_tolower($1_failmode)" ; then
              m4_tolower($1_failmode)="$m4_tolower($1_failmode), link with header"
            else
              m4_tolower($1_failmode)="link with header"
            fi])])
       LIBS=$ac_save_LIBS
       CPPFLAGS=$ac_save_CPPFLAGS
     fi])

dnl If we are still showing requested, then we can say yes. We have passed all
dnl the tests requested by the user (which might be none, but that is on the
dnl head of the user).

  if test $m4_tolower(coin_has_$1) = requested ; then
    m4_tolower(coin_has_$1)=yes
  fi

dnl The final value of coin_has_prim will be yes, no, or skipping. Skipping
dnl means that we defaulted to no or the user said no. No means that some
dnl check failed. Yes means that we defaulted to yes or the user made a
dnl specific request, and no check failed.

  # Define BUILDTOOLS_DEBUG to enable debugging output
  if test "$BUILDTOOLS_DEBUG" = 1 ; then
    AC_MSG_NOTICE([FIND_PRIM_LIB result for $1: "$m4_tolower(coin_has_$1)"])
    AC_MSG_NOTICE([Collected values for package '$1'])
    AC_MSG_NOTICE([m4_tolower($1_lflags) is "$m4_tolower($1_lflags)"])
    AC_MSG_NOTICE([m4_tolower($1_cflags) is "$m4_tolower($1_cflags)"])
    AC_MSG_NOTICE([m4_tolower($1_data) is "$m4_tolower($1_data)"])
    AC_MSG_NOTICE([m4_tolower($1_pcfiles) is "$m4_tolower($1_pcfiles)"])
  fi
])  # COIN_FIND_PRIM_LIBHDR


#                         COIN_CHK_LIBHDR                                 #

# COIN_CHK_LIBHDR([prim],[clients],[lflgs],[cflgs],[dflgs],
#              [function-body],[includes],
#              [dfltaction],[cmdopts])

# Determine whether we can use primary library [prim] ($1) and assemble
# information on the required linker flags (prim_lflags), compiler flags
# (prim_cflags), and data directories (prim_data). The call to DEF_PRIM_ARGS
# takes care of defining configure command line parameters for [prim]. The
# call to FIND_PRIM_LIBHDR sets variables for [prim] and does compile and
# link checks if requested. This macro controls the flow and sets up variables
# for use during the build.

# A compile check will be performed if [includes] ($7) is specified.  A link
# check will be performed if [function-body] ($6) is specified. The default is
# to use a program composed by concatenating [includes] and [function-body]
# with link flags specified by [lflgs] ($3) and compiler flags specified by
# [cflgs] ($4). If you need to override this, have a look at LIBHDR_DFLT_HACK
# and set [dfltaction] accordingly.

# The configure command line options offered to the user are controlled
# by [cmdopts] ($9). 'nodata' offers --with-prim, --with-prim-lflags, and
# --with-prim-cflags; 'dataonly' offers --with-prim and --with-prim-data;
# 'all' offers all four. Anything else defaults to 'nodata'. DEF_PRIM_ARGS and
# FIND_PRIM_LIB are tailored accordingly. The (hardwired) default is 'nodata'.

# Macro parameters [lflgs] ($3), [cflgs] ($4), and [dflgs] ($5) are used
# for --with-prim-lflags, --with-prim-cflags, and --with-prim-data if and
# only if there are no user-supplied values on the command line. A command
# line value will override the parameter value.

# [dfltaction] ($8) (default_skip, default_use) is used as the default for
# library usage if the user offers no guidance via command line parameters. The
# (hardwired) default is default_use.

# If you really wanted to use the obsolete 'default_build', use COIN_CHK_PKG
# instead. All COIN ThirdParty packages produce .pc files.

# Define an automake conditional COIN_HAS_PRIM to record the result. If we
# decide to use prim, also define a preprocessor symbol COIN_HAS_PRIM.

# Linker and compiler flag information will be propagated to the space-
# separated list of client packages [clients] ($2) using the _LFLAGS and
# _CFLAGS variables of the clients. These variables match Libs.private and
# Cflags.private, respectively, in a .pc file.

# Data directory information is used differently. Typically what's wanted is
# individual variables specifying the data directory for each primitive. Hence
# the macro defines PRIM_DATA for the primitive.

AC_DEFUN([AC_COIN_CHK_LIBHDR],
[ 
dnl Trap build or default_build right up front and fail in run_autotools.
  AC_COIN_DEPRECATE_BUILD_OPTION([COIN_CHK_LIBHDR],[$8])

  AC_MSG_CHECKING(
    m4_ifnblank([$6],
      m4_ifnblank([$7],
        m4_normalize([for library $1 with
        m4_if(
          AC_COIN_LIBHDR_DFLT_HACK([link],[$8],[unused],[together]),
          [together],
          [combined link and compile check],
          [separate link and compile checks])]),
        [for library $1 with link check]),
      m4_ifnblank([$7],[for library $1 with compile check],
                       [for library $1 (setting flags only)])))

dnl Make sure the necessary variables exist for each client package.

  m4_foreach_w([myvar],[$2],
    [AC_SUBST(m4_toupper(myvar)_LFLAGS)
     AC_SUBST(m4_toupper(myvar)_CFLAGS)
    ])

dnl Check to see if the user has overridden configure parameters from the
dnl environment. Caution! Quadrigraphs.

  m4_tolower(coin_has_$1)=noInfo
  if test x"$COIN_SKIP_PROJECTS" != x ; then
    for pkg in `echo $COIN_SKIP_PROJECTS | tr '@<:@:upper:@:>@' '@<:@:lower:@:>@'` ; do
      if test "$pkg" = "m4_tolower($1)" ; then
        m4_tolower(coin_has_$1)=skipping
      fi
    done
  fi

dnl If we are not skipping this project, define and process the command line
dnl options according to the cmdopts parameter. Then invoke FIND_PRIM_PKG
dnl to do the heavy lifting.

  if test "$m4_tolower(coin_has_$1)" != skipping ; then
    m4_case(m4_default($9,nodata),
      nodata,  [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,no,
                  AC_COIN_LIBHDR_DFLT_HACK([usage],[$8],[default_use]))],
      dataonly,[AC_COIN_DEF_PRIM_ARGS([$1],yes,no,no,yes,
                  AC_COIN_LIBHDR_DFLT_HACK([usage],[$8],[default_use]))],
               [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,yes,
                  AC_COIN_LIBHDR_DFLT_HACK([usage],[$8],[default_use]))],
      all,     [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,yes,
                  AC_COIN_LIBHDR_DFLT_HACK([usage],[$8],[default_use]))],
      [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,no,
         AC_COIN_LIBHDR_DFLT_HACK([usage],[$8],[default_use]))])
    AC_COIN_FIND_PRIM_LIBHDR(m4_tolower($1),
      [$3],[$4],[$5],[$6],[$7],
      m4_if([$8],[build],[default_build],[$8]),m4_default([$9],nodata))
    if test -n "$m4_tolower($1_failmode)" ; then
      AC_MSG_RESULT([$m4_tolower(coin_has_$1) ($m4_tolower($1_failmode))])
    else
      AC_MSG_RESULT([$m4_tolower(coin_has_$1)])
    fi
  else
    AC_MSG_RESULT([$m4_tolower(coin_has_$1) (COIN_SKIP_PROJECTS)])
  fi

dnl Possible outcomes are 'yes', 'no', or 'skipping'. 'Skipping' implies we
dnl decided to skip the package for some reason. 'No' means we wanted the
dnl package but failed a test. 'Yes' means we wanted the package and did not
dnl fail any tests.

dnl Try to offer some helpful advice in the event of failure.
  if test "$m4_tolower(coin_has_$1)" = 'no' ; then
    if expr "$m4_tolower($1_failmode)" : '.*header.*' &>/dev/null ; then
      AC_MSG_NOTICE([Compiler flags were "$m4_tolower($1)_cflags". Use --with-m4_tolower($1)-cflags to overwrite. Check config.log for details of failed compile attempt.])
    fi
    if expr "$m4_tolower($1_failmode)" : '.*link.*' &>/dev/null ; then
      AC_MSG_NOTICE([Linker flags are "$m4_tolower($1)_lflags". Use --with-m4_tolower($1)-lflags to overwrite. Check config.log for details of failed link attempt.])
    fi
    if test "$m4_tolower($1_userflags)" = 'yes' ; then
      AC_MSG_ERROR([user-specified flags for $1 do not work.])
    fi
  fi

dnl Normalise to yes or no for the remainder.

  if test "$m4_tolower(coin_has_$1)" != yes ; then
    m4_tolower(coin_has_$1)=no
  fi

dnl Create an automake conditional COIN_HAS_PRIM.

  AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),[test $m4_tolower(coin_has_$1) = yes])

dnl If we have located the package, define preprocessor symbol PACKAGE_HAS_PRIM
dnl and augment the necessary variables for the client packages.

  if test $m4_tolower(coin_has_$1) = yes ; then
    AC_DEFINE(m4_toupper(AC_PACKAGE_NAME)_HAS_[]m4_toupper($1),[1],
      [Define to 1 if the $1 package is available])
    m4_foreach_w([myvar],[$2],
      [m4_toupper(myvar)_LFLAGS="$m4_tolower($1_lflags) $m4_toupper(myvar)_LFLAGS"
       m4_toupper(myvar)_CFLAGS="$m4_tolower($1_cflags) $m4_toupper(myvar)_CFLAGS"
      ])

dnl Finally, set up PRIM_DATA, unless the user specified nodata.

    m4_if(m4_default([$9],nodata),nodata,[],
      [AC_SUBST(m4_toupper($1)_DATA)
       m4_toupper($1)_DATA=$m4_tolower($1_data)])
  fi

])dnl   # COIN_CHK_LIBHDR


#                          COIN_CHK_LIB                                   #

# COIN_CHK_LIB([prim],[clients],[lflgs],[cflgs],[dflgs],
#              [func],[header],[dfltaction],[cmdopts])

# This is a wrapper for COIN_CHK_LIBHDR (which see) that makes two simplifying
# assumptions:
#  * func is a function name and can be redeclared as
#      extern 'C' void func()
#    for the purpose of a link check
#  * header is a file name and can be used as
#      #include "header"
#    for the purpose of a compile check.
# It goes without saying that the the header file should not be included as
# part of the link check, lest the proper declaration of func in the header
# conflict with the simplified declaration used for the link check. If you
# want the real thing, use COIN_CHK_LIBHDR directly.

# The call to LIBHDR_DFLT_HACK picks up usage (yes or no) from configure.ac,
# then we force separate compilation.
# The calls to m4_ifnblank must expand here --- we want to pass the result
# of expansion to CHK_LIBHDR. Hence the double quoting for the code strings.

AC_DEFUN([AC_COIN_CHK_LIB],
[ 
dnl Trap build or default_build right up front and fail in run_autotools.
  AC_COIN_DEPRECATE_BUILD_OPTION([COIN_CHK_LIB],[$8])

dnl Fabricate the call to COIN_LIBHDR
  AC_COIN_CHK_LIBHDR([$1],[$2],[$3],[$4],[$5],
    m4_ifnblank([$6],
      [[
#ifdef __cplusplus
  extern "C"
#endif
  void $6() ;
  int main () { $6() ; return (0) ; }]],[]),
    m4_ifnblank([$7],[[#include "$7"]],[]),
    AC_COIN_LIBHDR_DFLT_HACK([usage],
      m4_if([$8],[build],[default_build],[$8]),[default_use])[ separate],
    [$9])
])dnl   # COIN_CHK_LIB



#                     COIN_FIND_PRIM_PKG                                  #
# COIN_FIND_PRIM_PKG([prim],[pcfile(s)],[default action],[cmdlineopts])

# Determine whether we can use primary package prim ($1) by passing [pcfile(s)]
# to pkgconf. From the result, assemble information on the required pkgconfig
# (.pc) files, library flags (prim_lflags), compiler flags (prim_cflags),
# and data directories (prim_data).

# cmdlineopts ($4) specifies the set of configure command line options
# defined and processed: 'nodata' produces --with-prim, --with-prim-libs, and
# --with-prim-cflags; 'dataonly' produces --with-prim and --with-prim-data;
# 'all' produces all four. Anything else defaults to 'nodata'.  Shell code
# produced by the macro is tailored based on cmdlineopts.

# --with-prim is interpreted as follows:
#   * --with-prim=no or --without-prim
#     prim status is set to skipping
#   * --with-prim or --with-prim=yes
#     prim status is set to requested
#   * --with-prim=build
#     prim status is set to requested but we look first for coinxxx.pc,
#     then xxx.pc
#   * Any other value is taken as equivalent to
#       --with-prim-data=value (dataonly) or
#       --with-prim-lflags=value (anything else)
#     prim status is set to requested

# The algorithm first checks for a user-specified value of --with-prim;
# values are interpreted as above.  Next, it looks for user specified values
# given with command line parameters --with-prim-lflags, --with-prim-cflags,
# and --with-prim-data. If any of these are specified, they override the .pc
# file(s) specified in [pcfiles(s)] and are accepted with no futher testing.
# Be careful!

# The usage portion of [dfltaction] ($3) (default_skip, default_use,
# default_build) is used as the default value of --with-prim if the user
# offers no guidance via command line parameters. Since FIND_PRIM_PKG does
# no checks, it's useless to specify linkage in [dfltaction].

# To repeat the warning with CHK_PKG, `default_build' ('build' from the
# command line) hasn't a hope of working for anything except a COIN
# Third-party package, where the convention is to name the .pc file
# coinprim.pc. The code below will take whatever is in [pcfile(s)] and prepend
# the string 'coin'. This will almost surely fail to do what the user expects
# if [pcfile(s)] actually contains more than one .pc file name.

# The macro doesn't test that the specified values actually work. This is
# deliberate.  There's no guarantee that user-specified libraries and/or
# directories actually exist yet. The same possibility exists for values
# returned when pkgconf checks the .pc file.

AC_DEFUN([AC_COIN_FIND_PRIM_PKG],
[
  AC_REQUIRE([AC_COIN_HAS_PKGCONFIG])

dnl Issue a warning that 'default_build' is deprecated.
  AC_COIN_DEPRECATE_BUILD_OPTION([COIN_FIND_PRIM_PKG],[$3])

dnl Set the default action and status from the macro parameters.

  m4_case(AC_COIN_LIBHDR_DFLT_HACK([usage],[$3],[default_use]),
    default_use,
      [m4_tolower(coin_has_$1)=requested],
    default_build,
      [m4_tolower(coin_has_$1)=requested
       m4_tolower($1_build)=yes],
    default_skip,
      [m4_tolower(coin_has_$1)=skipping
       m4_tolower($1_failmode)='default'],
    [m4_tolower(coin_has_$1)=requested])

dnl Set default values for flags, data, pcfiles
  m4_tolower($1_lflags)=
  m4_tolower($1_cflags)=
  m4_tolower($1_data)=
  m4_ifnblank([$2],
    [m4_tolower($1_pcfiles)="$2"],
    [m4_tolower($1_pcfiles)=m4_tolower($1)])

dnl See if the user specified --with-prim.  If the value is something other
dnl than 'yes', 'build', or 'no' and the client specified dataonly, the
dnl value is assigned to prim_data, otherwise to prim_lflags. But in the
dnl absence of an explicit 'no', allow non-null lflags, cflags, or data to
dnl be equivalent to 'yes'.

  withval="$m4_tolower(with_$1)"
  if test -n "$withval" ; then
    case "$withval" in
      no )
        m4_tolower(coin_has_$1)=skipping
        m4_tolower($1_failmode)='command line'
        ;;
      yes )
        m4_tolower(coin_has_$1)=requested
        m4_tolower($1_failmode)=''
        ;;
      build )
        m4_tolower(coin_has_$1)=requested
        m4_tolower($1_build)=yes
        m4_tolower($1_failmode)=''
        ;;
      * )
        m4_tolower(coin_has_$1)=requested
        m4_tolower($1_failmode)=''
        m4_if(m4_default($8,nodata),dataonly,
          [m4_tolower($1_data)="$withval"],
          [m4_tolower($1_lflags)="$withval"])
        m4_tolower($1_pcfiles)=''
        ;;
    esac
  else
    if test -n "$m4_tolower(with_$1_lflags)" ||
       test -n "$m4_tolower(with_$1_cflags)" ||
       test -n "$m4_tolower(with_$1_data)" ; then
      m4_tolower(coin_has_$1)=requested
      m4_tolower($1_failmode)=''
      m4_tolower($1_pcfiles)=''
    fi
  fi

dnl As long as we're not dataonly and we're not skipping prim, check for
dnl --with-prim-lflags and --with-prim-cflags. If we see any of no, yes,
dnl build, declare them invalid. Otherwise, we're not going to check so
dnl accept them on faith.

  m4_if(m4_default($4,nodata),dataonly,[],
    [if test "$m4_tolower(coin_has_$1)" != skipping ; then
       withval="$m4_tolower(with_$1_lflags)"
       if test -n "$withval" ; then
         case "$withval" in
           build | no | yes )
             AC_MSG_ERROR(m4_normalize(
               ["$withval" is not useful here; please specify link flags]
               [appropriate for your environment.]))
             ;;
           * )
             m4_tolower(coin_has_$1)=yes
             m4_tolower($1_lflags)="$withval"
             m4_tolower($1_pcfiles)=''
             ;;
         esac
       fi

       withval="$m4_tolower(with_$1_cflags)"
       if test -n "$withval" ; then
         case "$withval" in
           build | no | yes )
             AC_MSG_ERROR(m4_normalize(
               ["$withval" is not useful here; please specify compile flags]
               [appropriate for your environment.]))
             ;;
           * )
             m4_tolower(coin_has_$1)=yes
             m4_tolower($1_cflags)="$withval"
             m4_tolower($1_pcfiles)=''
             ;;
         esac
       fi
     fi])

dnl If we're not nodata and we're not skipping prim, check for
dnl --with-prim-data. A value will override the parameter value. Here we
dnl really can't do anything except accept on faith.

  m4_if(m4_default($4,nodata),nodata,[],
    [if test "$m4_tolower(coin_has_$1)" != skipping ; then
       withval="$m4_tolower(with_$1_data)"
       if test -n "$withval" ; then
         m4_tolower(coin_has_$1)=yes
         m4_tolower($1_data)="$withval"
         m4_tolower($1_pcfiles)=''
       fi
     fi])

dnl At this point, coin_has_prim can be one of skipping (user said no, or
dnl default was no without override), requested (user said yes or build,
dnl or default was yes or build without override), or yes (user provided
dnl --with-prim-lflags, --with-prim-cflags, or --with-prim-data on the
dnl command line).  If we're still looking, but pkg-config is not available,
dnl print a big warning and say the package isn't available.  Otherwise, look
dnl for the specified .pc files. If the user has requested build, prepend
dnl 'coin' to the pcfiles string and try that first. This will do the right
dnl thing if the .pc file name is the first word in the pcfile string

  if test $m4_tolower(coin_has_$1) = requested ; then
    if test -n "$PKG_CONFIG" ; then
      if test x"$m4_tolower($1_build)" = xyes ; then
        pcfile="coin$m4_tolower($1)_pcfiles"
        AC_COIN_CHK_MOD_EXISTS([$1],[$pcfile],
         [m4_tolower(coin_has_$1)=yes],
         [m4_tolower(coin_has_$1)=no])
      fi
      if ! test $m4_tolower(coin_has_$1) = yes ; then
        pcfile="$m4_tolower($1)_pcfiles"
        AC_COIN_CHK_MOD_EXISTS([$1],[$pcfile],
          [m4_tolower(coin_has_$1)=yes],
          [m4_tolower(coin_has_$1)=no])
      fi
      if test $m4_tolower(coin_has_$1) = yes ; then
        m4_tolower($1_data)=`PKG_CONFIG_PATH="$COIN_PKG_CONFIG_PATH" $PKG_CONFIG --variable=datadir "$pcfile" 2>/dev/null`
        m4_tolower($1_pcfiles)="$pcfile"
      fi
    else
      AC_MSG_WARN(m4_normalize(
        [Check for $1 via pkg-config could not be performed as there is no pkg-config available.]
        [Consider installing pkg-config or provide appropriate values for --with-m4_tolower($1)-lflags and --with-m4_tolower($1)-cflags.]))
      m4_tolower(coin_has_$1)=no
    fi
  fi

dnl The final value of coin_has_prim will be yes, no, or skipping.  No means
dnl we looked (with pkgconfig) and didn't find anything.  Skipping means
dnl the user said no or no was the default and not overridden.  Yes means we
dnl have something, from the user or from pkgconfig.  Note that we haven't
dnl run a useability test!

dnl Define BUILDTOOLS_DEBUG to enable debugging output
  if test "$BUILDTOOLS_DEBUG" = 1 ; then
    AC_MSG_NOTICE([FIND_PRIM_PKG result for $1: "$m4_tolower(coin_has_$1)"])
    AC_MSG_NOTICE([Collected values for package '$1'])
    AC_MSG_NOTICE([m4_tolower($1_lflags) is "$m4_tolower($1_lflags)"])
    AC_MSG_NOTICE([m4_tolower($1_cflags) is "$m4_tolower($1_cflags)"])
    AC_MSG_NOTICE([m4_tolower($1_data) is "$m4_tolower($1_data)"])
    AC_MSG_NOTICE([m4_tolower($1_pcfiles) is "$m4_tolower($1_pcfiles)"])
  fi

])  # COIN_FIND_PRIM_PKG


#                          COIN_CHK_PKG                                   #

# COIN_CHK_PKG([prim],[clients],[.pc file(s)],
#              [dfltaction],[cmdopts])

# Determine whether we can use primary package prim ($1) and assemble
# information on the required pkgconfig files (prim_pcfiles), linker
# flags (prim_lflags), compiler flags (prim_cflags), and data directories
# (prim_data).

# The configure command line options offered to the user are controlled
# by cmdopts ($5). 'nodata' offers --with-prim, --with-prim-lflags, and
# --with-prim-cflags. 'dataonly' offers --with-prim and --with-prim-data.
# 'all' offers all four. DEF_PRIM_ARGS and FIND_PRIM_PKG are tailored
# accordingly. The (hardwired) default is 'nodata'.

# Default action ($4) (no, yes, build) is the default action if the user
# offers no guidance via command line parameters. 'build' is deprecated
# and has no hope of working except for COIN-OR ThirdParty packages. Don't
# use it.

# If no .pc file names are specified, the macro will look for prim.pc. If
# there's exactly one .pc file name xxx.pc (either explicit or default) and
# [dfltaction] is build, the macro searches first for coinxxx.pc, then xxx.pc
# (this happens in CHK_PRIM_PKG).

# Define an automake conditional COIN_HAS_PRIM to record the result. If we
# decide to use prim, also define a preprocessor symbol COIN_HAS_PRIM.

# Linker and compiler flag information will be propagated to the space-
# separated list of client packages [clients] ($2) using the _PCFILES variable
# if a .pc file is used, otherwise by the _LFLAGS and _CFLAGS variables of
# client packages. These variables match Requires.private, Libs.private,
# and Cflags.private, respectively, in a .pc file.

# Data directory information is used differently. Typically what's wanted is
# individual variables specifying the data directory for each primitive. Hence
# the macro defines PRIM_DATA for the primitive.

# The macro doesn't test that the specified values actually work. This is
# deliberate.  There's no guarantee that user-specified libraries and/or
# directories actually exist yet. The same possibility exists for values
# returned when pkgconf checks the .pc file.

AC_DEFUN([AC_COIN_CHK_PKG],
[
dnl Issue a warning that 'default_build' is deprecated.
  AC_COIN_DEPRECATE_BUILD_OPTION([COIN_CHK_PKG],[$4])

  AC_REQUIRE([AC_COIN_HAS_PKGCONFIG])

  AC_MSG_CHECKING([for package $1])

dnl Make sure the necessary variables exist for each client package.

  m4_foreach_w([myvar],[$2],
    [AC_SUBST(m4_toupper(myvar)_LFLAGS)
     AC_SUBST(m4_toupper(myvar)_CFLAGS)
     AC_SUBST(m4_toupper(myvar)_PCFILES)
    ])

dnl Check if the user has overridden configure parameters from the
dnl environment. Caution! Quadrigraphs.

  m4_tolower(coin_has_$1)=noInfo
  if test x"$COIN_SKIP_PROJECTS" != x ; then
    for pkg in `echo $COIN_SKIP_PROJECTS | tr '@<:@:upper:@:>@' '@<:@:lower:@:>@'` ; do
      if test "$pkg" = "m4_tolower($1)" ; then
        m4_tolower(coin_has_$1)=skipping
      fi
    done
  fi

dnl If we are not skipping this primary, define and process the
dnl command line options according to the cmdopts parameter. Then invoke
dnl FIND_PRIM_PKG to do the heavy lifting.
dnl Note that we're grandfathering 'build' as an alternative keyword for
dnl 'default_build'. This should go away eventually. -lh,210118-

  if test "$m4_tolower(coin_has_$1)" != skipping ; then
    m4_case(m4_default($5,nodata),
      nodata,  [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,no,
                  AC_COIN_LIBHDR_DFLT_HACK([usage],[$4],[default_use]))],
      dataonly,[AC_COIN_DEF_PRIM_ARGS([$1],yes,no,no,yes,
                  AC_COIN_LIBHDR_DFLT_HACK([usage],[$4],[default_use]))],
      all,     [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,yes,
                  AC_COIN_LIBHDR_DFLT_HACK([usage],[$4],[default_use]))],
      [AC_COIN_DEF_PRIM_ARGS([$1],yes,yes,yes,no,
         AC_COIN_LIBHDR_DFLT_HACK([usage],[$4],[default_use]))])
    AC_COIN_FIND_PRIM_PKG(m4_tolower($1),m4_default([$3],m4_tolower($1)),
                          m4_if([$4],[build],[default_build],[$4]),
                          m4_default([$5],nodata))
    AC_MSG_RESULT([$m4_tolower(coin_has_$1)])
  else
    AC_MSG_RESULT([$m4_tolower(coin_has_$1) (COIN_SKIP_PROJECTS)])
  fi

dnl Possibilities are `yes', `no', or `skipping'. 'Skipping' implies we
dnl decided to skip the package for some reason. 'No' means we wanted the
dnl package but pkgconfig couldn't find the .pc files. 'Yes' means pkgconfig
dnl found the .pc files.  Normalise to yes or no for the remainder.

  if test "$m4_tolower(coin_has_$1)" != yes ; then
    m4_tolower(coin_has_$1)=no
  fi

dnl Create an automake conditional COIN_HAS_PRIM.

  AM_CONDITIONAL(m4_toupper(COIN_HAS_$1),[test $m4_tolower(coin_has_$1) = yes])

dnl If we've located the package, define preprocessor symbol PACKAGE_HAS_PRIM
dnl and augment the necessary variables for the client packages.

  if test $m4_tolower(coin_has_$1) = yes ; then
    AC_DEFINE(m4_toupper(AC_PACKAGE_NAME)_HAS_[]m4_toupper($1),[1],
      [Define to 1 if $1 is available.])
    m4_foreach_w([myvar],[$2],
      [m4_toupper(myvar)_PCFILES="$m4_tolower($1_pcfiles) $m4_toupper(myvar)_PCFILES"
       m4_toupper(myvar)_LFLAGS="$m4_tolower($1_lflags) $m4_toupper(myvar)_LFLAGS"
       m4_toupper(myvar)_CFLAGS="$m4_tolower($1_cflags) $m4_toupper(myvar)_CFLAGS"
      ])

dnl Finally, set up PRIM_DATA, unless the user specified nodata.
    m4_if(m4_default([$5],nodata),nodata,[],
      [AC_SUBST(m4_toupper($1)_DATA)
       m4_toupper($1)_DATA=$m4_tolower($1_data)])
  fi
])   # COIN_CHK_PKG



# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX(VERSION, [ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the specified
#   version of the C++ standard.  If necessary, add switches to CXX and
#   CXXCPP to enable support.  VERSION may be '11' (for the C++11 standard)
#   or '14' (for the C++14 standard).
#
#   The second argument, if specified, indicates whether you insist on an
#   extended mode (e.g. -std=gnu++11) or a strict conformance mode (e.g.
#   -std=c++11).  If neither is specified, you get whatever works, with
#   preference for no added switch, and then for an extended mode.
#
#   The third argument, if specified 'mandatory' or if left unspecified,
#   indicates that baseline support for the specified C++ standard is
#   required and that the macro should error out if no mode with that
#   support is found.  If specified 'optional', then configuration proceeds
#   regardless, after defining HAVE_CXX${VERSION} if and only if a
#   supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014, 2015 Google Inc.; contributed by Alexey Sokolov <sokolov@google.com>
#   Copyright (c) 2015 Paul Norman <penorman@mac.com>
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#   Copyright (c) 2016, 2018 Krzesimir Nowak <qdlacz@gmail.com>
#   Copyright (c) 2019 Enji Cooper <yaneurabeya@gmail.com>
#   Copyright (c) 2020 Jason Merrill <jason@redhat.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 12

dnl  This macro is based on the code from the AX_CXX_COMPILE_STDCXX_11 macro
dnl  (serial version number 13).

AC_DEFUN([AX_CXX_COMPILE_STDCXX], [dnl
  m4_if([$1], [11], [ax_cxx_compile_alternatives="11 0x"],
        [$1], [14], [ax_cxx_compile_alternatives="14 1y"],
        [$1], [17], [ax_cxx_compile_alternatives="17 1z"],
        [m4_fatal([invalid first argument `$1' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$2], [], [],
        [$2], [ext], [],
        [$2], [noext], [],
        [m4_fatal([invalid second argument `$2' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$3], [], [ax_cxx_compile_cxx$1_required=true],
        [$3], [mandatory], [ax_cxx_compile_cxx$1_required=true],
        [$3], [optional], [ax_cxx_compile_cxx$1_required=false],
        [m4_fatal([invalid third argument `$3' to AX_CXX_COMPILE_STDCXX])])
  AC_LANG_PUSH([C++])dnl
  ac_success=no

  m4_if([$2], [], [dnl
    AC_CACHE_CHECK(whether $CXX supports C++$1 features by default,
		   ax_cv_cxx_compile_cxx$1,
      [AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
        [ax_cv_cxx_compile_cxx$1=yes],
        [ax_cv_cxx_compile_cxx$1=no])])
    if test x$ax_cv_cxx_compile_cxx$1 = xyes; then
      ac_success=yes
    fi])

  m4_if([$2], [noext], [], [dnl
  if test x$ac_success = xno; then
    for alternative in ${ax_cxx_compile_alternatives}; do
      switch="-std=gnu++${alternative}"
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                     $cachevar,
        [ac_save_CXX="$CXX"
         CXX="$CXX $switch"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXX="$ac_save_CXX"])
      if eval test x\$$cachevar = xyes; then
        CXX="$CXX $switch"
        if test -n "$CXXCPP" ; then
          CXXCPP="$CXXCPP $switch"
        fi
        ac_success=yes
        break
      fi
    done
  fi])

  m4_if([$2], [ext], [], [dnl
  if test x$ac_success = xno; then
    dnl HP's aCC needs +std=c++11 according to:
    dnl http://h21007.www2.hp.com/portal/download/files/unprot/aCxx/PDF_Release_Notes/769149-001.pdf
    dnl Cray's crayCC needs "-h std=c++11"
    for alternative in ${ax_cxx_compile_alternatives}; do
      for switch in -std=c++${alternative} +std=c++${alternative} "-h std=c++${alternative}"; do
        cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
        AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                       $cachevar,
          [ac_save_CXX="$CXX"
           CXX="$CXX $switch"
           AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
            [eval $cachevar=yes],
            [eval $cachevar=no])
           CXX="$ac_save_CXX"])
        if eval test x\$$cachevar = xyes; then
          CXX="$CXX $switch"
          if test -n "$CXXCPP" ; then
            CXXCPP="$CXXCPP $switch"
          fi
          ac_success=yes
          break
        fi
      done
      if test x$ac_success = xyes; then
        break
      fi
    done
  fi])
  AC_LANG_POP([C++])
  if test x$ax_cxx_compile_cxx$1_required = xtrue; then
    if test x$ac_success = xno; then
      AC_MSG_ERROR([*** A compiler with support for C++$1 language features is required.])
    fi
  fi
  if test x$ac_success = xno; then
    HAVE_CXX$1=0
    AC_MSG_NOTICE([No compiler with C++$1 support was found])
  else
    HAVE_CXX$1=1
    AC_DEFINE(HAVE_CXX$1,1,
              [define if the compiler supports basic C++$1 syntax])
  fi
  AC_SUBST(HAVE_CXX$1)
])


dnl  Test body for checking C++11 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_11],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
)


dnl  Test body for checking C++14 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_14],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_14
)

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_17],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_14
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_17
)

dnl  Tests for new features in C++11

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_11], [[

// If the compiler admits that it is not ready for C++11, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201103L

#error "This is not a C++11 compiler"

#else

namespace cxx11
{

  namespace test_static_assert
  {

    template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

  }

  namespace test_final_override
  {

    struct Base
    {
      virtual ~Base() {}
      virtual void f() {}
    };

    struct Derived : public Base
    {
      virtual ~Derived() override {}
      virtual void f() override {}
    };

  }

  namespace test_double_right_angle_brackets
  {

    template < typename T >
    struct check {};

    typedef check<void> single_type;
    typedef check<check<void>> double_type;
    typedef check<check<check<void>>> triple_type;
    typedef check<check<check<check<void>>>> quadruple_type;

  }

  namespace test_decltype
  {

    int
    f()
    {
      int a = 1;
      decltype(a) b = 2;
      return a + b;
    }

  }

  namespace test_type_deduction
  {

    template < typename T1, typename T2 >
    struct is_same
    {
      static const bool value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static const bool value = true;
    };

    template < typename T1, typename T2 >
    auto
    add(T1 a1, T2 a2) -> decltype(a1 + a2)
    {
      return a1 + a2;
    }

    int
    test(const int c, volatile int v)
    {
      static_assert(is_same<int, decltype(0)>::value == true, "");
      static_assert(is_same<int, decltype(c)>::value == false, "");
      static_assert(is_same<int, decltype(v)>::value == false, "");
      auto ac = c;
      auto av = v;
      auto sumi = ac + av + 'x';
      auto sumf = ac + av + 1.0;
      static_assert(is_same<int, decltype(ac)>::value == true, "");
      static_assert(is_same<int, decltype(av)>::value == true, "");
      static_assert(is_same<int, decltype(sumi)>::value == true, "");
      static_assert(is_same<int, decltype(sumf)>::value == false, "");
      static_assert(is_same<int, decltype(add(c, v))>::value == true, "");
      return (sumf > 0.0) ? sumi : add(c, v);
    }

  }

  namespace test_noexcept
  {

    int f() { return 0; }
    int g() noexcept { return 0; }

    static_assert(noexcept(f()) == false, "");
    static_assert(noexcept(g()) == true, "");

  }

  namespace test_constexpr
  {

    template < typename CharT >
    unsigned long constexpr
    strlen_c_r(const CharT *const s, const unsigned long acc) noexcept
    {
      return *s ? strlen_c_r(s + 1, acc + 1) : acc;
    }

    template < typename CharT >
    unsigned long constexpr
    strlen_c(const CharT *const s) noexcept
    {
      return strlen_c_r(s, 0UL);
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("1") == 1UL, "");
    static_assert(strlen_c("example") == 7UL, "");
    static_assert(strlen_c("another\0example") == 7UL, "");

  }

  namespace test_rvalue_references
  {

    template < int N >
    struct answer
    {
      static constexpr int value = N;
    };

    answer<1> f(int&)       { return answer<1>(); }
    answer<2> f(const int&) { return answer<2>(); }
    answer<3> f(int&&)      { return answer<3>(); }

    void
    test()
    {
      int i = 0;
      const int c = 0;
      static_assert(decltype(f(i))::value == 1, "");
      static_assert(decltype(f(c))::value == 2, "");
      static_assert(decltype(f(0))::value == 3, "");
    }

  }

  namespace test_uniform_initialization
  {

    struct test
    {
      static const int zero {};
      static const int one {1};
    };

    static_assert(test::zero == 0, "");
    static_assert(test::one == 1, "");

  }

  namespace test_lambdas
  {

    void
    test1()
    {
      auto lambda1 = [](){};
      auto lambda2 = lambda1;
      lambda1();
      lambda2();
    }

    int
    test2()
    {
      auto a = [](int i, int j){ return i + j; }(1, 2);
      auto b = []() -> int { return '0'; }();
      auto c = [=](){ return a + b; }();
      auto d = [&](){ return c; }();
      auto e = [a, &b](int x) mutable {
        const auto identity = [](int y){ return y; };
        for (auto i = 0; i < a; ++i)
          a += b--;
        return x + identity(a + b);
      }(0);
      return a + b + c + d + e;
    }

    int
    test3()
    {
      const auto nullary = [](){ return 0; };
      const auto unary = [](int x){ return x; };
      using nullary_t = decltype(nullary);
      using unary_t = decltype(unary);
      const auto higher1st = [](nullary_t f){ return f(); };
      const auto higher2nd = [unary](nullary_t f1){
        return [unary, f1](unary_t f2){ return f2(unary(f1())); };
      };
      return higher1st(nullary) + higher2nd(nullary)(unary);
    }

  }

  namespace test_variadic_templates
  {

    template <int...>
    struct sum;

    template <int N0, int... N1toN>
    struct sum<N0, N1toN...>
    {
      static constexpr auto value = N0 + sum<N1toN...>::value;
    };

    template <>
    struct sum<>
    {
      static constexpr auto value = 0;
    };

    static_assert(sum<>::value == 0, "");
    static_assert(sum<1>::value == 1, "");
    static_assert(sum<23>::value == 23, "");
    static_assert(sum<1, 2>::value == 3, "");
    static_assert(sum<5, 5, 11>::value == 21, "");
    static_assert(sum<2, 3, 5, 7, 11, 13>::value == 41, "");

  }

  // http://stackoverflow.com/questions/13728184/template-aliases-and-sfinae
  // Clang 3.1 fails with headers of libstd++ 4.8.3 when using std::function
  // because of this.
  namespace test_template_alias_sfinae
  {

    struct foo {};

    template<typename T>
    using member = typename T::member_type;

    template<typename T>
    void func(...) {}

    template<typename T>
    void func(member<T>*) {}

    void test();

    void test() { func<foo>(0); }

  }

}  // namespace cxx11

#endif  // __cplusplus >= 201103L

]])


dnl  Tests for new features in C++14

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_14], [[

// If the compiler admits that it is not ready for C++14, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201402L

#error "This is not a C++14 compiler"

#else

namespace cxx14
{

  namespace test_polymorphic_lambdas
  {

    int
    test()
    {
      const auto lambda = [](auto&&... args){
        const auto istiny = [](auto x){
          return (sizeof(x) == 1UL) ? 1 : 0;
        };
        const int aretiny[] = { istiny(args)... };
        return aretiny[0];
      };
      return lambda(1, 1L, 1.0f, '1');
    }

  }

  namespace test_binary_literals
  {

    constexpr auto ivii = 0b0000000000101010;
    static_assert(ivii == 42, "wrong value");

  }

  namespace test_generalized_constexpr
  {

    template < typename CharT >
    constexpr unsigned long
    strlen_c(const CharT *const s) noexcept
    {
      auto length = 0UL;
      for (auto p = s; *p; ++p)
        ++length;
      return length;
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("x") == 1UL, "");
    static_assert(strlen_c("test") == 4UL, "");
    static_assert(strlen_c("another\0test") == 7UL, "");

  }

  namespace test_lambda_init_capture
  {

    int
    test()
    {
      auto x = 0;
      const auto lambda1 = [a = x](int b){ return a + b; };
      const auto lambda2 = [a = lambda1(x)](){ return a; };
      return lambda2();
    }

  }

  namespace test_digit_separators
  {

    constexpr auto ten_million = 100'000'000;
    static_assert(ten_million == 100000000, "");

  }

  namespace test_return_type_deduction
  {

    auto f(int& x) { return x; }
    decltype(auto) g(int& x) { return x; }

    template < typename T1, typename T2 >
    struct is_same
    {
      static constexpr auto value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static constexpr auto value = true;
    };

    int
    test()
    {
      auto x = 0;
      static_assert(is_same<int, decltype(f(x))>::value, "");
      static_assert(is_same<int&, decltype(g(x))>::value, "");
      return x;
    }

  }

}  // namespace cxx14

#endif  // __cplusplus >= 201402L

]])


dnl  Tests for new features in C++17

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_17], [[

// If the compiler admits that it is not ready for C++17, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201703L

#error "This is not a C++17 compiler"

#else

#include <initializer_list>
#include <utility>
#include <type_traits>

namespace cxx17
{

  namespace test_constexpr_lambdas
  {

    constexpr int foo = [](){return 42;}();

  }

  namespace test::nested_namespace::definitions
  {

  }

  namespace test_fold_expression
  {

    template<typename... Args>
    int multiply(Args... args)
    {
      return (args * ... * 1);
    }

    template<typename... Args>
    bool all(Args... args)
    {
      return (args && ...);
    }

  }

  namespace test_extended_static_assert
  {

    static_assert (true);

  }

  namespace test_auto_brace_init_list
  {

    auto foo = {5};
    auto bar {5};

    static_assert(std::is_same<std::initializer_list<int>, decltype(foo)>::value);
    static_assert(std::is_same<int, decltype(bar)>::value);
  }

  namespace test_typename_in_template_template_parameter
  {

    template<template<typename> typename X> struct D;

  }

  namespace test_fallthrough_nodiscard_maybe_unused_attributes
  {

    int f1()
    {
      return 42;
    }

    [[nodiscard]] int f2()
    {
      [[maybe_unused]] auto unused = f1();

      switch (f1())
      {
      case 17:
        f1();
        [[fallthrough]];
      case 42:
        f1();
      }
      return f1();
    }

  }

  namespace test_extended_aggregate_initialization
  {

    struct base1
    {
      int b1, b2 = 42;
    };

    struct base2
    {
      base2() {
        b3 = 42;
      }
      int b3;
    };

    struct derived : base1, base2
    {
        int d;
    };

    derived d1 {{1, 2}, {}, 4};  // full initialization
    derived d2 {{}, {}, 4};      // value-initialized bases

  }

  namespace test_general_range_based_for_loop
  {

    struct iter
    {
      int i;

      int& operator* ()
      {
        return i;
      }

      const int& operator* () const
      {
        return i;
      }

      iter& operator++()
      {
        ++i;
        return *this;
      }
    };

    struct sentinel
    {
      int i;
    };

    bool operator== (const iter& i, const sentinel& s)
    {
      return i.i == s.i;
    }

    bool operator!= (const iter& i, const sentinel& s)
    {
      return !(i == s);
    }

    struct range
    {
      iter begin() const
      {
        return {0};
      }

      sentinel end() const
      {
        return {5};
      }
    };

    void f()
    {
      range r {};

      for (auto i : r)
      {
        [[maybe_unused]] auto v = i;
      }
    }

  }

  namespace test_lambda_capture_asterisk_this_by_value
  {

    struct t
    {
      int i;
      int foo()
      {
        return [*this]()
        {
          return i;
        }();
      }
    };

  }

  namespace test_enum_class_construction
  {

    enum class byte : unsigned char
    {};

    byte foo {42};

  }

  namespace test_constexpr_if
  {

    template <bool cond>
    int f ()
    {
      if constexpr(cond)
      {
        return 13;
      }
      else
      {
        return 42;
      }
    }

  }

  namespace test_selection_statement_with_initializer
  {

    int f()
    {
      return 13;
    }

    int f2()
    {
      if (auto i = f(); i > 0)
      {
        return 3;
      }

      switch (auto i = f(); i + 4)
      {
      case 17:
        return 2;

      default:
        return 1;
      }
    }

  }

  namespace test_template_argument_deduction_for_class_templates
  {

    template <typename T1, typename T2>
    struct pair
    {
      pair (T1 p1, T2 p2)
        : m1 {p1},
          m2 {p2}
      {}

      T1 m1;
      T2 m2;
    };

    void f()
    {
      [[maybe_unused]] auto p = pair{13, 42u};
    }

  }

  namespace test_non_type_auto_template_parameters
  {

    template <auto n>
    struct B
    {};

    B<5> b1;
    B<'a'> b2;

  }

  namespace test_structured_bindings
  {

    int arr[2] = { 1, 2 };
    std::pair<int, int> pr = { 1, 2 };

    auto f1() -> int(&)[2]
    {
      return arr;
    }

    auto f2() -> std::pair<int, int>&
    {
      return pr;
    }

    struct S
    {
      int x1 : 2;
      volatile double y1;
    };

    S f3()
    {
      return {};
    }

    auto [ x1, y1 ] = f1();
    auto& [ xr1, yr1 ] = f1();
    auto [ x2, y2 ] = f2();
    auto& [ xr2, yr2 ] = f2();
    const auto [ x3, y3 ] = f3();

  }

  namespace test_exception_spec_type_system
  {

    struct Good {};
    struct Bad {};

    void g1() noexcept;
    void g2();

    template<typename T>
    Bad
    f(T*, T*);

    template<typename T1, typename T2>
    Good
    f(T1*, T2*);

    static_assert (std::is_same_v<Good, decltype(f(g1, g2))>);

  }

  namespace test_inline_variables
  {

    template<class T> void f(T)
    {}

    template<class T> inline T g(T)
    {
      return T{};
    }

    template<> inline void f<>(int)
    {}

    template<> int g<>(int)
    {
      return 5;
    }

  }

}  // namespace cxx17

#endif  // __cplusplus < 201703L

]])

# libtool.m4 - Configure libtool for the host system. -*-Autoconf-*-
#
#   Copyright (C) 1996-2001, 2003-2015 Free Software Foundation, Inc.
#   Written by Gordon Matzigkeit, 1996
#
# This file is free software; the Free Software Foundation gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.

m4_define([_LT_COPYING], [dnl
# Copyright (C) 2014 Free Software Foundation, Inc.
# This is free software; see the source for copying conditions.  There is NO
# warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# GNU Libtool is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of of the License, or
# (at your option) any later version.
#
# As a special exception to the GNU General Public License, if you
# distribute this file as part of a program or library that is built
# using GNU Libtool, you may include this file under the  same
# distribution terms that you use for the rest of that program.
#
# GNU Libtool is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
])

# serial 58 LT_INIT


# LT_PREREQ(VERSION)
# ------------------
# Complain and exit if this libtool version is less that VERSION.
m4_defun([LT_PREREQ],
[m4_if(m4_version_compare(m4_defn([LT_PACKAGE_VERSION]), [$1]), -1,
       [m4_default([$3],
		   [m4_fatal([Libtool version $1 or higher is required],
		             63)])],
       [$2])])


# _LT_CHECK_BUILDDIR
# ------------------
# Complain if the absolute build directory name contains unusual characters
m4_defun([_LT_CHECK_BUILDDIR],
[case `pwd` in
  *\ * | *\	*)
    AC_MSG_WARN([Libtool does not cope well with whitespace in `pwd`]) ;;
esac
])


# LT_INIT([OPTIONS])
# ------------------
AC_DEFUN([LT_INIT],
[AC_PREREQ([2.62])dnl We use AC_PATH_PROGS_FEATURE_CHECK
AC_REQUIRE([AC_CONFIG_AUX_DIR_DEFAULT])dnl
AC_BEFORE([$0], [LT_LANG])dnl
AC_BEFORE([$0], [LT_OUTPUT])dnl
AC_BEFORE([$0], [LTDL_INIT])dnl
m4_require([_LT_CHECK_BUILDDIR])dnl

dnl Autoconf doesn't catch unexpanded LT_ macros by default:
m4_pattern_forbid([^_?LT_[A-Z_]+$])dnl
m4_pattern_allow([^(_LT_EOF|LT_DLGLOBAL|LT_DLLAZY_OR_NOW|LT_MULTI_MODULE)$])dnl
dnl aclocal doesn't pull ltoptions.m4, ltsugar.m4, or ltversion.m4
dnl unless we require an AC_DEFUNed macro:
AC_REQUIRE([LTOPTIONS_VERSION])dnl
AC_REQUIRE([LTSUGAR_VERSION])dnl
AC_REQUIRE([LTVERSION_VERSION])dnl
AC_REQUIRE([LTOBSOLETE_VERSION])dnl
m4_require([_LT_PROG_LTMAIN])dnl

_LT_SHELL_INIT([SHELL=${CONFIG_SHELL-/bin/sh}])

dnl Parse OPTIONS
_LT_SET_OPTIONS([$0], [$1])

# This can be used to rebuild libtool when needed
LIBTOOL_DEPS=$ltmain

# Always use our own libtool.
LIBTOOL='$(SHELL) $(top_builddir)/libtool'
AC_SUBST(LIBTOOL)dnl

_LT_SETUP

# Only expand once:
m4_define([LT_INIT])
])# LT_INIT

# Old names:
AU_ALIAS([AC_PROG_LIBTOOL], [LT_INIT])
AU_ALIAS([AM_PROG_LIBTOOL], [LT_INIT])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_PROG_LIBTOOL], [])
dnl AC_DEFUN([AM_PROG_LIBTOOL], [])


# _LT_PREPARE_CC_BASENAME
# -----------------------
m4_defun([_LT_PREPARE_CC_BASENAME], [
# Calculate cc_basename.  Skip known compiler wrappers and cross-prefix.
func_cc_basename ()
{
    for cc_temp in @S|@*""; do
      case $cc_temp in
        compile | *[[\\/]]compile | ccache | *[[\\/]]ccache ) ;;
        distcc | *[[\\/]]distcc | purify | *[[\\/]]purify ) ;;
        \-*) ;;
        *) break;;
      esac
    done
    func_cc_basename_result=`$ECHO "$cc_temp" | $SED "s%.*/%%; s%^$host_alias-%%"`
}
])# _LT_PREPARE_CC_BASENAME


# _LT_CC_BASENAME(CC)
# -------------------
# It would be clearer to call AC_REQUIREs from _LT_PREPARE_CC_BASENAME,
# but that macro is also expanded into generated libtool script, which
# arranges for $SED and $ECHO to be set by different means.
m4_defun([_LT_CC_BASENAME],
[m4_require([_LT_PREPARE_CC_BASENAME])dnl
AC_REQUIRE([_LT_DECL_SED])dnl
AC_REQUIRE([_LT_PROG_ECHO_BACKSLASH])dnl
func_cc_basename $1
cc_basename=$func_cc_basename_result
])


# _LT_FILEUTILS_DEFAULTS
# ----------------------
# It is okay to use these file commands and assume they have been set
# sensibly after 'm4_require([_LT_FILEUTILS_DEFAULTS])'.
m4_defun([_LT_FILEUTILS_DEFAULTS],
[: ${CP="cp -f"}
: ${MV="mv -f"}
: ${RM="rm -f"}
])# _LT_FILEUTILS_DEFAULTS


# _LT_SETUP
# ---------
m4_defun([_LT_SETUP],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
AC_REQUIRE([_LT_PREPARE_SED_QUOTE_VARS])dnl
AC_REQUIRE([_LT_PROG_ECHO_BACKSLASH])dnl

_LT_DECL([], [PATH_SEPARATOR], [1], [The PATH separator for the build system])dnl
dnl
_LT_DECL([], [host_alias], [0], [The host system])dnl
_LT_DECL([], [host], [0])dnl
_LT_DECL([], [host_os], [0])dnl
dnl
_LT_DECL([], [build_alias], [0], [The build system])dnl
_LT_DECL([], [build], [0])dnl
_LT_DECL([], [build_os], [0])dnl
dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([LT_PATH_LD])dnl
AC_REQUIRE([LT_PATH_NM])dnl
dnl
AC_REQUIRE([AC_PROG_LN_S])dnl
test -z "$LN_S" && LN_S="ln -s"
_LT_DECL([], [LN_S], [1], [Whether we need soft or hard links])dnl
dnl
AC_REQUIRE([LT_CMD_MAX_LEN])dnl
_LT_DECL([objext], [ac_objext], [0], [Object file suffix (normally "o")])dnl
_LT_DECL([], [exeext], [0], [Executable file suffix (normally "")])dnl
dnl
m4_require([_LT_FILEUTILS_DEFAULTS])dnl
m4_require([_LT_CHECK_SHELL_FEATURES])dnl
m4_require([_LT_PATH_CONVERSION_FUNCTIONS])dnl
m4_require([_LT_CMD_RELOAD])dnl
m4_require([_LT_CHECK_MAGIC_METHOD])dnl
m4_require([_LT_CHECK_SHAREDLIB_FROM_LINKLIB])dnl
m4_require([_LT_CMD_OLD_ARCHIVE])dnl
m4_require([_LT_CMD_GLOBAL_SYMBOLS])dnl
m4_require([_LT_WITH_SYSROOT])dnl
m4_require([_LT_CMD_TRUNCATE])dnl

_LT_CONFIG_LIBTOOL_INIT([
# See if we are running on zsh, and set the options that allow our
# commands through without removal of \ escapes INIT.
if test -n "\${ZSH_VERSION+set}"; then
   setopt NO_GLOB_SUBST
fi
])
if test -n "${ZSH_VERSION+set}"; then
   setopt NO_GLOB_SUBST
fi

_LT_CHECK_OBJDIR

m4_require([_LT_TAG_COMPILER])dnl

case $host_os in
aix3*)
  # AIX sometimes has problems with the GCC collect2 program.  For some
  # reason, if we set the COLLECT_NAMES environment variable, the problems
  # vanish in a puff of smoke.
  if test set != "${COLLECT_NAMES+set}"; then
    COLLECT_NAMES=
    export COLLECT_NAMES
  fi
  ;;
esac

# Global variables:
ofile=libtool
can_build_shared=yes

# All known linkers require a '.a' archive for static linking (except MSVC,
# which needs '.lib').
libext=a

with_gnu_ld=$lt_cv_prog_gnu_ld

old_CC=$CC
old_CFLAGS=$CFLAGS

# Set sane defaults for various variables
test -z "$CC" && CC=cc
test -z "$LTCC" && LTCC=$CC
test -z "$LTCFLAGS" && LTCFLAGS=$CFLAGS
test -z "$LD" && LD=ld
test -z "$ac_objext" && ac_objext=o

_LT_CC_BASENAME([$compiler])

# Only perform the check for file, if the check method requires it
test -z "$MAGIC_CMD" && MAGIC_CMD=file
case $deplibs_check_method in
file_magic*)
  if test "$file_magic_cmd" = '$MAGIC_CMD'; then
    _LT_PATH_MAGIC
  fi
  ;;
esac

# Use C for the default configuration in the libtool script
LT_SUPPORTED_TAG([CC])
_LT_LANG_C_CONFIG
_LT_LANG_DEFAULT_CONFIG
_LT_CONFIG_COMMANDS
])# _LT_SETUP


# _LT_PREPARE_SED_QUOTE_VARS
# --------------------------
# Define a few sed substitution that help us do robust quoting.
m4_defun([_LT_PREPARE_SED_QUOTE_VARS],
[# Backslashify metacharacters that are still active within
# double-quoted strings.
sed_quote_subst='s/\([["`$\\]]\)/\\\1/g'

# Same as above, but do not quote variable references.
double_quote_subst='s/\([["`\\]]\)/\\\1/g'

# Sed substitution to delay expansion of an escaped shell variable in a
# double_quote_subst'ed string.
delay_variable_subst='s/\\\\\\\\\\\$/\\\\\\$/g'

# Sed substitution to delay expansion of an escaped single quote.
delay_single_quote_subst='s/'\''/'\'\\\\\\\'\''/g'

# Sed substitution to avoid accidental globbing in evaled expressions
no_glob_subst='s/\*/\\\*/g'
])

# _LT_PROG_LTMAIN
# ---------------
# Note that this code is called both from 'configure', and 'config.status'
# now that we use AC_CONFIG_COMMANDS to generate libtool.  Notably,
# 'config.status' has no value for ac_aux_dir unless we are using Automake,
# so we pass a copy along to make sure it has a sensible value anyway.
m4_defun([_LT_PROG_LTMAIN],
[m4_ifdef([AC_REQUIRE_AUX_FILE], [AC_REQUIRE_AUX_FILE([ltmain.sh])])dnl
_LT_CONFIG_LIBTOOL_INIT([ac_aux_dir='$ac_aux_dir'])
ltmain=$ac_aux_dir/ltmain.sh
])# _LT_PROG_LTMAIN



# So that we can recreate a full libtool script including additional
# tags, we accumulate the chunks of code to send to AC_CONFIG_COMMANDS
# in macros and then make a single call at the end using the 'libtool'
# label.


# _LT_CONFIG_LIBTOOL_INIT([INIT-COMMANDS])
# ----------------------------------------
# Register INIT-COMMANDS to be passed to AC_CONFIG_COMMANDS later.
m4_define([_LT_CONFIG_LIBTOOL_INIT],
[m4_ifval([$1],
          [m4_append([_LT_OUTPUT_LIBTOOL_INIT],
                     [$1
])])])

# Initialize.
m4_define([_LT_OUTPUT_LIBTOOL_INIT])


# _LT_CONFIG_LIBTOOL([COMMANDS])
# ------------------------------
# Register COMMANDS to be passed to AC_CONFIG_COMMANDS later.
m4_define([_LT_CONFIG_LIBTOOL],
[m4_ifval([$1],
          [m4_append([_LT_OUTPUT_LIBTOOL_COMMANDS],
                     [$1
])])])

# Initialize.
m4_define([_LT_OUTPUT_LIBTOOL_COMMANDS])


# _LT_CONFIG_SAVE_COMMANDS([COMMANDS], [INIT_COMMANDS])
# -----------------------------------------------------
m4_defun([_LT_CONFIG_SAVE_COMMANDS],
[_LT_CONFIG_LIBTOOL([$1])
_LT_CONFIG_LIBTOOL_INIT([$2])
])


# _LT_FORMAT_COMMENT([COMMENT])
# -----------------------------
# Add leading comment marks to the start of each line, and a trailing
# full-stop to the whole comment if one is not present already.
m4_define([_LT_FORMAT_COMMENT],
[m4_ifval([$1], [
m4_bpatsubst([m4_bpatsubst([$1], [^ *], [# ])],
              [['`$\]], [\\\&])]m4_bmatch([$1], [[!?.]$], [], [.])
)])





# _LT_DECL([CONFIGNAME], VARNAME, VALUE, [DESCRIPTION], [IS-TAGGED?])
# -------------------------------------------------------------------
# CONFIGNAME is the name given to the value in the libtool script.
# VARNAME is the (base) name used in the configure script.
# VALUE may be 0, 1 or 2 for a computed quote escaped value based on
# VARNAME.  Any other value will be used directly.
m4_define([_LT_DECL],
[lt_if_append_uniq([lt_decl_varnames], [$2], [, ],
    [lt_dict_add_subkey([lt_decl_dict], [$2], [libtool_name],
	[m4_ifval([$1], [$1], [$2])])
    lt_dict_add_subkey([lt_decl_dict], [$2], [value], [$3])
    m4_ifval([$4],
	[lt_dict_add_subkey([lt_decl_dict], [$2], [description], [$4])])
    lt_dict_add_subkey([lt_decl_dict], [$2],
	[tagged?], [m4_ifval([$5], [yes], [no])])])
])


# _LT_TAGDECL([CONFIGNAME], VARNAME, VALUE, [DESCRIPTION])
# --------------------------------------------------------
m4_define([_LT_TAGDECL], [_LT_DECL([$1], [$2], [$3], [$4], [yes])])


# lt_decl_tag_varnames([SEPARATOR], [VARNAME1...])
# ------------------------------------------------
m4_define([lt_decl_tag_varnames],
[_lt_decl_filter([tagged?], [yes], $@)])


# _lt_decl_filter(SUBKEY, VALUE, [SEPARATOR], [VARNAME1..])
# ---------------------------------------------------------
m4_define([_lt_decl_filter],
[m4_case([$#],
  [0], [m4_fatal([$0: too few arguments: $#])],
  [1], [m4_fatal([$0: too few arguments: $#: $1])],
  [2], [lt_dict_filter([lt_decl_dict], [$1], [$2], [], lt_decl_varnames)],
  [3], [lt_dict_filter([lt_decl_dict], [$1], [$2], [$3], lt_decl_varnames)],
  [lt_dict_filter([lt_decl_dict], $@)])[]dnl
])


# lt_decl_quote_varnames([SEPARATOR], [VARNAME1...])
# --------------------------------------------------
m4_define([lt_decl_quote_varnames],
[_lt_decl_filter([value], [1], $@)])


# lt_decl_dquote_varnames([SEPARATOR], [VARNAME1...])
# ---------------------------------------------------
m4_define([lt_decl_dquote_varnames],
[_lt_decl_filter([value], [2], $@)])


# lt_decl_varnames_tagged([SEPARATOR], [VARNAME1...])
# ---------------------------------------------------
m4_define([lt_decl_varnames_tagged],
[m4_assert([$# <= 2])dnl
_$0(m4_quote(m4_default([$1], [[, ]])),
    m4_ifval([$2], [[$2]], [m4_dquote(lt_decl_tag_varnames)]),
    m4_split(m4_normalize(m4_quote(_LT_TAGS)), [ ]))])
m4_define([_lt_decl_varnames_tagged],
[m4_ifval([$3], [lt_combine([$1], [$2], [_], $3)])])


# lt_decl_all_varnames([SEPARATOR], [VARNAME1...])
# ------------------------------------------------
m4_define([lt_decl_all_varnames],
[_$0(m4_quote(m4_default([$1], [[, ]])),
     m4_if([$2], [],
	   m4_quote(lt_decl_varnames),
	m4_quote(m4_shift($@))))[]dnl
])
m4_define([_lt_decl_all_varnames],
[lt_join($@, lt_decl_varnames_tagged([$1],
			lt_decl_tag_varnames([[, ]], m4_shift($@))))dnl
])


# _LT_CONFIG_STATUS_DECLARE([VARNAME])
# ------------------------------------
# Quote a variable value, and forward it to 'config.status' so that its
# declaration there will have the same value as in 'configure'.  VARNAME
# must have a single quote delimited value for this to work.
m4_define([_LT_CONFIG_STATUS_DECLARE],
[$1='`$ECHO "$][$1" | $SED "$delay_single_quote_subst"`'])


# _LT_CONFIG_STATUS_DECLARATIONS
# ------------------------------
# We delimit libtool config variables with single quotes, so when
# we write them to config.status, we have to be sure to quote all
# embedded single quotes properly.  In configure, this macro expands
# each variable declared with _LT_DECL (and _LT_TAGDECL) into:
#
#    <var>='`$ECHO "$<var>" | $SED "$delay_single_quote_subst"`'
m4_defun([_LT_CONFIG_STATUS_DECLARATIONS],
[m4_foreach([_lt_var], m4_quote(lt_decl_all_varnames),
    [m4_n([_LT_CONFIG_STATUS_DECLARE(_lt_var)])])])


# _LT_LIBTOOL_TAGS
# ----------------
# Output comment and list of tags supported by the script
m4_defun([_LT_LIBTOOL_TAGS],
[_LT_FORMAT_COMMENT([The names of the tagged configurations supported by this script])dnl
available_tags='_LT_TAGS'dnl
])


# _LT_LIBTOOL_DECLARE(VARNAME, [TAG])
# -----------------------------------
# Extract the dictionary values for VARNAME (optionally with TAG) and
# expand to a commented shell variable setting:
#
#    # Some comment about what VAR is for.
#    visible_name=$lt_internal_name
m4_define([_LT_LIBTOOL_DECLARE],
[_LT_FORMAT_COMMENT(m4_quote(lt_dict_fetch([lt_decl_dict], [$1],
					   [description])))[]dnl
m4_pushdef([_libtool_name],
    m4_quote(lt_dict_fetch([lt_decl_dict], [$1], [libtool_name])))[]dnl
m4_case(m4_quote(lt_dict_fetch([lt_decl_dict], [$1], [value])),
    [0], [_libtool_name=[$]$1],
    [1], [_libtool_name=$lt_[]$1],
    [2], [_libtool_name=$lt_[]$1],
    [_libtool_name=lt_dict_fetch([lt_decl_dict], [$1], [value])])[]dnl
m4_ifval([$2], [_$2])[]m4_popdef([_libtool_name])[]dnl
])


# _LT_LIBTOOL_CONFIG_VARS
# -----------------------
# Produce commented declarations of non-tagged libtool config variables
# suitable for insertion in the LIBTOOL CONFIG section of the 'libtool'
# script.  Tagged libtool config variables (even for the LIBTOOL CONFIG
# section) are produced by _LT_LIBTOOL_TAG_VARS.
m4_defun([_LT_LIBTOOL_CONFIG_VARS],
[m4_foreach([_lt_var],
    m4_quote(_lt_decl_filter([tagged?], [no], [], lt_decl_varnames)),
    [m4_n([_LT_LIBTOOL_DECLARE(_lt_var)])])])


# _LT_LIBTOOL_TAG_VARS(TAG)
# -------------------------
m4_define([_LT_LIBTOOL_TAG_VARS],
[m4_foreach([_lt_var], m4_quote(lt_decl_tag_varnames),
    [m4_n([_LT_LIBTOOL_DECLARE(_lt_var, [$1])])])])


# _LT_TAGVAR(VARNAME, [TAGNAME])
# ------------------------------
m4_define([_LT_TAGVAR], [m4_ifval([$2], [$1_$2], [$1])])


# _LT_CONFIG_COMMANDS
# -------------------
# Send accumulated output to $CONFIG_STATUS.  Thanks to the lists of
# variables for single and double quote escaping we saved from calls
# to _LT_DECL, we can put quote escaped variables declarations
# into 'config.status', and then the shell code to quote escape them in
# for loops in 'config.status'.  Finally, any additional code accumulated
# from calls to _LT_CONFIG_LIBTOOL_INIT is expanded.
m4_defun([_LT_CONFIG_COMMANDS],
[AC_PROVIDE_IFELSE([LT_OUTPUT],
	dnl If the libtool generation code has been placed in $CONFIG_LT,
	dnl instead of duplicating it all over again into config.status,
	dnl then we will have config.status run $CONFIG_LT later, so it
	dnl needs to know what name is stored there:
        [AC_CONFIG_COMMANDS([libtool],
            [$SHELL $CONFIG_LT || AS_EXIT(1)], [CONFIG_LT='$CONFIG_LT'])],
    dnl If the libtool generation code is destined for config.status,
    dnl expand the accumulated commands and init code now:
    [AC_CONFIG_COMMANDS([libtool],
        [_LT_OUTPUT_LIBTOOL_COMMANDS], [_LT_OUTPUT_LIBTOOL_COMMANDS_INIT])])
])#_LT_CONFIG_COMMANDS


# Initialize.
m4_define([_LT_OUTPUT_LIBTOOL_COMMANDS_INIT],
[

# The HP-UX ksh and POSIX shell print the target directory to stdout
# if CDPATH is set.
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH

sed_quote_subst='$sed_quote_subst'
double_quote_subst='$double_quote_subst'
delay_variable_subst='$delay_variable_subst'
_LT_CONFIG_STATUS_DECLARATIONS
LTCC='$LTCC'
LTCFLAGS='$LTCFLAGS'
compiler='$compiler_DEFAULT'

# A function that is used when there is no print builtin or printf.
func_fallback_echo ()
{
  eval 'cat <<_LTECHO_EOF
\$[]1
_LTECHO_EOF'
}

# Quote evaled strings.
for var in lt_decl_all_varnames([[ \
]], lt_decl_quote_varnames); do
    case \`eval \\\\\$ECHO \\\\""\\\\\$\$var"\\\\"\` in
    *[[\\\\\\\`\\"\\\$]]*)
      eval "lt_\$var=\\\\\\"\\\`\\\$ECHO \\"\\\$\$var\\" | \\\$SED \\"\\\$sed_quote_subst\\"\\\`\\\\\\"" ## exclude from sc_prohibit_nested_quotes
      ;;
    *)
      eval "lt_\$var=\\\\\\"\\\$\$var\\\\\\""
      ;;
    esac
done

# Double-quote double-evaled strings.
for var in lt_decl_all_varnames([[ \
]], lt_decl_dquote_varnames); do
    case \`eval \\\\\$ECHO \\\\""\\\\\$\$var"\\\\"\` in
    *[[\\\\\\\`\\"\\\$]]*)
      eval "lt_\$var=\\\\\\"\\\`\\\$ECHO \\"\\\$\$var\\" | \\\$SED -e \\"\\\$double_quote_subst\\" -e \\"\\\$sed_quote_subst\\" -e \\"\\\$delay_variable_subst\\"\\\`\\\\\\"" ## exclude from sc_prohibit_nested_quotes
      ;;
    *)
      eval "lt_\$var=\\\\\\"\\\$\$var\\\\\\""
      ;;
    esac
done

_LT_OUTPUT_LIBTOOL_INIT
])

# _LT_GENERATED_FILE_INIT(FILE, [COMMENT])
# ------------------------------------
# Generate a child script FILE with all initialization necessary to
# reuse the environment learned by the parent script, and make the
# file executable.  If COMMENT is supplied, it is inserted after the
# '#!' sequence but before initialization text begins.  After this
# macro, additional text can be appended to FILE to form the body of
# the child script.  The macro ends with non-zero status if the
# file could not be fully written (such as if the disk is full).
m4_ifdef([AS_INIT_GENERATED],
[m4_defun([_LT_GENERATED_FILE_INIT],[AS_INIT_GENERATED($@)])],
[m4_defun([_LT_GENERATED_FILE_INIT],
[m4_require([AS_PREPARE])]dnl
[m4_pushdef([AS_MESSAGE_LOG_FD])]dnl
[lt_write_fail=0
cat >$1 <<_ASEOF || lt_write_fail=1
#! $SHELL
# Generated by $as_me.
$2
SHELL=\${CONFIG_SHELL-$SHELL}
export SHELL
_ASEOF
cat >>$1 <<\_ASEOF || lt_write_fail=1
AS_SHELL_SANITIZE
_AS_PREPARE
exec AS_MESSAGE_FD>&1
_ASEOF
test 0 = "$lt_write_fail" && chmod +x $1[]dnl
m4_popdef([AS_MESSAGE_LOG_FD])])])# _LT_GENERATED_FILE_INIT

# LT_OUTPUT
# ---------
# This macro allows early generation of the libtool script (before
# AC_OUTPUT is called), incase it is used in configure for compilation
# tests.
AC_DEFUN([LT_OUTPUT],
[: ${CONFIG_LT=./config.lt}
AC_MSG_NOTICE([creating $CONFIG_LT])
_LT_GENERATED_FILE_INIT(["$CONFIG_LT"],
[# Run this file to recreate a libtool stub with the current configuration.])

cat >>"$CONFIG_LT" <<\_LTEOF
lt_cl_silent=false
exec AS_MESSAGE_LOG_FD>>config.log
{
  echo
  AS_BOX([Running $as_me.])
} >&AS_MESSAGE_LOG_FD

lt_cl_help="\
'$as_me' creates a local libtool stub from the current configuration,
for use in further configure time tests before the real libtool is
generated.

Usage: $[0] [[OPTIONS]]

  -h, --help      print this help, then exit
  -V, --version   print version number, then exit
  -q, --quiet     do not print progress messages
  -d, --debug     don't remove temporary files

Report bugs to <bug-libtool@gnu.org>."

lt_cl_version="\
m4_ifset([AC_PACKAGE_NAME], [AC_PACKAGE_NAME ])config.lt[]dnl
m4_ifset([AC_PACKAGE_VERSION], [ AC_PACKAGE_VERSION])
configured by $[0], generated by m4_PACKAGE_STRING.

Copyright (C) 2011 Free Software Foundation, Inc.
This config.lt script is free software; the Free Software Foundation
gives unlimited permision to copy, distribute and modify it."

while test 0 != $[#]
do
  case $[1] in
    --version | --v* | -V )
      echo "$lt_cl_version"; exit 0 ;;
    --help | --h* | -h )
      echo "$lt_cl_help"; exit 0 ;;
    --debug | --d* | -d )
      debug=: ;;
    --quiet | --q* | --silent | --s* | -q )
      lt_cl_silent=: ;;

    -*) AC_MSG_ERROR([unrecognized option: $[1]
Try '$[0] --help' for more information.]) ;;

    *) AC_MSG_ERROR([unrecognized argument: $[1]
Try '$[0] --help' for more information.]) ;;
  esac
  shift
done

if $lt_cl_silent; then
  exec AS_MESSAGE_FD>/dev/null
fi
_LTEOF

cat >>"$CONFIG_LT" <<_LTEOF
_LT_OUTPUT_LIBTOOL_COMMANDS_INIT
_LTEOF

cat >>"$CONFIG_LT" <<\_LTEOF
AC_MSG_NOTICE([creating $ofile])
_LT_OUTPUT_LIBTOOL_COMMANDS
AS_EXIT(0)
_LTEOF
chmod +x "$CONFIG_LT"

# configure is writing to config.log, but config.lt does its own redirection,
# appending to config.log, which fails on DOS, as config.log is still kept
# open by configure.  Here we exec the FD to /dev/null, effectively closing
# config.log, so it can be properly (re)opened and appended to by config.lt.
lt_cl_success=:
test yes = "$silent" &&
  lt_config_lt_args="$lt_config_lt_args --quiet"
exec AS_MESSAGE_LOG_FD>/dev/null
$SHELL "$CONFIG_LT" $lt_config_lt_args || lt_cl_success=false
exec AS_MESSAGE_LOG_FD>>config.log
$lt_cl_success || AS_EXIT(1)
])# LT_OUTPUT


# _LT_CONFIG(TAG)
# ---------------
# If TAG is the built-in tag, create an initial libtool script with a
# default configuration from the untagged config vars.  Otherwise add code
# to config.status for appending the configuration named by TAG from the
# matching tagged config vars.
m4_defun([_LT_CONFIG],
[m4_require([_LT_FILEUTILS_DEFAULTS])dnl
_LT_CONFIG_SAVE_COMMANDS([
  m4_define([_LT_TAG], m4_if([$1], [], [C], [$1]))dnl
  m4_if(_LT_TAG, [C], [
    # See if we are running on zsh, and set the options that allow our
    # commands through without removal of \ escapes.
    if test -n "${ZSH_VERSION+set}"; then
      setopt NO_GLOB_SUBST
    fi

    cfgfile=${ofile}T
    trap "$RM \"$cfgfile\"; exit 1" 1 2 15
    $RM "$cfgfile"

    cat <<_LT_EOF >> "$cfgfile"
#! $SHELL
# Generated automatically by $as_me ($PACKAGE) $VERSION
# Libtool was configured on host `(hostname || uname -n) 2>/dev/null | sed 1q`:
# NOTE: Changes made to this file will be lost: look at ltmain.sh.

# Provide generalized library-building support services.
# Written by Gordon Matzigkeit, 1996

_LT_COPYING
_LT_LIBTOOL_TAGS

# Configured defaults for sys_lib_dlsearch_path munging.
: \${LT_SYS_LIBRARY_PATH="$configure_time_lt_sys_library_path"}

# ### BEGIN LIBTOOL CONFIG
_LT_LIBTOOL_CONFIG_VARS
_LT_LIBTOOL_TAG_VARS
# ### END LIBTOOL CONFIG

_LT_EOF

    cat <<'_LT_EOF' >> "$cfgfile"

# ### BEGIN FUNCTIONS SHARED WITH CONFIGURE

_LT_PREPARE_MUNGE_PATH_LIST
_LT_PREPARE_CC_BASENAME

# ### END FUNCTIONS SHARED WITH CONFIGURE

_LT_EOF

  case $host_os in
  aix3*)
    cat <<\_LT_EOF >> "$cfgfile"
# AIX sometimes has problems with the GCC collect2 program.  For some
# reason, if we set the COLLECT_NAMES environment variable, the problems
# vanish in a puff of smoke.
if test set != "${COLLECT_NAMES+set}"; then
  COLLECT_NAMES=
  export COLLECT_NAMES
fi
_LT_EOF
    ;;
  esac

  _LT_PROG_LTMAIN

  # We use sed instead of cat because bash on DJGPP gets confused if
  # if finds mixed CR/LF and LF-only lines.  Since sed operates in
  # text mode, it properly converts lines to CR/LF.  This bash problem
  # is reportedly fixed, but why not run on old versions too?
  sed '$q' "$ltmain" >> "$cfgfile" \
     || (rm -f "$cfgfile"; exit 1)

   mv -f "$cfgfile" "$ofile" ||
    (rm -f "$ofile" && cp "$cfgfile" "$ofile" && rm -f "$cfgfile")
  chmod +x "$ofile"
],
[cat <<_LT_EOF >> "$ofile"

dnl Unfortunately we have to use $1 here, since _LT_TAG is not expanded
dnl in a comment (ie after a #).
# ### BEGIN LIBTOOL TAG CONFIG: $1
_LT_LIBTOOL_TAG_VARS(_LT_TAG)
# ### END LIBTOOL TAG CONFIG: $1
_LT_EOF
])dnl /m4_if
],
[m4_if([$1], [], [
    PACKAGE='$PACKAGE'
    VERSION='$VERSION'
    RM='$RM'
    ofile='$ofile'], [])
])dnl /_LT_CONFIG_SAVE_COMMANDS
])# _LT_CONFIG


# LT_SUPPORTED_TAG(TAG)
# ---------------------
# Trace this macro to discover what tags are supported by the libtool
# --tag option, using:
#    autoconf --trace 'LT_SUPPORTED_TAG:$1'
AC_DEFUN([LT_SUPPORTED_TAG], [])


# C support is built-in for now
m4_define([_LT_LANG_C_enabled], [])
m4_define([_LT_TAGS], [])


# LT_LANG(LANG)
# -------------
# Enable libtool support for the given language if not already enabled.
AC_DEFUN([LT_LANG],
[AC_BEFORE([$0], [LT_OUTPUT])dnl
m4_case([$1],
  [C],			[_LT_LANG(C)],
  [C++],		[_LT_LANG(CXX)],
  [Go],			[_LT_LANG(GO)],
  [Java],		[_LT_LANG(GCJ)],
  [Fortran 77],		[_LT_LANG(F77)],
  [Fortran],		[_LT_LANG(FC)],
  [Windows Resource],	[_LT_LANG(RC)],
  [m4_ifdef([_LT_LANG_]$1[_CONFIG],
    [_LT_LANG($1)],
    [m4_fatal([$0: unsupported language: "$1"])])])dnl
])# LT_LANG


# _LT_LANG(LANGNAME)
# ------------------
m4_defun([_LT_LANG],
[m4_ifdef([_LT_LANG_]$1[_enabled], [],
  [LT_SUPPORTED_TAG([$1])dnl
  m4_append([_LT_TAGS], [$1 ])dnl
  m4_define([_LT_LANG_]$1[_enabled], [])dnl
  _LT_LANG_$1_CONFIG($1)])dnl
])# _LT_LANG


m4_ifndef([AC_PROG_GO], [
# NOTE: This macro has been submitted for inclusion into   #
#  GNU Autoconf as AC_PROG_GO.  When it is available in    #
#  a released version of Autoconf we should remove this    #
#  macro and use it instead.                               #
m4_defun([AC_PROG_GO],
[AC_LANG_PUSH(Go)dnl
AC_ARG_VAR([GOC],     [Go compiler command])dnl
AC_ARG_VAR([GOFLAGS], [Go compiler flags])dnl
_AC_ARG_VAR_LDFLAGS()dnl
AC_CHECK_TOOL(GOC, gccgo)
if test -z "$GOC"; then
  if test -n "$ac_tool_prefix"; then
    AC_CHECK_PROG(GOC, [${ac_tool_prefix}gccgo], [${ac_tool_prefix}gccgo])
  fi
fi
if test -z "$GOC"; then
  AC_CHECK_PROG(GOC, gccgo, gccgo, false)
fi
])#m4_defun
])#m4_ifndef


# _LT_LANG_DEFAULT_CONFIG
# -----------------------
m4_defun([_LT_LANG_DEFAULT_CONFIG],
[AC_PROVIDE_IFELSE([AC_PROG_CXX],
  [LT_LANG(CXX)],
  [m4_define([AC_PROG_CXX], defn([AC_PROG_CXX])[LT_LANG(CXX)])])

AC_PROVIDE_IFELSE([AC_PROG_F77],
  [LT_LANG(F77)],
  [m4_define([AC_PROG_F77], defn([AC_PROG_F77])[LT_LANG(F77)])])

AC_PROVIDE_IFELSE([AC_PROG_FC],
  [LT_LANG(FC)],
  [m4_define([AC_PROG_FC], defn([AC_PROG_FC])[LT_LANG(FC)])])

dnl The call to [A][M_PROG_GCJ] is quoted like that to stop aclocal
dnl pulling things in needlessly.
AC_PROVIDE_IFELSE([AC_PROG_GCJ],
  [LT_LANG(GCJ)],
  [AC_PROVIDE_IFELSE([A][M_PROG_GCJ],
    [LT_LANG(GCJ)],
    [AC_PROVIDE_IFELSE([LT_PROG_GCJ],
      [LT_LANG(GCJ)],
      [m4_ifdef([AC_PROG_GCJ],
	[m4_define([AC_PROG_GCJ], defn([AC_PROG_GCJ])[LT_LANG(GCJ)])])
       m4_ifdef([A][M_PROG_GCJ],
	[m4_define([A][M_PROG_GCJ], defn([A][M_PROG_GCJ])[LT_LANG(GCJ)])])
       m4_ifdef([LT_PROG_GCJ],
	[m4_define([LT_PROG_GCJ], defn([LT_PROG_GCJ])[LT_LANG(GCJ)])])])])])

AC_PROVIDE_IFELSE([AC_PROG_GO],
  [LT_LANG(GO)],
  [m4_define([AC_PROG_GO], defn([AC_PROG_GO])[LT_LANG(GO)])])

AC_PROVIDE_IFELSE([LT_PROG_RC],
  [LT_LANG(RC)],
  [m4_define([LT_PROG_RC], defn([LT_PROG_RC])[LT_LANG(RC)])])
])# _LT_LANG_DEFAULT_CONFIG

# Obsolete macros:
AU_DEFUN([AC_LIBTOOL_CXX], [LT_LANG(C++)])
AU_DEFUN([AC_LIBTOOL_F77], [LT_LANG(Fortran 77)])
AU_DEFUN([AC_LIBTOOL_FC], [LT_LANG(Fortran)])
AU_DEFUN([AC_LIBTOOL_GCJ], [LT_LANG(Java)])
AU_DEFUN([AC_LIBTOOL_RC], [LT_LANG(Windows Resource)])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_CXX], [])
dnl AC_DEFUN([AC_LIBTOOL_F77], [])
dnl AC_DEFUN([AC_LIBTOOL_FC], [])
dnl AC_DEFUN([AC_LIBTOOL_GCJ], [])
dnl AC_DEFUN([AC_LIBTOOL_RC], [])


# _LT_TAG_COMPILER
# ----------------
m4_defun([_LT_TAG_COMPILER],
[AC_REQUIRE([AC_PROG_CC])dnl

_LT_DECL([LTCC], [CC], [1], [A C compiler])dnl
_LT_DECL([LTCFLAGS], [CFLAGS], [1], [LTCC compiler flags])dnl
_LT_TAGDECL([CC], [compiler], [1], [A language specific compiler])dnl
_LT_TAGDECL([with_gcc], [GCC], [0], [Is the compiler the GNU compiler?])dnl

# If no C compiler was specified, use CC.
LTCC=${LTCC-"$CC"}

# If no C compiler flags were specified, use CFLAGS.
LTCFLAGS=${LTCFLAGS-"$CFLAGS"}

# Allow CC to be a program name with arguments.
compiler=$CC
])# _LT_TAG_COMPILER


# _LT_COMPILER_BOILERPLATE
# ------------------------
# Check for compiler boilerplate output or warnings with
# the simple compiler test code.
m4_defun([_LT_COMPILER_BOILERPLATE],
[m4_require([_LT_DECL_SED])dnl
ac_outfile=conftest.$ac_objext
echo "$lt_simple_compile_test_code" >conftest.$ac_ext
eval "$ac_compile" 2>&1 >/dev/null | $SED '/^$/d; /^ *+/d' >conftest.err
_lt_compiler_boilerplate=`cat conftest.err`
$RM conftest*
])# _LT_COMPILER_BOILERPLATE


# _LT_LINKER_BOILERPLATE
# ----------------------
# Check for linker boilerplate output or warnings with
# the simple link test code.
m4_defun([_LT_LINKER_BOILERPLATE],
[m4_require([_LT_DECL_SED])dnl
ac_outfile=conftest.$ac_objext
echo "$lt_simple_link_test_code" >conftest.$ac_ext
eval "$ac_link" 2>&1 >/dev/null | $SED '/^$/d; /^ *+/d' >conftest.err
_lt_linker_boilerplate=`cat conftest.err`
$RM -r conftest*
])# _LT_LINKER_BOILERPLATE

# _LT_REQUIRED_DARWIN_CHECKS
# -------------------------
m4_defun_once([_LT_REQUIRED_DARWIN_CHECKS],[
  case $host_os in
    rhapsody* | darwin*)
    AC_CHECK_TOOL([DSYMUTIL], [dsymutil], [:])
    AC_CHECK_TOOL([NMEDIT], [nmedit], [:])
    AC_CHECK_TOOL([LIPO], [lipo], [:])
    AC_CHECK_TOOL([OTOOL], [otool], [:])
    AC_CHECK_TOOL([OTOOL64], [otool64], [:])
    _LT_DECL([], [DSYMUTIL], [1],
      [Tool to manipulate archived DWARF debug symbol files on Mac OS X])
    _LT_DECL([], [NMEDIT], [1],
      [Tool to change global to local symbols on Mac OS X])
    _LT_DECL([], [LIPO], [1],
      [Tool to manipulate fat objects and archives on Mac OS X])
    _LT_DECL([], [OTOOL], [1],
      [ldd/readelf like tool for Mach-O binaries on Mac OS X])
    _LT_DECL([], [OTOOL64], [1],
      [ldd/readelf like tool for 64 bit Mach-O binaries on Mac OS X 10.4])

    AC_CACHE_CHECK([for -single_module linker flag],[lt_cv_apple_cc_single_mod],
      [lt_cv_apple_cc_single_mod=no
      if test -z "$LT_MULTI_MODULE"; then
	# By default we will add the -single_module flag. You can override
	# by either setting the environment variable LT_MULTI_MODULE
	# non-empty at configure time, or by adding -multi_module to the
	# link flags.
	rm -rf libconftest.dylib*
	echo "int foo(void){return 1;}" > conftest.c
	echo "$LTCC $LTCFLAGS $LDFLAGS -o libconftest.dylib \
-dynamiclib -Wl,-single_module conftest.c" >&AS_MESSAGE_LOG_FD
	$LTCC $LTCFLAGS $LDFLAGS -o libconftest.dylib \
	  -dynamiclib -Wl,-single_module conftest.c 2>conftest.err
        _lt_result=$?
	# If there is a non-empty error log, and "single_module"
	# appears in it, assume the flag caused a linker warning
        if test -s conftest.err && $GREP single_module conftest.err; then
	  cat conftest.err >&AS_MESSAGE_LOG_FD
	# Otherwise, if the output was created with a 0 exit code from
	# the compiler, it worked.
	elif test -f libconftest.dylib && test 0 = "$_lt_result"; then
	  lt_cv_apple_cc_single_mod=yes
	else
	  cat conftest.err >&AS_MESSAGE_LOG_FD
	fi
	rm -rf libconftest.dylib*
	rm -f conftest.*
      fi])

    AC_CACHE_CHECK([for -exported_symbols_list linker flag],
      [lt_cv_ld_exported_symbols_list],
      [lt_cv_ld_exported_symbols_list=no
      save_LDFLAGS=$LDFLAGS
      echo "_main" > conftest.sym
      LDFLAGS="$LDFLAGS -Wl,-exported_symbols_list,conftest.sym"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],[])],
	[lt_cv_ld_exported_symbols_list=yes],
	[lt_cv_ld_exported_symbols_list=no])
	LDFLAGS=$save_LDFLAGS
    ])

    AC_CACHE_CHECK([for -force_load linker flag],[lt_cv_ld_force_load],
      [lt_cv_ld_force_load=no
      cat > conftest.c << _LT_EOF
int forced_loaded() { return 2;}
_LT_EOF
      echo "$LTCC $LTCFLAGS -c -o conftest.o conftest.c" >&AS_MESSAGE_LOG_FD
      $LTCC $LTCFLAGS -c -o conftest.o conftest.c 2>&AS_MESSAGE_LOG_FD
      echo "$AR cru libconftest.a conftest.o" >&AS_MESSAGE_LOG_FD
      $AR cru libconftest.a conftest.o 2>&AS_MESSAGE_LOG_FD
      echo "$RANLIB libconftest.a" >&AS_MESSAGE_LOG_FD
      $RANLIB libconftest.a 2>&AS_MESSAGE_LOG_FD
      cat > conftest.c << _LT_EOF
int main() { return 0;}
_LT_EOF
      echo "$LTCC $LTCFLAGS $LDFLAGS -o conftest conftest.c -Wl,-force_load,./libconftest.a" >&AS_MESSAGE_LOG_FD
      $LTCC $LTCFLAGS $LDFLAGS -o conftest conftest.c -Wl,-force_load,./libconftest.a 2>conftest.err
      _lt_result=$?
      if test -s conftest.err && $GREP force_load conftest.err; then
	cat conftest.err >&AS_MESSAGE_LOG_FD
      elif test -f conftest && test 0 = "$_lt_result" && $GREP forced_load conftest >/dev/null 2>&1; then
	lt_cv_ld_force_load=yes
      else
	cat conftest.err >&AS_MESSAGE_LOG_FD
      fi
        rm -f conftest.err libconftest.a conftest conftest.c
        rm -rf conftest.dSYM
    ])
    case $host_os in
    rhapsody* | darwin1.[[012]])
      _lt_dar_allow_undefined='$wl-undefined ${wl}suppress' ;;
    darwin1.*)
      _lt_dar_allow_undefined='$wl-flat_namespace $wl-undefined ${wl}suppress' ;;
    darwin*) # darwin 5.x on
      # if running on 10.5 or later, the deployment target defaults
      # to the OS version, if on x86, and 10.4, the deployment
      # target defaults to 10.4. Don't you love it?
      case ${MACOSX_DEPLOYMENT_TARGET-10.0},$host in
	10.0,*86*-darwin8*|10.0,*-darwin[[91]]*)
	  _lt_dar_allow_undefined='$wl-undefined ${wl}dynamic_lookup' ;;
	10.[[012]][[,.]]*)
	  _lt_dar_allow_undefined='$wl-flat_namespace $wl-undefined ${wl}suppress' ;;
	10.*)
	  _lt_dar_allow_undefined='$wl-undefined ${wl}dynamic_lookup' ;;
      esac
    ;;
  esac
    if test yes = "$lt_cv_apple_cc_single_mod"; then
      _lt_dar_single_mod='$single_module'
    fi
    if test yes = "$lt_cv_ld_exported_symbols_list"; then
      _lt_dar_export_syms=' $wl-exported_symbols_list,$output_objdir/$libname-symbols.expsym'
    else
      _lt_dar_export_syms='~$NMEDIT -s $output_objdir/$libname-symbols.expsym $lib'
    fi
    if test : != "$DSYMUTIL" && test no = "$lt_cv_ld_force_load"; then
      _lt_dsymutil='~$DSYMUTIL $lib || :'
    else
      _lt_dsymutil=
    fi
    ;;
  esac
])


# _LT_DARWIN_LINKER_FEATURES([TAG])
# ---------------------------------
# Checks for linker and compiler features on darwin
m4_defun([_LT_DARWIN_LINKER_FEATURES],
[
  m4_require([_LT_REQUIRED_DARWIN_CHECKS])
  _LT_TAGVAR(archive_cmds_need_lc, $1)=no
  _LT_TAGVAR(hardcode_direct, $1)=no
  _LT_TAGVAR(hardcode_automatic, $1)=yes
  _LT_TAGVAR(hardcode_shlibpath_var, $1)=unsupported
  if test yes = "$lt_cv_ld_force_load"; then
    _LT_TAGVAR(whole_archive_flag_spec, $1)='`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience $wl-force_load,$conv\"; done; func_echo_all \"$new_convenience\"`'
    m4_case([$1], [F77], [_LT_TAGVAR(compiler_needs_object, $1)=yes],
                  [FC],  [_LT_TAGVAR(compiler_needs_object, $1)=yes])
  else
    _LT_TAGVAR(whole_archive_flag_spec, $1)=''
  fi
  _LT_TAGVAR(link_all_deplibs, $1)=yes
  _LT_TAGVAR(allow_undefined_flag, $1)=$_lt_dar_allow_undefined
  case $cc_basename in
     ifort*|nagfor*) _lt_dar_can_shared=yes ;;
     *) _lt_dar_can_shared=$GCC ;;
  esac
  if test yes = "$_lt_dar_can_shared"; then
    output_verbose_link_cmd=func_echo_all
    _LT_TAGVAR(archive_cmds, $1)="\$CC -dynamiclib \$allow_undefined_flag -o \$lib \$libobjs \$deplibs \$compiler_flags -install_name \$rpath/\$soname \$verstring $_lt_dar_single_mod$_lt_dsymutil"
    _LT_TAGVAR(module_cmds, $1)="\$CC \$allow_undefined_flag -o \$lib -bundle \$libobjs \$deplibs \$compiler_flags$_lt_dsymutil"
    _LT_TAGVAR(archive_expsym_cmds, $1)="sed 's|^|_|' < \$export_symbols > \$output_objdir/\$libname-symbols.expsym~\$CC -dynamiclib \$allow_undefined_flag -o \$lib \$libobjs \$deplibs \$compiler_flags -install_name \$rpath/\$soname \$verstring $_lt_dar_single_mod$_lt_dar_export_syms$_lt_dsymutil"
    _LT_TAGVAR(module_expsym_cmds, $1)="sed -e 's|^|_|' < \$export_symbols > \$output_objdir/\$libname-symbols.expsym~\$CC \$allow_undefined_flag -o \$lib -bundle \$libobjs \$deplibs \$compiler_flags$_lt_dar_export_syms$_lt_dsymutil"
    m4_if([$1], [CXX],
[   if test yes != "$lt_cv_apple_cc_single_mod"; then
      _LT_TAGVAR(archive_cmds, $1)="\$CC -r -keep_private_externs -nostdlib -o \$lib-master.o \$libobjs~\$CC -dynamiclib \$allow_undefined_flag -o \$lib \$lib-master.o \$deplibs \$compiler_flags -install_name \$rpath/\$soname \$verstring$_lt_dsymutil"
      _LT_TAGVAR(archive_expsym_cmds, $1)="sed 's|^|_|' < \$export_symbols > \$output_objdir/\$libname-symbols.expsym~\$CC -r -keep_private_externs -nostdlib -o \$lib-master.o \$libobjs~\$CC -dynamiclib \$allow_undefined_flag -o \$lib \$lib-master.o \$deplibs \$compiler_flags -install_name \$rpath/\$soname \$verstring$_lt_dar_export_syms$_lt_dsymutil"
    fi
],[])
  else
  _LT_TAGVAR(ld_shlibs, $1)=no
  fi
])

# _LT_SYS_MODULE_PATH_AIX([TAGNAME])
# ----------------------------------
# Links a minimal program and checks the executable
# for the system default hardcoded library path. In most cases,
# this is /usr/lib:/lib, but when the MPI compilers are used
# the location of the communication and MPI libs are included too.
# If we don't find anything, use the default library path according
# to the aix ld manual.
# Store the results from the different compilers for each TAGNAME.
# Allow to override them for all tags through lt_cv_aix_libpath.
m4_defun([_LT_SYS_MODULE_PATH_AIX],
[m4_require([_LT_DECL_SED])dnl
if test set = "${lt_cv_aix_libpath+set}"; then
  aix_libpath=$lt_cv_aix_libpath
else
  AC_CACHE_VAL([_LT_TAGVAR([lt_cv_aix_libpath_], [$1])],
  [AC_LINK_IFELSE([AC_LANG_PROGRAM],[
  lt_aix_libpath_sed='[
      /Import File Strings/,/^$/ {
	  /^0/ {
	      s/^0  *\([^ ]*\) *$/\1/
	      p
	  }
      }]'
  _LT_TAGVAR([lt_cv_aix_libpath_], [$1])=`dump -H conftest$ac_exeext 2>/dev/null | $SED -n -e "$lt_aix_libpath_sed"`
  # Check for a 64-bit object if we didn't find anything.
  if test -z "$_LT_TAGVAR([lt_cv_aix_libpath_], [$1])"; then
    _LT_TAGVAR([lt_cv_aix_libpath_], [$1])=`dump -HX64 conftest$ac_exeext 2>/dev/null | $SED -n -e "$lt_aix_libpath_sed"`
  fi],[])
  if test -z "$_LT_TAGVAR([lt_cv_aix_libpath_], [$1])"; then
    _LT_TAGVAR([lt_cv_aix_libpath_], [$1])=/usr/lib:/lib
  fi
  ])
  aix_libpath=$_LT_TAGVAR([lt_cv_aix_libpath_], [$1])
fi
])# _LT_SYS_MODULE_PATH_AIX


# _LT_SHELL_INIT(ARG)
# -------------------
m4_define([_LT_SHELL_INIT],
[m4_divert_text([M4SH-INIT], [$1
])])# _LT_SHELL_INIT



# _LT_PROG_ECHO_BACKSLASH
# -----------------------
# Find how we can fake an echo command that does not interpret backslash.
# In particular, with Autoconf 2.60 or later we add some code to the start
# of the generated configure script that will find a shell with a builtin
# printf (that we can use as an echo command).
m4_defun([_LT_PROG_ECHO_BACKSLASH],
[ECHO='\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
ECHO=$ECHO$ECHO$ECHO$ECHO$ECHO
ECHO=$ECHO$ECHO$ECHO$ECHO$ECHO$ECHO

AC_MSG_CHECKING([how to print strings])
# Test print first, because it will be a builtin if present.
if test "X`( print -r -- -n ) 2>/dev/null`" = X-n && \
   test "X`print -r -- $ECHO 2>/dev/null`" = "X$ECHO"; then
  ECHO='print -r --'
elif test "X`printf %s $ECHO 2>/dev/null`" = "X$ECHO"; then
  ECHO='printf %s\n'
else
  # Use this function as a fallback that always works.
  func_fallback_echo ()
  {
    eval 'cat <<_LTECHO_EOF
$[]1
_LTECHO_EOF'
  }
  ECHO='func_fallback_echo'
fi

# func_echo_all arg...
# Invoke $ECHO with all args, space-separated.
func_echo_all ()
{
    $ECHO "$*"
}

case $ECHO in
  printf*) AC_MSG_RESULT([printf]) ;;
  print*) AC_MSG_RESULT([print -r]) ;;
  *) AC_MSG_RESULT([cat]) ;;
esac

m4_ifdef([_AS_DETECT_SUGGESTED],
[_AS_DETECT_SUGGESTED([
  test -n "${ZSH_VERSION+set}${BASH_VERSION+set}" || (
    ECHO='\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
    ECHO=$ECHO$ECHO$ECHO$ECHO$ECHO
    ECHO=$ECHO$ECHO$ECHO$ECHO$ECHO$ECHO
    PATH=/empty FPATH=/empty; export PATH FPATH
    test "X`printf %s $ECHO`" = "X$ECHO" \
      || test "X`print -r -- $ECHO`" = "X$ECHO" )])])

_LT_DECL([], [SHELL], [1], [Shell to use when invoking shell scripts])
_LT_DECL([], [ECHO], [1], [An echo program that protects backslashes])
])# _LT_PROG_ECHO_BACKSLASH


# _LT_WITH_SYSROOT
# ----------------
AC_DEFUN([_LT_WITH_SYSROOT],
[AC_MSG_CHECKING([for sysroot])
AC_ARG_WITH([sysroot],
[AS_HELP_STRING([--with-sysroot@<:@=DIR@:>@],
  [Search for dependent libraries within DIR (or the compiler's sysroot
   if not specified).])],
[], [with_sysroot=no])

dnl lt_sysroot will always be passed unquoted.  We quote it here
dnl in case the user passed a directory name.
lt_sysroot=
case $with_sysroot in #(
 yes)
   if test yes = "$GCC"; then
     lt_sysroot=`$CC --print-sysroot 2>/dev/null`
   fi
   ;; #(
 /*)
   lt_sysroot=`echo "$with_sysroot" | sed -e "$sed_quote_subst"`
   ;; #(
 no|'')
   ;; #(
 *)
   AC_MSG_RESULT([$with_sysroot])
   AC_MSG_ERROR([The sysroot must be an absolute path.])
   ;;
esac

 AC_MSG_RESULT([${lt_sysroot:-no}])
_LT_DECL([], [lt_sysroot], [0], [The root where to search for ]dnl
[dependent libraries, and where our libraries should be installed.])])

# _LT_ENABLE_LOCK
# ---------------
m4_defun([_LT_ENABLE_LOCK],
[AC_ARG_ENABLE([libtool-lock],
  [AS_HELP_STRING([--disable-libtool-lock],
    [avoid locking (might break parallel builds)])])
test no = "$enable_libtool_lock" || enable_libtool_lock=yes

# Some flags need to be propagated to the compiler or linker for good
# libtool support.
case $host in
ia64-*-hpux*)
  # Find out what ABI is being produced by ac_compile, and set mode
  # options accordingly.
  echo 'int i;' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.$ac_objext` in
      *ELF-32*)
	HPUX_IA64_MODE=32
	;;
      *ELF-64*)
	HPUX_IA64_MODE=64
	;;
    esac
  fi
  rm -rf conftest*
  ;;
*-*-irix6*)
  # Find out what ABI is being produced by ac_compile, and set linker
  # options accordingly.
  echo '[#]line '$LINENO' "configure"' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    if test yes = "$lt_cv_prog_gnu_ld"; then
      case `/usr/bin/file conftest.$ac_objext` in
	*32-bit*)
	  LD="${LD-ld} -melf32bsmip"
	  ;;
	*N32*)
	  LD="${LD-ld} -melf32bmipn32"
	  ;;
	*64-bit*)
	  LD="${LD-ld} -melf64bmip"
	;;
      esac
    else
      case `/usr/bin/file conftest.$ac_objext` in
	*32-bit*)
	  LD="${LD-ld} -32"
	  ;;
	*N32*)
	  LD="${LD-ld} -n32"
	  ;;
	*64-bit*)
	  LD="${LD-ld} -64"
	  ;;
      esac
    fi
  fi
  rm -rf conftest*
  ;;

mips64*-*linux*)
  # Find out what ABI is being produced by ac_compile, and set linker
  # options accordingly.
  echo '[#]line '$LINENO' "configure"' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    emul=elf
    case `/usr/bin/file conftest.$ac_objext` in
      *32-bit*)
	emul="${emul}32"
	;;
      *64-bit*)
	emul="${emul}64"
	;;
    esac
    case `/usr/bin/file conftest.$ac_objext` in
      *MSB*)
	emul="${emul}btsmip"
	;;
      *LSB*)
	emul="${emul}ltsmip"
	;;
    esac
    case `/usr/bin/file conftest.$ac_objext` in
      *N32*)
	emul="${emul}n32"
	;;
    esac
    LD="${LD-ld} -m $emul"
  fi
  rm -rf conftest*
  ;;

x86_64-*kfreebsd*-gnu|x86_64-*linux*|powerpc*-*linux*| \
s390*-*linux*|s390*-*tpf*|sparc*-*linux*)
  # Find out what ABI is being produced by ac_compile, and set linker
  # options accordingly.  Note that the listed cases only cover the
  # situations where additional linker options are needed (such as when
  # doing 32-bit compilation for a host where ld defaults to 64-bit, or
  # vice versa); the common cases where no linker options are needed do
  # not appear in the list.
  echo 'int i;' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.o` in
      *32-bit*)
	case $host in
	  x86_64-*kfreebsd*-gnu)
	    LD="${LD-ld} -m elf_i386_fbsd"
	    ;;
	  x86_64-*linux*)
	    case `/usr/bin/file conftest.o` in
	      *x86-64*)
		LD="${LD-ld} -m elf32_x86_64"
		;;
	      *)
		LD="${LD-ld} -m elf_i386"
		;;
	    esac
	    ;;
	  powerpc64le-*linux*)
	    LD="${LD-ld} -m elf32lppclinux"
	    ;;
	  powerpc64-*linux*)
	    LD="${LD-ld} -m elf32ppclinux"
	    ;;
	  s390x-*linux*)
	    LD="${LD-ld} -m elf_s390"
	    ;;
	  sparc64-*linux*)
	    LD="${LD-ld} -m elf32_sparc"
	    ;;
	esac
	;;
      *64-bit*)
	case $host in
	  x86_64-*kfreebsd*-gnu)
	    LD="${LD-ld} -m elf_x86_64_fbsd"
	    ;;
	  x86_64-*linux*)
	    LD="${LD-ld} -m elf_x86_64"
	    ;;
	  powerpcle-*linux*)
	    LD="${LD-ld} -m elf64lppc"
	    ;;
	  powerpc-*linux*)
	    LD="${LD-ld} -m elf64ppc"
	    ;;
	  s390*-*linux*|s390*-*tpf*)
	    LD="${LD-ld} -m elf64_s390"
	    ;;
	  sparc*-*linux*)
	    LD="${LD-ld} -m elf64_sparc"
	    ;;
	esac
	;;
    esac
  fi
  rm -rf conftest*
  ;;

*-*-sco3.2v5*)
  # On SCO OpenServer 5, we need -belf to get full-featured binaries.
  SAVE_CFLAGS=$CFLAGS
  CFLAGS="$CFLAGS -belf"
  AC_CACHE_CHECK([whether the C compiler needs -belf], lt_cv_cc_needs_belf,
    [AC_LANG_PUSH(C)
     AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],[[]])],[lt_cv_cc_needs_belf=yes],[lt_cv_cc_needs_belf=no])
     AC_LANG_POP])
  if test yes != "$lt_cv_cc_needs_belf"; then
    # this is probably gcc 2.8.0, egcs 1.0 or newer; no need for -belf
    CFLAGS=$SAVE_CFLAGS
  fi
  ;;
*-*solaris*)
  # Find out what ABI is being produced by ac_compile, and set linker
  # options accordingly.
  echo 'int i;' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.o` in
    *64-bit*)
      case $lt_cv_prog_gnu_ld in
      yes*)
        case $host in
        i?86-*-solaris*|x86_64-*-solaris*)
          LD="${LD-ld} -m elf_x86_64"
          ;;
        sparc*-*-solaris*)
          LD="${LD-ld} -m elf64_sparc"
          ;;
        esac
        # GNU ld 2.21 introduced _sol2 emulations.  Use them if available.
        if ${LD-ld} -V | grep _sol2 >/dev/null 2>&1; then
          LD=${LD-ld}_sol2
        fi
        ;;
      *)
	if ${LD-ld} -64 -r -o conftest2.o conftest.o >/dev/null 2>&1; then
	  LD="${LD-ld} -64"
	fi
	;;
      esac
      ;;
    esac
  fi
  rm -rf conftest*
  ;;
esac

need_locks=$enable_libtool_lock
])# _LT_ENABLE_LOCK


# _LT_PROG_AR
# -----------
m4_defun([_LT_PROG_AR],
[AC_CHECK_TOOLS(AR, [ar], false)
: ${AR=ar}
: ${AR_FLAGS=cru}
_LT_DECL([], [AR], [1], [The archiver])
_LT_DECL([], [AR_FLAGS], [1], [Flags to create an archive])

AC_CACHE_CHECK([for archiver @FILE support], [lt_cv_ar_at_file],
  [lt_cv_ar_at_file=no
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
     [echo conftest.$ac_objext > conftest.lst
      lt_ar_try='$AR $AR_FLAGS libconftest.a @conftest.lst >&AS_MESSAGE_LOG_FD'
      AC_TRY_EVAL([lt_ar_try])
      if test 0 -eq "$ac_status"; then
	# Ensure the archiver fails upon bogus file names.
	rm -f conftest.$ac_objext libconftest.a
	AC_TRY_EVAL([lt_ar_try])
	if test 0 -ne "$ac_status"; then
          lt_cv_ar_at_file=@
        fi
      fi
      rm -f conftest.* libconftest.a
     ])
  ])

if test no = "$lt_cv_ar_at_file"; then
  archiver_list_spec=
else
  archiver_list_spec=$lt_cv_ar_at_file
fi
_LT_DECL([], [archiver_list_spec], [1],
  [How to feed a file listing to the archiver])
])# _LT_PROG_AR


# _LT_CMD_OLD_ARCHIVE
# -------------------
m4_defun([_LT_CMD_OLD_ARCHIVE],
[_LT_PROG_AR

AC_CHECK_TOOL(STRIP, strip, :)
test -z "$STRIP" && STRIP=:
_LT_DECL([], [STRIP], [1], [A symbol stripping program])

AC_CHECK_TOOL(RANLIB, ranlib, :)
test -z "$RANLIB" && RANLIB=:
_LT_DECL([], [RANLIB], [1],
    [Commands used to install an old-style archive])

# Determine commands to create old-style static archives.
old_archive_cmds='$AR $AR_FLAGS $oldlib$oldobjs'
old_postinstall_cmds='chmod 644 $oldlib'
old_postuninstall_cmds=

if test -n "$RANLIB"; then
  case $host_os in
  bitrig* | openbsd*)
    old_postinstall_cmds="$old_postinstall_cmds~\$RANLIB -t \$tool_oldlib"
    ;;
  *)
    old_postinstall_cmds="$old_postinstall_cmds~\$RANLIB \$tool_oldlib"
    ;;
  esac
  old_archive_cmds="$old_archive_cmds~\$RANLIB \$tool_oldlib"
fi

case $host_os in
  darwin*)
    lock_old_archive_extraction=yes ;;
  *)
    lock_old_archive_extraction=no ;;
esac
_LT_DECL([], [old_postinstall_cmds], [2])
_LT_DECL([], [old_postuninstall_cmds], [2])
_LT_TAGDECL([], [old_archive_cmds], [2],
    [Commands used to build an old-style archive])
_LT_DECL([], [lock_old_archive_extraction], [0],
    [Whether to use a lock for old archive extraction])
])# _LT_CMD_OLD_ARCHIVE


# _LT_COMPILER_OPTION(MESSAGE, VARIABLE-NAME, FLAGS,
#		[OUTPUT-FILE], [ACTION-SUCCESS], [ACTION-FAILURE])
# ----------------------------------------------------------------
# Check whether the given compiler option works
AC_DEFUN([_LT_COMPILER_OPTION],
[m4_require([_LT_FILEUTILS_DEFAULTS])dnl
m4_require([_LT_DECL_SED])dnl
AC_CACHE_CHECK([$1], [$2],
  [$2=no
   m4_if([$4], , [ac_outfile=conftest.$ac_objext], [ac_outfile=$4])
   echo "$lt_simple_compile_test_code" > conftest.$ac_ext
   lt_compiler_flag="$3"  ## exclude from sc_useless_quotes_in_assignment
   # Insert the option either (1) after the last *FLAGS variable, or
   # (2) before a word containing "conftest.", or (3) at the end.
   # Note that $ac_compile itself does not contain backslashes and begins
   # with a dollar sign (not a hyphen), so the echo should work correctly.
   # The option is referenced via a variable to avoid confusing sed.
   lt_compile=`echo "$ac_compile" | $SED \
   -e 's:.*FLAGS}\{0,1\} :&$lt_compiler_flag :; t' \
   -e 's: [[^ ]]*conftest\.: $lt_compiler_flag&:; t' \
   -e 's:$: $lt_compiler_flag:'`
   (eval echo "\"\$as_me:$LINENO: $lt_compile\"" >&AS_MESSAGE_LOG_FD)
   (eval "$lt_compile" 2>conftest.err)
   ac_status=$?
   cat conftest.err >&AS_MESSAGE_LOG_FD
   echo "$as_me:$LINENO: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   if (exit $ac_status) && test -s "$ac_outfile"; then
     # The compiler can only warn and ignore the option if not recognized
     # So say no if there are warnings other than the usual output.
     $ECHO "$_lt_compiler_boilerplate" | $SED '/^$/d' >conftest.exp
     $SED '/^$/d; /^ *+/d' conftest.err >conftest.er2
     if test ! -s conftest.er2 || diff conftest.exp conftest.er2 >/dev/null; then
       $2=yes
     fi
   fi
   $RM conftest*
])

if test yes = "[$]$2"; then
    m4_if([$5], , :, [$5])
else
    m4_if([$6], , :, [$6])
fi
])# _LT_COMPILER_OPTION

# Old name:
AU_ALIAS([AC_LIBTOOL_COMPILER_OPTION], [_LT_COMPILER_OPTION])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_COMPILER_OPTION], [])


# _LT_LINKER_OPTION(MESSAGE, VARIABLE-NAME, FLAGS,
#                  [ACTION-SUCCESS], [ACTION-FAILURE])
# ----------------------------------------------------
# Check whether the given linker option works
AC_DEFUN([_LT_LINKER_OPTION],
[m4_require([_LT_FILEUTILS_DEFAULTS])dnl
m4_require([_LT_DECL_SED])dnl
AC_CACHE_CHECK([$1], [$2],
  [$2=no
   save_LDFLAGS=$LDFLAGS
   LDFLAGS="$LDFLAGS $3"
   echo "$lt_simple_link_test_code" > conftest.$ac_ext
   if (eval $ac_link 2>conftest.err) && test -s conftest$ac_exeext; then
     # The linker can only warn and ignore the option if not recognized
     # So say no if there are warnings
     if test -s conftest.err; then
       # Append any errors to the config.log.
       cat conftest.err 1>&AS_MESSAGE_LOG_FD
       $ECHO "$_lt_linker_boilerplate" | $SED '/^$/d' > conftest.exp
       $SED '/^$/d; /^ *+/d' conftest.err >conftest.er2
       if diff conftest.exp conftest.er2 >/dev/null; then
         $2=yes
       fi
     else
       $2=yes
     fi
   fi
   $RM -r conftest*
   LDFLAGS=$save_LDFLAGS
])

if test yes = "[$]$2"; then
    m4_if([$4], , :, [$4])
else
    m4_if([$5], , :, [$5])
fi
])# _LT_LINKER_OPTION

# Old name:
AU_ALIAS([AC_LIBTOOL_LINKER_OPTION], [_LT_LINKER_OPTION])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_LINKER_OPTION], [])


# LT_CMD_MAX_LEN
#---------------
AC_DEFUN([LT_CMD_MAX_LEN],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
# find the maximum length of command line arguments
AC_MSG_CHECKING([the maximum length of command line arguments])
AC_CACHE_VAL([lt_cv_sys_max_cmd_len], [dnl
  i=0
  teststring=ABCD

  case $build_os in
  msdosdjgpp*)
    # On DJGPP, this test can blow up pretty badly due to problems in libc
    # (any single argument exceeding 2000 bytes causes a buffer overrun
    # during glob expansion).  Even if it were fixed, the result of this
    # check would be larger than it should be.
    lt_cv_sys_max_cmd_len=12288;    # 12K is about right
    ;;

  gnu*)
    # Under GNU Hurd, this test is not required because there is
    # no limit to the length of command line arguments.
    # Libtool will interpret -1 as no limit whatsoever
    lt_cv_sys_max_cmd_len=-1;
    ;;

  cygwin* | mingw* | cegcc*)
    # On Win9x/ME, this test blows up -- it succeeds, but takes
    # about 5 minutes as the teststring grows exponentially.
    # Worse, since 9x/ME are not pre-emptively multitasking,
    # you end up with a "frozen" computer, even though with patience
    # the test eventually succeeds (with a max line length of 256k).
    # Instead, let's just punt: use the minimum linelength reported by
    # all of the supported platforms: 8192 (on NT/2K/XP).
    lt_cv_sys_max_cmd_len=8192;
    ;;

  mint*)
    # On MiNT this can take a long time and run out of memory.
    lt_cv_sys_max_cmd_len=8192;
    ;;

  amigaos*)
    # On AmigaOS with pdksh, this test takes hours, literally.
    # So we just punt and use a minimum line length of 8192.
    lt_cv_sys_max_cmd_len=8192;
    ;;

  bitrig* | darwin* | dragonfly* | freebsd* | netbsd* | openbsd*)
    # This has been around since 386BSD, at least.  Likely further.
    if test -x /sbin/sysctl; then
      lt_cv_sys_max_cmd_len=`/sbin/sysctl -n kern.argmax`
    elif test -x /usr/sbin/sysctl; then
      lt_cv_sys_max_cmd_len=`/usr/sbin/sysctl -n kern.argmax`
    else
      lt_cv_sys_max_cmd_len=65536	# usable default for all BSDs
    fi
    # And add a safety zone
    lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \/ 4`
    lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \* 3`
    ;;

  interix*)
    # We know the value 262144 and hardcode it with a safety zone (like BSD)
    lt_cv_sys_max_cmd_len=196608
    ;;

  os2*)
    # The test takes a long time on OS/2.
    lt_cv_sys_max_cmd_len=8192
    ;;

  osf*)
    # Dr. Hans Ekkehard Plesser reports seeing a kernel panic running configure
    # due to this test when exec_disable_arg_limit is 1 on Tru64. It is not
    # nice to cause kernel panics so lets avoid the loop below.
    # First set a reasonable default.
    lt_cv_sys_max_cmd_len=16384
    #
    if test -x /sbin/sysconfig; then
      case `/sbin/sysconfig -q proc exec_disable_arg_limit` in
        *1*) lt_cv_sys_max_cmd_len=-1 ;;
      esac
    fi
    ;;
  sco3.2v5*)
    lt_cv_sys_max_cmd_len=102400
    ;;
  sysv5* | sco5v6* | sysv4.2uw2*)
    kargmax=`grep ARG_MAX /etc/conf/cf.d/stune 2>/dev/null`
    if test -n "$kargmax"; then
      lt_cv_sys_max_cmd_len=`echo $kargmax | sed 's/.*[[	 ]]//'`
    else
      lt_cv_sys_max_cmd_len=32768
    fi
    ;;
  *)
    lt_cv_sys_max_cmd_len=`(getconf ARG_MAX) 2> /dev/null`
    if test -n "$lt_cv_sys_max_cmd_len" && \
       test undefined != "$lt_cv_sys_max_cmd_len"; then
      lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \/ 4`
      lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \* 3`
    else
      # Make teststring a little bigger before we do anything with it.
      # a 1K string should be a reasonable start.
      for i in 1 2 3 4 5 6 7 8; do
        teststring=$teststring$teststring
      done
      SHELL=${SHELL-${CONFIG_SHELL-/bin/sh}}
      # If test is not a shell built-in, we'll probably end up computing a
      # maximum length that is only half of the actual maximum length, but
      # we can't tell.
      while { test X`env echo "$teststring$teststring" 2>/dev/null` \
	         = "X$teststring$teststring"; } >/dev/null 2>&1 &&
	      test 17 != "$i" # 1/2 MB should be enough
      do
        i=`expr $i + 1`
        teststring=$teststring$teststring
      done
      # Only check the string length outside the loop.
      lt_cv_sys_max_cmd_len=`expr "X$teststring" : ".*" 2>&1`
      teststring=
      # Add a significant safety factor because C++ compilers can tack on
      # massive amounts of additional arguments before passing them to the
      # linker.  It appears as though 1/2 is a usable value.
      lt_cv_sys_max_cmd_len=`expr $lt_cv_sys_max_cmd_len \/ 2`
    fi
    ;;
  esac
])
if test -n "$lt_cv_sys_max_cmd_len"; then
  AC_MSG_RESULT($lt_cv_sys_max_cmd_len)
else
  AC_MSG_RESULT(none)
fi
max_cmd_len=$lt_cv_sys_max_cmd_len
_LT_DECL([], [max_cmd_len], [0],
    [What is the maximum length of a command?])
])# LT_CMD_MAX_LEN

# Old name:
AU_ALIAS([AC_LIBTOOL_SYS_MAX_CMD_LEN], [LT_CMD_MAX_LEN])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_SYS_MAX_CMD_LEN], [])


# _LT_HEADER_DLFCN
# ----------------
m4_defun([_LT_HEADER_DLFCN],
[AC_CHECK_HEADERS([dlfcn.h], [], [], [AC_INCLUDES_DEFAULT])dnl
])# _LT_HEADER_DLFCN


# _LT_TRY_DLOPEN_SELF (ACTION-IF-TRUE, ACTION-IF-TRUE-W-USCORE,
#                      ACTION-IF-FALSE, ACTION-IF-CROSS-COMPILING)
# ----------------------------------------------------------------
m4_defun([_LT_TRY_DLOPEN_SELF],
[m4_require([_LT_HEADER_DLFCN])dnl
if test yes = "$cross_compiling"; then :
  [$4]
else
  lt_dlunknown=0; lt_dlno_uscore=1; lt_dlneed_uscore=2
  lt_status=$lt_dlunknown
  cat > conftest.$ac_ext <<_LT_EOF
[#line $LINENO "configure"
#include "confdefs.h"

#if HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#include <stdio.h>

#ifdef RTLD_GLOBAL
#  define LT_DLGLOBAL		RTLD_GLOBAL
#else
#  ifdef DL_GLOBAL
#    define LT_DLGLOBAL		DL_GLOBAL
#  else
#    define LT_DLGLOBAL		0
#  endif
#endif

/* We may have to define LT_DLLAZY_OR_NOW in the command line if we
   find out it does not work in some platform. */
#ifndef LT_DLLAZY_OR_NOW
#  ifdef RTLD_LAZY
#    define LT_DLLAZY_OR_NOW		RTLD_LAZY
#  else
#    ifdef DL_LAZY
#      define LT_DLLAZY_OR_NOW		DL_LAZY
#    else
#      ifdef RTLD_NOW
#        define LT_DLLAZY_OR_NOW	RTLD_NOW
#      else
#        ifdef DL_NOW
#          define LT_DLLAZY_OR_NOW	DL_NOW
#        else
#          define LT_DLLAZY_OR_NOW	0
#        endif
#      endif
#    endif
#  endif
#endif

/* When -fvisibility=hidden is used, assume the code has been annotated
   correspondingly for the symbols needed.  */
#if defined __GNUC__ && (((__GNUC__ == 3) && (__GNUC_MINOR__ >= 3)) || (__GNUC__ > 3))
int fnord () __attribute__((visibility("default")));
#endif

int fnord () { return 42; }
int main ()
{
  void *self = dlopen (0, LT_DLGLOBAL|LT_DLLAZY_OR_NOW);
  int status = $lt_dlunknown;

  if (self)
    {
      if (dlsym (self,"fnord"))       status = $lt_dlno_uscore;
      else
        {
	  if (dlsym( self,"_fnord"))  status = $lt_dlneed_uscore;
          else puts (dlerror ());
	}
      /* dlclose (self); */
    }
  else
    puts (dlerror ());

  return status;
}]
_LT_EOF
  if AC_TRY_EVAL(ac_link) && test -s "conftest$ac_exeext" 2>/dev/null; then
    (./conftest; exit; ) >&AS_MESSAGE_LOG_FD 2>/dev/null
    lt_status=$?
    case x$lt_status in
      x$lt_dlno_uscore) $1 ;;
      x$lt_dlneed_uscore) $2 ;;
      x$lt_dlunknown|x*) $3 ;;
    esac
  else :
    # compilation failed
    $3
  fi
fi
rm -fr conftest*
])# _LT_TRY_DLOPEN_SELF


# LT_SYS_DLOPEN_SELF
# ------------------
AC_DEFUN([LT_SYS_DLOPEN_SELF],
[m4_require([_LT_HEADER_DLFCN])dnl
if test yes != "$enable_dlopen"; then
  enable_dlopen=unknown
  enable_dlopen_self=unknown
  enable_dlopen_self_static=unknown
else
  lt_cv_dlopen=no
  lt_cv_dlopen_libs=

  case $host_os in
  beos*)
    lt_cv_dlopen=load_add_on
    lt_cv_dlopen_libs=
    lt_cv_dlopen_self=yes
    ;;

  mingw* | pw32* | cegcc*)
    lt_cv_dlopen=LoadLibrary
    lt_cv_dlopen_libs=
    ;;

  cygwin*)
    lt_cv_dlopen=dlopen
    lt_cv_dlopen_libs=
    ;;

  darwin*)
    # if libdl is installed we need to link against it
    AC_CHECK_LIB([dl], [dlopen],
		[lt_cv_dlopen=dlopen lt_cv_dlopen_libs=-ldl],[
    lt_cv_dlopen=dyld
    lt_cv_dlopen_libs=
    lt_cv_dlopen_self=yes
    ])
    ;;

  tpf*)
    # Don't try to run any link tests for TPF.  We know it's impossible
    # because TPF is a cross-compiler, and we know how we open DSOs.
    lt_cv_dlopen=dlopen
    lt_cv_dlopen_libs=
    lt_cv_dlopen_self=no
    ;;

  *)
    AC_CHECK_FUNC([shl_load],
	  [lt_cv_dlopen=shl_load],
      [AC_CHECK_LIB([dld], [shl_load],
	    [lt_cv_dlopen=shl_load lt_cv_dlopen_libs=-ldld],
	[AC_CHECK_FUNC([dlopen],
	      [lt_cv_dlopen=dlopen],
	  [AC_CHECK_LIB([dl], [dlopen],
		[lt_cv_dlopen=dlopen lt_cv_dlopen_libs=-ldl],
	    [AC_CHECK_LIB([svld], [dlopen],
		  [lt_cv_dlopen=dlopen lt_cv_dlopen_libs=-lsvld],
	      [AC_CHECK_LIB([dld], [dld_link],
		    [lt_cv_dlopen=dld_link lt_cv_dlopen_libs=-ldld])
	      ])
	    ])
	  ])
	])
      ])
    ;;
  esac

  if test no = "$lt_cv_dlopen"; then
    enable_dlopen=no
  else
    enable_dlopen=yes
  fi

  case $lt_cv_dlopen in
  dlopen)
    save_CPPFLAGS=$CPPFLAGS
    test yes = "$ac_cv_header_dlfcn_h" && CPPFLAGS="$CPPFLAGS -DHAVE_DLFCN_H"

    save_LDFLAGS=$LDFLAGS
    wl=$lt_prog_compiler_wl eval LDFLAGS=\"\$LDFLAGS $export_dynamic_flag_spec\"

    save_LIBS=$LIBS
    LIBS="$lt_cv_dlopen_libs $LIBS"

    AC_CACHE_CHECK([whether a program can dlopen itself],
	  lt_cv_dlopen_self, [dnl
	  _LT_TRY_DLOPEN_SELF(
	    lt_cv_dlopen_self=yes, lt_cv_dlopen_self=yes,
	    lt_cv_dlopen_self=no, lt_cv_dlopen_self=cross)
    ])

    if test yes = "$lt_cv_dlopen_self"; then
      wl=$lt_prog_compiler_wl eval LDFLAGS=\"\$LDFLAGS $lt_prog_compiler_static\"
      AC_CACHE_CHECK([whether a statically linked program can dlopen itself],
	  lt_cv_dlopen_self_static, [dnl
	  _LT_TRY_DLOPEN_SELF(
	    lt_cv_dlopen_self_static=yes, lt_cv_dlopen_self_static=yes,
	    lt_cv_dlopen_self_static=no,  lt_cv_dlopen_self_static=cross)
      ])
    fi

    CPPFLAGS=$save_CPPFLAGS
    LDFLAGS=$save_LDFLAGS
    LIBS=$save_LIBS
    ;;
  esac

  case $lt_cv_dlopen_self in
  yes|no) enable_dlopen_self=$lt_cv_dlopen_self ;;
  *) enable_dlopen_self=unknown ;;
  esac

  case $lt_cv_dlopen_self_static in
  yes|no) enable_dlopen_self_static=$lt_cv_dlopen_self_static ;;
  *) enable_dlopen_self_static=unknown ;;
  esac
fi
_LT_DECL([dlopen_support], [enable_dlopen], [0],
	 [Whether dlopen is supported])
_LT_DECL([dlopen_self], [enable_dlopen_self], [0],
	 [Whether dlopen of programs is supported])
_LT_DECL([dlopen_self_static], [enable_dlopen_self_static], [0],
	 [Whether dlopen of statically linked programs is supported])
])# LT_SYS_DLOPEN_SELF

# Old name:
AU_ALIAS([AC_LIBTOOL_DLOPEN_SELF], [LT_SYS_DLOPEN_SELF])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_DLOPEN_SELF], [])


# _LT_COMPILER_C_O([TAGNAME])
# ---------------------------
# Check to see if options -c and -o are simultaneously supported by compiler.
# This macro does not hard code the compiler like AC_PROG_CC_C_O.
m4_defun([_LT_COMPILER_C_O],
[m4_require([_LT_DECL_SED])dnl
m4_require([_LT_FILEUTILS_DEFAULTS])dnl
m4_require([_LT_TAG_COMPILER])dnl
AC_CACHE_CHECK([if $compiler supports -c -o file.$ac_objext],
  [_LT_TAGVAR(lt_cv_prog_compiler_c_o, $1)],
  [_LT_TAGVAR(lt_cv_prog_compiler_c_o, $1)=no
   $RM -r conftest 2>/dev/null
   mkdir conftest
   cd conftest
   mkdir out
   echo "$lt_simple_compile_test_code" > conftest.$ac_ext

   lt_compiler_flag="-o out/conftest2.$ac_objext"
   # Insert the option either (1) after the last *FLAGS variable, or
   # (2) before a word containing "conftest.", or (3) at the end.
   # Note that $ac_compile itself does not contain backslashes and begins
   # with a dollar sign (not a hyphen), so the echo should work correctly.
   lt_compile=`echo "$ac_compile" | $SED \
   -e 's:.*FLAGS}\{0,1\} :&$lt_compiler_flag :; t' \
   -e 's: [[^ ]]*conftest\.: $lt_compiler_flag&:; t' \
   -e 's:$: $lt_compiler_flag:'`
   (eval echo "\"\$as_me:$LINENO: $lt_compile\"" >&AS_MESSAGE_LOG_FD)
   (eval "$lt_compile" 2>out/conftest.err)
   ac_status=$?
   cat out/conftest.err >&AS_MESSAGE_LOG_FD
   echo "$as_me:$LINENO: \$? = $ac_status" >&AS_MESSAGE_LOG_FD
   if (exit $ac_status) && test -s out/conftest2.$ac_objext
   then
     # The compiler can only warn and ignore the option if not recognized
     # So say no if there are warnings
     $ECHO "$_lt_compiler_boilerplate" | $SED '/^$/d' > out/conftest.exp
     $SED '/^$/d; /^ *+/d' out/conftest.err >out/conftest.er2
     if test ! -s out/conftest.er2 || diff out/conftest.exp out/conftest.er2 >/dev/null; then
       _LT_TAGVAR(lt_cv_prog_compiler_c_o, $1)=yes
     fi
   fi
   chmod u+w . 2>&AS_MESSAGE_LOG_FD
   $RM conftest*
   # SGI C++ compiler will create directory out/ii_files/ for
   # template instantiation
   test -d out/ii_files && $RM out/ii_files/* && rmdir out/ii_files
   $RM out/* && rmdir out
   cd ..
   $RM -r conftest
   $RM conftest*
])
_LT_TAGDECL([compiler_c_o], [lt_cv_prog_compiler_c_o], [1],
	[Does compiler simultaneously support -c and -o options?])
])# _LT_COMPILER_C_O


# _LT_COMPILER_FILE_LOCKS([TAGNAME])
# ----------------------------------
# Check to see if we can do hard links to lock some files if needed
m4_defun([_LT_COMPILER_FILE_LOCKS],
[m4_require([_LT_ENABLE_LOCK])dnl
m4_require([_LT_FILEUTILS_DEFAULTS])dnl
_LT_COMPILER_C_O([$1])

hard_links=nottested
if test no = "$_LT_TAGVAR(lt_cv_prog_compiler_c_o, $1)" && test no != "$need_locks"; then
  # do not overwrite the value of need_locks provided by the user
  AC_MSG_CHECKING([if we can lock with hard links])
  hard_links=yes
  $RM conftest*
  ln conftest.a conftest.b 2>/dev/null && hard_links=no
  touch conftest.a
  ln conftest.a conftest.b 2>&5 || hard_links=no
  ln conftest.a conftest.b 2>/dev/null && hard_links=no
  AC_MSG_RESULT([$hard_links])
  if test no = "$hard_links"; then
    AC_MSG_WARN(['$CC' does not support '-c -o', so 'make -j' may be unsafe])
    need_locks=warn
  fi
else
  need_locks=no
fi
_LT_DECL([], [need_locks], [1], [Must we lock files when doing compilation?])
])# _LT_COMPILER_FILE_LOCKS


# _LT_CHECK_OBJDIR
# ----------------
m4_defun([_LT_CHECK_OBJDIR],
[AC_CACHE_CHECK([for objdir], [lt_cv_objdir],
[rm -f .libs 2>/dev/null
mkdir .libs 2>/dev/null
if test -d .libs; then
  lt_cv_objdir=.libs
else
  # MS-DOS does not allow filenames that begin with a dot.
  lt_cv_objdir=_libs
fi
rmdir .libs 2>/dev/null])
objdir=$lt_cv_objdir
_LT_DECL([], [objdir], [0],
         [The name of the directory that contains temporary libtool files])dnl
m4_pattern_allow([LT_OBJDIR])dnl
AC_DEFINE_UNQUOTED([LT_OBJDIR], "$lt_cv_objdir/",
  [Define to the sub-directory where libtool stores uninstalled libraries.])
])# _LT_CHECK_OBJDIR


# _LT_LINKER_HARDCODE_LIBPATH([TAGNAME])
# --------------------------------------
# Check hardcoding attributes.
m4_defun([_LT_LINKER_HARDCODE_LIBPATH],
[AC_MSG_CHECKING([how to hardcode library paths into programs])
_LT_TAGVAR(hardcode_action, $1)=
if test -n "$_LT_TAGVAR(hardcode_libdir_flag_spec, $1)" ||
   test -n "$_LT_TAGVAR(runpath_var, $1)" ||
   test yes = "$_LT_TAGVAR(hardcode_automatic, $1)"; then

  # We can hardcode non-existent directories.
  if test no != "$_LT_TAGVAR(hardcode_direct, $1)" &&
     # If the only mechanism to avoid hardcoding is shlibpath_var, we
     # have to relink, otherwise we might link with an installed library
     # when we should be linking with a yet-to-be-installed one
     ## test no != "$_LT_TAGVAR(hardcode_shlibpath_var, $1)" &&
     test no != "$_LT_TAGVAR(hardcode_minus_L, $1)"; then
    # Linking always hardcodes the temporary library directory.
    _LT_TAGVAR(hardcode_action, $1)=relink
  else
    # We can link without hardcoding, and we can hardcode nonexisting dirs.
    _LT_TAGVAR(hardcode_action, $1)=immediate
  fi
else
  # We cannot hardcode anything, or else we can only hardcode existing
  # directories.
  _LT_TAGVAR(hardcode_action, $1)=unsupported
fi
AC_MSG_RESULT([$_LT_TAGVAR(hardcode_action, $1)])

if test relink = "$_LT_TAGVAR(hardcode_action, $1)" ||
   test yes = "$_LT_TAGVAR(inherit_rpath, $1)"; then
  # Fast installation is not supported
  enable_fast_install=no
elif test yes = "$shlibpath_overrides_runpath" ||
     test no = "$enable_shared"; then
  # Fast installation is not necessary
  enable_fast_install=needless
fi
_LT_TAGDECL([], [hardcode_action], [0],
    [How to hardcode a shared library path into an executable])
])# _LT_LINKER_HARDCODE_LIBPATH


# _LT_CMD_STRIPLIB
# ----------------
m4_defun([_LT_CMD_STRIPLIB],
[m4_require([_LT_DECL_EGREP])
striplib=
old_striplib=
AC_MSG_CHECKING([whether stripping libraries is possible])
if test -n "$STRIP" && $STRIP -V 2>&1 | $GREP "GNU strip" >/dev/null; then
  test -z "$old_striplib" && old_striplib="$STRIP --strip-debug"
  test -z "$striplib" && striplib="$STRIP --strip-unneeded"
  AC_MSG_RESULT([yes])
else
# FIXME - insert some real tests, host_os isn't really good enough
  case $host_os in
  darwin*)
    if test -n "$STRIP"; then
      striplib="$STRIP -x"
      old_striplib="$STRIP -S"
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
    fi
    ;;
  *)
    AC_MSG_RESULT([no])
    ;;
  esac
fi
_LT_DECL([], [old_striplib], [1], [Commands to strip libraries])
_LT_DECL([], [striplib], [1])
])# _LT_CMD_STRIPLIB


# _LT_PREPARE_MUNGE_PATH_LIST
# ---------------------------
# Make sure func_munge_path_list() is defined correctly.
m4_defun([_LT_PREPARE_MUNGE_PATH_LIST],
[[# func_munge_path_list VARIABLE PATH
# -----------------------------------
# VARIABLE is name of variable containing _space_ separated list of
# directories to be munged by the contents of PATH, which is string
# having a format:
# "DIR[:DIR]:"
#       string "DIR[ DIR]" will be prepended to VARIABLE
# ":DIR[:DIR]"
#       string "DIR[ DIR]" will be appended to VARIABLE
# "DIRP[:DIRP]::[DIRA:]DIRA"
#       string "DIRP[ DIRP]" will be prepended to VARIABLE and string
#       "DIRA[ DIRA]" will be appended to VARIABLE
# "DIR[:DIR]"
#       VARIABLE will be replaced by "DIR[ DIR]"
func_munge_path_list ()
{
    case x@S|@2 in
    x)
        ;;
    *:)
        eval @S|@1=\"`$ECHO @S|@2 | $SED 's/:/ /g'` \@S|@@S|@1\"
        ;;
    x:*)
        eval @S|@1=\"\@S|@@S|@1 `$ECHO @S|@2 | $SED 's/:/ /g'`\"
        ;;
    *::*)
        eval @S|@1=\"\@S|@@S|@1\ `$ECHO @S|@2 | $SED -e 's/.*:://' -e 's/:/ /g'`\"
        eval @S|@1=\"`$ECHO @S|@2 | $SED -e 's/::.*//' -e 's/:/ /g'`\ \@S|@@S|@1\"
        ;;
    *)
        eval @S|@1=\"`$ECHO @S|@2 | $SED 's/:/ /g'`\"
        ;;
    esac
}
]])# _LT_PREPARE_PATH_LIST


# _LT_SYS_DYNAMIC_LINKER([TAG])
# -----------------------------
# PORTME Fill in your ld.so characteristics
m4_defun([_LT_SYS_DYNAMIC_LINKER],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
m4_require([_LT_DECL_EGREP])dnl
m4_require([_LT_FILEUTILS_DEFAULTS])dnl
m4_require([_LT_DECL_OBJDUMP])dnl
m4_require([_LT_DECL_SED])dnl
m4_require([_LT_CHECK_SHELL_FEATURES])dnl
m4_require([_LT_PREPARE_MUNGE_PATH_LIST])dnl
AC_MSG_CHECKING([dynamic linker characteristics])
m4_if([$1],
	[], [
if test yes = "$GCC"; then
  case $host_os in
    darwin*) lt_awk_arg='/^libraries:/,/LR/' ;;
    *) lt_awk_arg='/^libraries:/' ;;
  esac
  case $host_os in
    mingw* | cegcc*) lt_sed_strip_eq='s|=\([[A-Za-z]]:\)|\1|g' ;;
    *) lt_sed_strip_eq='s|=/|/|g' ;;
  esac
  lt_search_path_spec=`$CC -print-search-dirs | awk $lt_awk_arg | $SED -e "s/^libraries://" -e $lt_sed_strip_eq`
  case $lt_search_path_spec in
  *\;*)
    # if the path contains ";" then we assume it to be the separator
    # otherwise default to the standard path separator (i.e. ":") - it is
    # assumed that no part of a normal pathname contains ";" but that should
    # okay in the real world where ";" in dirpaths is itself problematic.
    lt_search_path_spec=`$ECHO "$lt_search_path_spec" | $SED 's/;/ /g'`
    ;;
  *)
    lt_search_path_spec=`$ECHO "$lt_search_path_spec" | $SED "s/$PATH_SEPARATOR/ /g"`
    ;;
  esac
  # Ok, now we have the path, separated by spaces, we can step through it
  # and add multilib dir if necessary...
  lt_tmp_lt_search_path_spec=
  lt_multi_os_dir=/`$CC $CPPFLAGS $CFLAGS $LDFLAGS -print-multi-os-directory 2>/dev/null`
  # ...but if some path component already ends with the multilib dir we assume
  # that all is fine and trust -print-search-dirs as is (GCC 4.2? or newer).
  case "$lt_multi_os_dir; $lt_search_path_spec " in
  "/; "* | "/.; "* | "/./; "* | *"$lt_multi_os_dir "* | *"$lt_multi_os_dir/ "*)
    lt_multi_os_dir=
    ;;
  esac
  for lt_sys_path in $lt_search_path_spec; do
    if test -d "$lt_sys_path$lt_multi_os_dir"; then
      lt_tmp_lt_search_path_spec="$lt_tmp_lt_search_path_spec $lt_sys_path$lt_multi_os_dir"
    elif test -n "$lt_multi_os_dir"; then
      test -d "$lt_sys_path" && \
	lt_tmp_lt_search_path_spec="$lt_tmp_lt_search_path_spec $lt_sys_path"
    fi
  done
  lt_search_path_spec=`$ECHO "$lt_tmp_lt_search_path_spec" | awk '
BEGIN {RS = " "; FS = "/|\n";} {
  lt_foo = "";
  lt_count = 0;
  for (lt_i = NF; lt_i > 0; lt_i--) {
    if ($lt_i != "" && $lt_i != ".") {
      if ($lt_i == "..") {
        lt_count++;
      } else {
        if (lt_count == 0) {
          lt_foo = "/" $lt_i lt_foo;
        } else {
          lt_count--;
        }
      }
    }
  }
  if (lt_foo != "") { lt_freq[[lt_foo]]++; }
  if (lt_freq[[lt_foo]] == 1) { print lt_foo; }
}'`
  # AWK program above erroneously prepends '/' to C:/dos/paths
  # for these hosts.
  case $host_os in
    mingw* | cegcc*) lt_search_path_spec=`$ECHO "$lt_search_path_spec" |\
      $SED 's|/\([[A-Za-z]]:\)|\1|g'` ;;
  esac
  sys_lib_search_path_spec=`$ECHO "$lt_search_path_spec" | $lt_NL2SP`
else
  sys_lib_search_path_spec="/lib /usr/lib /usr/local/lib"
fi])
library_names_spec=
libname_spec='lib$name'
soname_spec=
shrext_cmds=.so
postinstall_cmds=
postuninstall_cmds=
finish_cmds=
finish_eval=
shlibpath_var=
shlibpath_overrides_runpath=unknown
version_type=none
dynamic_linker="$host_os ld.so"
sys_lib_dlsearch_path_spec="/lib /usr/lib"
need_lib_prefix=unknown
hardcode_into_libs=no

# when you set need_version to no, make sure it does not cause -set_version
# flags to be left without arguments
need_version=unknown

AC_ARG_VAR([LT_SYS_LIBRARY_PATH],
[User-defined run-time library search path.])

case $host_os in
aix3*)
  version_type=linux # correct to gnu/linux during the next big refactor
  library_names_spec='$libname$release$shared_ext$versuffix $libname.a'
  shlibpath_var=LIBPATH

  # AIX 3 has no versioning support, so we append a major version to the name.
  soname_spec='$libname$release$shared_ext$major'
  ;;

aix[[4-9]]*)
  version_type=linux # correct to gnu/linux during the next big refactor
  need_lib_prefix=no
  need_version=no
  hardcode_into_libs=yes
  if test ia64 = "$host_cpu"; then
    # AIX 5 supports IA64
    library_names_spec='$libname$release$shared_ext$major $libname$release$shared_ext$versuffix $libname$shared_ext'
    shlibpath_var=LD_LIBRARY_PATH
  else
    # With GCC up to 2.95.x, collect2 would create an import file
    # for dependence libraries.  The import file would start with
    # the line '#! .'.  This would cause the generated library to
    # depend on '.', always an invalid library.  This was fixed in
    # development snapshots of GCC prior to 3.0.
    case $host_os in
      aix4 | aix4.[[01]] | aix4.[[01]].*)
      if { echo '#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 97)'
	   echo ' yes '
	   echo '#endif'; } | $CC -E - | $GREP yes > /dev/null; then
	:
      else
	can_build_shared=no
      fi
      ;;
    esac
    # Using Import Files as archive members, it is possible to support
    # filename-based versioning of shared library archives on AIX. While
    # this would work for both with and without runtime linking, it will
    # prevent static linking of such archives. So we do filename-based
    # shared library versioning with .so extension only, which is used
    # when both runtime linking and shared linking is enabled.
    # Unfortunately, runtime linking may impact performance, so we do
    # not want this to be the default eventually. Also, we use the
    # versioned .so libs for executables only if there is the -brtl
    # linker flag in LDFLAGS as well, or --with-aix-soname=svr4 only.
    # To allow for filename-based versioning support, we need to create
    # libNAME.so.V as an archive file, containing:
    # *) an Import File, referring to the versioned filename of the
    #    archive as well as the shared archive member, telling the
    #    bitwidth (32 or 64) of that shared object, and providing the
    #    list of exported symbols of that shared object, eventually
    #    decorated with the 'weak' keyword
    # *) the shared object with the F_LOADONLY flag set, to really avoid
    #    it being seen by the linker.
    # At run time we better use the real file rather than another symlink,
    # but for link time we create the symlink libNAME.so -> libNAME.so.V

    case $with_aix_soname,$aix_use_runtimelinking in
    # AIX (on Power*) has no versioning support, so currently we cannot hardcode correct
    # soname into executable. Probably we can add versioning support to
    # collect2, so additional links can be useful in future.
    aix,yes) # traditional libtool
      dynamic_linker='AIX unversionable lib.so'
      # If using run time linking (on AIX 4.2 or later) use lib<name>.so
      # instead of lib<name>.a to let people know that these are not
      # typical AIX shared libraries.
      library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
      ;;
    aix,no) # traditional AIX only
      dynamic_linker='AIX lib.a[(]lib.so.V[)]'
      # We preserve .a as extension for shared libraries through AIX4.2
      # and later when we are not doing run time linking.
      library_names_spec='$libname$release.a $libname.a'
      soname_spec='$libname$release$shared_ext$major'
      ;;
    svr4,*) # full svr4 only
      dynamic_linker="AIX lib.so.V[(]$shared_archive_member_spec.o[)]"
      library_names_spec='$libname$release$shared_ext$major $libname$shared_ext'
      # We do not specify a path in Import Files, so LIBPATH fires.
      shlibpath_overrides_runpath=yes
      ;;
    *,yes) # both, prefer svr4
      dynamic_linker="AIX lib.so.V[(]$shared_archive_member_spec.o[)], lib.a[(]lib.so.V[)]"
      library_names_spec='$libname$release$shared_ext$major $libname$shared_ext'
      # unpreferred sharedlib libNAME.a needs extra handling
      postinstall_cmds='test -n "$linkname" || linkname="$realname"~func_stripname "" ".so" "$linkname"~$install_shared_prog "$dir/$func_stripname_result.$libext" "$destdir/$func_stripname_result.$libext"~test -z "$tstripme" || test -z "$striplib" || $striplib "$destdir/$func_stripname_result.$libext"'
      postuninstall_cmds='for n in $library_names $old_library; do :; done~func_stripname "" ".so" "$n"~test "$func_stripname_result" = "$n" || func_append rmfiles " $odir/$func_stripname_result.$libext"'
      # We do not specify a path in Import Files, so LIBPATH fires.
      shlibpath_overrides_runpath=yes
      ;;
    *,no) # both, prefer aix
      dynamic_linker="AIX lib.a[(]lib.so.V[)], lib.so.V[(]$shared_archive_member_spec.o[)]"
      library_names_spec='$libname$release.a $libname.a'
      soname_spec='$libname$release$shared_ext$major'
      # unpreferred sharedlib libNAME.so.V and symlink libNAME.so need extra handling
      postinstall_cmds='test -z "$dlname" || $install_shared_prog $dir/$dlname $destdir/$dlname~test -z "$tstripme" || test -z "$striplib" || $striplib $destdir/$dlname~test -n "$linkname" || linkname=$realname~func_stripname "" ".a" "$linkname"~(cd "$destdir" && $LN_S -f $dlname $func_stripname_result.so)'
      postuninstall_cmds='test -z "$dlname" || func_append rmfiles " $odir/$dlname"~for n in $old_library $library_names; do :; done~func_stripname "" ".a" "$n"~func_append rmfiles " $odir/$func_stripname_result.so"'
      ;;
    esac
    shlibpath_var=LIBPATH
  fi
  ;;

amigaos*)
  case $host_cpu in
  powerpc)
    # Since July 2007 AmigaOS4 officially supports .so libraries.
    # When compiling the executable, add -use-dynld -Lsobjs: to the compileline.
    library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
    ;;
  m68k)
    library_names_spec='$libname.ixlibrary $libname.a'
    # Create ${libname}_ixlibrary.a entries in /sys/libs.
    finish_eval='for lib in `ls $libdir/*.ixlibrary 2>/dev/null`; do libname=`func_echo_all "$lib" | $SED '\''s%^.*/\([[^/]]*\)\.ixlibrary$%\1%'\''`; $RM /sys/libs/${libname}_ixlibrary.a; $show "cd /sys/libs && $LN_S $lib ${libname}_ixlibrary.a"; cd /sys/libs && $LN_S $lib ${libname}_ixlibrary.a || exit 1; done'
    ;;
  esac
  ;;

beos*)
  library_names_spec='$libname$shared_ext'
  dynamic_linker="$host_os ld.so"
  shlibpath_var=LIBRARY_PATH
  ;;

bsdi[[45]]*)
  version_type=linux # correct to gnu/linux during the next big refactor
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  sys_lib_search_path_spec="/shlib /usr/lib /usr/X11/lib /usr/contrib/lib /lib /usr/local/lib"
  sys_lib_dlsearch_path_spec="/shlib /usr/lib /usr/local/lib"
  # the default ld.so.conf also contains /usr/contrib/lib and
  # /usr/X11R6/lib (/usr/X11 is a link to /usr/X11R6), but let us allow
  # libtool to hard-code these into programs
  ;;

cygwin* | mingw* | pw32* | cegcc*)
  version_type=windows
  shrext_cmds=.dll
  need_version=no
  need_lib_prefix=no

  case $GCC,$cc_basename in
  yes,*)
    # gcc
    library_names_spec='$libname.dll.a'
    # DLL is installed to $(libdir)/../bin by postinstall_cmds
    postinstall_cmds='base_file=`basename \$file`~
      dlpath=`$SHELL 2>&1 -c '\''. $dir/'\''\$base_file'\''i; echo \$dlname'\''`~
      dldir=$destdir/`dirname \$dlpath`~
      test -d \$dldir || mkdir -p \$dldir~
      $install_prog $dir/$dlname \$dldir/$dlname~
      chmod a+x \$dldir/$dlname~
      if test -n '\''$stripme'\'' && test -n '\''$striplib'\''; then
        eval '\''$striplib \$dldir/$dlname'\'' || exit \$?;
      fi'
    postuninstall_cmds='dldll=`$SHELL 2>&1 -c '\''. $file; echo \$dlname'\''`~
      dlpath=$dir/\$dldll~
       $RM \$dlpath'
    shlibpath_overrides_runpath=yes

    case $host_os in
    cygwin*)
      # Cygwin DLLs use 'cyg' prefix rather than 'lib'
      soname_spec='`echo $libname | sed -e 's/^lib/cyg/'``echo $release | $SED -e 's/[[.]]/-/g'`$versuffix$shared_ext'
m4_if([$1], [],[
      sys_lib_search_path_spec="$sys_lib_search_path_spec /usr/lib/w32api"])
      ;;
    mingw* | cegcc*)
      # MinGW DLLs use traditional 'lib' prefix
      soname_spec='$libname`echo $release | $SED -e 's/[[.]]/-/g'`$versuffix$shared_ext'
      ;;
    pw32*)
      # pw32 DLLs use 'pw' prefix rather than 'lib'
      library_names_spec='`echo $libname | sed -e 's/^lib/pw/'``echo $release | $SED -e 's/[[.]]/-/g'`$versuffix$shared_ext'
      ;;
    esac
    dynamic_linker='Win32 ld.exe'
    ;;

  *,cl* | *,icl* | *,ifort*)
    # Native MSVC or ICC or IFORT
    libname_spec='$name'
    soname_spec='$libname`echo $release | $SED -e 's/[[.]]/-/g'`$versuffix$shared_ext'
    library_names_spec='$libname.dll.lib'

    case $build_os in
    mingw*)
      sys_lib_search_path_spec=
      lt_save_ifs=$IFS
      IFS=';'
      for lt_path in $LIB
      do
        IFS=$lt_save_ifs
        # Let DOS variable expansion print the short 8.3 style file name.
        lt_path=`cd "$lt_path" 2>/dev/null && cmd //C "for %i in (".") do @echo %~si"`
        sys_lib_search_path_spec="$sys_lib_search_path_spec $lt_path"
      done
      IFS=$lt_save_ifs
      # Convert to MSYS style.
      sys_lib_search_path_spec=`$ECHO "$sys_lib_search_path_spec" | sed -e 's|\\\\|/|g' -e 's| \\([[a-zA-Z]]\\):| /\\1|g' -e 's|^ ||'`
      ;;
    cygwin*)
      # Convert to unix form, then to dos form, then back to unix form
      # but this time dos style (no spaces!) so that the unix form looks
      # like /cygdrive/c/PROGRA~1:/cygdr...
      sys_lib_search_path_spec=`cygpath --path --unix "$LIB"`
      sys_lib_search_path_spec=`cygpath --path --dos "$sys_lib_search_path_spec" 2>/dev/null`
      sys_lib_search_path_spec=`cygpath --path --unix "$sys_lib_search_path_spec" | $SED -e "s/$PATH_SEPARATOR/ /g"`
      ;;
    *)
      sys_lib_search_path_spec=$LIB
      if $ECHO "$sys_lib_search_path_spec" | [$GREP ';[c-zC-Z]:/' >/dev/null]; then
        # It is most probably a Windows format PATH.
        sys_lib_search_path_spec=`$ECHO "$sys_lib_search_path_spec" | $SED -e 's/;/ /g'`
      else
        sys_lib_search_path_spec=`$ECHO "$sys_lib_search_path_spec" | $SED -e "s/$PATH_SEPARATOR/ /g"`
      fi
      # FIXME: find the short name or the path components, as spaces are
      # common. (e.g. "Program Files" -> "PROGRA~1")
      ;;
    esac

    # DLL is installed to $(libdir)/../bin by postinstall_cmds
    postinstall_cmds='base_file=`basename \$file`~
      dlpath=`$SHELL 2>&1 -c '\''. $dir/'\''\$base_file'\''i; echo \$dlname'\''`~
      dldir=$destdir/`dirname \$dlpath`~
      test -d \$dldir || mkdir -p \$dldir~
      $install_prog $dir/$dlname \$dldir/$dlname'
    postuninstall_cmds='dldll=`$SHELL 2>&1 -c '\''. $file; echo \$dlname'\''`~
      dlpath=$dir/\$dldll~
       $RM \$dlpath'
    shlibpath_overrides_runpath=yes
    dynamic_linker='Win32 link.exe'
    ;;

  *)
    # Assume MSVC wrapper
    library_names_spec='$libname`echo $release | $SED -e 's/[[.]]/-/g'`$versuffix$shared_ext $libname.lib'
    dynamic_linker='Win32 ld.exe'
    ;;
  esac
  # FIXME: first we should search . and the directory the executable is in
  shlibpath_var=PATH
  ;;

darwin* | rhapsody*)
  dynamic_linker="$host_os dyld"
  version_type=darwin
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$major$shared_ext $libname$shared_ext'
  soname_spec='$libname$release$major$shared_ext'
  shlibpath_overrides_runpath=yes
  shlibpath_var=DYLD_LIBRARY_PATH
  shrext_cmds='`test .$module = .yes && echo .so || echo .dylib`'
m4_if([$1], [],[
  sys_lib_search_path_spec="$sys_lib_search_path_spec /usr/local/lib"])
  sys_lib_dlsearch_path_spec='/usr/local/lib /lib /usr/lib'
  ;;

dgux*)
  version_type=linux # correct to gnu/linux during the next big refactor
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

freebsd* | dragonfly*)
  # DragonFly does not have aout.  When/if they implement a new
  # versioning mechanism, adjust this.
  if test -x /usr/bin/objformat; then
    objformat=`/usr/bin/objformat`
  else
    case $host_os in
    freebsd[[23]].*) objformat=aout ;;
    *) objformat=elf ;;
    esac
  fi
  version_type=freebsd-$objformat
  case $version_type in
    freebsd-elf*)
      library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
      soname_spec='$libname$release$shared_ext$major'
      need_version=no
      need_lib_prefix=no
      ;;
    freebsd-*)
      library_names_spec='$libname$release$shared_ext$versuffix $libname$shared_ext$versuffix'
      need_version=yes
      ;;
  esac
  shlibpath_var=LD_LIBRARY_PATH
  case $host_os in
  freebsd2.*)
    shlibpath_overrides_runpath=yes
    ;;
  freebsd3.[[01]]* | freebsdelf3.[[01]]*)
    shlibpath_overrides_runpath=yes
    hardcode_into_libs=yes
    ;;
  freebsd3.[[2-9]]* | freebsdelf3.[[2-9]]* | \
  freebsd4.[[0-5]] | freebsdelf4.[[0-5]] | freebsd4.1.1 | freebsdelf4.1.1)
    shlibpath_overrides_runpath=no
    hardcode_into_libs=yes
    ;;
  *) # from 4.6 on, and DragonFly
    shlibpath_overrides_runpath=yes
    hardcode_into_libs=yes
    ;;
  esac
  ;;

haiku*)
  version_type=linux # correct to gnu/linux during the next big refactor
  need_lib_prefix=no
  need_version=no
  dynamic_linker="$host_os runtime_loader"
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  shlibpath_var=LIBRARY_PATH
  shlibpath_overrides_runpath=no
  sys_lib_dlsearch_path_spec='/boot/home/config/lib /boot/common/lib /boot/system/lib'
  hardcode_into_libs=yes
  ;;

hpux9* | hpux10* | hpux11*)
  # Give a soname corresponding to the major version so that dld.sl refuses to
  # link against other versions.
  version_type=sunos
  need_lib_prefix=no
  need_version=no
  case $host_cpu in
  ia64*)
    shrext_cmds='.so'
    hardcode_into_libs=yes
    dynamic_linker="$host_os dld.so"
    shlibpath_var=LD_LIBRARY_PATH
    shlibpath_overrides_runpath=yes # Unless +noenvvar is specified.
    library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
    soname_spec='$libname$release$shared_ext$major'
    if test 32 = "$HPUX_IA64_MODE"; then
      sys_lib_search_path_spec="/usr/lib/hpux32 /usr/local/lib/hpux32 /usr/local/lib"
      sys_lib_dlsearch_path_spec=/usr/lib/hpux32
    else
      sys_lib_search_path_spec="/usr/lib/hpux64 /usr/local/lib/hpux64"
      sys_lib_dlsearch_path_spec=/usr/lib/hpux64
    fi
    ;;
  hppa*64*)
    shrext_cmds='.sl'
    hardcode_into_libs=yes
    dynamic_linker="$host_os dld.sl"
    shlibpath_var=LD_LIBRARY_PATH # How should we handle SHLIB_PATH
    shlibpath_overrides_runpath=yes # Unless +noenvvar is specified.
    library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
    soname_spec='$libname$release$shared_ext$major'
    sys_lib_search_path_spec="/usr/lib/pa20_64 /usr/ccs/lib/pa20_64"
    sys_lib_dlsearch_path_spec=$sys_lib_search_path_spec
    ;;
  *)
    shrext_cmds='.sl'
    dynamic_linker="$host_os dld.sl"
    shlibpath_var=SHLIB_PATH
    shlibpath_overrides_runpath=no # +s is required to enable SHLIB_PATH
    library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
    soname_spec='$libname$release$shared_ext$major'
    ;;
  esac
  # HP-UX runs *really* slowly unless shared libraries are mode 555, ...
  postinstall_cmds='chmod 555 $lib'
  # or fails outright, so override atomically:
  install_override_mode=555
  ;;

interix[[3-9]]*)
  version_type=linux # correct to gnu/linux during the next big refactor
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  dynamic_linker='Interix 3.x ld.so.1 (PE, like ELF)'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  hardcode_into_libs=yes
  ;;

irix5* | irix6* | nonstopux*)
  case $host_os in
    nonstopux*) version_type=nonstopux ;;
    *)
	if test yes = "$lt_cv_prog_gnu_ld"; then
		version_type=linux # correct to gnu/linux during the next big refactor
	else
		version_type=irix
	fi ;;
  esac
  need_lib_prefix=no
  need_version=no
  soname_spec='$libname$release$shared_ext$major'
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$release$shared_ext $libname$shared_ext'
  case $host_os in
  irix5* | nonstopux*)
    libsuff= shlibsuff=
    ;;
  *)
    case $LD in # libtool.m4 will add one of these switches to LD
    *-32|*"-32 "|*-melf32bsmip|*"-melf32bsmip ")
      libsuff= shlibsuff= libmagic=32-bit;;
    *-n32|*"-n32 "|*-melf32bmipn32|*"-melf32bmipn32 ")
      libsuff=32 shlibsuff=N32 libmagic=N32;;
    *-64|*"-64 "|*-melf64bmip|*"-melf64bmip ")
      libsuff=64 shlibsuff=64 libmagic=64-bit;;
    *) libsuff= shlibsuff= libmagic=never-match;;
    esac
    ;;
  esac
  shlibpath_var=LD_LIBRARY${shlibsuff}_PATH
  shlibpath_overrides_runpath=no
  sys_lib_search_path_spec="/usr/lib$libsuff /lib$libsuff /usr/local/lib$libsuff"
  sys_lib_dlsearch_path_spec="/usr/lib$libsuff /lib$libsuff"
  hardcode_into_libs=yes
  ;;

# No shared lib support for Linux oldld, aout, or coff.
linux*oldld* | linux*aout* | linux*coff*)
  dynamic_linker=no
  ;;

linux*android*)
  version_type=none # Android doesn't support versioned libraries.
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext'
  soname_spec='$libname$release$shared_ext'
  finish_cmds=
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes

  # This implies no fast_install, which is unacceptable.
  # Some rework will be needed to allow for fast_install
  # before this can be enabled.
  hardcode_into_libs=yes

  dynamic_linker='Android linker'
  # Don't embed -rpath directories since the linker doesn't support them.
  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
  ;;

# This must be glibc/ELF.
linux* | k*bsd*-gnu | kopensolaris*-gnu | gnu*)
  version_type=linux # correct to gnu/linux during the next big refactor
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig -n $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no

  # Some binutils ld are patched to set DT_RUNPATH
  AC_CACHE_VAL([lt_cv_shlibpath_overrides_runpath],
    [lt_cv_shlibpath_overrides_runpath=no
    save_LDFLAGS=$LDFLAGS
    save_libdir=$libdir
    eval "libdir=/foo; wl=\"$_LT_TAGVAR(lt_prog_compiler_wl, $1)\"; \
	 LDFLAGS=\"\$LDFLAGS $_LT_TAGVAR(hardcode_libdir_flag_spec, $1)\""
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],[])],
      [AS_IF([ ($OBJDUMP -p conftest$ac_exeext) 2>/dev/null | grep "RUNPATH.*$libdir" >/dev/null],
	 [lt_cv_shlibpath_overrides_runpath=yes])])
    LDFLAGS=$save_LDFLAGS
    libdir=$save_libdir
    ])
  shlibpath_overrides_runpath=$lt_cv_shlibpath_overrides_runpath

  # This implies no fast_install, which is unacceptable.
  # Some rework will be needed to allow for fast_install
  # before this can be enabled.
  hardcode_into_libs=yes

  # Ideally, we could use ldconfig to report *all* directores which are
  # searched for libraries, however this is still not possible.  Aside from not
  # being certain /sbin/ldconfig is available, command
  # 'ldconfig -N -X -v | grep ^/' on 64bit Fedora does not report /usr/lib64,
  # even though it is searched at run-time.  Try to do the best guess by
  # appending ld.so.conf contents (and includes) to the search path.
  if test -f /etc/ld.so.conf; then
    lt_ld_extra=`awk '/^include / { system(sprintf("cd /etc; cat %s 2>/dev/null", \[$]2)); skip = 1; } { if (!skip) print \[$]0; skip = 0; }' < /etc/ld.so.conf | $SED -e 's/#.*//;/^[	 ]*hwcap[	 ]/d;s/[:,	]/ /g;s/=[^=]*$//;s/=[^= ]* / /g;s/"//g;/^$/d' | tr '\n' ' '`
    sys_lib_dlsearch_path_spec="/lib /usr/lib $lt_ld_extra"
  fi

  # We used to test for /lib/ld.so.1 and disable shared libraries on
  # powerpc, because MkLinux only supported shared libraries with the
  # GNU dynamic linker.  Since this was broken with cross compilers,
  # most powerpc-linux boxes support dynamic linking these days and
  # people can always --disable-shared, the test was removed, and we
  # assume the GNU/Linux dynamic linker is in use.
  dynamic_linker='GNU/Linux ld.so'
  ;;

netbsd*)
  version_type=sunos
  need_lib_prefix=no
  need_version=no
  if echo __ELF__ | $CC -E - | $GREP __ELF__ >/dev/null; then
    library_names_spec='$libname$release$shared_ext$versuffix $libname$shared_ext$versuffix'
    finish_cmds='PATH="\$PATH:/sbin" ldconfig -m $libdir'
    dynamic_linker='NetBSD (a.out) ld.so'
  else
    library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
    soname_spec='$libname$release$shared_ext$major'
    dynamic_linker='NetBSD ld.elf_so'
  fi
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  hardcode_into_libs=yes
  ;;

newsos6)
  version_type=linux # correct to gnu/linux during the next big refactor
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  ;;

*nto* | *qnx*)
  version_type=qnx
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  hardcode_into_libs=yes
  dynamic_linker='ldqnx.so'
  ;;

openbsd* | bitrig*)
  version_type=sunos
  sys_lib_dlsearch_path_spec=/usr/lib
  need_lib_prefix=no
  if test -z "`echo __ELF__ | $CC -E - | $GREP __ELF__`"; then
    need_version=no
  else
    need_version=yes
  fi
  library_names_spec='$libname$release$shared_ext$versuffix $libname$shared_ext$versuffix'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig -m $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  ;;

os2*)
  libname_spec='$name'
  version_type=windows
  shrext_cmds=.dll
  need_version=no
  need_lib_prefix=no
  # OS/2 can only load a DLL with a base name of 8 characters or less.
  soname_spec='`test -n "$os2dllname" && libname="$os2dllname";
    v=$($ECHO $release$versuffix | tr -d .-);
    n=$($ECHO $libname | cut -b -$((8 - ${#v})) | tr . _);
    $ECHO $n$v`$shared_ext'
  library_names_spec='${libname}_dll.$libext'
  dynamic_linker='OS/2 ld.exe'
  shlibpath_var=BEGINLIBPATH
  sys_lib_search_path_spec="/lib /usr/lib /usr/local/lib"
  sys_lib_dlsearch_path_spec=$sys_lib_search_path_spec
  postinstall_cmds='base_file=`basename \$file`~
    dlpath=`$SHELL 2>&1 -c '\''. $dir/'\''\$base_file'\''i; $ECHO \$dlname'\''`~
    dldir=$destdir/`dirname \$dlpath`~
    test -d \$dldir || mkdir -p \$dldir~
    $install_prog $dir/$dlname \$dldir/$dlname~
    chmod a+x \$dldir/$dlname~
    if test -n '\''$stripme'\'' && test -n '\''$striplib'\''; then
      eval '\''$striplib \$dldir/$dlname'\'' || exit \$?;
    fi'
  postuninstall_cmds='dldll=`$SHELL 2>&1 -c '\''. $file; $ECHO \$dlname'\''`~
    dlpath=$dir/\$dldll~
    $RM \$dlpath'
  ;;

osf3* | osf4* | osf5*)
  version_type=osf
  need_lib_prefix=no
  need_version=no
  soname_spec='$libname$release$shared_ext$major'
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  shlibpath_var=LD_LIBRARY_PATH
  sys_lib_search_path_spec="/usr/shlib /usr/ccs/lib /usr/lib/cmplrs/cc /usr/lib /usr/local/lib /var/shlib"
  sys_lib_dlsearch_path_spec=$sys_lib_search_path_spec
  ;;

rdos*)
  dynamic_linker=no
  ;;

solaris*)
  version_type=linux # correct to gnu/linux during the next big refactor
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  hardcode_into_libs=yes
  # ldd complains unless libraries are executable
  postinstall_cmds='chmod +x $lib'
  ;;

sunos4*)
  version_type=sunos
  library_names_spec='$libname$release$shared_ext$versuffix $libname$shared_ext$versuffix'
  finish_cmds='PATH="\$PATH:/usr/etc" ldconfig $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  if test yes = "$with_gnu_ld"; then
    need_lib_prefix=no
  fi
  need_version=yes
  ;;

sysv4 | sysv4.3*)
  version_type=linux # correct to gnu/linux during the next big refactor
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  shlibpath_var=LD_LIBRARY_PATH
  case $host_vendor in
    sni)
      shlibpath_overrides_runpath=no
      need_lib_prefix=no
      runpath_var=LD_RUN_PATH
      ;;
    siemens)
      need_lib_prefix=no
      ;;
    motorola)
      need_lib_prefix=no
      need_version=no
      shlibpath_overrides_runpath=no
      sys_lib_search_path_spec='/lib /usr/lib /usr/ccs/lib'
      ;;
  esac
  ;;

sysv4*MP*)
  if test -d /usr/nec; then
    version_type=linux # correct to gnu/linux during the next big refactor
    library_names_spec='$libname$shared_ext.$versuffix $libname$shared_ext.$major $libname$shared_ext'
    soname_spec='$libname$shared_ext.$major'
    shlibpath_var=LD_LIBRARY_PATH
  fi
  ;;

sysv5* | sco3.2v5* | sco5v6* | unixware* | OpenUNIX* | sysv4*uw2*)
  version_type=sco
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  hardcode_into_libs=yes
  if test yes = "$with_gnu_ld"; then
    sys_lib_search_path_spec='/usr/local/lib /usr/gnu/lib /usr/ccs/lib /usr/lib /lib'
  else
    sys_lib_search_path_spec='/usr/ccs/lib /usr/lib'
    case $host_os in
      sco3.2v5*)
        sys_lib_search_path_spec="$sys_lib_search_path_spec /lib"
	;;
    esac
  fi
  sys_lib_dlsearch_path_spec='/usr/lib'
  ;;

tpf*)
  # TPF is a cross-target only.  Preferred cross-host = GNU/Linux.
  version_type=linux # correct to gnu/linux during the next big refactor
  need_lib_prefix=no
  need_version=no
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  hardcode_into_libs=yes
  ;;

uts4*)
  version_type=linux # correct to gnu/linux during the next big refactor
  library_names_spec='$libname$release$shared_ext$versuffix $libname$release$shared_ext$major $libname$shared_ext'
  soname_spec='$libname$release$shared_ext$major'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

*)
  dynamic_linker=no
  ;;
esac
AC_MSG_RESULT([$dynamic_linker])
test no = "$dynamic_linker" && can_build_shared=no

variables_saved_for_relink="PATH $shlibpath_var $runpath_var"
if test yes = "$GCC"; then
  variables_saved_for_relink="$variables_saved_for_relink GCC_EXEC_PREFIX COMPILER_PATH LIBRARY_PATH"
fi

if test set = "${lt_cv_sys_lib_search_path_spec+set}"; then
  sys_lib_search_path_spec=$lt_cv_sys_lib_search_path_spec
fi

if test set = "${lt_cv_sys_lib_dlsearch_path_spec+set}"; then
  sys_lib_dlsearch_path_spec=$lt_cv_sys_lib_dlsearch_path_spec
fi

# remember unaugmented sys_lib_dlsearch_path content for libtool script decls...
configure_time_dlsearch_path=$sys_lib_dlsearch_path_spec

# ... but it needs LT_SYS_LIBRARY_PATH munging for other configure-time code
func_munge_path_list sys_lib_dlsearch_path_spec "$LT_SYS_LIBRARY_PATH"

# to be used as default LT_SYS_LIBRARY_PATH value in generated libtool
configure_time_lt_sys_library_path=$LT_SYS_LIBRARY_PATH

_LT_DECL([], [variables_saved_for_relink], [1],
    [Variables whose values should be saved in libtool wrapper scripts and
    restored at link time])
_LT_DECL([], [need_lib_prefix], [0],
    [Do we need the "lib" prefix for modules?])
_LT_DECL([], [need_version], [0], [Do we need a version for libraries?])
_LT_DECL([], [version_type], [0], [Library versioning type])
_LT_DECL([], [runpath_var], [0],  [Shared library runtime path variable])
_LT_DECL([], [shlibpath_var], [0],[Shared library path variable])
_LT_DECL([], [shlibpath_overrides_runpath], [0],
    [Is shlibpath searched before the hard-coded library search path?])
_LT_DECL([], [libname_spec], [1], [Format of library name prefix])
_LT_DECL([], [library_names_spec], [1],
    [[List of archive names.  First name is the real one, the rest are links.
    The last name is the one that the linker finds with -lNAME]])
_LT_DECL([], [soname_spec], [1],
    [[The coded name of the library, if different from the real name]])
_LT_DECL([], [install_override_mode], [1],
    [Permission mode override for installation of shared libraries])
_LT_DECL([], [postinstall_cmds], [2],
    [Command to use after installation of a shared archive])
_LT_DECL([], [postuninstall_cmds], [2],
    [Command to use after uninstallation of a shared archive])
_LT_DECL([], [finish_cmds], [2],
    [Commands used to finish a libtool library installation in a directory])
_LT_DECL([], [finish_eval], [1],
    [[As "finish_cmds", except a single script fragment to be evaled but
    not shown]])
_LT_DECL([], [hardcode_into_libs], [0],
    [Whether we should hardcode library paths into libraries])
_LT_DECL([], [sys_lib_search_path_spec], [2],
    [Compile-time system search path for libraries])
_LT_DECL([sys_lib_dlsearch_path_spec], [configure_time_dlsearch_path], [2],
    [Detected run-time system search path for libraries])
_LT_DECL([], [configure_time_lt_sys_library_path], [2],
    [Explicit LT_SYS_LIBRARY_PATH set during ./configure time])
])# _LT_SYS_DYNAMIC_LINKER


# _LT_PATH_TOOL_PREFIX(TOOL)
# --------------------------
# find a file program that can recognize shared library
AC_DEFUN([_LT_PATH_TOOL_PREFIX],
[m4_require([_LT_DECL_EGREP])dnl
AC_MSG_CHECKING([for $1])
AC_CACHE_VAL(lt_cv_path_MAGIC_CMD,
[case $MAGIC_CMD in
[[\\/*] |  ?:[\\/]*])
  lt_cv_path_MAGIC_CMD=$MAGIC_CMD # Let the user override the test with a path.
  ;;
*)
  lt_save_MAGIC_CMD=$MAGIC_CMD
  lt_save_ifs=$IFS; IFS=$PATH_SEPARATOR
dnl $ac_dummy forces splitting on constant user-supplied paths.
dnl POSIX.2 word splitting is done only on the output of word expansions,
dnl not every word.  This closes a longstanding sh security hole.
  ac_dummy="m4_if([$2], , $PATH, [$2])"
  for ac_dir in $ac_dummy; do
    IFS=$lt_save_ifs
    test -z "$ac_dir" && ac_dir=.
    if test -f "$ac_dir/$1"; then
      lt_cv_path_MAGIC_CMD=$ac_dir/"$1"
      if test -n "$file_magic_test_file"; then
	case $deplibs_check_method in
	"file_magic "*)
	  file_magic_regex=`expr "$deplibs_check_method" : "file_magic \(.*\)"`
	  MAGIC_CMD=$lt_cv_path_MAGIC_CMD
	  if eval $file_magic_cmd \$file_magic_test_file 2> /dev/null |
	    $EGREP "$file_magic_regex" > /dev/null; then
	    :
	  else
	    cat <<_LT_EOF 1>&2

*** Warning: the command libtool uses to detect shared libraries,
*** $file_magic_cmd, produces output that libtool cannot recognize.
*** The result is that libtool may fail to recognize shared libraries
*** as such.  This will affect the creation of libtool libraries that
*** depend on shared libraries, but programs linked with such libtool
*** libraries will work regardless of this problem.  Nevertheless, you
*** may want to report the problem to your system manager and/or to
*** bug-libtool@gnu.org

_LT_EOF
	  fi ;;
	esac
      fi
      break
    fi
  done
  IFS=$lt_save_ifs
  MAGIC_CMD=$lt_save_MAGIC_CMD
  ;;
esac])
MAGIC_CMD=$lt_cv_path_MAGIC_CMD
if test -n "$MAGIC_CMD"; then
  AC_MSG_RESULT($MAGIC_CMD)
else
  AC_MSG_RESULT(no)
fi
_LT_DECL([], [MAGIC_CMD], [0],
	 [Used to examine libraries when file_magic_cmd begins with "file"])dnl
])# _LT_PATH_TOOL_PREFIX

# Old name:
AU_ALIAS([AC_PATH_TOOL_PREFIX], [_LT_PATH_TOOL_PREFIX])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_PATH_TOOL_PREFIX], [])


# _LT_PATH_MAGIC
# --------------
# find a file program that can recognize a shared library
m4_defun([_LT_PATH_MAGIC],
[_LT_PATH_TOOL_PREFIX(${ac_tool_prefix}file, /usr/bin$PATH_SEPARATOR$PATH)
if test -z "$lt_cv_path_MAGIC_CMD"; then
  if test -n "$ac_tool_prefix"; then
    _LT_PATH_TOOL_PREFIX(file, /usr/bin$PATH_SEPARATOR$PATH)
  else
    MAGIC_CMD=:
  fi
fi
])# _LT_PATH_MAGIC


# LT_PATH_LD
# ----------
# find the pathname to the GNU or non-GNU linker
AC_DEFUN([LT_PATH_LD],
[AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
m4_require([_LT_DECL_SED])dnl
m4_require([_LT_DECL_EGREP])dnl
m4_require([_LT_PROG_ECHO_BACKSLASH])dnl

AC_ARG_WITH([gnu-ld],
    [AS_HELP_STRING([--with-gnu-ld],
	[assume the C compiler uses GNU ld @<:@default=no@:>@])],
    [test no = "$withval" || with_gnu_ld=yes],
    [with_gnu_ld=no])dnl

ac_prog=ld
if test yes = "$GCC"; then
  # Check if gcc -print-prog-name=ld gives a path.
  AC_MSG_CHECKING([for ld used by $CC])
  case $host in
  *-*-mingw*)
    # gcc leaves a trailing carriage return, which upsets mingw
    ac_prog=`($CC -print-prog-name=ld) 2>&5 | tr -d '\015'` ;;
  *)
    ac_prog=`($CC -print-prog-name=ld) 2>&5` ;;
  esac
  case $ac_prog in
    # Accept absolute paths.
    [[\\/]]* | ?:[[\\/]]*)
      re_direlt='/[[^/]][[^/]]*/\.\./'
      # Canonicalize the pathname of ld
      ac_prog=`$ECHO "$ac_prog"| $SED 's%\\\\%/%g'`
      while $ECHO "$ac_prog" | $GREP "$re_direlt" > /dev/null 2>&1; do
	ac_prog=`$ECHO $ac_prog| $SED "s%$re_direlt%/%"`
      done
      test -z "$LD" && LD=$ac_prog
      ;;
  "")
    # If it fails, then pretend we aren't using GCC.
    ac_prog=ld
    ;;
  *)
    # If it is relative, then search for the first ld in PATH.
    with_gnu_ld=unknown
    ;;
  esac
elif test yes = "$with_gnu_ld"; then
  AC_MSG_CHECKING([for GNU ld])
else
  AC_MSG_CHECKING([for non-GNU ld])
fi
AC_CACHE_VAL(lt_cv_path_LD,
[if test -z "$LD"; then
  lt_save_ifs=$IFS; IFS=$PATH_SEPARATOR
  for ac_dir in $PATH; do
    IFS=$lt_save_ifs
    test -z "$ac_dir" && ac_dir=.
    if test -f "$ac_dir/$ac_prog" || test -f "$ac_dir/$ac_prog$ac_exeext"; then
      lt_cv_path_LD=$ac_dir/$ac_prog
      # Check to see if the program is GNU ld.  I'd rather use --version,
      # but apparently some variants of GNU ld only accept -v.
      # Break only if it was the GNU/non-GNU ld that we prefer.
      case `"$lt_cv_path_LD" -v 2>&1 </dev/null` in
      *GNU* | *'with BFD'*)
	test no != "$with_gnu_ld" && break
	;;
      *)
	test yes != "$with_gnu_ld" && break
	;;
      esac
    fi
  done
  IFS=$lt_save_ifs
else
  lt_cv_path_LD=$LD # Let the user override the test with a path.
fi])
LD=$lt_cv_path_LD
if test -n "$LD"; then
  AC_MSG_RESULT($LD)
else
  AC_MSG_RESULT(no)
fi
test -z "$LD" && AC_MSG_ERROR([no acceptable ld found in \$PATH])
_LT_PATH_LD_GNU
AC_SUBST([LD])

_LT_TAGDECL([], [LD], [1], [The linker used to build libraries])
])# LT_PATH_LD

# Old names:
AU_ALIAS([AM_PROG_LD], [LT_PATH_LD])
AU_ALIAS([AC_PROG_LD], [LT_PATH_LD])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AM_PROG_LD], [])
dnl AC_DEFUN([AC_PROG_LD], [])


# _LT_PATH_LD_GNU
#- --------------
m4_defun([_LT_PATH_LD_GNU],
[AC_CACHE_CHECK([if the linker ($LD) is GNU ld], lt_cv_prog_gnu_ld,
[# I'd rather use --version here, but apparently some GNU lds only accept -v.
case `$LD -v 2>&1 </dev/null` in
*GNU* | *'with BFD'*)
  lt_cv_prog_gnu_ld=yes
  ;;
*)
  lt_cv_prog_gnu_ld=no
  ;;
esac])
with_gnu_ld=$lt_cv_prog_gnu_ld
])# _LT_PATH_LD_GNU


# _LT_CMD_RELOAD
# --------------
# find reload flag for linker
#   -- PORTME Some linkers may need a different reload flag.
m4_defun([_LT_CMD_RELOAD],
[AC_CACHE_CHECK([for $LD option to reload object files],
  lt_cv_ld_reload_flag,
  [lt_cv_ld_reload_flag='-r'])
reload_flag=$lt_cv_ld_reload_flag
case $reload_flag in
"" | " "*) ;;
*) reload_flag=" $reload_flag" ;;
esac
reload_cmds='$LD$reload_flag -o $output$reload_objs'
case $host_os in
  cygwin* | mingw* | pw32* | cegcc*)
    if test yes != "$GCC"; then
      reload_cmds=false
    fi
    ;;
  darwin*)
    if test yes = "$GCC"; then
      reload_cmds='$LTCC $LTCFLAGS -nostdlib $wl-r -o $output$reload_objs'
    else
      reload_cmds='$LD$reload_flag -o $output$reload_objs'
    fi
    ;;
esac
_LT_TAGDECL([], [reload_flag], [1], [How to create reloadable object files])dnl
_LT_TAGDECL([], [reload_cmds], [2])dnl
])# _LT_CMD_RELOAD


# _LT_PATH_DD
# -----------
# find a working dd
m4_defun([_LT_PATH_DD],
[AC_CACHE_CHECK([for a working dd], [ac_cv_path_lt_DD],
[printf 0123456789abcdef0123456789abcdef >conftest.i
cat conftest.i conftest.i >conftest2.i
: ${lt_DD:=$DD}
AC_PATH_PROGS_FEATURE_CHECK([lt_DD], [dd],
[if "$ac_path_lt_DD" bs=32 count=1 <conftest2.i >conftest.out 2>/dev/null; then
  cmp -s conftest.i conftest.out \
  && ac_cv_path_lt_DD="$ac_path_lt_DD" ac_path_lt_DD_found=:
fi])
rm -f conftest.i conftest2.i conftest.out])
])# _LT_PATH_DD


# _LT_CMD_TRUNCATE
# ----------------
# find command to truncate a binary pipe
m4_defun([_LT_CMD_TRUNCATE],
[m4_require([_LT_PATH_DD])
AC_CACHE_CHECK([how to truncate binary pipes], [lt_cv_truncate_bin],
[printf 0123456789abcdef0123456789abcdef >conftest.i
cat conftest.i conftest.i >conftest2.i
lt_cv_truncate_bin=
if "$ac_cv_path_lt_DD" bs=32 count=1 <conftest2.i >conftest.out 2>/dev/null; then
  cmp -s conftest.i conftest.out \
  && lt_cv_truncate_bin="$ac_cv_path_lt_DD bs=4096 count=1"
fi
rm -f conftest.i conftest2.i conftest.out
test -z "$lt_cv_truncate_bin" && lt_cv_truncate_bin="$SED -e 4q"])
_LT_DECL([lt_truncate_bin], [lt_cv_truncate_bin], [1],
  [Command to truncate a binary pipe])
])# _LT_CMD_TRUNCATE


# _LT_CHECK_MAGIC_METHOD
# ----------------------
# how to check for library dependencies
#  -- PORTME fill in with the dynamic library characteristics
m4_defun([_LT_CHECK_MAGIC_METHOD],
[m4_require([_LT_DECL_EGREP])
m4_require([_LT_DECL_OBJDUMP])
AC_CACHE_CHECK([how to recognize dependent libraries],
lt_cv_deplibs_check_method,
[lt_cv_file_magic_cmd='$MAGIC_CMD'
lt_cv_file_magic_test_file=
lt_cv_deplibs_check_method='unknown'
# Need to set the preceding variable on all platforms that support
# interlibrary dependencies.
# 'none' -- dependencies not supported.
# 'unknown' -- same as none, but documents that we really don't know.
# 'pass_all' -- all dependencies passed with no checks.
# 'test_compile' -- check by making test program.
# 'file_magic [[regex]]' -- check by looking for files in library path
# that responds to the $file_magic_cmd with a given extended regex.
# If you have 'file' or equivalent on your system and you're not sure
# whether 'pass_all' will *always* work, you probably want this one.

case $host_os in
aix[[4-9]]*)
  lt_cv_deplibs_check_method=pass_all
  ;;

beos*)
  lt_cv_deplibs_check_method=pass_all
  ;;

bsdi[[45]]*)
  lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (shared object|dynamic lib)'
  lt_cv_file_magic_cmd='/usr/bin/file -L'
  lt_cv_file_magic_test_file=/shlib/libc.so
  ;;

cygwin*)
  # func_win32_libid is a shell function defined in ltmain.sh
  lt_cv_deplibs_check_method='file_magic ^x86 archive import|^x86 DLL'
  lt_cv_file_magic_cmd='func_win32_libid'
  ;;

mingw* | pw32*)
  # Base MSYS/MinGW do not provide the 'file' command needed by
  # func_win32_libid shell function, so use a weaker test based on 'objdump',
  # unless we find 'file', for example because we are cross-compiling.
  if ( file / ) >/dev/null 2>&1; then
    lt_cv_deplibs_check_method='file_magic ^x86 archive import|^x86 DLL'
    lt_cv_file_magic_cmd='func_win32_libid'
  else
    # Keep this pattern in sync with the one in func_win32_libid.
    lt_cv_deplibs_check_method='file_magic file format (pei*-i386(.*architecture: i386)?|pe-arm-wince|pe-x86-64)'
    lt_cv_file_magic_cmd='$OBJDUMP -f'
  fi
  ;;

cegcc*)
  # use the weaker test based on 'objdump'. See mingw*.
  lt_cv_deplibs_check_method='file_magic file format pe-arm-.*little(.*architecture: arm)?'
  lt_cv_file_magic_cmd='$OBJDUMP -f'
  ;;

darwin* | rhapsody*)
  lt_cv_deplibs_check_method=pass_all
  ;;

freebsd* | dragonfly*)
  if echo __ELF__ | $CC -E - | $GREP __ELF__ > /dev/null; then
    case $host_cpu in
    i*86 )
      # Not sure whether the presence of OpenBSD here was a mistake.
      # Let's accept both of them until this is cleared up.
      lt_cv_deplibs_check_method='file_magic (FreeBSD|OpenBSD|DragonFly)/i[[3-9]]86 (compact )?demand paged shared library'
      lt_cv_file_magic_cmd=/usr/bin/file
      lt_cv_file_magic_test_file=`echo /usr/lib/libc.so.*`
      ;;
    esac
  else
    lt_cv_deplibs_check_method=pass_all
  fi
  ;;

haiku*)
  lt_cv_deplibs_check_method=pass_all
  ;;

hpux10.20* | hpux11*)
  lt_cv_file_magic_cmd=/usr/bin/file
  case $host_cpu in
  ia64*)
    lt_cv_deplibs_check_method='file_magic (s[[0-9]][[0-9]][[0-9]]|ELF-[[0-9]][[0-9]]) shared object file - IA64'
    lt_cv_file_magic_test_file=/usr/lib/hpux32/libc.so
    ;;
  hppa*64*)
    [lt_cv_deplibs_check_method='file_magic (s[0-9][0-9][0-9]|ELF[ -][0-9][0-9])(-bit)?( [LM]SB)? shared object( file)?[, -]* PA-RISC [0-9]\.[0-9]']
    lt_cv_file_magic_test_file=/usr/lib/pa20_64/libc.sl
    ;;
  *)
    lt_cv_deplibs_check_method='file_magic (s[[0-9]][[0-9]][[0-9]]|PA-RISC[[0-9]]\.[[0-9]]) shared library'
    lt_cv_file_magic_test_file=/usr/lib/libc.sl
    ;;
  esac
  ;;

interix[[3-9]]*)
  # PIC code is broken on Interix 3.x, that's why |\.a not |_pic\.a here
  lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so|\.a)$'
  ;;

irix5* | irix6* | nonstopux*)
  case $LD in
  *-32|*"-32 ") libmagic=32-bit;;
  *-n32|*"-n32 ") libmagic=N32;;
  *-64|*"-64 ") libmagic=64-bit;;
  *) libmagic=never-match;;
  esac
  lt_cv_deplibs_check_method=pass_all
  ;;

# This must be glibc/ELF.
linux* | k*bsd*-gnu | kopensolaris*-gnu | gnu*)
  lt_cv_deplibs_check_method=pass_all
  ;;

netbsd*)
  if echo __ELF__ | $CC -E - | $GREP __ELF__ > /dev/null; then
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so\.[[0-9]]+\.[[0-9]]+|_pic\.a)$'
  else
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so|_pic\.a)$'
  fi
  ;;

newos6*)
  lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (executable|dynamic lib)'
  lt_cv_file_magic_cmd=/usr/bin/file
  lt_cv_file_magic_test_file=/usr/lib/libnls.so
  ;;

*nto* | *qnx*)
  lt_cv_deplibs_check_method=pass_all
  ;;

openbsd* | bitrig*)
  if test -z "`echo __ELF__ | $CC -E - | $GREP __ELF__`"; then
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so\.[[0-9]]+\.[[0-9]]+|\.so|_pic\.a)$'
  else
    lt_cv_deplibs_check_method='match_pattern /lib[[^/]]+(\.so\.[[0-9]]+\.[[0-9]]+|_pic\.a)$'
  fi
  ;;

osf3* | osf4* | osf5*)
  lt_cv_deplibs_check_method=pass_all
  ;;

rdos*)
  lt_cv_deplibs_check_method=pass_all
  ;;

solaris*)
  lt_cv_deplibs_check_method=pass_all
  ;;

sysv5* | sco3.2v5* | sco5v6* | unixware* | OpenUNIX* | sysv4*uw2*)
  lt_cv_deplibs_check_method=pass_all
  ;;

sysv4 | sysv4.3*)
  case $host_vendor in
  motorola)
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (shared object|dynamic lib) M[[0-9]][[0-9]]* Version [[0-9]]'
    lt_cv_file_magic_test_file=`echo /usr/lib/libc.so*`
    ;;
  ncr)
    lt_cv_deplibs_check_method=pass_all
    ;;
  sequent)
    lt_cv_file_magic_cmd='/bin/file'
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB (shared object|dynamic lib )'
    ;;
  sni)
    lt_cv_file_magic_cmd='/bin/file'
    lt_cv_deplibs_check_method="file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB dynamic lib"
    lt_cv_file_magic_test_file=/lib/libc.so
    ;;
  siemens)
    lt_cv_deplibs_check_method=pass_all
    ;;
  pc)
    lt_cv_deplibs_check_method=pass_all
    ;;
  esac
  ;;

tpf*)
  lt_cv_deplibs_check_method=pass_all
  ;;
os2*)
  lt_cv_deplibs_check_method=pass_all
  ;;
esac
])

file_magic_glob=
want_nocaseglob=no
if test "$build" = "$host"; then
  case $host_os in
  mingw* | pw32*)
    if ( shopt | grep nocaseglob ) >/dev/null 2>&1; then
      want_nocaseglob=yes
    else
      file_magic_glob=`echo aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ | $SED -e "s/\(..\)/s\/[[\1]]\/[[\1]]\/g;/g"`
    fi
    ;;
  esac
fi

file_magic_cmd=$lt_cv_file_magic_cmd
deplibs_check_method=$lt_cv_deplibs_check_method
test -z "$deplibs_check_method" && deplibs_check_method=unknown

_LT_DECL([], [deplibs_check_method], [1],
    [Method to check whether dependent libraries are shared objects])
_LT_DECL([], [file_magic_cmd], [1],
    [Command to use when deplibs_check_method = "file_magic"])
_LT_DECL([], [file_magic_glob], [1],
    [How to find potential files when deplibs_check_method = "file_magic"])
_LT_DECL([], [want_nocaseglob], [1],
    [Find potential files using nocaseglob when deplibs_check_method = "file_magic"])
])# _LT_CHECK_MAGIC_METHOD


# LT_PATH_NM
# ----------
# find the pathname to a BSD- or MS-compatible name lister
AC_DEFUN([LT_PATH_NM],
[AC_REQUIRE([AC_PROG_CC])dnl
AC_CACHE_CHECK([for BSD- or MS-compatible name lister (nm)], lt_cv_path_NM,
[if test -n "$NM"; then
  # Let the user override the test.
  lt_cv_path_NM=$NM
else
  lt_nm_to_check=${ac_tool_prefix}nm
  if test -n "$ac_tool_prefix" && test "$build" = "$host"; then
    lt_nm_to_check="$lt_nm_to_check nm"
  fi
  for lt_tmp_nm in $lt_nm_to_check; do
    lt_save_ifs=$IFS; IFS=$PATH_SEPARATOR
    for ac_dir in $PATH /usr/ccs/bin/elf /usr/ccs/bin /usr/ucb /bin; do
      IFS=$lt_save_ifs
      test -z "$ac_dir" && ac_dir=.
      tmp_nm=$ac_dir/$lt_tmp_nm
      if test -f "$tmp_nm" || test -f "$tmp_nm$ac_exeext"; then
	# Check to see if the nm accepts a BSD-compat flag.
	# Adding the 'sed 1q' prevents false positives on HP-UX, which says:
	#   nm: unknown option "B" ignored
	# Tru64's nm complains that /dev/null is an invalid object file
	# MSYS converts /dev/null to NUL, MinGW nm treats NUL as empty
	case $build_os in
	mingw*) lt_bad_file=conftest.nm/nofile ;;
	*) lt_bad_file=/dev/null ;;
	esac
	case `"$tmp_nm" -B $lt_bad_file 2>&1 | sed '1q'` in
	*$lt_bad_file* | *'Invalid file or object type'*)
	  lt_cv_path_NM="$tmp_nm -B"
	  break 2
	  ;;
	*)
	  case `"$tmp_nm" -p /dev/null 2>&1 | sed '1q'` in
	  */dev/null*)
	    lt_cv_path_NM="$tmp_nm -p"
	    break 2
	    ;;
	  *)
	    lt_cv_path_NM=${lt_cv_path_NM="$tmp_nm"} # keep the first match, but
	    continue # so that we can try to find one that supports BSD flags
	    ;;
	  esac
	  ;;
	esac
      fi
    done
    IFS=$lt_save_ifs
  done
  : ${lt_cv_path_NM=no}
fi])
if test no != "$lt_cv_path_NM"; then
  NM=$lt_cv_path_NM
else
  # Didn't find any BSD compatible name lister, look for dumpbin.
  if test -n "$DUMPBIN"; then :
    # Let the user override the test.
  else
    AC_CHECK_TOOLS(DUMPBIN, [dumpbin "link -dump"], :)
    case `$DUMPBIN -symbols -headers /dev/null 2>&1 | sed '1q'` in
    *COFF*)
      DUMPBIN="$DUMPBIN -symbols -headers"
      ;;
    *)
      DUMPBIN=:
      ;;
    esac
  fi
  AC_SUBST([DUMPBIN])
  if test : != "$DUMPBIN"; then
    NM=$DUMPBIN
  fi
fi
test -z "$NM" && NM=nm
AC_SUBST([NM])
_LT_DECL([], [NM], [1], [A BSD- or MS-compatible name lister])dnl

AC_CACHE_CHECK([the name lister ($NM) interface], [lt_cv_nm_interface],
  [lt_cv_nm_interface="BSD nm"
  echo "int some_variable = 0;" > conftest.$ac_ext
  (eval echo "\"\$as_me:$LINENO: $ac_compile\"" >&AS_MESSAGE_LOG_FD)
  (eval "$ac_compile" 2>conftest.err)
  cat conftest.err >&AS_MESSAGE_LOG_FD
  (eval echo "\"\$as_me:$LINENO: $NM \\\"conftest.$ac_objext\\\"\"" >&AS_MESSAGE_LOG_FD)
  (eval "$NM \"conftest.$ac_objext\"" 2>conftest.err > conftest.out)
  cat conftest.err >&AS_MESSAGE_LOG_FD
  (eval echo "\"\$as_me:$LINENO: output\"" >&AS_MESSAGE_LOG_FD)
  cat conftest.out >&AS_MESSAGE_LOG_FD
  if $GREP 'External.*some_variable' conftest.out > /dev/null; then
    lt_cv_nm_interface="MS dumpbin"
  fi
  rm -f conftest*])
])# LT_PATH_NM

# Old names:
AU_ALIAS([AM_PROG_NM], [LT_PATH_NM])
AU_ALIAS([AC_PROG_NM], [LT_PATH_NM])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AM_PROG_NM], [])
dnl AC_DEFUN([AC_PROG_NM], [])

# _LT_CHECK_SHAREDLIB_FROM_LINKLIB
# --------------------------------
# how to determine the name of the shared library
# associated with a specific link library.
#  -- PORTME fill in with the dynamic library characteristics
m4_defun([_LT_CHECK_SHAREDLIB_FROM_LINKLIB],
[m4_require([_LT_DECL_EGREP])
m4_require([_LT_DECL_OBJDUMP])
m4_require([_LT_DECL_DLLTOOL])
AC_CACHE_CHECK([how to associate runtime and link libraries],
lt_cv_sharedlib_from_linklib_cmd,
[lt_cv_sharedlib_from_linklib_cmd='unknown'

case $host_os in
cygwin* | mingw* | pw32* | cegcc*)
  # two different shell functions defined in ltmain.sh;
  # decide which one to use based on capabilities of $DLLTOOL
  case `$DLLTOOL --help 2>&1` in
  *--identify-strict*)
    lt_cv_sharedlib_from_linklib_cmd=func_cygming_dll_for_implib
    ;;
  *)
    lt_cv_sharedlib_from_linklib_cmd=func_cygming_dll_for_implib_fallback
    ;;
  esac
  ;;
*)
  # fallback: assume linklib IS sharedlib
  lt_cv_sharedlib_from_linklib_cmd=$ECHO
  ;;
esac
])
sharedlib_from_linklib_cmd=$lt_cv_sharedlib_from_linklib_cmd
test -z "$sharedlib_from_linklib_cmd" && sharedlib_from_linklib_cmd=$ECHO

_LT_DECL([], [sharedlib_from_linklib_cmd], [1],
    [Command to associate shared and link libraries])
])# _LT_CHECK_SHAREDLIB_FROM_LINKLIB


# _LT_PATH_MANIFEST_TOOL
# ----------------------
# locate the manifest tool
m4_defun([_LT_PATH_MANIFEST_TOOL],
[AC_CHECK_TOOL(MANIFEST_TOOL, mt, :)
test -z "$MANIFEST_TOOL" && MANIFEST_TOOL=mt
AC_CACHE_CHECK([if $MANIFEST_TOOL is a manifest tool], [lt_cv_path_mainfest_tool],
  [lt_cv_path_mainfest_tool=no
  echo "$as_me:$LINENO: $MANIFEST_TOOL '-?'" >&AS_MESSAGE_LOG_FD
  $MANIFEST_TOOL '-?' 2>conftest.err > conftest.out
  cat conftest.err >&AS_MESSAGE_LOG_FD
  if $GREP 'Manifest Tool' conftest.out > /dev/null; then
    lt_cv_path_mainfest_tool=yes
  fi
  rm -f conftest*])
if test yes != "$lt_cv_path_mainfest_tool"; then
  MANIFEST_TOOL=:
fi
_LT_DECL([], [MANIFEST_TOOL], [1], [Manifest tool])dnl
])# _LT_PATH_MANIFEST_TOOL


# _LT_DLL_DEF_P([FILE])
# ---------------------
# True iff FILE is a Windows DLL '.def' file.
# Keep in sync with func_dll_def_p in the libtool script
AC_DEFUN([_LT_DLL_DEF_P],
[dnl
  test DEF = "`$SED -n dnl
    -e '\''s/^[[	 ]]*//'\'' dnl Strip leading whitespace
    -e '\''/^\(;.*\)*$/d'\'' dnl      Delete empty lines and comments
    -e '\''s/^\(EXPORTS\|LIBRARY\)\([[	 ]].*\)*$/DEF/p'\'' dnl
    -e q dnl                          Only consider the first "real" line
    $1`" dnl
])# _LT_DLL_DEF_P


# LT_LIB_M
# --------
# check for math library
AC_DEFUN([LT_LIB_M],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
LIBM=
case $host in
*-*-beos* | *-*-cegcc* | *-*-cygwin* | *-*-haiku* | *-*-pw32* | *-*-darwin*)
  # These system don't have libm, or don't need it
  ;;
*-ncr-sysv4.3*)
  AC_CHECK_LIB(mw, _mwvalidcheckl, LIBM=-lmw)
  AC_CHECK_LIB(m, cos, LIBM="$LIBM -lm")
  ;;
*)
  AC_CHECK_LIB(m, cos, LIBM=-lm)
  ;;
esac
AC_SUBST([LIBM])
])# LT_LIB_M

# Old name:
AU_ALIAS([AC_CHECK_LIBM], [LT_LIB_M])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_CHECK_LIBM], [])


# _LT_COMPILER_NO_RTTI([TAGNAME])
# -------------------------------
m4_defun([_LT_COMPILER_NO_RTTI],
[m4_require([_LT_TAG_COMPILER])dnl

_LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=

if test yes = "$GCC"; then
  case $cc_basename in
  nvcc*)
    _LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=' -Xcompiler -fno-builtin' ;;
  *)
    _LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=' -fno-builtin' ;;
  esac

  _LT_COMPILER_OPTION([if $compiler supports -fno-rtti -fno-exceptions],
    lt_cv_prog_compiler_rtti_exceptions,
    [-fno-rtti -fno-exceptions], [],
    [_LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)="$_LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1) -fno-rtti -fno-exceptions"])
fi
_LT_TAGDECL([no_builtin_flag], [lt_prog_compiler_no_builtin_flag], [1],
	[Compiler flag to turn off builtin functions])
])# _LT_COMPILER_NO_RTTI


# _LT_CMD_GLOBAL_SYMBOLS
# ----------------------
m4_defun([_LT_CMD_GLOBAL_SYMBOLS],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_AWK])dnl
AC_REQUIRE([LT_PATH_NM])dnl
AC_REQUIRE([LT_PATH_LD])dnl
m4_require([_LT_DECL_SED])dnl
m4_require([_LT_DECL_EGREP])dnl
m4_require([_LT_TAG_COMPILER])dnl

# Check for command to grab the raw symbol name followed by C symbol from nm.
AC_MSG_CHECKING([command to parse $NM output from $compiler object])
AC_CACHE_VAL([lt_cv_sys_global_symbol_pipe],
[
# These are sane defaults that work on at least a few old systems.
# [They come from Ultrix.  What could be older than Ultrix?!! ;)]

# Character class describing NM global symbol codes.
symcode='[[BCDEGRST]]'

# Regexp to match symbols that can be accessed directly from C.
sympat='\([[_A-Za-z]][[_A-Za-z0-9]]*\)'

# Define system-specific variables.
case $host_os in
aix*)
  symcode='[[BCDT]]'
  ;;
cygwin* | mingw* | pw32* | cegcc*)
  symcode='[[ABCDGISTW]]'
  ;;
hpux*)
  if test ia64 = "$host_cpu"; then
    symcode='[[ABCDEGRST]]'
  fi
  ;;
irix* | nonstopux*)
  symcode='[[BCDEGRST]]'
  ;;
osf*)
  symcode='[[BCDEGQRST]]'
  ;;
solaris*)
  symcode='[[BDRT]]'
  ;;
sco3.2v5*)
  symcode='[[DT]]'
  ;;
sysv4.2uw2*)
  symcode='[[DT]]'
  ;;
sysv5* | sco5v6* | unixware* | OpenUNIX*)
  symcode='[[ABDT]]'
  ;;
sysv4)
  symcode='[[DFNSTU]]'
  ;;
esac

# If we're using GNU nm, then use its standard symbol codes.
case `$NM -V 2>&1` in
*GNU* | *'with BFD'*)
  symcode='[[ABCDGIRSTW]]' ;;
esac

if test "$lt_cv_nm_interface" = "MS dumpbin"; then
  # Gets list of data symbols to import.
  lt_cv_sys_global_symbol_to_import="sed -n -e 's/^I .* \(.*\)$/\1/p'"
  # Adjust the below global symbol transforms to fixup imported variables.
  lt_cdecl_hook=" -e 's/^I .* \(.*\)$/extern __declspec(dllimport) char \1;/p'"
  lt_c_name_hook=" -e 's/^I .* \(.*\)$/  {\"\1\", (void *) 0},/p'"
  lt_c_name_lib_hook="\
  -e 's/^I .* \(lib.*\)$/  {\"\1\", (void *) 0},/p'\
  -e 's/^I .* \(.*\)$/  {\"lib\1\", (void *) 0},/p'"
else
  # Disable hooks by default.
  lt_cv_sys_global_symbol_to_import=
  lt_cdecl_hook=
  lt_c_name_hook=
  lt_c_name_lib_hook=
fi

# Transform an extracted symbol line into a proper C declaration.
# Some systems (esp. on ia64) link data and code symbols differently,
# so use this general approach.
lt_cv_sys_global_symbol_to_cdecl="sed -n"\
$lt_cdecl_hook\
" -e 's/^T .* \(.*\)$/extern int \1();/p'"\
" -e 's/^$symcode$symcode* .* \(.*\)$/extern char \1;/p'"

# Transform an extracted symbol line into symbol name and symbol address
lt_cv_sys_global_symbol_to_c_name_address="sed -n"\
$lt_c_name_hook\
" -e 's/^: \(.*\) .*$/  {\"\1\", (void *) 0},/p'"\
" -e 's/^$symcode$symcode* .* \(.*\)$/  {\"\1\", (void *) \&\1},/p'"

# Transform an extracted symbol line into symbol name with lib prefix and
# symbol address.
lt_cv_sys_global_symbol_to_c_name_address_lib_prefix="sed -n"\
$lt_c_name_lib_hook\
" -e 's/^: \(.*\) .*$/  {\"\1\", (void *) 0},/p'"\
" -e 's/^$symcode$symcode* .* \(lib.*\)$/  {\"\1\", (void *) \&\1},/p'"\
" -e 's/^$symcode$symcode* .* \(.*\)$/  {\"lib\1\", (void *) \&\1},/p'"

# Handle CRLF in mingw tool chain
opt_cr=
case $build_os in
mingw*)
  opt_cr=`$ECHO 'x\{0,1\}' | tr x '\015'` # option cr in regexp
  ;;
esac

# Try without a prefix underscore, then with it.
for ac_symprfx in "" "_"; do

  # Transform symcode, sympat, and symprfx into a raw symbol and a C symbol.
  symxfrm="\\1 $ac_symprfx\\2 \\2"

  # Write the raw and C identifiers.
  if test "$lt_cv_nm_interface" = "MS dumpbin"; then
    # Fake it for dumpbin and say T for any non-static function,
    # D for any global variable and I for any imported variable.
    # Also find C++ and __fastcall symbols from MSVC++,
    # which start with @ or ?.
    lt_cv_sys_global_symbol_pipe="$AWK ['"\
"     {last_section=section; section=\$ 3};"\
"     /^COFF SYMBOL TABLE/{for(i in hide) delete hide[i]};"\
"     /Section length .*#relocs.*(pick any)/{hide[last_section]=1};"\
"     /^ *Symbol name *: /{split(\$ 0,sn,\":\"); si=substr(sn[2],2)};"\
"     /^ *Type *: code/{print \"T\",si,substr(si,length(prfx))};"\
"     /^ *Type *: data/{print \"I\",si,substr(si,length(prfx))};"\
"     \$ 0!~/External *\|/{next};"\
"     / 0+ UNDEF /{next}; / UNDEF \([^|]\)*()/{next};"\
"     {if(hide[section]) next};"\
"     {f=\"D\"}; \$ 0~/\(\).*\|/{f=\"T\"};"\
"     {split(\$ 0,a,/\||\r/); split(a[2],s)};"\
"     s[1]~/^[@?]/{print f,s[1],s[1]; next};"\
"     s[1]~prfx {split(s[1],t,\"@\"); print f,t[1],substr(t[1],length(prfx))}"\
"     ' prfx=^$ac_symprfx]"
  else
    lt_cv_sys_global_symbol_pipe="sed -n -e 's/^.*[[	 ]]\($symcode$symcode*\)[[	 ]][[	 ]]*$ac_symprfx$sympat$opt_cr$/$symxfrm/p'"
  fi
  lt_cv_sys_global_symbol_pipe="$lt_cv_sys_global_symbol_pipe | sed '/ __gnu_lto/d'"

  # Check to see that the pipe works correctly.
  pipe_works=no

  rm -f conftest*
  cat > conftest.$ac_ext <<_LT_EOF
#ifdef __cplusplus
extern "C" {
#endif
char nm_test_var;
void nm_test_func(void);
void nm_test_func(void){}
#ifdef __cplusplus
}
#endif
int main(){nm_test_var='a';nm_test_func();return(0);}
_LT_EOF

  if AC_TRY_EVAL(ac_compile); then
    # Now try to grab the symbols.
    nlist=conftest.nm
    if AC_TRY_EVAL(NM conftest.$ac_objext \| "$lt_cv_sys_global_symbol_pipe" \> $nlist) && test -s "$nlist"; then
      # Try sorting and uniquifying the output.
      if sort "$nlist" | uniq > "$nlist"T; then
	mv -f "$nlist"T "$nlist"
      else
	rm -f "$nlist"T
      fi

      # Make sure that we snagged all the symbols we need.
      if $GREP ' nm_test_var$' "$nlist" >/dev/null; then
	if $GREP ' nm_test_func$' "$nlist" >/dev/null; then
	  cat <<_LT_EOF > conftest.$ac_ext
/* Keep this code in sync between libtool.m4, ltmain, lt_system.h, and tests.  */
#if defined _WIN32 || defined __CYGWIN__ || defined _WIN32_WCE
/* DATA imports from DLLs on WIN32 can't be const, because runtime
   relocations are performed -- see ld's documentation on pseudo-relocs.  */
# define LT@&t@_DLSYM_CONST
#elif defined __osf__
/* This system does not cope well with relocations in const data.  */
# define LT@&t@_DLSYM_CONST
#else
# define LT@&t@_DLSYM_CONST const
#endif

#ifdef __cplusplus
extern "C" {
#endif

_LT_EOF
	  # Now generate the symbol file.
	  eval "$lt_cv_sys_global_symbol_to_cdecl"' < "$nlist" | $GREP -v main >> conftest.$ac_ext'

	  cat <<_LT_EOF >> conftest.$ac_ext

/* The mapping between symbol names and symbols.  */
LT@&t@_DLSYM_CONST struct {
  const char *name;
  void       *address;
}
lt__PROGRAM__LTX_preloaded_symbols[[]] =
{
  { "@PROGRAM@", (void *) 0 },
_LT_EOF
	  $SED "s/^$symcode$symcode* .* \(.*\)$/  {\"\1\", (void *) \&\1},/" < "$nlist" | $GREP -v main >> conftest.$ac_ext
	  cat <<\_LT_EOF >> conftest.$ac_ext
  {0, (void *) 0}
};

/* This works around a problem in FreeBSD linker */
#ifdef FREEBSD_WORKAROUND
static const void *lt_preloaded_setup() {
  return lt__PROGRAM__LTX_preloaded_symbols;
}
#endif

#ifdef __cplusplus
}
#endif
_LT_EOF
	  # Now try linking the two files.
	  mv conftest.$ac_objext conftstm.$ac_objext
	  lt_globsym_save_LIBS=$LIBS
	  lt_globsym_save_CFLAGS=$CFLAGS
	  LIBS=conftstm.$ac_objext
	  CFLAGS="$CFLAGS$_LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)"
	  if AC_TRY_EVAL(ac_link) && test -s conftest$ac_exeext; then
	    pipe_works=yes
	  fi
	  LIBS=$lt_globsym_save_LIBS
	  CFLAGS=$lt_globsym_save_CFLAGS
	else
	  echo "cannot find nm_test_func in $nlist" >&AS_MESSAGE_LOG_FD
	fi
      else
	echo "cannot find nm_test_var in $nlist" >&AS_MESSAGE_LOG_FD
      fi
    else
      echo "cannot run $lt_cv_sys_global_symbol_pipe" >&AS_MESSAGE_LOG_FD
    fi
  else
    echo "$progname: failed program was:" >&AS_MESSAGE_LOG_FD
    cat conftest.$ac_ext >&5
  fi
  rm -rf conftest* conftst*

  # Do not use the global_symbol_pipe unless it works.
  if test yes = "$pipe_works"; then
    break
  else
    lt_cv_sys_global_symbol_pipe=
  fi
done
])
if test -z "$lt_cv_sys_global_symbol_pipe"; then
  lt_cv_sys_global_symbol_to_cdecl=
fi
if test -z "$lt_cv_sys_global_symbol_pipe$lt_cv_sys_global_symbol_to_cdecl"; then
  AC_MSG_RESULT(failed)
else
  AC_MSG_RESULT(ok)
fi

# Response file support.
if test "$lt_cv_nm_interface" = "MS dumpbin"; then
  nm_file_list_spec='@'
elif $NM --help 2>/dev/null | grep '[[@]]FILE' >/dev/null; then
  nm_file_list_spec='@'
fi

_LT_DECL([global_symbol_pipe], [lt_cv_sys_global_symbol_pipe], [1],
    [Take the output of nm and produce a listing of raw symbols and C names])
_LT_DECL([global_symbol_to_cdecl], [lt_cv_sys_global_symbol_to_cdecl], [1],
    [Transform the output of nm in a proper C declaration])
_LT_DECL([global_symbol_to_import], [lt_cv_sys_global_symbol_to_import], [1],
    [Transform the output of nm into a list of symbols to manually relocate])
_LT_DECL([global_symbol_to_c_name_address],
    [lt_cv_sys_global_symbol_to_c_name_address], [1],
    [Transform the output of nm in a C name address pair])
_LT_DECL([global_symbol_to_c_name_address_lib_prefix],
    [lt_cv_sys_global_symbol_to_c_name_address_lib_prefix], [1],
    [Transform the output of nm in a C name address pair when lib prefix is needed])
_LT_DECL([nm_interface], [lt_cv_nm_interface], [1],
    [The name lister interface])
_LT_DECL([], [nm_file_list_spec], [1],
    [Specify filename containing input files for $NM])
]) # _LT_CMD_GLOBAL_SYMBOLS


# _LT_COMPILER_PIC([TAGNAME])
# ---------------------------
m4_defun([_LT_COMPILER_PIC],
[m4_require([_LT_TAG_COMPILER])dnl
_LT_TAGVAR(lt_prog_compiler_wl, $1)=
_LT_TAGVAR(lt_prog_compiler_pic, $1)=
_LT_TAGVAR(lt_prog_compiler_static, $1)=

m4_if([$1], [CXX], [
  # C++ specific cases for pic, static, wl, etc.
  if test yes = "$GXX"; then
    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
    _LT_TAGVAR(lt_prog_compiler_static, $1)='-static'

    case $host_os in
    aix*)
      # All AIX code is PIC.
      if test ia64 = "$host_cpu"; then
	# AIX 5 now supports IA64 processor
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      fi
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
      ;;

    amigaos*)
      case $host_cpu in
      powerpc)
            # see comment about AmigaOS4 .so support
            _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
        ;;
      m68k)
            # FIXME: we need at least 68020 code to build shared libraries, but
            # adding the '-m68020' flag to GCC prevents building anything better,
            # like '-m68040'.
            _LT_TAGVAR(lt_prog_compiler_pic, $1)='-m68020 -resident32 -malways-restore-a4'
        ;;
      esac
      ;;

    beos* | irix5* | irix6* | nonstopux* | osf3* | osf4* | osf5*)
      # PIC is the default for these OSes.
      ;;
    mingw* | cygwin* | os2* | pw32* | cegcc*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      # Although the cygwin gcc ignores -fPIC, still need this for old-style
      # (--disable-auto-import) libraries
      m4_if([$1], [GCJ], [],
	[_LT_TAGVAR(lt_prog_compiler_pic, $1)='-DDLL_EXPORT'])
      case $host_os in
      os2*)
	_LT_TAGVAR(lt_prog_compiler_static, $1)='$wl-static'
	;;
      esac
      ;;
    darwin* | rhapsody*)
      # PIC is the default on this platform
      # Common symbols not allowed in MH_DYLIB files
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fno-common'
      ;;
    *djgpp*)
      # DJGPP does not support shared libraries at all
      _LT_TAGVAR(lt_prog_compiler_pic, $1)=
      ;;
    haiku*)
      # PIC is the default for Haiku.
      # The "-static" flag exists, but is broken.
      _LT_TAGVAR(lt_prog_compiler_static, $1)=
      ;;
    interix[[3-9]]*)
      # Interix 3.x gcc -fpic/-fPIC options generate broken code.
      # Instead, we relocate shared libraries at runtime.
      ;;
    sysv4*MP*)
      if test -d /usr/nec; then
	_LT_TAGVAR(lt_prog_compiler_pic, $1)=-Kconform_pic
      fi
      ;;
    hpux*)
      # PIC is the default for 64-bit PA HP-UX, but not for 32-bit
      # PA HP-UX.  On IA64 HP-UX, PIC is the default but the pic flag
      # sets the default TLS model and affects inlining.
      case $host_cpu in
      hppa*64*)
	;;
      *)
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	;;
      esac
      ;;
    *qnx* | *nto*)
      # QNX uses GNU C++, but need to define -shared option too, otherwise
      # it will coredump.
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC -shared'
      ;;
    *)
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
      ;;
    esac
  else
    case $host_os in
      aix[[4-9]]*)
	# All AIX code is PIC.
	if test ia64 = "$host_cpu"; then
	  # AIX 5 now supports IA64 processor
	  _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	else
	  _LT_TAGVAR(lt_prog_compiler_static, $1)='-bnso -bI:/lib/syscalls.exp'
	fi
	;;
      chorus*)
	case $cc_basename in
	cxch68*)
	  # Green Hills C++ Compiler
	  # _LT_TAGVAR(lt_prog_compiler_static, $1)="--no_auto_instantiation -u __main -u __premain -u _abort -r $COOL_DIR/lib/libOrb.a $MVME_DIR/lib/CC/libC.a $MVME_DIR/lib/classix/libcx.s.a"
	  ;;
	esac
	;;
      mingw* | cygwin* | os2* | pw32* | cegcc*)
	# This hack is so that the source file can tell whether it is being
	# built for inclusion in a dll (and should export symbols for example).
	m4_if([$1], [GCJ], [],
	  [_LT_TAGVAR(lt_prog_compiler_pic, $1)='-DDLL_EXPORT'])
	;;
      dgux*)
	case $cc_basename in
	  ec++*)
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    ;;
	  ghcx*)
	    # Green Hills C++ Compiler
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    ;;
	  *)
	    ;;
	esac
	;;
      freebsd* | dragonfly*)
	# FreeBSD uses GNU C++
	;;
      hpux9* | hpux10* | hpux11*)
	case $cc_basename in
	  CC*)
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='$wl-a ${wl}archive'
	    if test ia64 != "$host_cpu"; then
	      _LT_TAGVAR(lt_prog_compiler_pic, $1)='+Z'
	    fi
	    ;;
	  aCC*)
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='$wl-a ${wl}archive'
	    case $host_cpu in
	    hppa*64*|ia64*)
	      # +Z the default
	      ;;
	    *)
	      _LT_TAGVAR(lt_prog_compiler_pic, $1)='+Z'
	      ;;
	    esac
	    ;;
	  *)
	    ;;
	esac
	;;
      interix*)
	# This is c89, which is MS Visual C++ (no shared libs)
	# Anyone wants to do a port?
	;;
      irix5* | irix6* | nonstopux*)
	case $cc_basename in
	  CC*)
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
	    # CC pic flag -KPIC is the default.
	    ;;
	  *)
	    ;;
	esac
	;;
      linux* | k*bsd*-gnu | kopensolaris*-gnu | gnu*)
	case $cc_basename in
	  KCC*)
	    # KAI C++ Compiler
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='--backend -Wl,'
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	    ;;
	  ecpc* )
	    # old Intel C++ for x86_64, which still supported -KPIC.
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-static'
	    ;;
	  icpc* )
	    # Intel C++, used to be incompatible with GCC.
	    # ICC 10 doesn't accept -KPIC any more.
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-static'
	    ;;
	  pgCC* | pgcpp*)
	    # Portland Group C++ compiler
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fpic'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    ;;
	  cxx*)
	    # Compaq C++
	    # Make sure the PIC flag is empty.  It appears that all Alpha
	    # Linux and Compaq Tru64 Unix objects are PIC.
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)=
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
	    ;;
	  xlc* | xlC* | bgxl[[cC]]* | mpixl[[cC]]*)
	    # IBM XL 8.0, 9.0 on PPC and BlueGene
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-qpic'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-qstaticlink'
	    ;;
	  *)
	    case `$CC -V 2>&1 | sed 5q` in
	    *Sun\ C*)
	      # Sun C++ 5.9
	      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	      _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld '
	      ;;
	    esac
	    ;;
	esac
	;;
      lynxos*)
	;;
      m88k*)
	;;
      mvs*)
	case $cc_basename in
	  cxx*)
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-W c,exportall'
	    ;;
	  *)
	    ;;
	esac
	;;
      netbsd*)
	;;
      *qnx* | *nto*)
        # QNX uses GNU C++, but need to define -shared option too, otherwise
        # it will coredump.
        _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC -shared'
        ;;
      osf3* | osf4* | osf5*)
	case $cc_basename in
	  KCC*)
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='--backend -Wl,'
	    ;;
	  RCC*)
	    # Rational C++ 2.4.1
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    ;;
	  cxx*)
	    # Digital/Compaq C++
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    # Make sure the PIC flag is empty.  It appears that all Alpha
	    # Linux and Compaq Tru64 Unix objects are PIC.
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)=
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
	    ;;
	  *)
	    ;;
	esac
	;;
      psos*)
	;;
      solaris*)
	case $cc_basename in
	  CC* | sunCC*)
	    # Sun C++ 4.2, 5.x and Centerline C++
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld '
	    ;;
	  gcx*)
	    # Green Hills C++ Compiler
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-PIC'
	    ;;
	  *)
	    ;;
	esac
	;;
      sunos4*)
	case $cc_basename in
	  CC*)
	    # Sun C++ 4.x
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    ;;
	  lcc*)
	    # Lucid
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
	    ;;
	  *)
	    ;;
	esac
	;;
      sysv5* | unixware* | sco3.2v5* | sco5v6* | OpenUNIX*)
	case $cc_basename in
	  CC*)
	    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	    ;;
	esac
	;;
      tandem*)
	case $cc_basename in
	  NCC*)
	    # NonStop-UX NCC 3.20
	    _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	    ;;
	  *)
	    ;;
	esac
	;;
      vxworks*)
	;;
      *)
	_LT_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
	;;
    esac
  fi
],
[
  if test yes = "$GCC"; then
    _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
    _LT_TAGVAR(lt_prog_compiler_static, $1)='-static'

    case $host_os in
      aix*)
      # All AIX code is PIC.
      if test ia64 = "$host_cpu"; then
	# AIX 5 now supports IA64 processor
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      fi
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
      ;;

    amigaos*)
      case $host_cpu in
      powerpc)
            # see comment about AmigaOS4 .so support
            _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
        ;;
      m68k)
            # FIXME: we need at least 68020 code to build shared libraries, but
            # adding the '-m68020' flag to GCC prevents building anything better,
            # like '-m68040'.
            _LT_TAGVAR(lt_prog_compiler_pic, $1)='-m68020 -resident32 -malways-restore-a4'
        ;;
      esac
      ;;

    beos* | irix5* | irix6* | nonstopux* | osf3* | osf4* | osf5*)
      # PIC is the default for these OSes.
      ;;

    mingw* | cygwin* | pw32* | os2* | cegcc*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      # Although the cygwin gcc ignores -fPIC, still need this for old-style
      # (--disable-auto-import) libraries
      m4_if([$1], [GCJ], [],
	[_LT_TAGVAR(lt_prog_compiler_pic, $1)='-DDLL_EXPORT'])
      case $host_os in
      os2*)
	_LT_TAGVAR(lt_prog_compiler_static, $1)='$wl-static'
	;;
      esac
      ;;

    darwin* | rhapsody*)
      # PIC is the default on this platform
      # Common symbols not allowed in MH_DYLIB files
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fno-common'
      ;;

    haiku*)
      # PIC is the default for Haiku.
      # The "-static" flag exists, but is broken.
      _LT_TAGVAR(lt_prog_compiler_static, $1)=
      ;;

    hpux*)
      # PIC is the default for 64-bit PA HP-UX, but not for 32-bit
      # PA HP-UX.  On IA64 HP-UX, PIC is the default but the pic flag
      # sets the default TLS model and affects inlining.
      case $host_cpu in
      hppa*64*)
	# +Z the default
	;;
      *)
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	;;
      esac
      ;;

    interix[[3-9]]*)
      # Interix 3.x gcc -fpic/-fPIC options generate broken code.
      # Instead, we relocate shared libraries at runtime.
      ;;

    msdosdjgpp*)
      # Just because we use GCC doesn't mean we suddenly get shared libraries
      # on systems that don't support them.
      _LT_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
      enable_shared=no
      ;;

    *nto* | *qnx*)
      # QNX uses GNU C++, but need to define -shared option too, otherwise
      # it will coredump.
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC -shared'
      ;;

    sysv4*MP*)
      if test -d /usr/nec; then
	_LT_TAGVAR(lt_prog_compiler_pic, $1)=-Kconform_pic
      fi
      ;;

    *)
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
      ;;
    esac

    case $cc_basename in
    nvcc*) # Cuda Compiler Driver 2.2
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Xlinker '
      if test -n "$_LT_TAGVAR(lt_prog_compiler_pic, $1)"; then
        _LT_TAGVAR(lt_prog_compiler_pic, $1)="-Xcompiler $_LT_TAGVAR(lt_prog_compiler_pic, $1)"
      fi
      ;;
    esac
  else
    # PORTME Check for flag to pass linker flags through the system compiler.
    case $host_os in
    aix*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      if test ia64 = "$host_cpu"; then
	# AIX 5 now supports IA64 processor
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      else
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-bnso -bI:/lib/syscalls.exp'
      fi
      ;;

    darwin* | rhapsody*)
      # PIC is the default on this platform
      # Common symbols not allowed in MH_DYLIB files
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fno-common'
      case $cc_basename in
      nagfor*)
        # NAG Fortran compiler
        _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,-Wl,,'
        _LT_TAGVAR(lt_prog_compiler_pic, $1)='-PIC'
        _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
        ;;
      esac
      ;;

    mingw* | cygwin* | pw32* | os2* | cegcc*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      m4_if([$1], [GCJ], [],
	[_LT_TAGVAR(lt_prog_compiler_pic, $1)='-DDLL_EXPORT'])
      case $host_os in
      os2*)
	_LT_TAGVAR(lt_prog_compiler_static, $1)='$wl-static'
	;;
      esac
      ;;

    hpux9* | hpux10* | hpux11*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      # PIC is the default for IA64 HP-UX and 64-bit HP-UX, but
      # not for PA HP-UX.
      case $host_cpu in
      hppa*64*|ia64*)
	# +Z the default
	;;
      *)
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='+Z'
	;;
      esac
      # Is there a better lt_prog_compiler_static that works with the bundled CC?
      _LT_TAGVAR(lt_prog_compiler_static, $1)='$wl-a ${wl}archive'
      ;;

    irix5* | irix6* | nonstopux*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      # PIC (with -KPIC) is the default.
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
      ;;

    linux* | k*bsd*-gnu | kopensolaris*-gnu | gnu*)
      case $cc_basename in
      # old Intel for x86_64, which still supported -KPIC.
      ecc*)
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-static'
        ;;
      # icc used to be incompatible with GCC.
      # ICC 10 doesn't accept -KPIC any more.
      icc* | ifort*)
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-static'
        ;;
      # Lahey Fortran 8.1.
      lf95*)
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='--shared'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='--static'
	;;
      nagfor*)
	# NAG Fortran compiler
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,-Wl,,'
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-PIC'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	;;
      tcc*)
	# Fabrice Bellard et al's Tiny C Compiler
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-static'
	;;
      pgcc* | pgf77* | pgf90* | pgf95* | pgfortran*)
        # Portland Group compilers (*not* the Pentium gcc compiler,
	# which looks to be a dead project)
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-fpic'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
        ;;
      ccc*)
        _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
        # All Alpha code is PIC.
        _LT_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
        ;;
      xl* | bgxl* | bgf* | mpixl*)
	# IBM XL C 8.0/Fortran 10.1, 11.1 on PPC and BlueGene
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-qpic'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-qstaticlink'
	;;
      *)
	case `$CC -V 2>&1 | sed 5q` in
	*Sun\ Ceres\ Fortran* | *Sun*Fortran*\ [[1-7]].* | *Sun*Fortran*\ 8.[[0-3]]*)
	  # Sun Fortran 8.3 passes all unrecognized flags to the linker
	  _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	  _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	  _LT_TAGVAR(lt_prog_compiler_wl, $1)=''
	  ;;
	*Sun\ F* | *Sun*Fortran*)
	  _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	  _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	  _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld '
	  ;;
	*Sun\ C*)
	  # Sun C 5.9
	  _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
	  _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	  _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	  ;;
        *Intel*\ [[CF]]*Compiler*)
	  _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	  _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC'
	  _LT_TAGVAR(lt_prog_compiler_static, $1)='-static'
	  ;;
	*Portland\ Group*)
	  _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
	  _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fpic'
	  _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
	  ;;
	esac
	;;
      esac
      ;;

    newsos6)
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    *nto* | *qnx*)
      # QNX uses GNU C++, but need to define -shared option too, otherwise
      # it will coredump.
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-fPIC -shared'
      ;;

    osf3* | osf4* | osf5*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      # All OSF/1 code is PIC.
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
      ;;

    rdos*)
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-non_shared'
      ;;

    solaris*)
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      case $cc_basename in
      f77* | f90* | f95* | sunf77* | sunf90* | sunf95*)
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld ';;
      *)
	_LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,';;
      esac
      ;;

    sunos4*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Qoption ld '
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-PIC'
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    sysv4 | sysv4.2uw2* | sysv4.3*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    sysv4*MP*)
      if test -d /usr/nec; then
	_LT_TAGVAR(lt_prog_compiler_pic, $1)='-Kconform_pic'
	_LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      fi
      ;;

    sysv5* | unixware* | sco3.2v5* | sco5v6* | OpenUNIX*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-KPIC'
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    unicos*)
      _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Wl,'
      _LT_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
      ;;

    uts4*)
      _LT_TAGVAR(lt_prog_compiler_pic, $1)='-pic'
      _LT_TAGVAR(lt_prog_compiler_static, $1)='-Bstatic'
      ;;

    *)
      _LT_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no
      ;;
    esac
  fi
])
case $host_os in
  # For platforms that do not support PIC, -DPIC is meaningless:
  *djgpp*)
    _LT_TAGVAR(lt_prog_compiler_pic, $1)=
    ;;
  *)
    _LT_TAGVAR(lt_prog_compiler_pic, $1)="$_LT_TAGVAR(lt_prog_compiler_pic, $1)@&t@m4_if([$1],[],[ -DPIC],[m4_if([$1],[CXX],[ -DPIC],[])])"
    ;;
esac

AC_CACHE_CHECK([for $compiler option to produce PIC],
  [_LT_TAGVAR(lt_cv_prog_compiler_pic, $1)],
  [_LT_TAGVAR(lt_cv_prog_compiler_pic, $1)=$_LT_TAGVAR(lt_prog_compiler_pic, $1)])
_LT_TAGVAR(lt_prog_compiler_pic, $1)=$_LT_TAGVAR(lt_cv_prog_compiler_pic, $1)

#
# Check to make sure the PIC flag actually works.
#
if test -n "$_LT_TAGVAR(lt_prog_compiler_pic, $1)"; then
  _LT_COMPILER_OPTION([if $compiler PIC flag $_LT_TAGVAR(lt_prog_compiler_pic, $1) works],
    [_LT_TAGVAR(lt_cv_prog_compiler_pic_works, $1)],
    [$_LT_TAGVAR(lt_prog_compiler_pic, $1)@&t@m4_if([$1],[],[ -DPIC],[m4_if([$1],[CXX],[ -DPIC],[])])], [],
    [case $_LT_TAGVAR(lt_prog_compiler_pic, $1) in
     "" | " "*) ;;
     *) _LT_TAGVAR(lt_prog_compiler_pic, $1)=" $_LT_TAGVAR(lt_prog_compiler_pic, $1)" ;;
     esac],
    [_LT_TAGVAR(lt_prog_compiler_pic, $1)=
     _LT_TAGVAR(lt_prog_compiler_can_build_shared, $1)=no])
fi
_LT_TAGDECL([pic_flag], [lt_prog_compiler_pic], [1],
	[Additional compiler flags for building library objects])

_LT_TAGDECL([wl], [lt_prog_compiler_wl], [1],
	[How to pass a linker flag through the compiler])
#
# Check to make sure the static flag actually works.
#
wl=$_LT_TAGVAR(lt_prog_compiler_wl, $1) eval lt_tmp_static_flag=\"$_LT_TAGVAR(lt_prog_compiler_static, $1)\"
_LT_LINKER_OPTION([if $compiler static flag $lt_tmp_static_flag works],
  _LT_TAGVAR(lt_cv_prog_compiler_static_works, $1),
  $lt_tmp_static_flag,
  [],
  [_LT_TAGVAR(lt_prog_compiler_static, $1)=])
_LT_TAGDECL([link_static_flag], [lt_prog_compiler_static], [1],
	[Compiler flag to prevent dynamic linking])
])# _LT_COMPILER_PIC


# _LT_LINKER_SHLIBS([TAGNAME])
# ----------------------------
# See if the linker supports building shared libraries.
m4_defun([_LT_LINKER_SHLIBS],
[AC_REQUIRE([LT_PATH_LD])dnl
AC_REQUIRE([LT_PATH_NM])dnl
m4_require([_LT_PATH_MANIFEST_TOOL])dnl
m4_require([_LT_FILEUTILS_DEFAULTS])dnl
m4_require([_LT_DECL_EGREP])dnl
m4_require([_LT_DECL_SED])dnl
m4_require([_LT_CMD_GLOBAL_SYMBOLS])dnl
m4_require([_LT_TAG_COMPILER])dnl
AC_MSG_CHECKING([whether the $compiler linker ($LD) supports shared libraries])
m4_if([$1], [CXX], [
  _LT_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED '\''s/.* //'\'' | sort | uniq > $export_symbols'
  _LT_TAGVAR(exclude_expsyms, $1)=['_GLOBAL_OFFSET_TABLE_|_GLOBAL__F[ID]_.*']
  case $host_os in
  aix[[4-9]]*)
    # If we're using GNU nm, then we don't want the "-C" option.
    # -C means demangle to GNU nm, but means don't demangle to AIX nm.
    # Without the "-l" option, or with the "-B" option, AIX nm treats
    # weak defined symbols like other global defined symbols, whereas
    # GNU nm marks them as "W".
    # While the 'weak' keyword is ignored in the Export File, we need
    # it in the Import File for the 'aix-soname' feature, so we have
    # to replace the "-B" option with "-P" for AIX nm.
    if $NM -V 2>&1 | $GREP 'GNU' > /dev/null; then
      _LT_TAGVAR(export_symbols_cmds, $1)='$NM -Bpg $libobjs $convenience | awk '\''{ if (((\$ 2 == "T") || (\$ 2 == "D") || (\$ 2 == "B") || (\$ 2 == "W")) && ([substr](\$ 3,1,1) != ".")) { if (\$ 2 == "W") { print \$ 3 " weak" } else { print \$ 3 } } }'\'' | sort -u > $export_symbols'
    else
      _LT_TAGVAR(export_symbols_cmds, $1)='`func_echo_all $NM | $SED -e '\''s/B\([[^B]]*\)$/P\1/'\''` -PCpgl $libobjs $convenience | awk '\''{ if (((\$ 2 == "T") || (\$ 2 == "D") || (\$ 2 == "B") || (\$ 2 == "W") || (\$ 2 == "V") || (\$ 2 == "Z")) && ([substr](\$ 1,1,1) != ".")) { if ((\$ 2 == "W") || (\$ 2 == "V") || (\$ 2 == "Z")) { print \$ 1 " weak" } else { print \$ 1 } } }'\'' | sort -u > $export_symbols'
    fi
    ;;
  pw32*)
    _LT_TAGVAR(export_symbols_cmds, $1)=$ltdll_cmds
    ;;
  cygwin* | mingw* | cegcc*)
    case $cc_basename in
    cl* | icl* | ifort*)
      _LT_TAGVAR(exclude_expsyms, $1)='_NULL_IMPORT_DESCRIPTOR|_IMPORT_DESCRIPTOR_.*'
      ;;
    *)
      _LT_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED -e '\''/^[[BCDGRS]][[ ]]/s/.*[[ ]]\([[^ ]]*\)/\1 DATA/;s/^.*[[ ]]__nm__\([[^ ]]*\)[[ ]][[^ ]]*/\1 DATA/;/^I[[ ]]/d;/^[[AITW]][[ ]]/s/.* //'\'' | sort | uniq > $export_symbols'
      _LT_TAGVAR(exclude_expsyms, $1)=['[_]+GLOBAL_OFFSET_TABLE_|[_]+GLOBAL__[FID]_.*|[_]+head_[A-Za-z0-9_]+_dll|[A-Za-z0-9_]+_dll_iname']
      ;;
    esac
    ;;
  *)
    _LT_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED '\''s/.* //'\'' | sort | uniq > $export_symbols'
    ;;
  esac
], [
  runpath_var=
  _LT_TAGVAR(allow_undefined_flag, $1)=
  _LT_TAGVAR(always_export_symbols, $1)=no
  _LT_TAGVAR(archive_cmds, $1)=
  _LT_TAGVAR(archive_expsym_cmds, $1)=
  _LT_TAGVAR(compiler_needs_object, $1)=no
  _LT_TAGVAR(enable_shared_with_static_runtimes, $1)=no
  _LT_TAGVAR(export_dynamic_flag_spec, $1)=
  _LT_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED '\''s/.* //'\'' | sort | uniq > $export_symbols'
  _LT_TAGVAR(hardcode_automatic, $1)=no
  _LT_TAGVAR(hardcode_direct, $1)=no
  _LT_TAGVAR(hardcode_direct_absolute, $1)=no
  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)=
  _LT_TAGVAR(hardcode_libdir_separator, $1)=
  _LT_TAGVAR(hardcode_minus_L, $1)=no
  _LT_TAGVAR(hardcode_shlibpath_var, $1)=unsupported
  _LT_TAGVAR(inherit_rpath, $1)=no
  _LT_TAGVAR(link_all_deplibs, $1)=unknown
  _LT_TAGVAR(module_cmds, $1)=
  _LT_TAGVAR(module_expsym_cmds, $1)=
  _LT_TAGVAR(old_archive_from_new_cmds, $1)=
  _LT_TAGVAR(old_archive_from_expsyms_cmds, $1)=
  _LT_TAGVAR(thread_safe_flag_spec, $1)=
  _LT_TAGVAR(whole_archive_flag_spec, $1)=
  # include_expsyms should be a list of space-separated symbols to be *always*
  # included in the symbol list
  _LT_TAGVAR(include_expsyms, $1)=
  # exclude_expsyms can be an extended regexp of symbols to exclude
  # it will be wrapped by ' (' and ')$', so one must not match beginning or
  # end of line.  Example: 'a|bc|.*d.*' will exclude the symbols 'a' and 'bc',
  # as well as any symbol that contains 'd'.
  _LT_TAGVAR(exclude_expsyms, $1)=['_GLOBAL_OFFSET_TABLE_|_GLOBAL__F[ID]_.*']
  # Although _GLOBAL_OFFSET_TABLE_ is a valid symbol C name, most a.out
  # platforms (ab)use it in PIC code, but their linkers get confused if
  # the symbol is explicitly referenced.  Since portable code cannot
  # rely on this symbol name, it's probably fine to never include it in
  # preloaded symbol tables.
  # Exclude shared library initialization/finalization symbols.
dnl Note also adjust exclude_expsyms for C++ above.
  extract_expsyms_cmds=

  case $host_os in
  cygwin* | mingw* | pw32* | cegcc*)
    # FIXME: the MSVC++ port hasn't been tested in a loooong time
    # When not using gcc, we currently assume that we are using
    # Microsoft Visual C++.
    if test yes != "$GCC"; then
      with_gnu_ld=no
    fi
    ;;
  interix*)
    # we just hope/assume this is gcc and not c89 (= MSVC++)
    with_gnu_ld=yes
    ;;
  openbsd* | bitrig*)
    with_gnu_ld=no
    ;;
  esac

  _LT_TAGVAR(ld_shlibs, $1)=yes

  # On some targets, GNU ld is compatible enough with the native linker
  # that we're better off using the native interface for both.
  lt_use_gnu_ld_interface=no
  if test yes = "$with_gnu_ld"; then
    case $host_os in
      aix*)
	# The AIX port of GNU ld has always aspired to compatibility
	# with the native linker.  However, as the warning in the GNU ld
	# block says, versions before 2.19.5* couldn't really create working
	# shared libraries, regardless of the interface used.
	case `$LD -v 2>&1` in
	  *\ \(GNU\ Binutils\)\ 2.19.5*) ;;
	  *\ \(GNU\ Binutils\)\ 2.[[2-9]]*) ;;
	  *\ \(GNU\ Binutils\)\ [[3-9]]*) ;;
	  *)
	    lt_use_gnu_ld_interface=yes
	    ;;
	esac
	;;
      *)
	lt_use_gnu_ld_interface=yes
	;;
    esac
  fi

  if test yes = "$lt_use_gnu_ld_interface"; then
    # If archive_cmds runs LD, not CC, wlarc should be empty
    wlarc='$wl'

    # Set some defaults for GNU ld with shared library support. These
    # are reset later if shared libraries are not supported. Putting them
    # here allows them to be overridden if necessary.
    runpath_var=LD_RUN_PATH
    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
    _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-dynamic'
    # ancient GNU ld didn't support --whole-archive et. al.
    if $LD --help 2>&1 | $GREP 'no-whole-archive' > /dev/null; then
      _LT_TAGVAR(whole_archive_flag_spec, $1)=$wlarc'--whole-archive$convenience '$wlarc'--no-whole-archive'
    else
      _LT_TAGVAR(whole_archive_flag_spec, $1)=
    fi
    supports_anon_versioning=no
    case `$LD -v | $SED -e 's/([^)]\+)\s\+//' 2>&1` in
      *GNU\ gold*) supports_anon_versioning=yes ;;
      *\ [[01]].* | *\ 2.[[0-9]].* | *\ 2.10.*) ;; # catch versions < 2.11
      *\ 2.11.93.0.2\ *) supports_anon_versioning=yes ;; # RH7.3 ...
      *\ 2.11.92.0.12\ *) supports_anon_versioning=yes ;; # Mandrake 8.2 ...
      *\ 2.11.*) ;; # other 2.11 versions
      *) supports_anon_versioning=yes ;;
    esac

    # See if GNU ld supports shared libraries.
    case $host_os in
    aix[[3-9]]*)
      # On AIX/PPC, the GNU linker is very broken
      if test ia64 != "$host_cpu"; then
	_LT_TAGVAR(ld_shlibs, $1)=no
	cat <<_LT_EOF 1>&2

*** Warning: the GNU linker, at least up to release 2.19, is reported
*** to be unable to reliably create shared libraries on AIX.
*** Therefore, libtool is disabling shared libraries support.  If you
*** really care for shared libraries, you may want to install binutils
*** 2.20 or above, or modify your PATH so that a non-GNU linker is found.
*** You will then need to restart the configuration process.

_LT_EOF
      fi
      ;;

    amigaos*)
      case $host_cpu in
      powerpc)
            # see comment about AmigaOS4 .so support
            _LT_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
            _LT_TAGVAR(archive_expsym_cmds, $1)=''
        ;;
      m68k)
            _LT_TAGVAR(archive_cmds, $1)='$RM $output_objdir/a2ixlibrary.data~$ECHO "#define NAME $libname" > $output_objdir/a2ixlibrary.data~$ECHO "#define LIBRARY_ID 1" >> $output_objdir/a2ixlibrary.data~$ECHO "#define VERSION $major" >> $output_objdir/a2ixlibrary.data~$ECHO "#define REVISION $revision" >> $output_objdir/a2ixlibrary.data~$AR $AR_FLAGS $lib $libobjs~$RANLIB $lib~(cd $output_objdir && a2ixlibrary -32)'
            _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
            _LT_TAGVAR(hardcode_minus_L, $1)=yes
        ;;
      esac
      ;;

    beos*)
      if $LD --help 2>&1 | $GREP ': supported targets:.* elf' > /dev/null; then
	_LT_TAGVAR(allow_undefined_flag, $1)=unsupported
	# Joseph Beckenbach <jrb3@best.com> says some releases of gcc
	# support --undefined.  This deserves some investigation.  FIXME
	_LT_TAGVAR(archive_cmds, $1)='$CC -nostart $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
      else
	_LT_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    cygwin* | mingw* | pw32* | cegcc*)
      # _LT_TAGVAR(hardcode_libdir_flag_spec, $1) is actually meaningless,
      # as there is no search path for DLLs.
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-all-symbols'
      _LT_TAGVAR(allow_undefined_flag, $1)=unsupported
      _LT_TAGVAR(always_export_symbols, $1)=no
      _LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
      _LT_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED -e '\''/^[[BCDGRS]][[ ]]/s/.*[[ ]]\([[^ ]]*\)/\1 DATA/;s/^.*[[ ]]__nm__\([[^ ]]*\)[[ ]][[^ ]]*/\1 DATA/;/^I[[ ]]/d;/^[[AITW]][[ ]]/s/.* //'\'' | sort | uniq > $export_symbols'
      _LT_TAGVAR(exclude_expsyms, $1)=['[_]+GLOBAL_OFFSET_TABLE_|[_]+GLOBAL__[FID]_.*|[_]+head_[A-Za-z0-9_]+_dll|[A-Za-z0-9_]+_dll_iname']

      if $LD --help 2>&1 | $GREP 'auto-import' > /dev/null; then
        _LT_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags -o $output_objdir/$soname $wl--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
	# If the export-symbols file already is a .def file, use it as
	# is; otherwise, prepend EXPORTS...
	_LT_TAGVAR(archive_expsym_cmds, $1)='if _LT_DLL_DEF_P([$export_symbols]); then
          cp $export_symbols $output_objdir/$soname.def;
        else
          echo EXPORTS > $output_objdir/$soname.def;
          cat $export_symbols >> $output_objdir/$soname.def;
        fi~
        $CC -shared $output_objdir/$soname.def $libobjs $deplibs $compiler_flags -o $output_objdir/$soname $wl--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
      else
	_LT_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    haiku*)
      _LT_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
      _LT_TAGVAR(link_all_deplibs, $1)=yes
      ;;

    os2*)
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_TAGVAR(allow_undefined_flag, $1)=unsupported
      shrext_cmds=.dll
      _LT_TAGVAR(archive_cmds, $1)='$ECHO "LIBRARY ${soname%$shared_ext} INITINSTANCE TERMINSTANCE" > $output_objdir/$libname.def~
	$ECHO "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~
	$ECHO "DATA MULTIPLE NONSHARED" >> $output_objdir/$libname.def~
	$ECHO EXPORTS >> $output_objdir/$libname.def~
	emxexp $libobjs | $SED /"_DLL_InitTerm"/d >> $output_objdir/$libname.def~
	$CC -Zdll -Zcrtdll -o $output_objdir/$soname $libobjs $deplibs $compiler_flags $output_objdir/$libname.def~
	emximp -o $lib $output_objdir/$libname.def'
      _LT_TAGVAR(archive_expsym_cmds, $1)='$ECHO "LIBRARY ${soname%$shared_ext} INITINSTANCE TERMINSTANCE" > $output_objdir/$libname.def~
	$ECHO "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~
	$ECHO "DATA MULTIPLE NONSHARED" >> $output_objdir/$libname.def~
	$ECHO EXPORTS >> $output_objdir/$libname.def~
	prefix_cmds="$SED"~
	if test EXPORTS = "`$SED 1q $export_symbols`"; then
	  prefix_cmds="$prefix_cmds -e 1d";
	fi~
	prefix_cmds="$prefix_cmds -e \"s/^\(.*\)$/_\1/g\""~
	cat $export_symbols | $prefix_cmds >> $output_objdir/$libname.def~
	$CC -Zdll -Zcrtdll -o $output_objdir/$soname $libobjs $deplibs $compiler_flags $output_objdir/$libname.def~
	emximp -o $lib $output_objdir/$libname.def'
      _LT_TAGVAR(old_archive_From_new_cmds, $1)='emximp -o $output_objdir/${libname}_dll.a $output_objdir/$libname.def'
      _LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
      ;;

    interix[[3-9]]*)
      _LT_TAGVAR(hardcode_direct, $1)=no
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
      _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
      # Hack: On Interix 3.x, we cannot compile PIC because of a broken gcc.
      # Instead, shared libraries are loaded at an image base (0x10000000 by
      # default) and relocated if they conflict, which is a slow very memory
      # consuming and fragmenting process.  To avoid this, we pick a random,
      # 256 KiB-aligned image base between 0x50000000 and 0x6FFC0000 at link
      # time.  Moving up from 0x10000000 also allows more sbrk(2) space.
      _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-h,$soname $wl--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
      _LT_TAGVAR(archive_expsym_cmds, $1)='sed "s|^|_|" $export_symbols >$output_objdir/$soname.expsym~$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-h,$soname $wl--retain-symbols-file,$output_objdir/$soname.expsym $wl--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
      ;;

    gnu* | linux* | tpf* | k*bsd*-gnu | kopensolaris*-gnu)
      tmp_diet=no
      if test linux-dietlibc = "$host_os"; then
	case $cc_basename in
	  diet\ *) tmp_diet=yes;;	# linux-dietlibc with static linking (!diet-dyn)
	esac
      fi
      if $LD --help 2>&1 | $EGREP ': supported targets:.* elf' > /dev/null \
	 && test no = "$tmp_diet"
      then
	tmp_addflag=' $pic_flag'
	tmp_sharedflag='-shared'
	case $cc_basename,$host_cpu in
        pgcc*)				# Portland Group C compiler
	  _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; func_echo_all \"$new_convenience\"` $wl--no-whole-archive'
	  tmp_addflag=' $pic_flag'
	  ;;
	pgf77* | pgf90* | pgf95* | pgfortran*)
					# Portland Group f77 and f90 compilers
	  _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; func_echo_all \"$new_convenience\"` $wl--no-whole-archive'
	  tmp_addflag=' $pic_flag -Mnomain' ;;
	ecc*,ia64* | icc*,ia64*)	# Intel C compiler on ia64
	  tmp_addflag=' -i_dynamic' ;;
	efc*,ia64* | ifort*,ia64*)	# Intel Fortran compiler on ia64
	  tmp_addflag=' -i_dynamic -nofor_main' ;;
	ifc* | ifort*)			# Intel Fortran compiler
	  tmp_addflag=' -nofor_main' ;;
	lf95*)				# Lahey Fortran 8.1
	  _LT_TAGVAR(whole_archive_flag_spec, $1)=
	  tmp_sharedflag='--shared' ;;
        nagfor*)                        # NAGFOR 5.3
          tmp_sharedflag='-Wl,-shared' ;;
	xl[[cC]]* | bgxl[[cC]]* | mpixl[[cC]]*) # IBM XL C 8.0 on PPC (deal with xlf below)
	  tmp_sharedflag='-qmkshrobj'
	  tmp_addflag= ;;
	nvcc*)	# Cuda Compiler Driver 2.2
	  _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; func_echo_all \"$new_convenience\"` $wl--no-whole-archive'
	  _LT_TAGVAR(compiler_needs_object, $1)=yes
	  ;;
	esac
	case `$CC -V 2>&1 | sed 5q` in
	*Sun\ C*)			# Sun C 5.9
	  _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive`new_convenience=; for conv in $convenience\"\"; do test -z \"$conv\" || new_convenience=\"$new_convenience,$conv\"; done; func_echo_all \"$new_convenience\"` $wl--no-whole-archive'
	  _LT_TAGVAR(compiler_needs_object, $1)=yes
	  tmp_sharedflag='-G' ;;
	*Sun\ F*)			# Sun Fortran 8.3
	  tmp_sharedflag='-G' ;;
	esac
	_LT_TAGVAR(archive_cmds, $1)='$CC '"$tmp_sharedflag""$tmp_addflag"' $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'

        if test yes = "$supports_anon_versioning"; then
          _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $output_objdir/$libname.ver~
            cat $export_symbols | sed -e "s/\(.*\)/\1;/" >> $output_objdir/$libname.ver~
            echo "local: *; };" >> $output_objdir/$libname.ver~
            $CC '"$tmp_sharedflag""$tmp_addflag"' $libobjs $deplibs $compiler_flags $wl-soname $wl$soname $wl-version-script $wl$output_objdir/$libname.ver -o $lib'
        fi

	case $cc_basename in
	tcc*)
	  _LT_TAGVAR(export_dynamic_flag_spec, $1)='-rdynamic'
	  ;;
	xlf* | bgf* | bgxlf* | mpixlf*)
	  # IBM XL Fortran 10.1 on PPC cannot create shared libs itself
	  _LT_TAGVAR(whole_archive_flag_spec, $1)='--whole-archive$convenience --no-whole-archive'
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
	  _LT_TAGVAR(archive_cmds, $1)='$LD -shared $libobjs $deplibs $linker_flags -soname $soname -o $lib'
	  if test yes = "$supports_anon_versioning"; then
	    _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $output_objdir/$libname.ver~
              cat $export_symbols | sed -e "s/\(.*\)/\1;/" >> $output_objdir/$libname.ver~
              echo "local: *; };" >> $output_objdir/$libname.ver~
              $LD -shared $libobjs $deplibs $linker_flags -soname $soname -version-script $output_objdir/$libname.ver -o $lib'
	  fi
	  ;;
	esac
      else
        _LT_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    netbsd*)
      if echo __ELF__ | $CC -E - | $GREP __ELF__ >/dev/null; then
	_LT_TAGVAR(archive_cmds, $1)='$LD -Bshareable $libobjs $deplibs $linker_flags -o $lib'
	wlarc=
      else
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
      fi
      ;;

    solaris*)
      if $LD -v 2>&1 | $GREP 'BFD 2\.8' > /dev/null; then
	_LT_TAGVAR(ld_shlibs, $1)=no
	cat <<_LT_EOF 1>&2

*** Warning: The releases 2.8.* of the GNU linker cannot reliably
*** create shared libraries on Solaris systems.  Therefore, libtool
*** is disabling shared libraries support.  We urge you to upgrade GNU
*** binutils to release 2.9.1 or newer.  Another option is to modify
*** your PATH or compiler configuration so that the native linker is
*** used, and then restart.

_LT_EOF
      elif $LD --help 2>&1 | $GREP ': supported targets:.* elf' > /dev/null; then
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
      else
	_LT_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    sysv5* | sco3.2v5* | sco5v6* | unixware* | OpenUNIX*)
      case `$LD -v 2>&1` in
        *\ [[01]].* | *\ 2.[[0-9]].* | *\ 2.1[[0-5]].*)
	_LT_TAGVAR(ld_shlibs, $1)=no
	cat <<_LT_EOF 1>&2

*** Warning: Releases of the GNU linker prior to 2.16.91.0.3 cannot
*** reliably create shared libraries on SCO systems.  Therefore, libtool
*** is disabling shared libraries support.  We urge you to upgrade GNU
*** binutils to release 2.16.91.0.3 or newer.  Another option is to modify
*** your PATH or compiler configuration so that the native linker is
*** used, and then restart.

_LT_EOF
	;;
	*)
	  # For security reasons, it is highly recommended that you always
	  # use absolute paths for naming shared libraries, and exclude the
	  # DT_RUNPATH tag from executables and libraries.  But doing so
	  # requires that you compile everything twice, which is a pain.
	  if $LD --help 2>&1 | $GREP ': supported targets:.* elf' > /dev/null; then
	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
	    _LT_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
	    _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
	  else
	    _LT_TAGVAR(ld_shlibs, $1)=no
	  fi
	;;
      esac
      ;;

    sunos4*)
      _LT_TAGVAR(archive_cmds, $1)='$LD -assert pure-text -Bshareable -o $lib $libobjs $deplibs $linker_flags'
      wlarc=
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    *)
      if $LD --help 2>&1 | $GREP ': supported targets:.* elf' > /dev/null; then
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
      else
	_LT_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;
    esac

    if test no = "$_LT_TAGVAR(ld_shlibs, $1)"; then
      runpath_var=
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)=
      _LT_TAGVAR(export_dynamic_flag_spec, $1)=
      _LT_TAGVAR(whole_archive_flag_spec, $1)=
    fi
  else
    # PORTME fill in a description of your system's linker (not GNU ld)
    case $host_os in
    aix3*)
      _LT_TAGVAR(allow_undefined_flag, $1)=unsupported
      _LT_TAGVAR(always_export_symbols, $1)=yes
      _LT_TAGVAR(archive_expsym_cmds, $1)='$LD -o $output_objdir/$soname $libobjs $deplibs $linker_flags -bE:$export_symbols -T512 -H512 -bM:SRE~$AR $AR_FLAGS $lib $output_objdir/$soname'
      # Note: this linker hardcodes the directories in LIBPATH if there
      # are no directories specified by -L.
      _LT_TAGVAR(hardcode_minus_L, $1)=yes
      if test yes = "$GCC" && test -z "$lt_prog_compiler_static"; then
	# Neither direct hardcoding nor static linking is supported with a
	# broken collect2.
	_LT_TAGVAR(hardcode_direct, $1)=unsupported
      fi
      ;;

    aix[[4-9]]*)
      if test ia64 = "$host_cpu"; then
	# On IA64, the linker does run time linking by default, so we don't
	# have to do anything special.
	aix_use_runtimelinking=no
	exp_sym_flag='-Bexport'
	no_entry_flag=
      else
	# If we're using GNU nm, then we don't want the "-C" option.
	# -C means demangle to GNU nm, but means don't demangle to AIX nm.
	# Without the "-l" option, or with the "-B" option, AIX nm treats
	# weak defined symbols like other global defined symbols, whereas
	# GNU nm marks them as "W".
	# While the 'weak' keyword is ignored in the Export File, we need
	# it in the Import File for the 'aix-soname' feature, so we have
	# to replace the "-B" option with "-P" for AIX nm.
	if $NM -V 2>&1 | $GREP 'GNU' > /dev/null; then
	  _LT_TAGVAR(export_symbols_cmds, $1)='$NM -Bpg $libobjs $convenience | awk '\''{ if (((\$ 2 == "T") || (\$ 2 == "D") || (\$ 2 == "B") || (\$ 2 == "W")) && ([substr](\$ 3,1,1) != ".")) { if (\$ 2 == "W") { print \$ 3 " weak" } else { print \$ 3 } } }'\'' | sort -u > $export_symbols'
	else
	  _LT_TAGVAR(export_symbols_cmds, $1)='`func_echo_all $NM | $SED -e '\''s/B\([[^B]]*\)$/P\1/'\''` -PCpgl $libobjs $convenience | awk '\''{ if (((\$ 2 == "T") || (\$ 2 == "D") || (\$ 2 == "B") || (\$ 2 == "W") || (\$ 2 == "V") || (\$ 2 == "Z")) && ([substr](\$ 1,1,1) != ".")) { if ((\$ 2 == "W") || (\$ 2 == "V") || (\$ 2 == "Z")) { print \$ 1 " weak" } else { print \$ 1 } } }'\'' | sort -u > $export_symbols'
	fi
	aix_use_runtimelinking=no

	# Test if we are trying to use run time linking or normal
	# AIX style linking. If -brtl is somewhere in LDFLAGS, we
	# have runtime linking enabled, and use it for executables.
	# For shared libraries, we enable/disable runtime linking
	# depending on the kind of the shared library created -
	# when "with_aix_soname,aix_use_runtimelinking" is:
	# "aix,no"   lib.a(lib.so.V) shared, rtl:no,  for executables
	# "aix,yes"  lib.so          shared, rtl:yes, for executables
	#            lib.a           static archive
	# "both,no"  lib.so.V(shr.o) shared, rtl:yes
	#            lib.a(lib.so.V) shared, rtl:no,  for executables
	# "both,yes" lib.so.V(shr.o) shared, rtl:yes, for executables
	#            lib.a(lib.so.V) shared, rtl:no
	# "svr4,*"   lib.so.V(shr.o) shared, rtl:yes, for executables
	#            lib.a           static archive
	case $host_os in aix4.[[23]]|aix4.[[23]].*|aix[[5-9]]*)
	  for ld_flag in $LDFLAGS; do
	  if (test x-brtl = "x$ld_flag" || test x-Wl,-brtl = "x$ld_flag"); then
	    aix_use_runtimelinking=yes
	    break
	  fi
	  done
	  if test svr4,no = "$with_aix_soname,$aix_use_runtimelinking"; then
	    # With aix-soname=svr4, we create the lib.so.V shared archives only,
	    # so we don't have lib.a shared libs to link our executables.
	    # We have to force runtime linking in this case.
	    aix_use_runtimelinking=yes
	    LDFLAGS="$LDFLAGS -Wl,-brtl"
	  fi
	  ;;
	esac

	exp_sym_flag='-bexport'
	no_entry_flag='-bnoentry'
      fi

      # When large executables or shared objects are built, AIX ld can
      # have problems creating the table of contents.  If linking a library
      # or program results in "error TOC overflow" add -mminimal-toc to
      # CXXFLAGS/CFLAGS for g++/gcc.  In the cases where that is not
      # enough to fix the problem, add -Wl,-bbigtoc to LDFLAGS.

      _LT_TAGVAR(archive_cmds, $1)=''
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_direct_absolute, $1)=yes
      _LT_TAGVAR(hardcode_libdir_separator, $1)=':'
      _LT_TAGVAR(link_all_deplibs, $1)=yes
      _LT_TAGVAR(file_list_spec, $1)='$wl-f,'
      case $with_aix_soname,$aix_use_runtimelinking in
      aix,*) ;; # traditional, no import file
      svr4,* | *,yes) # use import file
	# The Import File defines what to hardcode.
	_LT_TAGVAR(hardcode_direct, $1)=no
	_LT_TAGVAR(hardcode_direct_absolute, $1)=no
	;;
      esac

      if test yes = "$GCC"; then
	case $host_os in aix4.[[012]]|aix4.[[012]].*)
	# We only want to do this on AIX 4.2 and lower, the check
	# below for broken collect2 doesn't work under 4.3+
	  collect2name=`$CC -print-prog-name=collect2`
	  if test -f "$collect2name" &&
	   strings "$collect2name" | $GREP resolve_lib_name >/dev/null
	  then
	  # We have reworked collect2
	  :
	  else
	  # We have old collect2
	  _LT_TAGVAR(hardcode_direct, $1)=unsupported
	  # It fails to find uninstalled libraries when the uninstalled
	  # path is not listed in the libpath.  Setting hardcode_minus_L
	  # to unsupported forces relinking
	  _LT_TAGVAR(hardcode_minus_L, $1)=yes
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
	  _LT_TAGVAR(hardcode_libdir_separator, $1)=
	  fi
	  ;;
	esac
	shared_flag='-shared'
	if test yes = "$aix_use_runtimelinking"; then
	  shared_flag="$shared_flag "'$wl-G'
	fi
	# Need to ensure runtime linking is disabled for the traditional
	# shared library, or the linker may eventually find shared libraries
	# /with/ Import File - we do not want to mix them.
	shared_flag_aix='-shared'
	shared_flag_svr4='-shared $wl-G'
      else
	# not using gcc
	if test ia64 = "$host_cpu"; then
	# VisualAge C++, Version 5.5 for AIX 5L for IA-64, Beta 3 Release
	# chokes on -Wl,-G. The following line is correct:
	  shared_flag='-G'
	else
	  if test yes = "$aix_use_runtimelinking"; then
	    shared_flag='$wl-G'
	  else
	    shared_flag='$wl-bM:SRE'
	  fi
	  shared_flag_aix='$wl-bM:SRE'
	  shared_flag_svr4='$wl-G'
	fi
      fi

      _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-bexpall'
      # It seems that -bexpall does not export symbols beginning with
      # underscore (_), so it is better to generate a list of symbols to export.
      _LT_TAGVAR(always_export_symbols, $1)=yes
      if test aix,yes = "$with_aix_soname,$aix_use_runtimelinking"; then
	# Warning - without using the other runtime loading flags (-brtl),
	# -berok will link without error, but may produce a broken library.
	_LT_TAGVAR(allow_undefined_flag, $1)='-berok'
        # Determine the default libpath from the value encoded in an
        # empty executable.
        _LT_SYS_MODULE_PATH_AIX([$1])
        _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-blibpath:$libdir:'"$aix_libpath"
        _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -o $output_objdir/$soname $libobjs $deplibs $wl'$no_entry_flag' $compiler_flags `if test -n "$allow_undefined_flag"; then func_echo_all "$wl$allow_undefined_flag"; else :; fi` $wl'$exp_sym_flag:\$export_symbols' '$shared_flag
      else
	if test ia64 = "$host_cpu"; then
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-R $libdir:/usr/lib:/lib'
	  _LT_TAGVAR(allow_undefined_flag, $1)="-z nodefs"
	  _LT_TAGVAR(archive_expsym_cmds, $1)="\$CC $shared_flag"' -o $output_objdir/$soname $libobjs $deplibs '"\$wl$no_entry_flag"' $compiler_flags $wl$allow_undefined_flag '"\$wl$exp_sym_flag:\$export_symbols"
	else
	 # Determine the default libpath from the value encoded in an
	 # empty executable.
	 _LT_SYS_MODULE_PATH_AIX([$1])
	 _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-blibpath:$libdir:'"$aix_libpath"
	  # Warning - without using the other run time loading flags,
	  # -berok will link without error, but may produce a broken library.
	  _LT_TAGVAR(no_undefined_flag, $1)=' $wl-bernotok'
	  _LT_TAGVAR(allow_undefined_flag, $1)=' $wl-berok'
	  if test yes = "$with_gnu_ld"; then
	    # We only use this code for GNU lds that support --whole-archive.
	    _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive$convenience $wl--no-whole-archive'
	  else
	    # Exported symbols can be pulled into shared objects from archives
	    _LT_TAGVAR(whole_archive_flag_spec, $1)='$convenience'
	  fi
	  _LT_TAGVAR(archive_cmds_need_lc, $1)=yes
	  _LT_TAGVAR(archive_expsym_cmds, $1)='$RM -r $output_objdir/$realname.d~$MKDIR $output_objdir/$realname.d'
	  # -brtl affects multiple linker settings, -berok does not and is overridden later
	  compiler_flags_filtered='`func_echo_all "$compiler_flags " | $SED -e "s%-brtl\\([[, ]]\\)%-berok\\1%g"`'
	  if test svr4 != "$with_aix_soname"; then
	    # This is similar to how AIX traditionally builds its shared libraries.
	    _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$CC '$shared_flag_aix' -o $output_objdir/$realname.d/$soname $libobjs $deplibs $wl-bnoentry '$compiler_flags_filtered'$wl-bE:$export_symbols$allow_undefined_flag~$AR $AR_FLAGS $output_objdir/$libname$release.a $output_objdir/$realname.d/$soname'
	  fi
	  if test aix != "$with_aix_soname"; then
	    _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$CC '$shared_flag_svr4' -o $output_objdir/$realname.d/$shared_archive_member_spec.o $libobjs $deplibs $wl-bnoentry '$compiler_flags_filtered'$wl-bE:$export_symbols$allow_undefined_flag~$STRIP -e $output_objdir/$realname.d/$shared_archive_member_spec.o~( func_echo_all "#! $soname($shared_archive_member_spec.o)"; if test shr_64 = "$shared_archive_member_spec"; then func_echo_all "# 64"; else func_echo_all "# 32"; fi; cat $export_symbols ) > $output_objdir/$realname.d/$shared_archive_member_spec.imp~$AR $AR_FLAGS $output_objdir/$soname $output_objdir/$realname.d/$shared_archive_member_spec.o $output_objdir/$realname.d/$shared_archive_member_spec.imp'
	  else
	    # used by -dlpreopen to get the symbols
	    _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$MV  $output_objdir/$realname.d/$soname $output_objdir'
	  fi
	  _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$RM -r $output_objdir/$realname.d'
	fi
      fi
      ;;

    amigaos*)
      case $host_cpu in
      powerpc)
            # see comment about AmigaOS4 .so support
            _LT_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
            _LT_TAGVAR(archive_expsym_cmds, $1)=''
        ;;
      m68k)
            _LT_TAGVAR(archive_cmds, $1)='$RM $output_objdir/a2ixlibrary.data~$ECHO "#define NAME $libname" > $output_objdir/a2ixlibrary.data~$ECHO "#define LIBRARY_ID 1" >> $output_objdir/a2ixlibrary.data~$ECHO "#define VERSION $major" >> $output_objdir/a2ixlibrary.data~$ECHO "#define REVISION $revision" >> $output_objdir/a2ixlibrary.data~$AR $AR_FLAGS $lib $libobjs~$RANLIB $lib~(cd $output_objdir && a2ixlibrary -32)'
            _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
            _LT_TAGVAR(hardcode_minus_L, $1)=yes
        ;;
      esac
      ;;

    bsdi[[45]]*)
      _LT_TAGVAR(export_dynamic_flag_spec, $1)=-rdynamic
      ;;

    cygwin* | mingw* | pw32* | cegcc*)
      # When not using gcc, we currently assume that we are using
      # Microsoft Visual C++ or Intel C++/Fortran Compiler.
      # hardcode_libdir_flag_spec is actually meaningless, as there is
      # no search path for DLLs.
      case $cc_basename in
      cl* | icl*| ifort*)
	# Native MSVC or ICC or IFORT
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)=' '
	_LT_TAGVAR(allow_undefined_flag, $1)=unsupported
	_LT_TAGVAR(always_export_symbols, $1)=yes
	_LT_TAGVAR(file_list_spec, $1)='@'
	# Tell ltmain to make .lib files, not .a files.
	libext=lib
	# Tell ltmain to make .dll files, not .so files.
	shrext_cmds=.dll
	# FIXME: Setting linknames here is a bad hack.
	_LT_TAGVAR(archive_cmds, $1)='$CC -o $output_objdir/$soname $libobjs $compiler_flags $deplibs -Wl,-DLL,-IMPLIB:"$tool_output_objdir$libname.dll.lib"~linknames='
	_LT_TAGVAR(archive_expsym_cmds, $1)='if _LT_DLL_DEF_P([$export_symbols]); then
            cp "$export_symbols" "$output_objdir/$soname.def";
            echo "$tool_output_objdir$soname.def" > "$output_objdir/$soname.exp";
          else
            $SED -e '\''s/^/-link -EXPORT:/'\'' < $export_symbols > $output_objdir/$soname.exp;
          fi~
          $CC -o $tool_output_objdir$soname $libobjs $compiler_flags $deplibs "@$tool_output_objdir$soname.exp" -Wl,-DLL,-IMPLIB:"$tool_output_objdir$libname.dll.lib"~
          linknames='
	# The linker will not automatically build a static lib if we build a DLL.
	# _LT_TAGVAR(old_archive_from_new_cmds, $1)='true'
	_LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
	_LT_TAGVAR(exclude_expsyms, $1)='_NULL_IMPORT_DESCRIPTOR|_IMPORT_DESCRIPTOR_.*'
	_LT_TAGVAR(export_symbols_cmds, $1)='$NM $libobjs $convenience | $global_symbol_pipe | $SED -e '\''/^[[BCDGRS]][[ ]]/s/.*[[ ]]\([[^ ]]*\)/\1,DATA/'\'' | $SED -e '\''/^[[AITW]][[ ]]/s/.*[[ ]]//'\'' | sort | uniq > $export_symbols'
	# Don't use ranlib
	_LT_TAGVAR(old_postinstall_cmds, $1)='chmod 644 $oldlib'
	_LT_TAGVAR(postlink_cmds, $1)='lt_outputfile="@OUTPUT@"~
          lt_tool_outputfile="@TOOL_OUTPUT@"~
          case $lt_outputfile in
            *.exe|*.EXE) ;;
            *)
              lt_outputfile=$lt_outputfile.exe
              lt_tool_outputfile=$lt_tool_outputfile.exe
              ;;
          esac~
          if test : != "$MANIFEST_TOOL" && test -f "$lt_outputfile.manifest"; then
            $MANIFEST_TOOL -manifest "$lt_tool_outputfile.manifest" -outputresource:"$lt_tool_outputfile" || exit 1;
            $RM "$lt_outputfile.manifest";
          fi'
	;;
      *)
	# Assume MSVC and ICC and IFORT wrapper
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)=' '
	_LT_TAGVAR(allow_undefined_flag, $1)=unsupported
	# Tell ltmain to make .lib files, not .a files.
	libext=lib
	# Tell ltmain to make .dll files, not .so files.
	shrext_cmds=.dll
	# FIXME: Setting linknames here is a bad hack.
	_LT_TAGVAR(archive_cmds, $1)='$CC -o $lib $libobjs $compiler_flags `func_echo_all "$deplibs" | $SED '\''s/ -lc$//'\''` -link -dll~linknames='
	# The linker will automatically build a .lib file if we build a DLL.
	_LT_TAGVAR(old_archive_from_new_cmds, $1)='true'
	# FIXME: Should let the user specify the lib program.
	_LT_TAGVAR(old_archive_cmds, $1)='lib -OUT:$oldlib$oldobjs$old_deplibs'
	_LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
	;;
      esac
      ;;

    darwin* | rhapsody*)
      _LT_DARWIN_LINKER_FEATURES($1)
      ;;

    dgux*)
      _LT_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    # FreeBSD 2.2.[012] allows us to include c++rt0.o to get C++ constructor
    # support.  Future versions do this automatically, but an explicit c++rt0.o
    # does not break anything, and helps significantly (at the cost of a little
    # extra space).
    freebsd2.2*)
      _LT_TAGVAR(archive_cmds, $1)='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags /usr/lib/c++rt0.o'
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    # Unfortunately, older versions of FreeBSD 2 do not have this feature.
    freebsd2.*)
      _LT_TAGVAR(archive_cmds, $1)='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    # FreeBSD 3 and greater uses gcc -shared to do shared libraries.
    freebsd* | dragonfly*)
      _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    hpux9*)
      if test yes = "$GCC"; then
	_LT_TAGVAR(archive_cmds, $1)='$RM $output_objdir/$soname~$CC -shared $pic_flag $wl+b $wl$install_libdir -o $output_objdir/$soname $libobjs $deplibs $compiler_flags~test "x$output_objdir/$soname" = "x$lib" || mv $output_objdir/$soname $lib'
      else
	_LT_TAGVAR(archive_cmds, $1)='$RM $output_objdir/$soname~$LD -b +b $install_libdir -o $output_objdir/$soname $libobjs $deplibs $linker_flags~test "x$output_objdir/$soname" = "x$lib" || mv $output_objdir/$soname $lib'
      fi
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl+b $wl$libdir'
      _LT_TAGVAR(hardcode_libdir_separator, $1)=:
      _LT_TAGVAR(hardcode_direct, $1)=yes

      # hardcode_minus_L: Not really in the search PATH,
      # but as the default location of the library.
      _LT_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
      ;;

    hpux10*)
      if test yes,no = "$GCC,$with_gnu_ld"; then
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $wl+h $wl$soname $wl+b $wl$install_libdir -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_TAGVAR(archive_cmds, $1)='$LD -b +h $soname +b $install_libdir -o $lib $libobjs $deplibs $linker_flags'
      fi
      if test no = "$with_gnu_ld"; then
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl+b $wl$libdir'
	_LT_TAGVAR(hardcode_libdir_separator, $1)=:
	_LT_TAGVAR(hardcode_direct, $1)=yes
	_LT_TAGVAR(hardcode_direct_absolute, $1)=yes
	_LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
	# hardcode_minus_L: Not really in the search PATH,
	# but as the default location of the library.
	_LT_TAGVAR(hardcode_minus_L, $1)=yes
      fi
      ;;

    hpux11*)
      if test yes,no = "$GCC,$with_gnu_ld"; then
	case $host_cpu in
	hppa*64*)
	  _LT_TAGVAR(archive_cmds, $1)='$CC -shared $wl+h $wl$soname -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	ia64*)
	  _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $wl+h $wl$soname $wl+nodefaultrpath -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	*)
	  _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $wl+h $wl$soname $wl+b $wl$install_libdir -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	esac
      else
	case $host_cpu in
	hppa*64*)
	  _LT_TAGVAR(archive_cmds, $1)='$CC -b $wl+h $wl$soname -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	ia64*)
	  _LT_TAGVAR(archive_cmds, $1)='$CC -b $wl+h $wl$soname $wl+nodefaultrpath -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	*)
	m4_if($1, [], [
	  # Older versions of the 11.00 compiler do not understand -b yet
	  # (HP92453-01 A.11.01.20 doesn't, HP92453-01 B.11.X.35175-35176.GP does)
	  _LT_LINKER_OPTION([if $CC understands -b],
	    _LT_TAGVAR(lt_cv_prog_compiler__b, $1), [-b],
	    [_LT_TAGVAR(archive_cmds, $1)='$CC -b $wl+h $wl$soname $wl+b $wl$install_libdir -o $lib $libobjs $deplibs $compiler_flags'],
	    [_LT_TAGVAR(archive_cmds, $1)='$LD -b +h $soname +b $install_libdir -o $lib $libobjs $deplibs $linker_flags'])],
	  [_LT_TAGVAR(archive_cmds, $1)='$CC -b $wl+h $wl$soname $wl+b $wl$install_libdir -o $lib $libobjs $deplibs $compiler_flags'])
	  ;;
	esac
      fi
      if test no = "$with_gnu_ld"; then
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl+b $wl$libdir'
	_LT_TAGVAR(hardcode_libdir_separator, $1)=:

	case $host_cpu in
	hppa*64*|ia64*)
	  _LT_TAGVAR(hardcode_direct, $1)=no
	  _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	  ;;
	*)
	  _LT_TAGVAR(hardcode_direct, $1)=yes
	  _LT_TAGVAR(hardcode_direct_absolute, $1)=yes
	  _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'

	  # hardcode_minus_L: Not really in the search PATH,
	  # but as the default location of the library.
	  _LT_TAGVAR(hardcode_minus_L, $1)=yes
	  ;;
	esac
      fi
      ;;

    irix5* | irix6* | nonstopux*)
      if test yes = "$GCC"; then
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` $wl-update_registry $wl$output_objdir/so_locations -o $lib'
	# Try to use the -exported_symbol ld option, if it does not
	# work, assume that -exports_file does not work either and
	# implicitly export all symbols.
	# This should be the same for all languages, so no per-tag cache variable.
	AC_CACHE_CHECK([whether the $host_os linker accepts -exported_symbol],
	  [lt_cv_irix_exported_symbol],
	  [save_LDFLAGS=$LDFLAGS
	   LDFLAGS="$LDFLAGS -shared $wl-exported_symbol ${wl}foo $wl-update_registry $wl/dev/null"
	   AC_LINK_IFELSE(
	     [AC_LANG_SOURCE(
	        [AC_LANG_CASE([C], [[int foo (void) { return 0; }]],
			      [C++], [[int foo (void) { return 0; }]],
			      [Fortran 77], [[
      subroutine foo
      end]],
			      [Fortran], [[
      subroutine foo
      end]])])],
	      [lt_cv_irix_exported_symbol=yes],
	      [lt_cv_irix_exported_symbol=no])
           LDFLAGS=$save_LDFLAGS])
	if test yes = "$lt_cv_irix_exported_symbol"; then
          _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` $wl-update_registry $wl$output_objdir/so_locations $wl-exports_file $wl$export_symbols -o $lib'
	fi
      else
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags -soname $soname `test -n "$verstring" && func_echo_all "-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags -soname $soname `test -n "$verstring" && func_echo_all "-set_version $verstring"` -update_registry $output_objdir/so_locations -exports_file $export_symbols -o $lib'
      fi
      _LT_TAGVAR(archive_cmds_need_lc, $1)='no'
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
      _LT_TAGVAR(hardcode_libdir_separator, $1)=:
      _LT_TAGVAR(inherit_rpath, $1)=yes
      _LT_TAGVAR(link_all_deplibs, $1)=yes
      ;;

    linux*)
      case $cc_basename in
      tcc*)
	# Fabrice Bellard et al's Tiny C Compiler
	_LT_TAGVAR(ld_shlibs, $1)=yes
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
	;;
      esac
      ;;

    netbsd*)
      if echo __ELF__ | $CC -E - | $GREP __ELF__ >/dev/null; then
	_LT_TAGVAR(archive_cmds, $1)='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'  # a.out
      else
	_LT_TAGVAR(archive_cmds, $1)='$LD -shared -o $lib $libobjs $deplibs $linker_flags'      # ELF
      fi
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    newsos6)
      _LT_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
      _LT_TAGVAR(hardcode_libdir_separator, $1)=:
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    *nto* | *qnx*)
      ;;

    openbsd* | bitrig*)
      if test -f /usr/libexec/ld.so; then
	_LT_TAGVAR(hardcode_direct, $1)=yes
	_LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	_LT_TAGVAR(hardcode_direct_absolute, $1)=yes
	if test -z "`echo __ELF__ | $CC -E - | $GREP __ELF__`"; then
	  _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
	  _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags $wl-retain-symbols-file,$export_symbols'
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
	  _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
	else
	  _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
	fi
      else
	_LT_TAGVAR(ld_shlibs, $1)=no
      fi
      ;;

    os2*)
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_TAGVAR(allow_undefined_flag, $1)=unsupported
      shrext_cmds=.dll
      _LT_TAGVAR(archive_cmds, $1)='$ECHO "LIBRARY ${soname%$shared_ext} INITINSTANCE TERMINSTANCE" > $output_objdir/$libname.def~
	$ECHO "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~
	$ECHO "DATA MULTIPLE NONSHARED" >> $output_objdir/$libname.def~
	$ECHO EXPORTS >> $output_objdir/$libname.def~
	emxexp $libobjs | $SED /"_DLL_InitTerm"/d >> $output_objdir/$libname.def~
	$CC -Zdll -Zcrtdll -o $output_objdir/$soname $libobjs $deplibs $compiler_flags $output_objdir/$libname.def~
	emximp -o $lib $output_objdir/$libname.def'
      _LT_TAGVAR(archive_expsym_cmds, $1)='$ECHO "LIBRARY ${soname%$shared_ext} INITINSTANCE TERMINSTANCE" > $output_objdir/$libname.def~
	$ECHO "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~
	$ECHO "DATA MULTIPLE NONSHARED" >> $output_objdir/$libname.def~
	$ECHO EXPORTS >> $output_objdir/$libname.def~
	prefix_cmds="$SED"~
	if test EXPORTS = "`$SED 1q $export_symbols`"; then
	  prefix_cmds="$prefix_cmds -e 1d";
	fi~
	prefix_cmds="$prefix_cmds -e \"s/^\(.*\)$/_\1/g\""~
	cat $export_symbols | $prefix_cmds >> $output_objdir/$libname.def~
	$CC -Zdll -Zcrtdll -o $output_objdir/$soname $libobjs $deplibs $compiler_flags $output_objdir/$libname.def~
	emximp -o $lib $output_objdir/$libname.def'
      _LT_TAGVAR(old_archive_From_new_cmds, $1)='emximp -o $output_objdir/${libname}_dll.a $output_objdir/$libname.def'
      _LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
      ;;

    osf3*)
      if test yes = "$GCC"; then
	_LT_TAGVAR(allow_undefined_flag, $1)=' $wl-expect_unresolved $wl\*'
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared$allow_undefined_flag $libobjs $deplibs $compiler_flags $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` $wl-update_registry $wl$output_objdir/so_locations -o $lib'
      else
	_LT_TAGVAR(allow_undefined_flag, $1)=' -expect_unresolved \*'
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared$allow_undefined_flag $libobjs $deplibs $compiler_flags -soname $soname `test -n "$verstring" && func_echo_all "-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib'
      fi
      _LT_TAGVAR(archive_cmds_need_lc, $1)='no'
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
      _LT_TAGVAR(hardcode_libdir_separator, $1)=:
      ;;

    osf4* | osf5*)	# as osf3* with the addition of -msym flag
      if test yes = "$GCC"; then
	_LT_TAGVAR(allow_undefined_flag, $1)=' $wl-expect_unresolved $wl\*'
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared$allow_undefined_flag $pic_flag $libobjs $deplibs $compiler_flags $wl-msym $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` $wl-update_registry $wl$output_objdir/so_locations -o $lib'
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
      else
	_LT_TAGVAR(allow_undefined_flag, $1)=' -expect_unresolved \*'
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared$allow_undefined_flag $libobjs $deplibs $compiler_flags -msym -soname $soname `test -n "$verstring" && func_echo_all "-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib'
	_LT_TAGVAR(archive_expsym_cmds, $1)='for i in `cat $export_symbols`; do printf "%s %s\\n" -exported_symbol "\$i" >> $lib.exp; done; printf "%s\\n" "-hidden">> $lib.exp~
          $CC -shared$allow_undefined_flag $wl-input $wl$lib.exp $compiler_flags $libobjs $deplibs -soname $soname `test -n "$verstring" && $ECHO "-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib~$RM $lib.exp'

	# Both c and cxx compiler support -rpath directly
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-rpath $libdir'
      fi
      _LT_TAGVAR(archive_cmds_need_lc, $1)='no'
      _LT_TAGVAR(hardcode_libdir_separator, $1)=:
      ;;

    solaris*)
      _LT_TAGVAR(no_undefined_flag, $1)=' -z defs'
      if test yes = "$GCC"; then
	wlarc='$wl'
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $wl-z ${wl}text $wl-h $wl$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~echo "local: *; };" >> $lib.exp~
          $CC -shared $pic_flag $wl-z ${wl}text $wl-M $wl$lib.exp $wl-h $wl$soname -o $lib $libobjs $deplibs $compiler_flags~$RM $lib.exp'
      else
	case `$CC -V 2>&1` in
	*"Compilers 5.0"*)
	  wlarc=''
	  _LT_TAGVAR(archive_cmds, $1)='$LD -G$allow_undefined_flag -h $soname -o $lib $libobjs $deplibs $linker_flags'
	  _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~echo "local: *; };" >> $lib.exp~
            $LD -G$allow_undefined_flag -M $lib.exp -h $soname -o $lib $libobjs $deplibs $linker_flags~$RM $lib.exp'
	  ;;
	*)
	  wlarc='$wl'
	  _LT_TAGVAR(archive_cmds, $1)='$CC -G$allow_undefined_flag -h $soname -o $lib $libobjs $deplibs $compiler_flags'
	  _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~echo "local: *; };" >> $lib.exp~
            $CC -G$allow_undefined_flag -M $lib.exp -h $soname -o $lib $libobjs $deplibs $compiler_flags~$RM $lib.exp'
	  ;;
	esac
      fi
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      case $host_os in
      solaris2.[[0-5]] | solaris2.[[0-5]].*) ;;
      *)
	# The compiler driver will combine and reorder linker options,
	# but understands '-z linker_flag'.  GCC discards it without '$wl',
	# but is careful enough not to reorder.
	# Supported since Solaris 2.6 (maybe 2.5.1?)
	if test yes = "$GCC"; then
	  _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl-z ${wl}allextract$convenience $wl-z ${wl}defaultextract'
	else
	  _LT_TAGVAR(whole_archive_flag_spec, $1)='-z allextract$convenience -z defaultextract'
	fi
	;;
      esac
      _LT_TAGVAR(link_all_deplibs, $1)=yes
      ;;

    sunos4*)
      if test sequent = "$host_vendor"; then
	# Use $CC to link under sequent, because it throws in some extra .o
	# files that make .init and .fini sections work.
	_LT_TAGVAR(archive_cmds, $1)='$CC -G $wl-h $soname -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_TAGVAR(archive_cmds, $1)='$LD -assert pure-text -Bstatic -o $lib $libobjs $deplibs $linker_flags'
      fi
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_TAGVAR(hardcode_direct, $1)=yes
      _LT_TAGVAR(hardcode_minus_L, $1)=yes
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    sysv4)
      case $host_vendor in
	sni)
	  _LT_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
	  _LT_TAGVAR(hardcode_direct, $1)=yes # is this really true???
	;;
	siemens)
	  ## LD is ld it makes a PLAMLIB
	  ## CC just makes a GrossModule.
	  _LT_TAGVAR(archive_cmds, $1)='$LD -G -o $lib $libobjs $deplibs $linker_flags'
	  _LT_TAGVAR(reload_cmds, $1)='$CC -r -o $output$reload_objs'
	  _LT_TAGVAR(hardcode_direct, $1)=no
        ;;
	motorola)
	  _LT_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
	  _LT_TAGVAR(hardcode_direct, $1)=no #Motorola manual says yes, but my tests say they lie
	;;
      esac
      runpath_var='LD_RUN_PATH'
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    sysv4.3*)
      _LT_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      _LT_TAGVAR(export_dynamic_flag_spec, $1)='-Bexport'
      ;;

    sysv4*MP*)
      if test -d /usr/nec; then
	_LT_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
	_LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	runpath_var=LD_RUN_PATH
	hardcode_runpath_var=yes
	_LT_TAGVAR(ld_shlibs, $1)=yes
      fi
      ;;

    sysv4*uw2* | sysv5OpenUNIX* | sysv5UnixWare7.[[01]].[[10]]* | unixware7* | sco3.2v5.0.[[024]]*)
      _LT_TAGVAR(no_undefined_flag, $1)='$wl-z,text'
      _LT_TAGVAR(archive_cmds_need_lc, $1)=no
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      runpath_var='LD_RUN_PATH'

      if test yes = "$GCC"; then
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_TAGVAR(archive_cmds, $1)='$CC -G $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -G $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
      fi
      ;;

    sysv5* | sco3.2v5* | sco5v6*)
      # Note: We CANNOT use -z defs as we might desire, because we do not
      # link with -lc, and that would cause any symbols used from libc to
      # always be unresolved, which means just about no library would
      # ever link correctly.  If we're not using GNU ld we use -z text
      # though, which does catch some bad symbols but isn't as heavy-handed
      # as -z defs.
      _LT_TAGVAR(no_undefined_flag, $1)='$wl-z,text'
      _LT_TAGVAR(allow_undefined_flag, $1)='$wl-z,nodefs'
      _LT_TAGVAR(archive_cmds_need_lc, $1)=no
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-R,$libdir'
      _LT_TAGVAR(hardcode_libdir_separator, $1)=':'
      _LT_TAGVAR(link_all_deplibs, $1)=yes
      _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-Bexport'
      runpath_var='LD_RUN_PATH'

      if test yes = "$GCC"; then
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
      else
	_LT_TAGVAR(archive_cmds, $1)='$CC -G $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -G $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
      fi
      ;;

    uts4*)
      _LT_TAGVAR(archive_cmds, $1)='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      ;;

    *)
      _LT_TAGVAR(ld_shlibs, $1)=no
      ;;
    esac

    if test sni = "$host_vendor"; then
      case $host in
      sysv4 | sysv4.2uw2* | sysv4.3* | sysv5*)
	_LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-Blargedynsym'
	;;
      esac
    fi
  fi
])
AC_MSG_RESULT([$_LT_TAGVAR(ld_shlibs, $1)])
test no = "$_LT_TAGVAR(ld_shlibs, $1)" && can_build_shared=no

_LT_TAGVAR(with_gnu_ld, $1)=$with_gnu_ld

_LT_DECL([], [libext], [0], [Old archive suffix (normally "a")])dnl
_LT_DECL([], [shrext_cmds], [1], [Shared library suffix (normally ".so")])dnl
_LT_DECL([], [extract_expsyms_cmds], [2],
    [The commands to extract the exported symbol list from a shared archive])

#
# Do we need to explicitly link libc?
#
case "x$_LT_TAGVAR(archive_cmds_need_lc, $1)" in
x|xyes)
  # Assume -lc should be added
  _LT_TAGVAR(archive_cmds_need_lc, $1)=yes

  if test yes,yes = "$GCC,$enable_shared"; then
    case $_LT_TAGVAR(archive_cmds, $1) in
    *'~'*)
      # FIXME: we may have to deal with multi-command sequences.
      ;;
    '$CC '*)
      # Test whether the compiler implicitly links with -lc since on some
      # systems, -lgcc has to come before -lc. If gcc already passes -lc
      # to ld, don't add -lc before -lgcc.
      AC_CACHE_CHECK([whether -lc should be explicitly linked in],
	[lt_cv_]_LT_TAGVAR(archive_cmds_need_lc, $1),
	[$RM conftest*
	echo "$lt_simple_compile_test_code" > conftest.$ac_ext

	if AC_TRY_EVAL(ac_compile) 2>conftest.err; then
	  soname=conftest
	  lib=conftest
	  libobjs=conftest.$ac_objext
	  deplibs=
	  wl=$_LT_TAGVAR(lt_prog_compiler_wl, $1)
	  pic_flag=$_LT_TAGVAR(lt_prog_compiler_pic, $1)
	  compiler_flags=-v
	  linker_flags=-v
	  verstring=
	  output_objdir=.
	  libname=conftest
	  lt_save_allow_undefined_flag=$_LT_TAGVAR(allow_undefined_flag, $1)
	  _LT_TAGVAR(allow_undefined_flag, $1)=
	  if AC_TRY_EVAL(_LT_TAGVAR(archive_cmds, $1) 2\>\&1 \| $GREP \" -lc \" \>/dev/null 2\>\&1)
	  then
	    lt_cv_[]_LT_TAGVAR(archive_cmds_need_lc, $1)=no
	  else
	    lt_cv_[]_LT_TAGVAR(archive_cmds_need_lc, $1)=yes
	  fi
	  _LT_TAGVAR(allow_undefined_flag, $1)=$lt_save_allow_undefined_flag
	else
	  cat conftest.err 1>&5
	fi
	$RM conftest*
	])
      _LT_TAGVAR(archive_cmds_need_lc, $1)=$lt_cv_[]_LT_TAGVAR(archive_cmds_need_lc, $1)
      ;;
    esac
  fi
  ;;
esac

_LT_TAGDECL([build_libtool_need_lc], [archive_cmds_need_lc], [0],
    [Whether or not to add -lc for building shared libraries])
_LT_TAGDECL([allow_libtool_libs_with_static_runtimes],
    [enable_shared_with_static_runtimes], [0],
    [Whether or not to disallow shared libs when runtime libs are static])
_LT_TAGDECL([], [export_dynamic_flag_spec], [1],
    [Compiler flag to allow reflexive dlopens])
_LT_TAGDECL([], [whole_archive_flag_spec], [1],
    [Compiler flag to generate shared objects directly from archives])
_LT_TAGDECL([], [compiler_needs_object], [1],
    [Whether the compiler copes with passing no objects directly])
_LT_TAGDECL([], [old_archive_from_new_cmds], [2],
    [Create an old-style archive from a shared archive])
_LT_TAGDECL([], [old_archive_from_expsyms_cmds], [2],
    [Create a temporary old-style archive to link instead of a shared archive])
_LT_TAGDECL([], [archive_cmds], [2], [Commands used to build a shared archive])
_LT_TAGDECL([], [archive_expsym_cmds], [2])
_LT_TAGDECL([], [module_cmds], [2],
    [Commands used to build a loadable module if different from building
    a shared archive.])
_LT_TAGDECL([], [module_expsym_cmds], [2])
_LT_TAGDECL([], [with_gnu_ld], [1],
    [Whether we are building with GNU ld or not])
_LT_TAGDECL([], [allow_undefined_flag], [1],
    [Flag that allows shared libraries with undefined symbols to be built])
_LT_TAGDECL([], [no_undefined_flag], [1],
    [Flag that enforces no undefined symbols])
_LT_TAGDECL([], [hardcode_libdir_flag_spec], [1],
    [Flag to hardcode $libdir into a binary during linking.
    This must work even if $libdir does not exist])
_LT_TAGDECL([], [hardcode_libdir_separator], [1],
    [Whether we need a single "-rpath" flag with a separated argument])
_LT_TAGDECL([], [hardcode_direct], [0],
    [Set to "yes" if using DIR/libNAME$shared_ext during linking hardcodes
    DIR into the resulting binary])
_LT_TAGDECL([], [hardcode_direct_absolute], [0],
    [Set to "yes" if using DIR/libNAME$shared_ext during linking hardcodes
    DIR into the resulting binary and the resulting library dependency is
    "absolute", i.e impossible to change by setting $shlibpath_var if the
    library is relocated])
_LT_TAGDECL([], [hardcode_minus_L], [0],
    [Set to "yes" if using the -LDIR flag during linking hardcodes DIR
    into the resulting binary])
_LT_TAGDECL([], [hardcode_shlibpath_var], [0],
    [Set to "yes" if using SHLIBPATH_VAR=DIR during linking hardcodes DIR
    into the resulting binary])
_LT_TAGDECL([], [hardcode_automatic], [0],
    [Set to "yes" if building a shared library automatically hardcodes DIR
    into the library and all subsequent libraries and executables linked
    against it])
_LT_TAGDECL([], [inherit_rpath], [0],
    [Set to yes if linker adds runtime paths of dependent libraries
    to runtime path list])
_LT_TAGDECL([], [link_all_deplibs], [0],
    [Whether libtool must link a program against all its dependency libraries])
_LT_TAGDECL([], [always_export_symbols], [0],
    [Set to "yes" if exported symbols are required])
_LT_TAGDECL([], [export_symbols_cmds], [2],
    [The commands to list exported symbols])
_LT_TAGDECL([], [exclude_expsyms], [1],
    [Symbols that should not be listed in the preloaded symbols])
_LT_TAGDECL([], [include_expsyms], [1],
    [Symbols that must always be exported])
_LT_TAGDECL([], [prelink_cmds], [2],
    [Commands necessary for linking programs (against libraries) with templates])
_LT_TAGDECL([], [postlink_cmds], [2],
    [Commands necessary for finishing linking programs])
_LT_TAGDECL([], [file_list_spec], [1],
    [Specify filename containing input files])
dnl FIXME: Not yet implemented
dnl _LT_TAGDECL([], [thread_safe_flag_spec], [1],
dnl    [Compiler flag to generate thread safe objects])
])# _LT_LINKER_SHLIBS


# _LT_LANG_C_CONFIG([TAG])
# ------------------------
# Ensure that the configuration variables for a C compiler are suitably
# defined.  These variables are subsequently used by _LT_CONFIG to write
# the compiler configuration to 'libtool'.
m4_defun([_LT_LANG_C_CONFIG],
[m4_require([_LT_DECL_EGREP])dnl
lt_save_CC=$CC
AC_LANG_PUSH(C)

# Source file extension for C test sources.
ac_ext=c

# Object file extension for compiled C test sources.
objext=o
_LT_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code="int some_variable = 0;"

# Code to be used in simple link tests
lt_simple_link_test_code='int main(){return(0);}'

_LT_TAG_COMPILER
# Save the default compiler, since it gets overwritten when the other
# tags are being tested, and _LT_TAGVAR(compiler, []) is a NOP.
compiler_DEFAULT=$CC

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

if test -n "$compiler"; then
  _LT_COMPILER_NO_RTTI($1)
  _LT_COMPILER_PIC($1)
  _LT_COMPILER_C_O($1)
  _LT_COMPILER_FILE_LOCKS($1)
  _LT_LINKER_SHLIBS($1)
  _LT_SYS_DYNAMIC_LINKER($1)
  _LT_LINKER_HARDCODE_LIBPATH($1)
  LT_SYS_DLOPEN_SELF
  _LT_CMD_STRIPLIB

  # Report what library types will actually be built
  AC_MSG_CHECKING([if libtool supports shared libraries])
  AC_MSG_RESULT([$can_build_shared])

  AC_MSG_CHECKING([whether to build shared libraries])
  test no = "$can_build_shared" && enable_shared=no

  # On AIX, shared libraries and static libraries use the same namespace, and
  # are all built from PIC.
  case $host_os in
  aix3*)
    test yes = "$enable_shared" && enable_static=no
    if test -n "$RANLIB"; then
      archive_cmds="$archive_cmds~\$RANLIB \$lib"
      postinstall_cmds='$RANLIB $lib'
    fi
    ;;

  aix[[4-9]]*)
    if test ia64 != "$host_cpu"; then
      case $enable_shared,$with_aix_soname,$aix_use_runtimelinking in
      yes,aix,yes) ;;			# shared object as lib.so file only
      yes,svr4,*) ;;			# shared object as lib.so archive member only
      yes,*) enable_static=no ;;	# shared object in lib.a archive as well
      esac
    fi
    ;;
  esac
  AC_MSG_RESULT([$enable_shared])

  AC_MSG_CHECKING([whether to build static libraries])
  # Make sure either enable_shared or enable_static is yes.
  test yes = "$enable_shared" || enable_static=yes
  AC_MSG_RESULT([$enable_static])

  _LT_CONFIG($1)
fi
AC_LANG_POP
CC=$lt_save_CC
])# _LT_LANG_C_CONFIG


# _LT_LANG_CXX_CONFIG([TAG])
# --------------------------
# Ensure that the configuration variables for a C++ compiler are suitably
# defined.  These variables are subsequently used by _LT_CONFIG to write
# the compiler configuration to 'libtool'.
m4_defun([_LT_LANG_CXX_CONFIG],
[m4_require([_LT_FILEUTILS_DEFAULTS])dnl
m4_require([_LT_DECL_EGREP])dnl
m4_require([_LT_PATH_MANIFEST_TOOL])dnl
if test -n "$CXX" && ( test no != "$CXX" &&
    ( (test g++ = "$CXX" && `g++ -v >/dev/null 2>&1` ) ||
    (test g++ != "$CXX"))); then
  AC_PROG_CXXCPP
else
  _lt_caught_CXX_error=yes
fi

AC_LANG_PUSH(C++)
_LT_TAGVAR(archive_cmds_need_lc, $1)=no
_LT_TAGVAR(allow_undefined_flag, $1)=
_LT_TAGVAR(always_export_symbols, $1)=no
_LT_TAGVAR(archive_expsym_cmds, $1)=
_LT_TAGVAR(compiler_needs_object, $1)=no
_LT_TAGVAR(export_dynamic_flag_spec, $1)=
_LT_TAGVAR(hardcode_direct, $1)=no
_LT_TAGVAR(hardcode_direct_absolute, $1)=no
_LT_TAGVAR(hardcode_libdir_flag_spec, $1)=
_LT_TAGVAR(hardcode_libdir_separator, $1)=
_LT_TAGVAR(hardcode_minus_L, $1)=no
_LT_TAGVAR(hardcode_shlibpath_var, $1)=unsupported
_LT_TAGVAR(hardcode_automatic, $1)=no
_LT_TAGVAR(inherit_rpath, $1)=no
_LT_TAGVAR(module_cmds, $1)=
_LT_TAGVAR(module_expsym_cmds, $1)=
_LT_TAGVAR(link_all_deplibs, $1)=unknown
_LT_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
_LT_TAGVAR(reload_flag, $1)=$reload_flag
_LT_TAGVAR(reload_cmds, $1)=$reload_cmds
_LT_TAGVAR(no_undefined_flag, $1)=
_LT_TAGVAR(whole_archive_flag_spec, $1)=
_LT_TAGVAR(enable_shared_with_static_runtimes, $1)=no

# Source file extension for C++ test sources.
ac_ext=cpp

# Object file extension for compiled C++ test sources.
objext=o
_LT_TAGVAR(objext, $1)=$objext

# No sense in running all these tests if we already determined that
# the CXX compiler isn't working.  Some variables (like enable_shared)
# are currently assumed to apply to all compilers on this platform,
# and will be corrupted by setting them based on a non-working compiler.
if test yes != "$_lt_caught_CXX_error"; then
  # Code to be used in simple compile tests
  lt_simple_compile_test_code="int some_variable = 0;"

  # Code to be used in simple link tests
  lt_simple_link_test_code='int main(int, char *[[]]) { return(0); }'

  # ltmain only uses $CC for tagged configurations so make sure $CC is set.
  _LT_TAG_COMPILER

  # save warnings/boilerplate of simple test code
  _LT_COMPILER_BOILERPLATE
  _LT_LINKER_BOILERPLATE

  # Allow CC to be a program name with arguments.
  lt_save_CC=$CC
  lt_save_CFLAGS=$CFLAGS
  lt_save_LD=$LD
  lt_save_GCC=$GCC
  GCC=$GXX
  lt_save_with_gnu_ld=$with_gnu_ld
  lt_save_path_LD=$lt_cv_path_LD
  if test -n "${lt_cv_prog_gnu_ldcxx+set}"; then
    lt_cv_prog_gnu_ld=$lt_cv_prog_gnu_ldcxx
  else
    $as_unset lt_cv_prog_gnu_ld
  fi
  if test -n "${lt_cv_path_LDCXX+set}"; then
    lt_cv_path_LD=$lt_cv_path_LDCXX
  else
    $as_unset lt_cv_path_LD
  fi
  test -z "${LDCXX+set}" || LD=$LDCXX
  CC=${CXX-"c++"}
  CFLAGS=$CXXFLAGS
  compiler=$CC
  _LT_TAGVAR(compiler, $1)=$CC
  _LT_CC_BASENAME([$compiler])

  if test -n "$compiler"; then
    # We don't want -fno-exception when compiling C++ code, so set the
    # no_builtin_flag separately
    if test yes = "$GXX"; then
      _LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=' -fno-builtin'
    else
      _LT_TAGVAR(lt_prog_compiler_no_builtin_flag, $1)=
    fi

    if test yes = "$GXX"; then
      # Set up default GNU C++ configuration

      LT_PATH_LD

      # Check if GNU C++ uses GNU ld as the underlying linker, since the
      # archiving commands below assume that GNU ld is being used.
      if test yes = "$with_gnu_ld"; then
        _LT_TAGVAR(archive_cmds, $1)='$CC $pic_flag -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname -o $lib'
        _LT_TAGVAR(archive_expsym_cmds, $1)='$CC $pic_flag -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'

        _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
        _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-dynamic'

        # If archive_cmds runs LD, not CC, wlarc should be empty
        # XXX I think wlarc can be eliminated in ltcf-cxx, but I need to
        #     investigate it a little bit more. (MM)
        wlarc='$wl'

        # ancient GNU ld didn't support --whole-archive et. al.
        if eval "`$CC -print-prog-name=ld` --help 2>&1" |
	  $GREP 'no-whole-archive' > /dev/null; then
          _LT_TAGVAR(whole_archive_flag_spec, $1)=$wlarc'--whole-archive$convenience '$wlarc'--no-whole-archive'
        else
          _LT_TAGVAR(whole_archive_flag_spec, $1)=
        fi
      else
        with_gnu_ld=no
        wlarc=

        # A generic and very simple default shared library creation
        # command for GNU C++ for the case where it uses the native
        # linker, instead of GNU ld.  If possible, this setting should
        # overridden to take advantage of the native linker features on
        # the platform it is being used on.
        _LT_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $lib'
      fi

      # Commands to make compiler produce verbose output that lists
      # what "hidden" libraries, object files and flags are used when
      # linking a shared library.
      output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | $GREP -v "^Configured with:" | $GREP "\-L"'

    else
      GXX=no
      with_gnu_ld=no
      wlarc=
    fi

    # PORTME: fill in a description of your system's C++ link characteristics
    AC_MSG_CHECKING([whether the $compiler linker ($LD) supports shared libraries])
    _LT_TAGVAR(ld_shlibs, $1)=yes
    case $host_os in
      aix3*)
        # FIXME: insert proper C++ library support
        _LT_TAGVAR(ld_shlibs, $1)=no
        ;;
      aix[[4-9]]*)
        if test ia64 = "$host_cpu"; then
          # On IA64, the linker does run time linking by default, so we don't
          # have to do anything special.
          aix_use_runtimelinking=no
          exp_sym_flag='-Bexport'
          no_entry_flag=
        else
          aix_use_runtimelinking=no

          # Test if we are trying to use run time linking or normal
          # AIX style linking. If -brtl is somewhere in LDFLAGS, we
          # have runtime linking enabled, and use it for executables.
          # For shared libraries, we enable/disable runtime linking
          # depending on the kind of the shared library created -
          # when "with_aix_soname,aix_use_runtimelinking" is:
          # "aix,no"   lib.a(lib.so.V) shared, rtl:no,  for executables
          # "aix,yes"  lib.so          shared, rtl:yes, for executables
          #            lib.a           static archive
          # "both,no"  lib.so.V(shr.o) shared, rtl:yes
          #            lib.a(lib.so.V) shared, rtl:no,  for executables
          # "both,yes" lib.so.V(shr.o) shared, rtl:yes, for executables
          #            lib.a(lib.so.V) shared, rtl:no
          # "svr4,*"   lib.so.V(shr.o) shared, rtl:yes, for executables
          #            lib.a           static archive
          case $host_os in aix4.[[23]]|aix4.[[23]].*|aix[[5-9]]*)
	    for ld_flag in $LDFLAGS; do
	      case $ld_flag in
	      *-brtl*)
	        aix_use_runtimelinking=yes
	        break
	        ;;
	      esac
	    done
	    if test svr4,no = "$with_aix_soname,$aix_use_runtimelinking"; then
	      # With aix-soname=svr4, we create the lib.so.V shared archives only,
	      # so we don't have lib.a shared libs to link our executables.
	      # We have to force runtime linking in this case.
	      aix_use_runtimelinking=yes
	      LDFLAGS="$LDFLAGS -Wl,-brtl"
	    fi
	    ;;
          esac

          exp_sym_flag='-bexport'
          no_entry_flag='-bnoentry'
        fi

        # When large executables or shared objects are built, AIX ld can
        # have problems creating the table of contents.  If linking a library
        # or program results in "error TOC overflow" add -mminimal-toc to
        # CXXFLAGS/CFLAGS for g++/gcc.  In the cases where that is not
        # enough to fix the problem, add -Wl,-bbigtoc to LDFLAGS.

        _LT_TAGVAR(archive_cmds, $1)=''
        _LT_TAGVAR(hardcode_direct, $1)=yes
        _LT_TAGVAR(hardcode_direct_absolute, $1)=yes
        _LT_TAGVAR(hardcode_libdir_separator, $1)=':'
        _LT_TAGVAR(link_all_deplibs, $1)=yes
        _LT_TAGVAR(file_list_spec, $1)='$wl-f,'
        case $with_aix_soname,$aix_use_runtimelinking in
        aix,*) ;;	# no import file
        svr4,* | *,yes) # use import file
          # The Import File defines what to hardcode.
          _LT_TAGVAR(hardcode_direct, $1)=no
          _LT_TAGVAR(hardcode_direct_absolute, $1)=no
          ;;
        esac

        if test yes = "$GXX"; then
          case $host_os in aix4.[[012]]|aix4.[[012]].*)
          # We only want to do this on AIX 4.2 and lower, the check
          # below for broken collect2 doesn't work under 4.3+
	  collect2name=`$CC -print-prog-name=collect2`
	  if test -f "$collect2name" &&
	     strings "$collect2name" | $GREP resolve_lib_name >/dev/null
	  then
	    # We have reworked collect2
	    :
	  else
	    # We have old collect2
	    _LT_TAGVAR(hardcode_direct, $1)=unsupported
	    # It fails to find uninstalled libraries when the uninstalled
	    # path is not listed in the libpath.  Setting hardcode_minus_L
	    # to unsupported forces relinking
	    _LT_TAGVAR(hardcode_minus_L, $1)=yes
	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
	    _LT_TAGVAR(hardcode_libdir_separator, $1)=
	  fi
          esac
          shared_flag='-shared'
	  if test yes = "$aix_use_runtimelinking"; then
	    shared_flag=$shared_flag' $wl-G'
	  fi
	  # Need to ensure runtime linking is disabled for the traditional
	  # shared library, or the linker may eventually find shared libraries
	  # /with/ Import File - we do not want to mix them.
	  shared_flag_aix='-shared'
	  shared_flag_svr4='-shared $wl-G'
        else
          # not using gcc
          if test ia64 = "$host_cpu"; then
	  # VisualAge C++, Version 5.5 for AIX 5L for IA-64, Beta 3 Release
	  # chokes on -Wl,-G. The following line is correct:
	  shared_flag='-G'
          else
	    if test yes = "$aix_use_runtimelinking"; then
	      shared_flag='$wl-G'
	    else
	      shared_flag='$wl-bM:SRE'
	    fi
	    shared_flag_aix='$wl-bM:SRE'
	    shared_flag_svr4='$wl-G'
          fi
        fi

        _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-bexpall'
        # It seems that -bexpall does not export symbols beginning with
        # underscore (_), so it is better to generate a list of symbols to
	# export.
        _LT_TAGVAR(always_export_symbols, $1)=yes
	if test aix,yes = "$with_aix_soname,$aix_use_runtimelinking"; then
          # Warning - without using the other runtime loading flags (-brtl),
          # -berok will link without error, but may produce a broken library.
          # The "-G" linker flag allows undefined symbols.
          _LT_TAGVAR(no_undefined_flag, $1)='-bernotok'
          # Determine the default libpath from the value encoded in an empty
          # executable.
          _LT_SYS_MODULE_PATH_AIX([$1])
          _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-blibpath:$libdir:'"$aix_libpath"

          _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -o $output_objdir/$soname $libobjs $deplibs $wl'$no_entry_flag' $compiler_flags `if test -n "$allow_undefined_flag"; then func_echo_all "$wl$allow_undefined_flag"; else :; fi` $wl'$exp_sym_flag:\$export_symbols' '$shared_flag
        else
          if test ia64 = "$host_cpu"; then
	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-R $libdir:/usr/lib:/lib'
	    _LT_TAGVAR(allow_undefined_flag, $1)="-z nodefs"
	    _LT_TAGVAR(archive_expsym_cmds, $1)="\$CC $shared_flag"' -o $output_objdir/$soname $libobjs $deplibs '"\$wl$no_entry_flag"' $compiler_flags $wl$allow_undefined_flag '"\$wl$exp_sym_flag:\$export_symbols"
          else
	    # Determine the default libpath from the value encoded in an
	    # empty executable.
	    _LT_SYS_MODULE_PATH_AIX([$1])
	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-blibpath:$libdir:'"$aix_libpath"
	    # Warning - without using the other run time loading flags,
	    # -berok will link without error, but may produce a broken library.
	    _LT_TAGVAR(no_undefined_flag, $1)=' $wl-bernotok'
	    _LT_TAGVAR(allow_undefined_flag, $1)=' $wl-berok'
	    if test yes = "$with_gnu_ld"; then
	      # We only use this code for GNU lds that support --whole-archive.
	      _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive$convenience $wl--no-whole-archive'
	    else
	      # Exported symbols can be pulled into shared objects from archives
	      _LT_TAGVAR(whole_archive_flag_spec, $1)='$convenience'
	    fi
	    _LT_TAGVAR(archive_cmds_need_lc, $1)=yes
	    _LT_TAGVAR(archive_expsym_cmds, $1)='$RM -r $output_objdir/$realname.d~$MKDIR $output_objdir/$realname.d'
	    # -brtl affects multiple linker settings, -berok does not and is overridden later
	    compiler_flags_filtered='`func_echo_all "$compiler_flags " | $SED -e "s%-brtl\\([[, ]]\\)%-berok\\1%g"`'
	    if test svr4 != "$with_aix_soname"; then
	      # This is similar to how AIX traditionally builds its shared
	      # libraries. Need -bnortl late, we may have -brtl in LDFLAGS.
	      _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$CC '$shared_flag_aix' -o $output_objdir/$realname.d/$soname $libobjs $deplibs $wl-bnoentry '$compiler_flags_filtered'$wl-bE:$export_symbols$allow_undefined_flag~$AR $AR_FLAGS $output_objdir/$libname$release.a $output_objdir/$realname.d/$soname'
	    fi
	    if test aix != "$with_aix_soname"; then
	      _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$CC '$shared_flag_svr4' -o $output_objdir/$realname.d/$shared_archive_member_spec.o $libobjs $deplibs $wl-bnoentry '$compiler_flags_filtered'$wl-bE:$export_symbols$allow_undefined_flag~$STRIP -e $output_objdir/$realname.d/$shared_archive_member_spec.o~( func_echo_all "#! $soname($shared_archive_member_spec.o)"; if test shr_64 = "$shared_archive_member_spec"; then func_echo_all "# 64"; else func_echo_all "# 32"; fi; cat $export_symbols ) > $output_objdir/$realname.d/$shared_archive_member_spec.imp~$AR $AR_FLAGS $output_objdir/$soname $output_objdir/$realname.d/$shared_archive_member_spec.o $output_objdir/$realname.d/$shared_archive_member_spec.imp'
	    else
	      # used by -dlpreopen to get the symbols
	      _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$MV  $output_objdir/$realname.d/$soname $output_objdir'
	    fi
	    _LT_TAGVAR(archive_expsym_cmds, $1)="$_LT_TAGVAR(archive_expsym_cmds, $1)"'~$RM -r $output_objdir/$realname.d'
          fi
        fi
        ;;

      beos*)
	if $LD --help 2>&1 | $GREP ': supported targets:.* elf' > /dev/null; then
	  _LT_TAGVAR(allow_undefined_flag, $1)=unsupported
	  # Joseph Beckenbach <jrb3@best.com> says some releases of gcc
	  # support --undefined.  This deserves some investigation.  FIXME
	  _LT_TAGVAR(archive_cmds, $1)='$CC -nostart $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
	else
	  _LT_TAGVAR(ld_shlibs, $1)=no
	fi
	;;

      chorus*)
        case $cc_basename in
          *)
	  # FIXME: insert proper C++ library support
	  _LT_TAGVAR(ld_shlibs, $1)=no
	  ;;
        esac
        ;;

      cygwin* | mingw* | pw32* | cegcc*)
	case $GXX,$cc_basename in
	,cl* | no,cl* | ,icl* | no,icl* | ,ifort* | no,ifort*)
	  # Native MSVC or ICC or IFORT
	  # hardcode_libdir_flag_spec is actually meaningless, as there is
	  # no search path for DLLs.
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)=' '
	  _LT_TAGVAR(allow_undefined_flag, $1)=unsupported
	  _LT_TAGVAR(always_export_symbols, $1)=yes
	  _LT_TAGVAR(file_list_spec, $1)='@'
	  # Tell ltmain to make .lib files, not .a files.
	  libext=lib
	  # Tell ltmain to make .dll files, not .so files.
	  shrext_cmds=.dll
	  # FIXME: Setting linknames here is a bad hack.
	  _LT_TAGVAR(archive_cmds, $1)='$CC -o $output_objdir/$soname $libobjs $compiler_flags $deplibs -Wl,-DLL,-IMPLIB:"$tool_output_objdir$libname.dll.lib"~linknames='
	  _LT_TAGVAR(archive_expsym_cmds, $1)='if _LT_DLL_DEF_P([$export_symbols]); then
              cp "$export_symbols" "$output_objdir/$soname.def";
              echo "$tool_output_objdir$soname.def" > "$output_objdir/$soname.exp";
            else
              $SED -e '\''s/^/-link -EXPORT:/'\'' < $export_symbols > $output_objdir/$soname.exp;
            fi~
            $CC -o $tool_output_objdir$soname $libobjs $compiler_flags $deplibs "@$tool_output_objdir$soname.exp" -Wl,-DLL,-IMPLIB:"$tool_output_objdir$libname.dll.lib"~
            linknames='
	  # The linker will not automatically build a static lib if we build a DLL.
	  # _LT_TAGVAR(old_archive_from_new_cmds, $1)='true'
	  _LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
	  # Don't use ranlib
	  _LT_TAGVAR(old_postinstall_cmds, $1)='chmod 644 $oldlib'
	  _LT_TAGVAR(postlink_cmds, $1)='lt_outputfile="@OUTPUT@"~
            lt_tool_outputfile="@TOOL_OUTPUT@"~
            case $lt_outputfile in
              *.exe|*.EXE) ;;
              *)
                lt_outputfile=$lt_outputfile.exe
                lt_tool_outputfile=$lt_tool_outputfile.exe
                ;;
            esac~
            func_to_tool_file "$lt_outputfile"~
            if test : != "$MANIFEST_TOOL" && test -f "$lt_outputfile.manifest"; then
              $MANIFEST_TOOL -manifest "$lt_tool_outputfile.manifest" -outputresource:"$lt_tool_outputfile" || exit 1;
              $RM "$lt_outputfile.manifest";
            fi'
	  ;;
	*)
	  # g++
	  # _LT_TAGVAR(hardcode_libdir_flag_spec, $1) is actually meaningless,
	  # as there is no search path for DLLs.
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
	  _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-all-symbols'
	  _LT_TAGVAR(allow_undefined_flag, $1)=unsupported
	  _LT_TAGVAR(always_export_symbols, $1)=no
	  _LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes

	  if $LD --help 2>&1 | $GREP 'auto-import' > /dev/null; then
	    _LT_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $output_objdir/$soname $wl--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
	    # If the export-symbols file already is a .def file, use it as
	    # is; otherwise, prepend EXPORTS...
	    _LT_TAGVAR(archive_expsym_cmds, $1)='if _LT_DLL_DEF_P([$export_symbols]); then
              cp $export_symbols $output_objdir/$soname.def;
            else
              echo EXPORTS > $output_objdir/$soname.def;
              cat $export_symbols >> $output_objdir/$soname.def;
            fi~
            $CC -shared -nostdlib $output_objdir/$soname.def $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $output_objdir/$soname $wl--enable-auto-image-base -Xlinker --out-implib -Xlinker $lib'
	  else
	    _LT_TAGVAR(ld_shlibs, $1)=no
	  fi
	  ;;
	esac
	;;
      darwin* | rhapsody*)
        _LT_DARWIN_LINKER_FEATURES($1)
	;;

      os2*)
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-L$libdir'
	_LT_TAGVAR(hardcode_minus_L, $1)=yes
	_LT_TAGVAR(allow_undefined_flag, $1)=unsupported
	shrext_cmds=.dll
	_LT_TAGVAR(archive_cmds, $1)='$ECHO "LIBRARY ${soname%$shared_ext} INITINSTANCE TERMINSTANCE" > $output_objdir/$libname.def~
	  $ECHO "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~
	  $ECHO "DATA MULTIPLE NONSHARED" >> $output_objdir/$libname.def~
	  $ECHO EXPORTS >> $output_objdir/$libname.def~
	  emxexp $libobjs | $SED /"_DLL_InitTerm"/d >> $output_objdir/$libname.def~
	  $CC -Zdll -Zcrtdll -o $output_objdir/$soname $libobjs $deplibs $compiler_flags $output_objdir/$libname.def~
	  emximp -o $lib $output_objdir/$libname.def'
	_LT_TAGVAR(archive_expsym_cmds, $1)='$ECHO "LIBRARY ${soname%$shared_ext} INITINSTANCE TERMINSTANCE" > $output_objdir/$libname.def~
	  $ECHO "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~
	  $ECHO "DATA MULTIPLE NONSHARED" >> $output_objdir/$libname.def~
	  $ECHO EXPORTS >> $output_objdir/$libname.def~
	  prefix_cmds="$SED"~
	  if test EXPORTS = "`$SED 1q $export_symbols`"; then
	    prefix_cmds="$prefix_cmds -e 1d";
	  fi~
	  prefix_cmds="$prefix_cmds -e \"s/^\(.*\)$/_\1/g\""~
	  cat $export_symbols | $prefix_cmds >> $output_objdir/$libname.def~
	  $CC -Zdll -Zcrtdll -o $output_objdir/$soname $libobjs $deplibs $compiler_flags $output_objdir/$libname.def~
	  emximp -o $lib $output_objdir/$libname.def'
	_LT_TAGVAR(old_archive_From_new_cmds, $1)='emximp -o $output_objdir/${libname}_dll.a $output_objdir/$libname.def'
	_LT_TAGVAR(enable_shared_with_static_runtimes, $1)=yes
	;;

      dgux*)
        case $cc_basename in
          ec++*)
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
          ghcx*)
	    # Green Hills C++ Compiler
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
          *)
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
        esac
        ;;

      freebsd2.*)
        # C++ shared libraries reported to be fairly broken before
	# switch to ELF
        _LT_TAGVAR(ld_shlibs, $1)=no
        ;;

      freebsd-elf*)
        _LT_TAGVAR(archive_cmds_need_lc, $1)=no
        ;;

      freebsd* | dragonfly*)
        # FreeBSD 3 and later use GNU C++ and GNU ld with standard ELF
        # conventions
        _LT_TAGVAR(ld_shlibs, $1)=yes
        ;;

      haiku*)
        _LT_TAGVAR(archive_cmds, $1)='$CC -shared $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
        _LT_TAGVAR(link_all_deplibs, $1)=yes
        ;;

      hpux9*)
        _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl+b $wl$libdir'
        _LT_TAGVAR(hardcode_libdir_separator, $1)=:
        _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
        _LT_TAGVAR(hardcode_direct, $1)=yes
        _LT_TAGVAR(hardcode_minus_L, $1)=yes # Not in the search PATH,
				             # but as the default
				             # location of the library.

        case $cc_basename in
          CC*)
            # FIXME: insert proper C++ library support
            _LT_TAGVAR(ld_shlibs, $1)=no
            ;;
          aCC*)
            _LT_TAGVAR(archive_cmds, $1)='$RM $output_objdir/$soname~$CC -b $wl+b $wl$install_libdir -o $output_objdir/$soname $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~test "x$output_objdir/$soname" = "x$lib" || mv $output_objdir/$soname $lib'
            # Commands to make compiler produce verbose output that lists
            # what "hidden" libraries, object files and flags are used when
            # linking a shared library.
            #
            # There doesn't appear to be a way to prevent this compiler from
            # explicitly linking system object files so we need to strip them
            # from the output so that they don't get included in the library
            # dependencies.
            output_verbose_link_cmd='templist=`($CC -b $CFLAGS -v conftest.$objext 2>&1) | $EGREP "\-L"`; list= ; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; func_echo_all "$list"'
            ;;
          *)
            if test yes = "$GXX"; then
              _LT_TAGVAR(archive_cmds, $1)='$RM $output_objdir/$soname~$CC -shared -nostdlib $pic_flag $wl+b $wl$install_libdir -o $output_objdir/$soname $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~test "x$output_objdir/$soname" = "x$lib" || mv $output_objdir/$soname $lib'
            else
              # FIXME: insert proper C++ library support
              _LT_TAGVAR(ld_shlibs, $1)=no
            fi
            ;;
        esac
        ;;

      hpux10*|hpux11*)
        if test no = "$with_gnu_ld"; then
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl+b $wl$libdir'
	  _LT_TAGVAR(hardcode_libdir_separator, $1)=:

          case $host_cpu in
            hppa*64*|ia64*)
              ;;
            *)
	      _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
              ;;
          esac
        fi
        case $host_cpu in
          hppa*64*|ia64*)
            _LT_TAGVAR(hardcode_direct, $1)=no
            _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
            ;;
          *)
            _LT_TAGVAR(hardcode_direct, $1)=yes
            _LT_TAGVAR(hardcode_direct_absolute, $1)=yes
            _LT_TAGVAR(hardcode_minus_L, $1)=yes # Not in the search PATH,
					         # but as the default
					         # location of the library.
            ;;
        esac

        case $cc_basename in
          CC*)
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
          aCC*)
	    case $host_cpu in
	      hppa*64*)
	        _LT_TAGVAR(archive_cmds, $1)='$CC -b $wl+h $wl$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	        ;;
	      ia64*)
	        _LT_TAGVAR(archive_cmds, $1)='$CC -b $wl+h $wl$soname $wl+nodefaultrpath -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	        ;;
	      *)
	        _LT_TAGVAR(archive_cmds, $1)='$CC -b $wl+h $wl$soname $wl+b $wl$install_libdir -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	        ;;
	    esac
	    # Commands to make compiler produce verbose output that lists
	    # what "hidden" libraries, object files and flags are used when
	    # linking a shared library.
	    #
	    # There doesn't appear to be a way to prevent this compiler from
	    # explicitly linking system object files so we need to strip them
	    # from the output so that they don't get included in the library
	    # dependencies.
	    output_verbose_link_cmd='templist=`($CC -b $CFLAGS -v conftest.$objext 2>&1) | $GREP "\-L"`; list= ; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; func_echo_all "$list"'
	    ;;
          *)
	    if test yes = "$GXX"; then
	      if test no = "$with_gnu_ld"; then
	        case $host_cpu in
	          hppa*64*)
	            _LT_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib -fPIC $wl+h $wl$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	            ;;
	          ia64*)
	            _LT_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $pic_flag $wl+h $wl$soname $wl+nodefaultrpath -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	            ;;
	          *)
	            _LT_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $pic_flag $wl+h $wl$soname $wl+b $wl$install_libdir -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	            ;;
	        esac
	      fi
	    else
	      # FIXME: insert proper C++ library support
	      _LT_TAGVAR(ld_shlibs, $1)=no
	    fi
	    ;;
        esac
        ;;

      interix[[3-9]]*)
	_LT_TAGVAR(hardcode_direct, $1)=no
	_LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
	_LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
	# Hack: On Interix 3.x, we cannot compile PIC because of a broken gcc.
	# Instead, shared libraries are loaded at an image base (0x10000000 by
	# default) and relocated if they conflict, which is a slow very memory
	# consuming and fragmenting process.  To avoid this, we pick a random,
	# 256 KiB-aligned image base between 0x50000000 and 0x6FFC0000 at link
	# time.  Moving up from 0x10000000 also allows more sbrk(2) space.
	_LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-h,$soname $wl--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
	_LT_TAGVAR(archive_expsym_cmds, $1)='sed "s|^|_|" $export_symbols >$output_objdir/$soname.expsym~$CC -shared $pic_flag $libobjs $deplibs $compiler_flags $wl-h,$soname $wl--retain-symbols-file,$output_objdir/$soname.expsym $wl--image-base,`expr ${RANDOM-$$} % 4096 / 2 \* 262144 + 1342177280` -o $lib'
	;;
      irix5* | irix6*)
        case $cc_basename in
          CC*)
	    # SGI C++
	    _LT_TAGVAR(archive_cmds, $1)='$CC -shared -all -multigot $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -soname $soname `test -n "$verstring" && func_echo_all "-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib'

	    # Archives containing C++ object files must be created using
	    # "CC -ar", where "CC" is the IRIX C++ compiler.  This is
	    # necessary to make sure instantiated templates are included
	    # in the archive.
	    _LT_TAGVAR(old_archive_cmds, $1)='$CC -ar -WR,-u -o $oldlib $oldobjs'
	    ;;
          *)
	    if test yes = "$GXX"; then
	      if test no = "$with_gnu_ld"; then
	        _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` $wl-update_registry $wl$output_objdir/so_locations -o $lib'
	      else
	        _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` -o $lib'
	      fi
	    fi
	    _LT_TAGVAR(link_all_deplibs, $1)=yes
	    ;;
        esac
        _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
        _LT_TAGVAR(hardcode_libdir_separator, $1)=:
        _LT_TAGVAR(inherit_rpath, $1)=yes
        ;;

      linux* | k*bsd*-gnu | kopensolaris*-gnu | gnu*)
        case $cc_basename in
          KCC*)
	    # Kuck and Associates, Inc. (KAI) C++ Compiler

	    # KCC will only create a shared library if the output file
	    # ends with ".so" (or ".sl" for HP-UX), so rename the library
	    # to its proper name (with version) after linking.
	    _LT_TAGVAR(archive_cmds, $1)='tempext=`echo $shared_ext | $SED -e '\''s/\([[^()0-9A-Za-z{}]]\)/\\\\\1/g'\''`; templib=`echo $lib | $SED -e "s/\$tempext\..*/.so/"`; $CC $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags --soname $soname -o \$templib; mv \$templib $lib'
	    _LT_TAGVAR(archive_expsym_cmds, $1)='tempext=`echo $shared_ext | $SED -e '\''s/\([[^()0-9A-Za-z{}]]\)/\\\\\1/g'\''`; templib=`echo $lib | $SED -e "s/\$tempext\..*/.so/"`; $CC $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags --soname $soname -o \$templib $wl-retain-symbols-file,$export_symbols; mv \$templib $lib'
	    # Commands to make compiler produce verbose output that lists
	    # what "hidden" libraries, object files and flags are used when
	    # linking a shared library.
	    #
	    # There doesn't appear to be a way to prevent this compiler from
	    # explicitly linking system object files so we need to strip them
	    # from the output so that they don't get included in the library
	    # dependencies.
	    output_verbose_link_cmd='templist=`$CC $CFLAGS -v conftest.$objext -o libconftest$shared_ext 2>&1 | $GREP "ld"`; rm -f libconftest$shared_ext; list= ; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; func_echo_all "$list"'

	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
	    _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-dynamic'

	    # Archives containing C++ object files must be created using
	    # "CC -Bstatic", where "CC" is the KAI C++ compiler.
	    _LT_TAGVAR(old_archive_cmds, $1)='$CC -Bstatic -o $oldlib $oldobjs'
	    ;;
	  icpc* | ecpc* )
	    # Intel C++
	    with_gnu_ld=yes
	    # version 8.0 and above of icpc choke on multiply defined symbols
	    # if we add $predep_objects and $postdep_objects, however 7.1 and
	    # earlier do not add the objects themselves.
	    case `$CC -V 2>&1` in
	      *"Version 7."*)
	        _LT_TAGVAR(archive_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname -o $lib'
		_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
		;;
	      *)  # Version 8.0 or newer
	        tmp_idyn=
	        case $host_cpu in
		  ia64*) tmp_idyn=' -i_dynamic';;
		esac
	        _LT_TAGVAR(archive_cmds, $1)='$CC -shared'"$tmp_idyn"' $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
		_LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared'"$tmp_idyn"' $libobjs $deplibs $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
		;;
	    esac
	    _LT_TAGVAR(archive_cmds_need_lc, $1)=no
	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
	    _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-dynamic'
	    _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive$convenience $wl--no-whole-archive'
	    ;;
          pgCC* | pgcpp*)
            # Portland Group C++ compiler
	    case `$CC -V` in
	    *pgCC\ [[1-5]].* | *pgcpp\ [[1-5]].*)
	      _LT_TAGVAR(prelink_cmds, $1)='tpldir=Template.dir~
               rm -rf $tpldir~
               $CC --prelink_objects --instantiation_dir $tpldir $objs $libobjs $compile_deplibs~
               compile_command="$compile_command `find $tpldir -name \*.o | sort | $NL2SP`"'
	      _LT_TAGVAR(old_archive_cmds, $1)='tpldir=Template.dir~
                rm -rf $tpldir~
                $CC --prelink_objects --instantiation_dir $tpldir $oldobjs$old_deplibs~
                $AR $AR_FLAGS $oldlib$oldobjs$old_deplibs `find $tpldir -name \*.o | sort | $NL2SP`~
                $RANLIB $oldlib'
	      _LT_TAGVAR(archive_cmds, $1)='tpldir=Template.dir~
                rm -rf $tpldir~
                $CC --prelink_objects --instantiation_dir $tpldir $predep_objects $libobjs $deplibs $convenience $postdep_objects~
                $CC -shared $pic_flag $predep_objects $libobjs $deplibs `find $tpldir -name \*.o | sort | $NL2SP` $postdep_objects $compiler_flags $wl-soname $wl$soname -o $lib'
	      _LT_TAGVAR(archive_expsym_cmds, $1)='tpldir=Template.dir~
                rm -rf $tpldir~
                $CC --prelink_objects --instantiation_dir $tpldir $predep_objects $libobjs $deplibs $convenience $postdep_objects~
                $CC -shared $pic_flag $predep_objects $libobjs $deplibs `find $tpldir -name \*.o | sort | $NL2SP` $postdep_objects $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
	      ;;
	    *) # Version 6 and above use weak symbols
	      _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname -o $lib'
	      _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname $wl-retain-symbols-file $wl$export_symbols -o $lib'
	      ;;
	    esac

	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl--rpath $wl$libdir'
	    _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-dynamic'
	    _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive`for conv in $convenience\"\"; do test  -n \"$conv\" && new_convenience=\"$new_convenience,$conv\"; done; func_echo_all \"$new_convenience\"` $wl--no-whole-archive'
            ;;
	  cxx*)
	    # Compaq C++
	    _LT_TAGVAR(archive_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname -o $lib'
	    _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname  -o $lib $wl-retain-symbols-file $wl$export_symbols'

	    runpath_var=LD_RUN_PATH
	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-rpath $libdir'
	    _LT_TAGVAR(hardcode_libdir_separator, $1)=:

	    # Commands to make compiler produce verbose output that lists
	    # what "hidden" libraries, object files and flags are used when
	    # linking a shared library.
	    #
	    # There doesn't appear to be a way to prevent this compiler from
	    # explicitly linking system object files so we need to strip them
	    # from the output so that they don't get included in the library
	    # dependencies.
	    output_verbose_link_cmd='templist=`$CC -shared $CFLAGS -v conftest.$objext 2>&1 | $GREP "ld"`; templist=`func_echo_all "$templist" | $SED "s/\(^.*ld.*\)\( .*ld .*$\)/\1/"`; list= ; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; func_echo_all "X$list" | $Xsed'
	    ;;
	  xl* | mpixl* | bgxl*)
	    # IBM XL 8.0 on PPC, with GNU ld
	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
	    _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl--export-dynamic'
	    _LT_TAGVAR(archive_cmds, $1)='$CC -qmkshrobj $libobjs $deplibs $compiler_flags $wl-soname $wl$soname -o $lib'
	    if test yes = "$supports_anon_versioning"; then
	      _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $output_objdir/$libname.ver~
                cat $export_symbols | sed -e "s/\(.*\)/\1;/" >> $output_objdir/$libname.ver~
                echo "local: *; };" >> $output_objdir/$libname.ver~
                $CC -qmkshrobj $libobjs $deplibs $compiler_flags $wl-soname $wl$soname $wl-version-script $wl$output_objdir/$libname.ver -o $lib'
	    fi
	    ;;
	  *)
	    case `$CC -V 2>&1 | sed 5q` in
	    *Sun\ C*)
	      # Sun C++ 5.9
	      _LT_TAGVAR(no_undefined_flag, $1)=' -zdefs'
	      _LT_TAGVAR(archive_cmds, $1)='$CC -G$allow_undefined_flag -h$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	      _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -G$allow_undefined_flag -h$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-retain-symbols-file $wl$export_symbols'
	      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
	      _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl--whole-archive`new_convenience=; for conv in $convenience\"\"; do test -z \"$conv\" || new_convenience=\"$new_convenience,$conv\"; done; func_echo_all \"$new_convenience\"` $wl--no-whole-archive'
	      _LT_TAGVAR(compiler_needs_object, $1)=yes

	      # Not sure whether something based on
	      # $CC $CFLAGS -v conftest.$objext -o libconftest$shared_ext 2>&1
	      # would be better.
	      output_verbose_link_cmd='func_echo_all'

	      # Archives containing C++ object files must be created using
	      # "CC -xar", where "CC" is the Sun C++ compiler.  This is
	      # necessary to make sure instantiated templates are included
	      # in the archive.
	      _LT_TAGVAR(old_archive_cmds, $1)='$CC -xar -o $oldlib $oldobjs'
	      ;;
	    esac
	    ;;
	esac
	;;

      lynxos*)
        # FIXME: insert proper C++ library support
	_LT_TAGVAR(ld_shlibs, $1)=no
	;;

      m88k*)
        # FIXME: insert proper C++ library support
        _LT_TAGVAR(ld_shlibs, $1)=no
	;;

      mvs*)
        case $cc_basename in
          cxx*)
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
	  *)
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
	esac
	;;

      netbsd*)
        if echo __ELF__ | $CC -E - | $GREP __ELF__ >/dev/null; then
	  _LT_TAGVAR(archive_cmds, $1)='$LD -Bshareable  -o $lib $predep_objects $libobjs $deplibs $postdep_objects $linker_flags'
	  wlarc=
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
	  _LT_TAGVAR(hardcode_direct, $1)=yes
	  _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	fi
	# Workaround some broken pre-1.5 toolchains
	output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | $GREP conftest.$objext | $SED -e "s:-lgcc -lc -lgcc::"'
	;;

      *nto* | *qnx*)
        _LT_TAGVAR(ld_shlibs, $1)=yes
	;;

      openbsd* | bitrig*)
	if test -f /usr/libexec/ld.so; then
	  _LT_TAGVAR(hardcode_direct, $1)=yes
	  _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	  _LT_TAGVAR(hardcode_direct_absolute, $1)=yes
	  _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -o $lib'
	  _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
	  if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`"; then
	    _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $pic_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-retain-symbols-file,$export_symbols -o $lib'
	    _LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-E'
	    _LT_TAGVAR(whole_archive_flag_spec, $1)=$wlarc'--whole-archive$convenience '$wlarc'--no-whole-archive'
	  fi
	  output_verbose_link_cmd=func_echo_all
	else
	  _LT_TAGVAR(ld_shlibs, $1)=no
	fi
	;;

      osf3* | osf4* | osf5*)
        case $cc_basename in
          KCC*)
	    # Kuck and Associates, Inc. (KAI) C++ Compiler

	    # KCC will only create a shared library if the output file
	    # ends with ".so" (or ".sl" for HP-UX), so rename the library
	    # to its proper name (with version) after linking.
	    _LT_TAGVAR(archive_cmds, $1)='tempext=`echo $shared_ext | $SED -e '\''s/\([[^()0-9A-Za-z{}]]\)/\\\\\1/g'\''`; templib=`echo "$lib" | $SED -e "s/\$tempext\..*/.so/"`; $CC $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags --soname $soname -o \$templib; mv \$templib $lib'

	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath,$libdir'
	    _LT_TAGVAR(hardcode_libdir_separator, $1)=:

	    # Archives containing C++ object files must be created using
	    # the KAI C++ compiler.
	    case $host in
	      osf3*) _LT_TAGVAR(old_archive_cmds, $1)='$CC -Bstatic -o $oldlib $oldobjs' ;;
	      *) _LT_TAGVAR(old_archive_cmds, $1)='$CC -o $oldlib $oldobjs' ;;
	    esac
	    ;;
          RCC*)
	    # Rational C++ 2.4.1
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
          cxx*)
	    case $host in
	      osf3*)
	        _LT_TAGVAR(allow_undefined_flag, $1)=' $wl-expect_unresolved $wl\*'
	        _LT_TAGVAR(archive_cmds, $1)='$CC -shared$allow_undefined_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $soname `test -n "$verstring" && func_echo_all "$wl-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib'
	        _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
		;;
	      *)
	        _LT_TAGVAR(allow_undefined_flag, $1)=' -expect_unresolved \*'
	        _LT_TAGVAR(archive_cmds, $1)='$CC -shared$allow_undefined_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -msym -soname $soname `test -n "$verstring" && func_echo_all "-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib'
	        _LT_TAGVAR(archive_expsym_cmds, $1)='for i in `cat $export_symbols`; do printf "%s %s\\n" -exported_symbol "\$i" >> $lib.exp; done~
                  echo "-hidden">> $lib.exp~
                  $CC -shared$allow_undefined_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags -msym -soname $soname $wl-input $wl$lib.exp  `test -n "$verstring" && $ECHO "-set_version $verstring"` -update_registry $output_objdir/so_locations -o $lib~
                  $RM $lib.exp'
	        _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-rpath $libdir'
		;;
	    esac

	    _LT_TAGVAR(hardcode_libdir_separator, $1)=:

	    # Commands to make compiler produce verbose output that lists
	    # what "hidden" libraries, object files and flags are used when
	    # linking a shared library.
	    #
	    # There doesn't appear to be a way to prevent this compiler from
	    # explicitly linking system object files so we need to strip them
	    # from the output so that they don't get included in the library
	    # dependencies.
	    output_verbose_link_cmd='templist=`$CC -shared $CFLAGS -v conftest.$objext 2>&1 | $GREP "ld" | $GREP -v "ld:"`; templist=`func_echo_all "$templist" | $SED "s/\(^.*ld.*\)\( .*ld.*$\)/\1/"`; list= ; for z in $templist; do case $z in conftest.$objext) list="$list $z";; *.$objext);; *) list="$list $z";;esac; done; func_echo_all "$list"'
	    ;;
	  *)
	    if test yes,no = "$GXX,$with_gnu_ld"; then
	      _LT_TAGVAR(allow_undefined_flag, $1)=' $wl-expect_unresolved $wl\*'
	      case $host in
	        osf3*)
	          _LT_TAGVAR(archive_cmds, $1)='$CC -shared -nostdlib $allow_undefined_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` $wl-update_registry $wl$output_objdir/so_locations -o $lib'
		  ;;
	        *)
	          _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -nostdlib $allow_undefined_flag $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-msym $wl-soname $wl$soname `test -n "$verstring" && func_echo_all "$wl-set_version $wl$verstring"` $wl-update_registry $wl$output_objdir/so_locations -o $lib'
		  ;;
	      esac

	      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-rpath $wl$libdir'
	      _LT_TAGVAR(hardcode_libdir_separator, $1)=:

	      # Commands to make compiler produce verbose output that lists
	      # what "hidden" libraries, object files and flags are used when
	      # linking a shared library.
	      output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | $GREP -v "^Configured with:" | $GREP "\-L"'

	    else
	      # FIXME: insert proper C++ library support
	      _LT_TAGVAR(ld_shlibs, $1)=no
	    fi
	    ;;
        esac
        ;;

      psos*)
        # FIXME: insert proper C++ library support
        _LT_TAGVAR(ld_shlibs, $1)=no
        ;;

      sunos4*)
        case $cc_basename in
          CC*)
	    # Sun C++ 4.x
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
          lcc*)
	    # Lucid
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
          *)
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
        esac
        ;;

      solaris*)
        case $cc_basename in
          CC* | sunCC*)
	    # Sun C++ 4.2, 5.x and Centerline C++
            _LT_TAGVAR(archive_cmds_need_lc,$1)=yes
	    _LT_TAGVAR(no_undefined_flag, $1)=' -zdefs'
	    _LT_TAGVAR(archive_cmds, $1)='$CC -G$allow_undefined_flag -h$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags'
	    _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~echo "local: *; };" >> $lib.exp~
              $CC -G$allow_undefined_flag $wl-M $wl$lib.exp -h$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~$RM $lib.exp'

	    _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='-R$libdir'
	    _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	    case $host_os in
	      solaris2.[[0-5]] | solaris2.[[0-5]].*) ;;
	      *)
		# The compiler driver will combine and reorder linker options,
		# but understands '-z linker_flag'.
	        # Supported since Solaris 2.6 (maybe 2.5.1?)
		_LT_TAGVAR(whole_archive_flag_spec, $1)='-z allextract$convenience -z defaultextract'
	        ;;
	    esac
	    _LT_TAGVAR(link_all_deplibs, $1)=yes

	    output_verbose_link_cmd='func_echo_all'

	    # Archives containing C++ object files must be created using
	    # "CC -xar", where "CC" is the Sun C++ compiler.  This is
	    # necessary to make sure instantiated templates are included
	    # in the archive.
	    _LT_TAGVAR(old_archive_cmds, $1)='$CC -xar -o $oldlib $oldobjs'
	    ;;
          gcx*)
	    # Green Hills C++ Compiler
	    _LT_TAGVAR(archive_cmds, $1)='$CC -shared $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-h $wl$soname -o $lib'

	    # The C++ compiler must be used to create the archive.
	    _LT_TAGVAR(old_archive_cmds, $1)='$CC $LDFLAGS -archive -o $oldlib $oldobjs'
	    ;;
          *)
	    # GNU C++ compiler with Solaris linker
	    if test yes,no = "$GXX,$with_gnu_ld"; then
	      _LT_TAGVAR(no_undefined_flag, $1)=' $wl-z ${wl}defs'
	      if $CC --version | $GREP -v '^2\.7' > /dev/null; then
	        _LT_TAGVAR(archive_cmds, $1)='$CC -shared $pic_flag -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-h $wl$soname -o $lib'
	        _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~echo "local: *; };" >> $lib.exp~
                  $CC -shared $pic_flag -nostdlib $wl-M $wl$lib.exp $wl-h $wl$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~$RM $lib.exp'

	        # Commands to make compiler produce verbose output that lists
	        # what "hidden" libraries, object files and flags are used when
	        # linking a shared library.
	        output_verbose_link_cmd='$CC -shared $CFLAGS -v conftest.$objext 2>&1 | $GREP -v "^Configured with:" | $GREP "\-L"'
	      else
	        # g++ 2.7 appears to require '-G' NOT '-shared' on this
	        # platform.
	        _LT_TAGVAR(archive_cmds, $1)='$CC -G -nostdlib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags $wl-h $wl$soname -o $lib'
	        _LT_TAGVAR(archive_expsym_cmds, $1)='echo "{ global:" > $lib.exp~cat $export_symbols | $SED -e "s/\(.*\)/\1;/" >> $lib.exp~echo "local: *; };" >> $lib.exp~
                  $CC -G -nostdlib $wl-M $wl$lib.exp $wl-h $wl$soname -o $lib $predep_objects $libobjs $deplibs $postdep_objects $compiler_flags~$RM $lib.exp'

	        # Commands to make compiler produce verbose output that lists
	        # what "hidden" libraries, object files and flags are used when
	        # linking a shared library.
	        output_verbose_link_cmd='$CC -G $CFLAGS -v conftest.$objext 2>&1 | $GREP -v "^Configured with:" | $GREP "\-L"'
	      fi

	      _LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-R $wl$libdir'
	      case $host_os in
		solaris2.[[0-5]] | solaris2.[[0-5]].*) ;;
		*)
		  _LT_TAGVAR(whole_archive_flag_spec, $1)='$wl-z ${wl}allextract$convenience $wl-z ${wl}defaultextract'
		  ;;
	      esac
	    fi
	    ;;
        esac
        ;;

    sysv4*uw2* | sysv5OpenUNIX* | sysv5UnixWare7.[[01]].[[10]]* | unixware7* | sco3.2v5.0.[[024]]*)
      _LT_TAGVAR(no_undefined_flag, $1)='$wl-z,text'
      _LT_TAGVAR(archive_cmds_need_lc, $1)=no
      _LT_TAGVAR(hardcode_shlibpath_var, $1)=no
      runpath_var='LD_RUN_PATH'

      case $cc_basename in
        CC*)
	  _LT_TAGVAR(archive_cmds, $1)='$CC -G $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	  _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -G $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
	*)
	  _LT_TAGVAR(archive_cmds, $1)='$CC -shared $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	  _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	  ;;
      esac
      ;;

      sysv5* | sco3.2v5* | sco5v6*)
	# Note: We CANNOT use -z defs as we might desire, because we do not
	# link with -lc, and that would cause any symbols used from libc to
	# always be unresolved, which means just about no library would
	# ever link correctly.  If we're not using GNU ld we use -z text
	# though, which does catch some bad symbols but isn't as heavy-handed
	# as -z defs.
	_LT_TAGVAR(no_undefined_flag, $1)='$wl-z,text'
	_LT_TAGVAR(allow_undefined_flag, $1)='$wl-z,nodefs'
	_LT_TAGVAR(archive_cmds_need_lc, $1)=no
	_LT_TAGVAR(hardcode_shlibpath_var, $1)=no
	_LT_TAGVAR(hardcode_libdir_flag_spec, $1)='$wl-R,$libdir'
	_LT_TAGVAR(hardcode_libdir_separator, $1)=':'
	_LT_TAGVAR(link_all_deplibs, $1)=yes
	_LT_TAGVAR(export_dynamic_flag_spec, $1)='$wl-Bexport'
	runpath_var='LD_RUN_PATH'

	case $cc_basename in
          CC*)
	    _LT_TAGVAR(archive_cmds, $1)='$CC -G $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	    _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -G $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	    _LT_TAGVAR(old_archive_cmds, $1)='$CC -Tprelink_objects $oldobjs~
              '"$_LT_TAGVAR(old_archive_cmds, $1)"
	    _LT_TAGVAR(reload_cmds, $1)='$CC -Tprelink_objects $reload_objs~
              '"$_LT_TAGVAR(reload_cmds, $1)"
	    ;;
	  *)
	    _LT_TAGVAR(archive_cmds, $1)='$CC -shared $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	    _LT_TAGVAR(archive_expsym_cmds, $1)='$CC -shared $wl-Bexport:$export_symbols $wl-h,$soname -o $lib $libobjs $deplibs $compiler_flags'
	    ;;
	esac
      ;;

      tandem*)
        case $cc_basename in
          NCC*)
	    # NonStop-UX NCC 3.20
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
          *)
	    # FIXME: insert proper C++ library support
	    _LT_TAGVAR(ld_shlibs, $1)=no
	    ;;
        esac
        ;;

      vxworks*)
        # FIXME: insert proper C++ library support
        _LT_TAGVAR(ld_shlibs, $1)=no
        ;;

      *)
        # FIXME: insert proper C++ library support
        _LT_TAGVAR(ld_shlibs, $1)=no
        ;;
    esac

    AC_MSG_RESULT([$_LT_TAGVAR(ld_shlibs, $1)])
    test no = "$_LT_TAGVAR(ld_shlibs, $1)" && can_build_shared=no

    _LT_TAGVAR(GCC, $1)=$GXX
    _LT_TAGVAR(LD, $1)=$LD

    ## CAVEAT EMPTOR:
    ## There is no encapsulation within the following macros, do not change
    ## the running order or otherwise move them around unless you know exactly
    ## what you are doing...
    _LT_SYS_HIDDEN_LIBDEPS($1)
    _LT_COMPILER_PIC($1)
    _LT_COMPILER_C_O($1)
    _LT_COMPILER_FILE_LOCKS($1)
    _LT_LINKER_SHLIBS($1)
    _LT_SYS_DYNAMIC_LINKER($1)
    _LT_LINKER_HARDCODE_LIBPATH($1)

    _LT_CONFIG($1)
  fi # test -n "$compiler"

  CC=$lt_save_CC
  CFLAGS=$lt_save_CFLAGS
  LDCXX=$LD
  LD=$lt_save_LD
  GCC=$lt_save_GCC
  with_gnu_ld=$lt_save_with_gnu_ld
  lt_cv_path_LDCXX=$lt_cv_path_LD
  lt_cv_path_LD=$lt_save_path_LD
  lt_cv_prog_gnu_ldcxx=$lt_cv_prog_gnu_ld
  lt_cv_prog_gnu_ld=$lt_save_with_gnu_ld
fi # test yes != "$_lt_caught_CXX_error"

AC_LANG_POP
])# _LT_LANG_CXX_CONFIG


# _LT_FUNC_STRIPNAME_CNF
# ----------------------
# func_stripname_cnf prefix suffix name
# strip PREFIX and SUFFIX off of NAME.
# PREFIX and SUFFIX must not contain globbing or regex special
# characters, hashes, percent signs, but SUFFIX may contain a leading
# dot (in which case that matches only a dot).
#
# This function is identical to the (non-XSI) version of func_stripname,
# except this one can be used by m4 code that may be executed by configure,
# rather than the libtool script.
m4_defun([_LT_FUNC_STRIPNAME_CNF],[dnl
AC_REQUIRE([_LT_DECL_SED])
AC_REQUIRE([_LT_PROG_ECHO_BACKSLASH])
func_stripname_cnf ()
{
  case @S|@2 in
  .*) func_stripname_result=`$ECHO "@S|@3" | $SED "s%^@S|@1%%; s%\\\\@S|@2\$%%"`;;
  *)  func_stripname_result=`$ECHO "@S|@3" | $SED "s%^@S|@1%%; s%@S|@2\$%%"`;;
  esac
} # func_stripname_cnf
])# _LT_FUNC_STRIPNAME_CNF


# _LT_SYS_HIDDEN_LIBDEPS([TAGNAME])
# ---------------------------------
# Figure out "hidden" library dependencies from verbose
# compiler output when linking a shared library.
# Parse the compiler output and extract the necessary
# objects, libraries and library flags.
m4_defun([_LT_SYS_HIDDEN_LIBDEPS],
[m4_require([_LT_FILEUTILS_DEFAULTS])dnl
AC_REQUIRE([_LT_FUNC_STRIPNAME_CNF])dnl
# Dependencies to place before and after the object being linked:
_LT_TAGVAR(predep_objects, $1)=
_LT_TAGVAR(postdep_objects, $1)=
_LT_TAGVAR(predeps, $1)=
_LT_TAGVAR(postdeps, $1)=
_LT_TAGVAR(compiler_lib_search_path, $1)=

dnl we can't use the lt_simple_compile_test_code here,
dnl because it contains code intended for an executable,
dnl not a library.  It's possible we should let each
dnl tag define a new lt_????_link_test_code variable,
dnl but it's only used here...
m4_if([$1], [], [cat > conftest.$ac_ext <<_LT_EOF
int a;
void foo (void) { a = 0; }
_LT_EOF
], [$1], [CXX], [cat > conftest.$ac_ext <<_LT_EOF
class Foo
{
public:
  Foo (void) { a = 0; }
private:
  int a;
};
_LT_EOF
], [$1], [F77], [cat > conftest.$ac_ext <<_LT_EOF
      subroutine foo
      implicit none
      integer*4 a
      a=0
      return
      end
_LT_EOF
], [$1], [FC], [cat > conftest.$ac_ext <<_LT_EOF
      subroutine foo
      implicit none
      integer a
      a=0
      return
      end
_LT_EOF
], [$1], [GCJ], [cat > conftest.$ac_ext <<_LT_EOF
public class foo {
  private int a;
  public void bar (void) {
    a = 0;
  }
};
_LT_EOF
], [$1], [GO], [cat > conftest.$ac_ext <<_LT_EOF
package foo
func foo() {
}
_LT_EOF
])

_lt_libdeps_save_CFLAGS=$CFLAGS
case "$CC $CFLAGS " in #(
*\ -flto*\ *) CFLAGS="$CFLAGS -fno-lto" ;;
*\ -fwhopr*\ *) CFLAGS="$CFLAGS -fno-whopr" ;;
*\ -fuse-linker-plugin*\ *) CFLAGS="$CFLAGS -fno-use-linker-plugin" ;;
esac

dnl Parse the compiler output and extract the necessary
dnl objects, libraries and library flags.
if AC_TRY_EVAL(ac_compile); then
  # Parse the compiler output and extract the necessary
  # objects, libraries and library flags.

  # Sentinel used to keep track of whether or not we are before
  # the conftest object file.
  pre_test_object_deps_done=no

  for p in `eval "$output_verbose_link_cmd"`; do
    case $prev$p in

    -L* | -R* | -l*)
       # Some compilers place space between "-{L,R}" and the path.
       # Remove the space.
       if test x-L = "$p" ||
          test x-R = "$p"; then
	 prev=$p
	 continue
       fi

       # Expand the sysroot to ease extracting the directories later.
       if test -z "$prev"; then
         case $p in
         -L*) func_stripname_cnf '-L' '' "$p"; prev=-L; p=$func_stripname_result ;;
         -R*) func_stripname_cnf '-R' '' "$p"; prev=-R; p=$func_stripname_result ;;
         -l*) func_stripname_cnf '-l' '' "$p"; prev=-l; p=$func_stripname_result ;;
         esac
       fi
       case $p in
       =*) func_stripname_cnf '=' '' "$p"; p=$lt_sysroot$func_stripname_result ;;
       esac
       if test no = "$pre_test_object_deps_done"; then
	 case $prev in
	 -L | -R)
	   # Internal compiler library paths should come after those
	   # provided the user.  The postdeps already come after the
	   # user supplied libs so there is no need to process them.
	   if test -z "$_LT_TAGVAR(compiler_lib_search_path, $1)"; then
	     _LT_TAGVAR(compiler_lib_search_path, $1)=$prev$p
	   else
	     _LT_TAGVAR(compiler_lib_search_path, $1)="${_LT_TAGVAR(compiler_lib_search_path, $1)} $prev$p"
	   fi
	   ;;
	 # The "-l" case would never come before the object being
	 # linked, so don't bother handling this case.
	 esac
       else
	 if test -z "$_LT_TAGVAR(postdeps, $1)"; then
	   _LT_TAGVAR(postdeps, $1)=$prev$p
	 else
	   _LT_TAGVAR(postdeps, $1)="${_LT_TAGVAR(postdeps, $1)} $prev$p"
	 fi
       fi
       prev=
       ;;

    *.lto.$objext) ;; # Ignore GCC LTO objects
    *.$objext)
       # This assumes that the test object file only shows up
       # once in the compiler output.
       if test "$p" = "conftest.$objext"; then
	 pre_test_object_deps_done=yes
	 continue
       fi

       if test no = "$pre_test_object_deps_done"; then
	 if test -z "$_LT_TAGVAR(predep_objects, $1)"; then
	   _LT_TAGVAR(predep_objects, $1)=$p
	 else
	   _LT_TAGVAR(predep_objects, $1)="$_LT_TAGVAR(predep_objects, $1) $p"
	 fi
       else
	 if test -z "$_LT_TAGVAR(postdep_objects, $1)"; then
	   _LT_TAGVAR(postdep_objects, $1)=$p
	 else
	   _LT_TAGVAR(postdep_objects, $1)="$_LT_TAGVAR(postdep_objects, $1) $p"
	 fi
       fi
       ;;

    *) ;; # Ignore the rest.

    esac
  done

  # Clean up.
  rm -f a.out a.exe
else
  echo "libtool.m4: error: problem compiling $1 test program"
fi

$RM -f confest.$objext
CFLAGS=$_lt_libdeps_save_CFLAGS

# PORTME: override above test on systems where it is broken
m4_if([$1], [CXX],
[case $host_os in
interix[[3-9]]*)
  # Interix 3.5 installs completely hosed .la files for C++, so rather than
  # hack all around it, let's just trust "g++" to DTRT.
  _LT_TAGVAR(predep_objects,$1)=
  _LT_TAGVAR(postdep_objects,$1)=
  _LT_TAGVAR(postdeps,$1)=
  ;;
esac
])

case " $_LT_TAGVAR(postdeps, $1) " in
*" -lc "*) _LT_TAGVAR(archive_cmds_need_lc, $1)=no ;;
esac
 _LT_TAGVAR(compiler_lib_search_dirs, $1)=
if test -n "${_LT_TAGVAR(compiler_lib_search_path, $1)}"; then
 _LT_TAGVAR(compiler_lib_search_dirs, $1)=`echo " ${_LT_TAGVAR(compiler_lib_search_path, $1)}" | $SED -e 's! -L! !g' -e 's!^ !!'`
fi
_LT_TAGDECL([], [compiler_lib_search_dirs], [1],
    [The directories searched by this compiler when creating a shared library])
_LT_TAGDECL([], [predep_objects], [1],
    [Dependencies to place before and after the objects being linked to
    create a shared library])
_LT_TAGDECL([], [postdep_objects], [1])
_LT_TAGDECL([], [predeps], [1])
_LT_TAGDECL([], [postdeps], [1])
_LT_TAGDECL([], [compiler_lib_search_path], [1],
    [The library search path used internally by the compiler when linking
    a shared library])
])# _LT_SYS_HIDDEN_LIBDEPS


# _LT_LANG_F77_CONFIG([TAG])
# --------------------------
# Ensure that the configuration variables for a Fortran 77 compiler are
# suitably defined.  These variables are subsequently used by _LT_CONFIG
# to write the compiler configuration to 'libtool'.
m4_defun([_LT_LANG_F77_CONFIG],
[AC_LANG_PUSH(Fortran 77)
if test -z "$F77" || test no = "$F77"; then
  _lt_disable_F77=yes
fi

_LT_TAGVAR(archive_cmds_need_lc, $1)=no
_LT_TAGVAR(allow_undefined_flag, $1)=
_LT_TAGVAR(always_export_symbols, $1)=no
_LT_TAGVAR(archive_expsym_cmds, $1)=
_LT_TAGVAR(export_dynamic_flag_spec, $1)=
_LT_TAGVAR(hardcode_direct, $1)=no
_LT_TAGVAR(hardcode_direct_absolute, $1)=no
_LT_TAGVAR(hardcode_libdir_flag_spec, $1)=
_LT_TAGVAR(hardcode_libdir_separator, $1)=
_LT_TAGVAR(hardcode_minus_L, $1)=no
_LT_TAGVAR(hardcode_automatic, $1)=no
_LT_TAGVAR(inherit_rpath, $1)=no
_LT_TAGVAR(module_cmds, $1)=
_LT_TAGVAR(module_expsym_cmds, $1)=
_LT_TAGVAR(link_all_deplibs, $1)=unknown
_LT_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
_LT_TAGVAR(reload_flag, $1)=$reload_flag
_LT_TAGVAR(reload_cmds, $1)=$reload_cmds
_LT_TAGVAR(no_undefined_flag, $1)=
_LT_TAGVAR(whole_archive_flag_spec, $1)=
_LT_TAGVAR(enable_shared_with_static_runtimes, $1)=no

# Source file extension for f77 test sources.
ac_ext=f

# Object file extension for compiled f77 test sources.
objext=o
_LT_TAGVAR(objext, $1)=$objext

# No sense in running all these tests if we already determined that
# the F77 compiler isn't working.  Some variables (like enable_shared)
# are currently assumed to apply to all compilers on this platform,
# and will be corrupted by setting them based on a non-working compiler.
if test yes != "$_lt_disable_F77"; then
  # Code to be used in simple compile tests
  lt_simple_compile_test_code="\
      subroutine t
      return
      end
"

  # Code to be used in simple link tests
  lt_simple_link_test_code="\
      program t
      end
"

  # ltmain only uses $CC for tagged configurations so make sure $CC is set.
  _LT_TAG_COMPILER

  # save warnings/boilerplate of simple test code
  _LT_COMPILER_BOILERPLATE
  _LT_LINKER_BOILERPLATE

  # Allow CC to be a program name with arguments.
  lt_save_CC=$CC
  lt_save_GCC=$GCC
  lt_save_CFLAGS=$CFLAGS
  CC=${F77-"f77"}
  CFLAGS=$FFLAGS
  compiler=$CC
  _LT_TAGVAR(compiler, $1)=$CC
  _LT_CC_BASENAME([$compiler])
  GCC=$G77
  if test -n "$compiler"; then
    AC_MSG_CHECKING([if libtool supports shared libraries])
    AC_MSG_RESULT([$can_build_shared])

    AC_MSG_CHECKING([whether to build shared libraries])
    test no = "$can_build_shared" && enable_shared=no

    # On AIX, shared libraries and static libraries use the same namespace, and
    # are all built from PIC.
    case $host_os in
      aix3*)
        test yes = "$enable_shared" && enable_static=no
        if test -n "$RANLIB"; then
          archive_cmds="$archive_cmds~\$RANLIB \$lib"
          postinstall_cmds='$RANLIB $lib'
        fi
        ;;
      aix[[4-9]]*)
	if test ia64 != "$host_cpu"; then
	  case $enable_shared,$with_aix_soname,$aix_use_runtimelinking in
	  yes,aix,yes) ;;		# shared object as lib.so file only
	  yes,svr4,*) ;;		# shared object as lib.so archive member only
	  yes,*) enable_static=no ;;	# shared object in lib.a archive as well
	  esac
	fi
        ;;
    esac
    AC_MSG_RESULT([$enable_shared])

    AC_MSG_CHECKING([whether to build static libraries])
    # Make sure either enable_shared or enable_static is yes.
    test yes = "$enable_shared" || enable_static=yes
    AC_MSG_RESULT([$enable_static])

    _LT_TAGVAR(GCC, $1)=$G77
    _LT_TAGVAR(LD, $1)=$LD

    ## CAVEAT EMPTOR:
    ## There is no encapsulation within the following macros, do not change
    ## the running order or otherwise move them around unless you know exactly
    ## what you are doing...
    _LT_COMPILER_PIC($1)
    _LT_COMPILER_C_O($1)
    _LT_COMPILER_FILE_LOCKS($1)
    _LT_LINKER_SHLIBS($1)
    _LT_SYS_DYNAMIC_LINKER($1)
    _LT_LINKER_HARDCODE_LIBPATH($1)

    _LT_CONFIG($1)
  fi # test -n "$compiler"

  GCC=$lt_save_GCC
  CC=$lt_save_CC
  CFLAGS=$lt_save_CFLAGS
fi # test yes != "$_lt_disable_F77"

AC_LANG_POP
])# _LT_LANG_F77_CONFIG


# _LT_LANG_FC_CONFIG([TAG])
# -------------------------
# Ensure that the configuration variables for a Fortran compiler are
# suitably defined.  These variables are subsequently used by _LT_CONFIG
# to write the compiler configuration to 'libtool'.
m4_defun([_LT_LANG_FC_CONFIG],
[AC_LANG_PUSH(Fortran)

if test -z "$FC" || test no = "$FC"; then
  _lt_disable_FC=yes
fi

_LT_TAGVAR(archive_cmds_need_lc, $1)=no
_LT_TAGVAR(allow_undefined_flag, $1)=
_LT_TAGVAR(always_export_symbols, $1)=no
_LT_TAGVAR(archive_expsym_cmds, $1)=
_LT_TAGVAR(export_dynamic_flag_spec, $1)=
_LT_TAGVAR(hardcode_direct, $1)=no
_LT_TAGVAR(hardcode_direct_absolute, $1)=no
_LT_TAGVAR(hardcode_libdir_flag_spec, $1)=
_LT_TAGVAR(hardcode_libdir_separator, $1)=
_LT_TAGVAR(hardcode_minus_L, $1)=no
_LT_TAGVAR(hardcode_automatic, $1)=no
_LT_TAGVAR(inherit_rpath, $1)=no
_LT_TAGVAR(module_cmds, $1)=
_LT_TAGVAR(module_expsym_cmds, $1)=
_LT_TAGVAR(link_all_deplibs, $1)=unknown
_LT_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
_LT_TAGVAR(reload_flag, $1)=$reload_flag
_LT_TAGVAR(reload_cmds, $1)=$reload_cmds
_LT_TAGVAR(no_undefined_flag, $1)=
_LT_TAGVAR(whole_archive_flag_spec, $1)=
_LT_TAGVAR(enable_shared_with_static_runtimes, $1)=no

# Source file extension for fc test sources.
ac_ext=${ac_fc_srcext-f}

# Object file extension for compiled fc test sources.
objext=o
_LT_TAGVAR(objext, $1)=$objext

# No sense in running all these tests if we already determined that
# the FC compiler isn't working.  Some variables (like enable_shared)
# are currently assumed to apply to all compilers on this platform,
# and will be corrupted by setting them based on a non-working compiler.
if test yes != "$_lt_disable_FC"; then
  # Code to be used in simple compile tests
  lt_simple_compile_test_code="\
      subroutine t
      return
      end
"

  # Code to be used in simple link tests
  lt_simple_link_test_code="\
      program t
      end
"

  # ltmain only uses $CC for tagged configurations so make sure $CC is set.
  _LT_TAG_COMPILER

  # save warnings/boilerplate of simple test code
  _LT_COMPILER_BOILERPLATE
  _LT_LINKER_BOILERPLATE

  # Allow CC to be a program name with arguments.
  lt_save_CC=$CC
  lt_save_GCC=$GCC
  lt_save_CFLAGS=$CFLAGS
  CC=${FC-"f95"}
  CFLAGS=$FCFLAGS
  compiler=$CC
  GCC=$ac_cv_fc_compiler_gnu

  _LT_TAGVAR(compiler, $1)=$CC
  _LT_CC_BASENAME([$compiler])

  if test -n "$compiler"; then
    AC_MSG_CHECKING([if libtool supports shared libraries])
    AC_MSG_RESULT([$can_build_shared])

    AC_MSG_CHECKING([whether to build shared libraries])
    test no = "$can_build_shared" && enable_shared=no

    # On AIX, shared libraries and static libraries use the same namespace, and
    # are all built from PIC.
    case $host_os in
      aix3*)
        test yes = "$enable_shared" && enable_static=no
        if test -n "$RANLIB"; then
          archive_cmds="$archive_cmds~\$RANLIB \$lib"
          postinstall_cmds='$RANLIB $lib'
        fi
        ;;
      aix[[4-9]]*)
	if test ia64 != "$host_cpu"; then
	  case $enable_shared,$with_aix_soname,$aix_use_runtimelinking in
	  yes,aix,yes) ;;		# shared object as lib.so file only
	  yes,svr4,*) ;;		# shared object as lib.so archive member only
	  yes,*) enable_static=no ;;	# shared object in lib.a archive as well
	  esac
	fi
        ;;
    esac
    AC_MSG_RESULT([$enable_shared])

    AC_MSG_CHECKING([whether to build static libraries])
    # Make sure either enable_shared or enable_static is yes.
    test yes = "$enable_shared" || enable_static=yes
    AC_MSG_RESULT([$enable_static])

    _LT_TAGVAR(GCC, $1)=$ac_cv_fc_compiler_gnu
    _LT_TAGVAR(LD, $1)=$LD

    ## CAVEAT EMPTOR:
    ## There is no encapsulation within the following macros, do not change
    ## the running order or otherwise move them around unless you know exactly
    ## what you are doing...
    _LT_SYS_HIDDEN_LIBDEPS($1)
    _LT_COMPILER_PIC($1)
    _LT_COMPILER_C_O($1)
    _LT_COMPILER_FILE_LOCKS($1)
    _LT_LINKER_SHLIBS($1)
    _LT_SYS_DYNAMIC_LINKER($1)
    _LT_LINKER_HARDCODE_LIBPATH($1)

    _LT_CONFIG($1)
  fi # test -n "$compiler"

  GCC=$lt_save_GCC
  CC=$lt_save_CC
  CFLAGS=$lt_save_CFLAGS
fi # test yes != "$_lt_disable_FC"

AC_LANG_POP
])# _LT_LANG_FC_CONFIG


# _LT_LANG_GCJ_CONFIG([TAG])
# --------------------------
# Ensure that the configuration variables for the GNU Java Compiler compiler
# are suitably defined.  These variables are subsequently used by _LT_CONFIG
# to write the compiler configuration to 'libtool'.
m4_defun([_LT_LANG_GCJ_CONFIG],
[AC_REQUIRE([LT_PROG_GCJ])dnl
AC_LANG_SAVE

# Source file extension for Java test sources.
ac_ext=java

# Object file extension for compiled Java test sources.
objext=o
_LT_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code="class foo {}"

# Code to be used in simple link tests
lt_simple_link_test_code='public class conftest { public static void main(String[[]] argv) {}; }'

# ltmain only uses $CC for tagged configurations so make sure $CC is set.
_LT_TAG_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

# Allow CC to be a program name with arguments.
lt_save_CC=$CC
lt_save_CFLAGS=$CFLAGS
lt_save_GCC=$GCC
GCC=yes
CC=${GCJ-"gcj"}
CFLAGS=$GCJFLAGS
compiler=$CC
_LT_TAGVAR(compiler, $1)=$CC
_LT_TAGVAR(LD, $1)=$LD
_LT_CC_BASENAME([$compiler])

# GCJ did not exist at the time GCC didn't implicitly link libc in.
_LT_TAGVAR(archive_cmds_need_lc, $1)=no

_LT_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
_LT_TAGVAR(reload_flag, $1)=$reload_flag
_LT_TAGVAR(reload_cmds, $1)=$reload_cmds

if test -n "$compiler"; then
  _LT_COMPILER_NO_RTTI($1)
  _LT_COMPILER_PIC($1)
  _LT_COMPILER_C_O($1)
  _LT_COMPILER_FILE_LOCKS($1)
  _LT_LINKER_SHLIBS($1)
  _LT_LINKER_HARDCODE_LIBPATH($1)

  _LT_CONFIG($1)
fi

AC_LANG_RESTORE

GCC=$lt_save_GCC
CC=$lt_save_CC
CFLAGS=$lt_save_CFLAGS
])# _LT_LANG_GCJ_CONFIG


# _LT_LANG_GO_CONFIG([TAG])
# --------------------------
# Ensure that the configuration variables for the GNU Go compiler
# are suitably defined.  These variables are subsequently used by _LT_CONFIG
# to write the compiler configuration to 'libtool'.
m4_defun([_LT_LANG_GO_CONFIG],
[AC_REQUIRE([LT_PROG_GO])dnl
AC_LANG_SAVE

# Source file extension for Go test sources.
ac_ext=go

# Object file extension for compiled Go test sources.
objext=o
_LT_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code="package main; func main() { }"

# Code to be used in simple link tests
lt_simple_link_test_code='package main; func main() { }'

# ltmain only uses $CC for tagged configurations so make sure $CC is set.
_LT_TAG_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

# Allow CC to be a program name with arguments.
lt_save_CC=$CC
lt_save_CFLAGS=$CFLAGS
lt_save_GCC=$GCC
GCC=yes
CC=${GOC-"gccgo"}
CFLAGS=$GOFLAGS
compiler=$CC
_LT_TAGVAR(compiler, $1)=$CC
_LT_TAGVAR(LD, $1)=$LD
_LT_CC_BASENAME([$compiler])

# Go did not exist at the time GCC didn't implicitly link libc in.
_LT_TAGVAR(archive_cmds_need_lc, $1)=no

_LT_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
_LT_TAGVAR(reload_flag, $1)=$reload_flag
_LT_TAGVAR(reload_cmds, $1)=$reload_cmds

if test -n "$compiler"; then
  _LT_COMPILER_NO_RTTI($1)
  _LT_COMPILER_PIC($1)
  _LT_COMPILER_C_O($1)
  _LT_COMPILER_FILE_LOCKS($1)
  _LT_LINKER_SHLIBS($1)
  _LT_LINKER_HARDCODE_LIBPATH($1)

  _LT_CONFIG($1)
fi

AC_LANG_RESTORE

GCC=$lt_save_GCC
CC=$lt_save_CC
CFLAGS=$lt_save_CFLAGS
])# _LT_LANG_GO_CONFIG


# _LT_LANG_RC_CONFIG([TAG])
# -------------------------
# Ensure that the configuration variables for the Windows resource compiler
# are suitably defined.  These variables are subsequently used by _LT_CONFIG
# to write the compiler configuration to 'libtool'.
m4_defun([_LT_LANG_RC_CONFIG],
[AC_REQUIRE([LT_PROG_RC])dnl
AC_LANG_SAVE

# Source file extension for RC test sources.
ac_ext=rc

# Object file extension for compiled RC test sources.
objext=o
_LT_TAGVAR(objext, $1)=$objext

# Code to be used in simple compile tests
lt_simple_compile_test_code='sample MENU { MENUITEM "&Soup", 100, CHECKED }'

# Code to be used in simple link tests
lt_simple_link_test_code=$lt_simple_compile_test_code

# ltmain only uses $CC for tagged configurations so make sure $CC is set.
_LT_TAG_COMPILER

# save warnings/boilerplate of simple test code
_LT_COMPILER_BOILERPLATE
_LT_LINKER_BOILERPLATE

# Allow CC to be a program name with arguments.
lt_save_CC=$CC
lt_save_CFLAGS=$CFLAGS
lt_save_GCC=$GCC
GCC=
CC=${RC-"windres"}
CFLAGS=
compiler=$CC
_LT_TAGVAR(compiler, $1)=$CC
_LT_CC_BASENAME([$compiler])
_LT_TAGVAR(lt_cv_prog_compiler_c_o, $1)=yes

if test -n "$compiler"; then
  :
  _LT_CONFIG($1)
fi

GCC=$lt_save_GCC
AC_LANG_RESTORE
CC=$lt_save_CC
CFLAGS=$lt_save_CFLAGS
])# _LT_LANG_RC_CONFIG


# LT_PROG_GCJ
# -----------
AC_DEFUN([LT_PROG_GCJ],
[m4_ifdef([AC_PROG_GCJ], [AC_PROG_GCJ],
  [m4_ifdef([A][M_PROG_GCJ], [A][M_PROG_GCJ],
    [AC_CHECK_TOOL(GCJ, gcj,)
      test set = "${GCJFLAGS+set}" || GCJFLAGS="-g -O2"
      AC_SUBST(GCJFLAGS)])])[]dnl
])

# Old name:
AU_ALIAS([LT_AC_PROG_GCJ], [LT_PROG_GCJ])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([LT_AC_PROG_GCJ], [])


# LT_PROG_GO
# ----------
AC_DEFUN([LT_PROG_GO],
[AC_CHECK_TOOL(GOC, gccgo,)
])


# LT_PROG_RC
# ----------
AC_DEFUN([LT_PROG_RC],
[AC_CHECK_TOOL(RC, windres,)
])

# Old name:
AU_ALIAS([LT_AC_PROG_RC], [LT_PROG_RC])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([LT_AC_PROG_RC], [])


# _LT_DECL_EGREP
# --------------
# If we don't have a new enough Autoconf to choose the best grep
# available, choose the one first in the user's PATH.
m4_defun([_LT_DECL_EGREP],
[AC_REQUIRE([AC_PROG_EGREP])dnl
AC_REQUIRE([AC_PROG_FGREP])dnl
test -z "$GREP" && GREP=grep
_LT_DECL([], [GREP], [1], [A grep program that handles long lines])
_LT_DECL([], [EGREP], [1], [An ERE matcher])
_LT_DECL([], [FGREP], [1], [A literal string matcher])
dnl Non-bleeding-edge autoconf doesn't subst GREP, so do it here too
AC_SUBST([GREP])
])


# _LT_DECL_OBJDUMP
# --------------
# If we don't have a new enough Autoconf to choose the best objdump
# available, choose the one first in the user's PATH.
m4_defun([_LT_DECL_OBJDUMP],
[AC_CHECK_TOOL(OBJDUMP, objdump, false)
test -z "$OBJDUMP" && OBJDUMP=objdump
_LT_DECL([], [OBJDUMP], [1], [An object symbol dumper])
AC_SUBST([OBJDUMP])
])

# _LT_DECL_DLLTOOL
# ----------------
# Ensure DLLTOOL variable is set.
m4_defun([_LT_DECL_DLLTOOL],
[AC_CHECK_TOOL(DLLTOOL, dlltool, false)
test -z "$DLLTOOL" && DLLTOOL=dlltool
_LT_DECL([], [DLLTOOL], [1], [DLL creation program])
AC_SUBST([DLLTOOL])
])

# _LT_DECL_SED
# ------------
# Check for a fully-functional sed program, that truncates
# as few characters as possible.  Prefer GNU sed if found.
m4_defun([_LT_DECL_SED],
[AC_PROG_SED
test -z "$SED" && SED=sed
Xsed="$SED -e 1s/^X//"
_LT_DECL([], [SED], [1], [A sed program that does not truncate output])
_LT_DECL([], [Xsed], ["\$SED -e 1s/^X//"],
    [Sed that helps us avoid accidentally triggering echo(1) options like -n])
])# _LT_DECL_SED

m4_ifndef([AC_PROG_SED], [
# NOTE: This macro has been submitted for inclusion into   #
#  GNU Autoconf as AC_PROG_SED.  When it is available in   #
#  a released version of Autoconf we should remove this    #
#  macro and use it instead.                               #

m4_defun([AC_PROG_SED],
[AC_MSG_CHECKING([for a sed that does not truncate output])
AC_CACHE_VAL(lt_cv_path_SED,
[# Loop through the user's path and test for sed and gsed.
# Then use that list of sed's as ones to test for truncation.
as_save_IFS=$IFS; IFS=$PATH_SEPARATOR
for as_dir in $PATH
do
  IFS=$as_save_IFS
  test -z "$as_dir" && as_dir=.
  for lt_ac_prog in sed gsed; do
    for ac_exec_ext in '' $ac_executable_extensions; do
      if $as_executable_p "$as_dir/$lt_ac_prog$ac_exec_ext"; then
        lt_ac_sed_list="$lt_ac_sed_list $as_dir/$lt_ac_prog$ac_exec_ext"
      fi
    done
  done
done
IFS=$as_save_IFS
lt_ac_max=0
lt_ac_count=0
# Add /usr/xpg4/bin/sed as it is typically found on Solaris
# along with /bin/sed that truncates output.
for lt_ac_sed in $lt_ac_sed_list /usr/xpg4/bin/sed; do
  test ! -f "$lt_ac_sed" && continue
  cat /dev/null > conftest.in
  lt_ac_count=0
  echo $ECHO_N "0123456789$ECHO_C" >conftest.in
  # Check for GNU sed and select it if it is found.
  if "$lt_ac_sed" --version 2>&1 < /dev/null | grep 'GNU' > /dev/null; then
    lt_cv_path_SED=$lt_ac_sed
    break
  fi
  while true; do
    cat conftest.in conftest.in >conftest.tmp
    mv conftest.tmp conftest.in
    cp conftest.in conftest.nl
    echo >>conftest.nl
    $lt_ac_sed -e 's/a$//' < conftest.nl >conftest.out || break
    cmp -s conftest.out conftest.nl || break
    # 10000 chars as input seems more than enough
    test 10 -lt "$lt_ac_count" && break
    lt_ac_count=`expr $lt_ac_count + 1`
    if test "$lt_ac_count" -gt "$lt_ac_max"; then
      lt_ac_max=$lt_ac_count
      lt_cv_path_SED=$lt_ac_sed
    fi
  done
done
])
SED=$lt_cv_path_SED
AC_SUBST([SED])
AC_MSG_RESULT([$SED])
])#AC_PROG_SED
])#m4_ifndef

# Old name:
AU_ALIAS([LT_AC_PROG_SED], [AC_PROG_SED])
dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([LT_AC_PROG_SED], [])


# _LT_CHECK_SHELL_FEATURES
# ------------------------
# Find out whether the shell is Bourne or XSI compatible,
# or has some other useful features.
m4_defun([_LT_CHECK_SHELL_FEATURES],
[if ( (MAIL=60; unset MAIL) || exit) >/dev/null 2>&1; then
  lt_unset=unset
else
  lt_unset=false
fi
_LT_DECL([], [lt_unset], [0], [whether the shell understands "unset"])dnl

# test EBCDIC or ASCII
case `echo X|tr X '\101'` in
 A) # ASCII based system
    # \n is not interpreted correctly by Solaris 8 /usr/ucb/tr
  lt_SP2NL='tr \040 \012'
  lt_NL2SP='tr \015\012 \040\040'
  ;;
 *) # EBCDIC based system
  lt_SP2NL='tr \100 \n'
  lt_NL2SP='tr \r\n \100\100'
  ;;
esac
_LT_DECL([SP2NL], [lt_SP2NL], [1], [turn spaces into newlines])dnl
_LT_DECL([NL2SP], [lt_NL2SP], [1], [turn newlines into spaces])dnl
])# _LT_CHECK_SHELL_FEATURES


# _LT_PATH_CONVERSION_FUNCTIONS
# -----------------------------
# Determine what file name conversion functions should be used by
# func_to_host_file (and, implicitly, by func_to_host_path).  These are needed
# for certain cross-compile configurations and native mingw.
m4_defun([_LT_PATH_CONVERSION_FUNCTIONS],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
AC_MSG_CHECKING([how to convert $build file names to $host format])
AC_CACHE_VAL(lt_cv_to_host_file_cmd,
[case $host in
  *-*-mingw* )
    case $build in
      *-*-mingw* ) # actually msys
        lt_cv_to_host_file_cmd=func_convert_file_msys_to_w32
        ;;
      *-*-cygwin* )
        lt_cv_to_host_file_cmd=func_convert_file_cygwin_to_w32
        ;;
      * ) # otherwise, assume *nix
        lt_cv_to_host_file_cmd=func_convert_file_nix_to_w32
        ;;
    esac
    ;;
  *-*-cygwin* )
    case $build in
      *-*-mingw* ) # actually msys
        lt_cv_to_host_file_cmd=func_convert_file_msys_to_cygwin
        ;;
      *-*-cygwin* )
        lt_cv_to_host_file_cmd=func_convert_file_noop
        ;;
      * ) # otherwise, assume *nix
        lt_cv_to_host_file_cmd=func_convert_file_nix_to_cygwin
        ;;
    esac
    ;;
  * ) # unhandled hosts (and "normal" native builds)
    lt_cv_to_host_file_cmd=func_convert_file_noop
    ;;
esac
])
to_host_file_cmd=$lt_cv_to_host_file_cmd
AC_MSG_RESULT([$lt_cv_to_host_file_cmd])
_LT_DECL([to_host_file_cmd], [lt_cv_to_host_file_cmd],
         [0], [convert $build file names to $host format])dnl

AC_MSG_CHECKING([how to convert $build file names to toolchain format])
AC_CACHE_VAL(lt_cv_to_tool_file_cmd,
[#assume ordinary cross tools, or native build.
lt_cv_to_tool_file_cmd=func_convert_file_noop
case $host in
  *-*-mingw* )
    case $build in
      *-*-mingw* ) # actually msys
        lt_cv_to_tool_file_cmd=func_convert_file_msys_to_w32
        ;;
    esac
    ;;
esac
])
to_tool_file_cmd=$lt_cv_to_tool_file_cmd
AC_MSG_RESULT([$lt_cv_to_tool_file_cmd])
_LT_DECL([to_tool_file_cmd], [lt_cv_to_tool_file_cmd],
         [0], [convert $build files to toolchain format])dnl
])# _LT_PATH_CONVERSION_FUNCTIONS

# Helper functions for option handling.                    -*- Autoconf -*-
#
#   Copyright (C) 2004-2005, 2007-2009, 2011-2015 Free Software
#   Foundation, Inc.
#   Written by Gary V. Vaughan, 2004
#
# This file is free software; the Free Software Foundation gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.

# serial 8 ltoptions.m4

# This is to help aclocal find these macros, as it can't see m4_define.
AC_DEFUN([LTOPTIONS_VERSION], [m4_if([1])])


# _LT_MANGLE_OPTION(MACRO-NAME, OPTION-NAME)
# ------------------------------------------
m4_define([_LT_MANGLE_OPTION],
[[_LT_OPTION_]m4_bpatsubst($1__$2, [[^a-zA-Z0-9_]], [_])])


# _LT_SET_OPTION(MACRO-NAME, OPTION-NAME)
# ---------------------------------------
# Set option OPTION-NAME for macro MACRO-NAME, and if there is a
# matching handler defined, dispatch to it.  Other OPTION-NAMEs are
# saved as a flag.
m4_define([_LT_SET_OPTION],
[m4_define(_LT_MANGLE_OPTION([$1], [$2]))dnl
m4_ifdef(_LT_MANGLE_DEFUN([$1], [$2]),
        _LT_MANGLE_DEFUN([$1], [$2]),
    [m4_warning([Unknown $1 option '$2'])])[]dnl
])


# _LT_IF_OPTION(MACRO-NAME, OPTION-NAME, IF-SET, [IF-NOT-SET])
# ------------------------------------------------------------
# Execute IF-SET if OPTION is set, IF-NOT-SET otherwise.
m4_define([_LT_IF_OPTION],
[m4_ifdef(_LT_MANGLE_OPTION([$1], [$2]), [$3], [$4])])


# _LT_UNLESS_OPTIONS(MACRO-NAME, OPTION-LIST, IF-NOT-SET)
# -------------------------------------------------------
# Execute IF-NOT-SET unless all options in OPTION-LIST for MACRO-NAME
# are set.
m4_define([_LT_UNLESS_OPTIONS],
[m4_foreach([_LT_Option], m4_split(m4_normalize([$2])),
	    [m4_ifdef(_LT_MANGLE_OPTION([$1], _LT_Option),
		      [m4_define([$0_found])])])[]dnl
m4_ifdef([$0_found], [m4_undefine([$0_found])], [$3
])[]dnl
])


# _LT_SET_OPTIONS(MACRO-NAME, OPTION-LIST)
# ----------------------------------------
# OPTION-LIST is a space-separated list of Libtool options associated
# with MACRO-NAME.  If any OPTION has a matching handler declared with
# LT_OPTION_DEFINE, dispatch to that macro; otherwise complain about
# the unknown option and exit.
m4_defun([_LT_SET_OPTIONS],
[# Set options
m4_foreach([_LT_Option], m4_split(m4_normalize([$2])),
    [_LT_SET_OPTION([$1], _LT_Option)])

m4_if([$1],[LT_INIT],[
  dnl
  dnl Simply set some default values (i.e off) if boolean options were not
  dnl specified:
  _LT_UNLESS_OPTIONS([LT_INIT], [dlopen], [enable_dlopen=no
  ])
  _LT_UNLESS_OPTIONS([LT_INIT], [win32-dll], [enable_win32_dll=no
  ])
  dnl
  dnl If no reference was made to various pairs of opposing options, then
  dnl we run the default mode handler for the pair.  For example, if neither
  dnl 'shared' nor 'disable-shared' was passed, we enable building of shared
  dnl archives by default:
  _LT_UNLESS_OPTIONS([LT_INIT], [shared disable-shared], [_LT_ENABLE_SHARED])
  _LT_UNLESS_OPTIONS([LT_INIT], [static disable-static], [_LT_ENABLE_STATIC])
  _LT_UNLESS_OPTIONS([LT_INIT], [pic-only no-pic], [_LT_WITH_PIC])
  _LT_UNLESS_OPTIONS([LT_INIT], [fast-install disable-fast-install],
		   [_LT_ENABLE_FAST_INSTALL])
  _LT_UNLESS_OPTIONS([LT_INIT], [aix-soname=aix aix-soname=both aix-soname=svr4],
		   [_LT_WITH_AIX_SONAME([aix])])
  ])
])# _LT_SET_OPTIONS



# _LT_MANGLE_DEFUN(MACRO-NAME, OPTION-NAME)
# -----------------------------------------
m4_define([_LT_MANGLE_DEFUN],
[[_LT_OPTION_DEFUN_]m4_bpatsubst(m4_toupper([$1__$2]), [[^A-Z0-9_]], [_])])


# LT_OPTION_DEFINE(MACRO-NAME, OPTION-NAME, CODE)
# -----------------------------------------------
m4_define([LT_OPTION_DEFINE],
[m4_define(_LT_MANGLE_DEFUN([$1], [$2]), [$3])[]dnl
])# LT_OPTION_DEFINE


# dlopen
# ------
LT_OPTION_DEFINE([LT_INIT], [dlopen], [enable_dlopen=yes
])

AU_DEFUN([AC_LIBTOOL_DLOPEN],
[_LT_SET_OPTION([LT_INIT], [dlopen])
AC_DIAGNOSE([obsolete],
[$0: Remove this warning and the call to _LT_SET_OPTION when you
put the 'dlopen' option into LT_INIT's first parameter.])
])

dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_DLOPEN], [])


# win32-dll
# ---------
# Declare package support for building win32 dll's.
LT_OPTION_DEFINE([LT_INIT], [win32-dll],
[enable_win32_dll=yes

case $host in
*-*-cygwin* | *-*-mingw* | *-*-pw32* | *-*-cegcc*)
  AC_CHECK_TOOL(AS, as, false)
  AC_CHECK_TOOL(DLLTOOL, dlltool, false)
  AC_CHECK_TOOL(OBJDUMP, objdump, false)
  ;;
esac

test -z "$AS" && AS=as
_LT_DECL([], [AS],      [1], [Assembler program])dnl

test -z "$DLLTOOL" && DLLTOOL=dlltool
_LT_DECL([], [DLLTOOL], [1], [DLL creation program])dnl

test -z "$OBJDUMP" && OBJDUMP=objdump
_LT_DECL([], [OBJDUMP], [1], [Object dumper program])dnl
])# win32-dll

AU_DEFUN([AC_LIBTOOL_WIN32_DLL],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
_LT_SET_OPTION([LT_INIT], [win32-dll])
AC_DIAGNOSE([obsolete],
[$0: Remove this warning and the call to _LT_SET_OPTION when you
put the 'win32-dll' option into LT_INIT's first parameter.])
])

dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_WIN32_DLL], [])


# _LT_ENABLE_SHARED([DEFAULT])
# ----------------------------
# implement the --enable-shared flag, and supports the 'shared' and
# 'disable-shared' LT_INIT options.
# DEFAULT is either 'yes' or 'no'.  If omitted, it defaults to 'yes'.
m4_define([_LT_ENABLE_SHARED],
[m4_define([_LT_ENABLE_SHARED_DEFAULT], [m4_if($1, no, no, yes)])dnl
AC_ARG_ENABLE([shared],
    [AS_HELP_STRING([--enable-shared@<:@=PKGS@:>@],
	[build shared libraries @<:@default=]_LT_ENABLE_SHARED_DEFAULT[@:>@])],
    [p=${PACKAGE-default}
    case $enableval in
    yes) enable_shared=yes ;;
    no) enable_shared=no ;;
    *)
      enable_shared=no
      # Look at the argument we got.  We use all the common list separators.
      lt_save_ifs=$IFS; IFS=$IFS$PATH_SEPARATOR,
      for pkg in $enableval; do
	IFS=$lt_save_ifs
	if test "X$pkg" = "X$p"; then
	  enable_shared=yes
	fi
      done
      IFS=$lt_save_ifs
      ;;
    esac],
    [enable_shared=]_LT_ENABLE_SHARED_DEFAULT)

    _LT_DECL([build_libtool_libs], [enable_shared], [0],
	[Whether or not to build shared libraries])
])# _LT_ENABLE_SHARED

LT_OPTION_DEFINE([LT_INIT], [shared], [_LT_ENABLE_SHARED([yes])])
LT_OPTION_DEFINE([LT_INIT], [disable-shared], [_LT_ENABLE_SHARED([no])])

# Old names:
AC_DEFUN([AC_ENABLE_SHARED],
[_LT_SET_OPTION([LT_INIT], m4_if([$1], [no], [disable-])[shared])
])

AC_DEFUN([AC_DISABLE_SHARED],
[_LT_SET_OPTION([LT_INIT], [disable-shared])
])

AU_DEFUN([AM_ENABLE_SHARED], [AC_ENABLE_SHARED($@)])
AU_DEFUN([AM_DISABLE_SHARED], [AC_DISABLE_SHARED($@)])

dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AM_ENABLE_SHARED], [])
dnl AC_DEFUN([AM_DISABLE_SHARED], [])



# _LT_ENABLE_STATIC([DEFAULT])
# ----------------------------
# implement the --enable-static flag, and support the 'static' and
# 'disable-static' LT_INIT options.
# DEFAULT is either 'yes' or 'no'.  If omitted, it defaults to 'yes'.
m4_define([_LT_ENABLE_STATIC],
[m4_define([_LT_ENABLE_STATIC_DEFAULT], [m4_if($1, no, no, yes)])dnl
AC_ARG_ENABLE([static],
    [AS_HELP_STRING([--enable-static@<:@=PKGS@:>@],
	[build static libraries @<:@default=]_LT_ENABLE_STATIC_DEFAULT[@:>@])],
    [p=${PACKAGE-default}
    case $enableval in
    yes) enable_static=yes ;;
    no) enable_static=no ;;
    *)
     enable_static=no
      # Look at the argument we got.  We use all the common list separators.
      lt_save_ifs=$IFS; IFS=$IFS$PATH_SEPARATOR,
      for pkg in $enableval; do
	IFS=$lt_save_ifs
	if test "X$pkg" = "X$p"; then
	  enable_static=yes
	fi
      done
      IFS=$lt_save_ifs
      ;;
    esac],
    [enable_static=]_LT_ENABLE_STATIC_DEFAULT)

    _LT_DECL([build_old_libs], [enable_static], [0],
	[Whether or not to build static libraries])
])# _LT_ENABLE_STATIC

LT_OPTION_DEFINE([LT_INIT], [static], [_LT_ENABLE_STATIC([yes])])
LT_OPTION_DEFINE([LT_INIT], [disable-static], [_LT_ENABLE_STATIC([no])])

# Old names:
AC_DEFUN([AC_ENABLE_STATIC],
[_LT_SET_OPTION([LT_INIT], m4_if([$1], [no], [disable-])[static])
])

AC_DEFUN([AC_DISABLE_STATIC],
[_LT_SET_OPTION([LT_INIT], [disable-static])
])

AU_DEFUN([AM_ENABLE_STATIC], [AC_ENABLE_STATIC($@)])
AU_DEFUN([AM_DISABLE_STATIC], [AC_DISABLE_STATIC($@)])

dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AM_ENABLE_STATIC], [])
dnl AC_DEFUN([AM_DISABLE_STATIC], [])



# _LT_ENABLE_FAST_INSTALL([DEFAULT])
# ----------------------------------
# implement the --enable-fast-install flag, and support the 'fast-install'
# and 'disable-fast-install' LT_INIT options.
# DEFAULT is either 'yes' or 'no'.  If omitted, it defaults to 'yes'.
m4_define([_LT_ENABLE_FAST_INSTALL],
[m4_define([_LT_ENABLE_FAST_INSTALL_DEFAULT], [m4_if($1, no, no, yes)])dnl
AC_ARG_ENABLE([fast-install],
    [AS_HELP_STRING([--enable-fast-install@<:@=PKGS@:>@],
    [optimize for fast installation @<:@default=]_LT_ENABLE_FAST_INSTALL_DEFAULT[@:>@])],
    [p=${PACKAGE-default}
    case $enableval in
    yes) enable_fast_install=yes ;;
    no) enable_fast_install=no ;;
    *)
      enable_fast_install=no
      # Look at the argument we got.  We use all the common list separators.
      lt_save_ifs=$IFS; IFS=$IFS$PATH_SEPARATOR,
      for pkg in $enableval; do
	IFS=$lt_save_ifs
	if test "X$pkg" = "X$p"; then
	  enable_fast_install=yes
	fi
      done
      IFS=$lt_save_ifs
      ;;
    esac],
    [enable_fast_install=]_LT_ENABLE_FAST_INSTALL_DEFAULT)

_LT_DECL([fast_install], [enable_fast_install], [0],
	 [Whether or not to optimize for fast installation])dnl
])# _LT_ENABLE_FAST_INSTALL

LT_OPTION_DEFINE([LT_INIT], [fast-install], [_LT_ENABLE_FAST_INSTALL([yes])])
LT_OPTION_DEFINE([LT_INIT], [disable-fast-install], [_LT_ENABLE_FAST_INSTALL([no])])

# Old names:
AU_DEFUN([AC_ENABLE_FAST_INSTALL],
[_LT_SET_OPTION([LT_INIT], m4_if([$1], [no], [disable-])[fast-install])
AC_DIAGNOSE([obsolete],
[$0: Remove this warning and the call to _LT_SET_OPTION when you put
the 'fast-install' option into LT_INIT's first parameter.])
])

AU_DEFUN([AC_DISABLE_FAST_INSTALL],
[_LT_SET_OPTION([LT_INIT], [disable-fast-install])
AC_DIAGNOSE([obsolete],
[$0: Remove this warning and the call to _LT_SET_OPTION when you put
the 'disable-fast-install' option into LT_INIT's first parameter.])
])

dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_ENABLE_FAST_INSTALL], [])
dnl AC_DEFUN([AM_DISABLE_FAST_INSTALL], [])


# _LT_WITH_AIX_SONAME([DEFAULT])
# ----------------------------------
# implement the --with-aix-soname flag, and support the `aix-soname=aix'
# and `aix-soname=both' and `aix-soname=svr4' LT_INIT options. DEFAULT
# is either `aix', `both' or `svr4'.  If omitted, it defaults to `aix'.
m4_define([_LT_WITH_AIX_SONAME],
[m4_define([_LT_WITH_AIX_SONAME_DEFAULT], [m4_if($1, svr4, svr4, m4_if($1, both, both, aix))])dnl
shared_archive_member_spec=
case $host,$enable_shared in
power*-*-aix[[5-9]]*,yes)
  AC_MSG_CHECKING([which variant of shared library versioning to provide])
  AC_ARG_WITH([aix-soname],
    [AS_HELP_STRING([--with-aix-soname=aix|svr4|both],
      [shared library versioning (aka "SONAME") variant to provide on AIX, @<:@default=]_LT_WITH_AIX_SONAME_DEFAULT[@:>@.])],
    [case $withval in
    aix|svr4|both)
      ;;
    *)
      AC_MSG_ERROR([Unknown argument to --with-aix-soname])
      ;;
    esac
    lt_cv_with_aix_soname=$with_aix_soname],
    [AC_CACHE_VAL([lt_cv_with_aix_soname],
      [lt_cv_with_aix_soname=]_LT_WITH_AIX_SONAME_DEFAULT)
    with_aix_soname=$lt_cv_with_aix_soname])
  AC_MSG_RESULT([$with_aix_soname])
  if test aix != "$with_aix_soname"; then
    # For the AIX way of multilib, we name the shared archive member
    # based on the bitwidth used, traditionally 'shr.o' or 'shr_64.o',
    # and 'shr.imp' or 'shr_64.imp', respectively, for the Import File.
    # Even when GNU compilers ignore OBJECT_MODE but need '-maix64' flag,
    # the AIX toolchain works better with OBJECT_MODE set (default 32).
    if test 64 = "${OBJECT_MODE-32}"; then
      shared_archive_member_spec=shr_64
    else
      shared_archive_member_spec=shr
    fi
  fi
  ;;
*)
  with_aix_soname=aix
  ;;
esac

_LT_DECL([], [shared_archive_member_spec], [0],
    [Shared archive member basename, for filename based shared library versioning on AIX])dnl
])# _LT_WITH_AIX_SONAME

LT_OPTION_DEFINE([LT_INIT], [aix-soname=aix], [_LT_WITH_AIX_SONAME([aix])])
LT_OPTION_DEFINE([LT_INIT], [aix-soname=both], [_LT_WITH_AIX_SONAME([both])])
LT_OPTION_DEFINE([LT_INIT], [aix-soname=svr4], [_LT_WITH_AIX_SONAME([svr4])])


# _LT_WITH_PIC([MODE])
# --------------------
# implement the --with-pic flag, and support the 'pic-only' and 'no-pic'
# LT_INIT options.
# MODE is either 'yes' or 'no'.  If omitted, it defaults to 'both'.
m4_define([_LT_WITH_PIC],
[AC_ARG_WITH([pic],
    [AS_HELP_STRING([--with-pic@<:@=PKGS@:>@],
	[try to use only PIC/non-PIC objects @<:@default=use both@:>@])],
    [lt_p=${PACKAGE-default}
    case $withval in
    yes|no) pic_mode=$withval ;;
    *)
      pic_mode=default
      # Look at the argument we got.  We use all the common list separators.
      lt_save_ifs=$IFS; IFS=$IFS$PATH_SEPARATOR,
      for lt_pkg in $withval; do
	IFS=$lt_save_ifs
	if test "X$lt_pkg" = "X$lt_p"; then
	  pic_mode=yes
	fi
      done
      IFS=$lt_save_ifs
      ;;
    esac],
    [pic_mode=m4_default([$1], [default])])

_LT_DECL([], [pic_mode], [0], [What type of objects to build])dnl
])# _LT_WITH_PIC

LT_OPTION_DEFINE([LT_INIT], [pic-only], [_LT_WITH_PIC([yes])])
LT_OPTION_DEFINE([LT_INIT], [no-pic], [_LT_WITH_PIC([no])])

# Old name:
AU_DEFUN([AC_LIBTOOL_PICMODE],
[_LT_SET_OPTION([LT_INIT], [pic-only])
AC_DIAGNOSE([obsolete],
[$0: Remove this warning and the call to _LT_SET_OPTION when you
put the 'pic-only' option into LT_INIT's first parameter.])
])

dnl aclocal-1.4 backwards compatibility:
dnl AC_DEFUN([AC_LIBTOOL_PICMODE], [])


m4_define([_LTDL_MODE], [])
LT_OPTION_DEFINE([LTDL_INIT], [nonrecursive],
		 [m4_define([_LTDL_MODE], [nonrecursive])])
LT_OPTION_DEFINE([LTDL_INIT], [recursive],
		 [m4_define([_LTDL_MODE], [recursive])])
LT_OPTION_DEFINE([LTDL_INIT], [subproject],
		 [m4_define([_LTDL_MODE], [subproject])])

m4_define([_LTDL_TYPE], [])
LT_OPTION_DEFINE([LTDL_INIT], [installable],
		 [m4_define([_LTDL_TYPE], [installable])])
LT_OPTION_DEFINE([LTDL_INIT], [convenience],
		 [m4_define([_LTDL_TYPE], [convenience])])

# ltsugar.m4 -- libtool m4 base layer.                         -*-Autoconf-*-
#
# Copyright (C) 2004-2005, 2007-2008, 2011-2015 Free Software
# Foundation, Inc.
# Written by Gary V. Vaughan, 2004
#
# This file is free software; the Free Software Foundation gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.

# serial 6 ltsugar.m4

# This is to help aclocal find these macros, as it can't see m4_define.
AC_DEFUN([LTSUGAR_VERSION], [m4_if([0.1])])


# lt_join(SEP, ARG1, [ARG2...])
# -----------------------------
# Produce ARG1SEPARG2...SEPARGn, omitting [] arguments and their
# associated separator.
# Needed until we can rely on m4_join from Autoconf 2.62, since all earlier
# versions in m4sugar had bugs.
m4_define([lt_join],
[m4_if([$#], [1], [],
       [$#], [2], [[$2]],
       [m4_if([$2], [], [], [[$2]_])$0([$1], m4_shift(m4_shift($@)))])])
m4_define([_lt_join],
[m4_if([$#$2], [2], [],
       [m4_if([$2], [], [], [[$1$2]])$0([$1], m4_shift(m4_shift($@)))])])


# lt_car(LIST)
# lt_cdr(LIST)
# ------------
# Manipulate m4 lists.
# These macros are necessary as long as will still need to support
# Autoconf-2.59, which quotes differently.
m4_define([lt_car], [[$1]])
m4_define([lt_cdr],
[m4_if([$#], 0, [m4_fatal([$0: cannot be called without arguments])],
       [$#], 1, [],
       [m4_dquote(m4_shift($@))])])
m4_define([lt_unquote], $1)


# lt_append(MACRO-NAME, STRING, [SEPARATOR])
# ------------------------------------------
# Redefine MACRO-NAME to hold its former content plus 'SEPARATOR''STRING'.
# Note that neither SEPARATOR nor STRING are expanded; they are appended
# to MACRO-NAME as is (leaving the expansion for when MACRO-NAME is invoked).
# No SEPARATOR is output if MACRO-NAME was previously undefined (different
# than defined and empty).
#
# This macro is needed until we can rely on Autoconf 2.62, since earlier
# versions of m4sugar mistakenly expanded SEPARATOR but not STRING.
m4_define([lt_append],
[m4_define([$1],
	   m4_ifdef([$1], [m4_defn([$1])[$3]])[$2])])



# lt_combine(SEP, PREFIX-LIST, INFIX, SUFFIX1, [SUFFIX2...])
# ----------------------------------------------------------
# Produce a SEP delimited list of all paired combinations of elements of
# PREFIX-LIST with SUFFIX1 through SUFFIXn.  Each element of the list
# has the form PREFIXmINFIXSUFFIXn.
# Needed until we can rely on m4_combine added in Autoconf 2.62.
m4_define([lt_combine],
[m4_if(m4_eval([$# > 3]), [1],
       [m4_pushdef([_Lt_sep], [m4_define([_Lt_sep], m4_defn([lt_car]))])]]dnl
[[m4_foreach([_Lt_prefix], [$2],
	     [m4_foreach([_Lt_suffix],
		]m4_dquote(m4_dquote(m4_shift(m4_shift(m4_shift($@)))))[,
	[_Lt_sep([$1])[]m4_defn([_Lt_prefix])[$3]m4_defn([_Lt_suffix])])])])])


# lt_if_append_uniq(MACRO-NAME, VARNAME, [SEPARATOR], [UNIQ], [NOT-UNIQ])
# -----------------------------------------------------------------------
# Iff MACRO-NAME does not yet contain VARNAME, then append it (delimited
# by SEPARATOR if supplied) and expand UNIQ, else NOT-UNIQ.
m4_define([lt_if_append_uniq],
[m4_ifdef([$1],
	  [m4_if(m4_index([$3]m4_defn([$1])[$3], [$3$2$3]), [-1],
		 [lt_append([$1], [$2], [$3])$4],
		 [$5])],
	  [lt_append([$1], [$2], [$3])$4])])


# lt_dict_add(DICT, KEY, VALUE)
# -----------------------------
m4_define([lt_dict_add],
[m4_define([$1($2)], [$3])])


# lt_dict_add_subkey(DICT, KEY, SUBKEY, VALUE)
# --------------------------------------------
m4_define([lt_dict_add_subkey],
[m4_define([$1($2:$3)], [$4])])


# lt_dict_fetch(DICT, KEY, [SUBKEY])
# ----------------------------------
m4_define([lt_dict_fetch],
[m4_ifval([$3],
	m4_ifdef([$1($2:$3)], [m4_defn([$1($2:$3)])]),
    m4_ifdef([$1($2)], [m4_defn([$1($2)])]))])


# lt_if_dict_fetch(DICT, KEY, [SUBKEY], VALUE, IF-TRUE, [IF-FALSE])
# -----------------------------------------------------------------
m4_define([lt_if_dict_fetch],
[m4_if(lt_dict_fetch([$1], [$2], [$3]), [$4],
	[$5],
    [$6])])


# lt_dict_filter(DICT, [SUBKEY], VALUE, [SEPARATOR], KEY, [...])
# --------------------------------------------------------------
m4_define([lt_dict_filter],
[m4_if([$5], [], [],
  [lt_join(m4_quote(m4_default([$4], [[, ]])),
           lt_unquote(m4_split(m4_normalize(m4_foreach(_Lt_key, lt_car([m4_shiftn(4, $@)]),
		      [lt_if_dict_fetch([$1], _Lt_key, [$2], [$3], [_Lt_key ])])))))])[]dnl
])

# ltversion.m4 -- version numbers			-*- Autoconf -*-
#
#   Copyright (C) 2004, 2011-2015 Free Software Foundation, Inc.
#   Written by Scott James Remnant, 2004
#
# This file is free software; the Free Software Foundation gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.

# @configure_input@

# serial 4179 ltversion.m4
# This file is part of GNU Libtool

m4_define([LT_PACKAGE_VERSION], [2.4.6])
m4_define([LT_PACKAGE_REVISION], [2.4.6])

AC_DEFUN([LTVERSION_VERSION],
[macro_version='2.4.6'
macro_revision='2.4.6'
_LT_DECL(, macro_version, 0, [Which release of libtool.m4 was used?])
_LT_DECL(, macro_revision, 0)
])

# lt~obsolete.m4 -- aclocal satisfying obsolete definitions.    -*-Autoconf-*-
#
#   Copyright (C) 2004-2005, 2007, 2009, 2011-2015 Free Software
#   Foundation, Inc.
#   Written by Scott James Remnant, 2004.
#
# This file is free software; the Free Software Foundation gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.

# serial 5 lt~obsolete.m4

# These exist entirely to fool aclocal when bootstrapping libtool.
#
# In the past libtool.m4 has provided macros via AC_DEFUN (or AU_DEFUN),
# which have later been changed to m4_define as they aren't part of the
# exported API, or moved to Autoconf or Automake where they belong.
#
# The trouble is, aclocal is a bit thick.  It'll see the old AC_DEFUN
# in /usr/share/aclocal/libtool.m4 and remember it, then when it sees us
# using a macro with the same name in our local m4/libtool.m4 it'll
# pull the old libtool.m4 in (it doesn't see our shiny new m4_define
# and doesn't know about Autoconf macros at all.)
#
# So we provide this file, which has a silly filename so it's always
# included after everything else.  This provides aclocal with the
# AC_DEFUNs it wants, but when m4 processes it, it doesn't do anything
# because those macros already exist, or will be overwritten later.
# We use AC_DEFUN over AU_DEFUN for compatibility with aclocal-1.6.
#
# Anytime we withdraw an AC_DEFUN or AU_DEFUN, remember to add it here.
# Yes, that means every name once taken will need to remain here until
# we give up compatibility with versions before 1.7, at which point
# we need to keep only those names which we still refer to.

# This is to help aclocal find these macros, as it can't see m4_define.
AC_DEFUN([LTOBSOLETE_VERSION], [m4_if([1])])

m4_ifndef([AC_LIBTOOL_LINKER_OPTION],	[AC_DEFUN([AC_LIBTOOL_LINKER_OPTION])])
m4_ifndef([AC_PROG_EGREP],		[AC_DEFUN([AC_PROG_EGREP])])
m4_ifndef([_LT_AC_PROG_ECHO_BACKSLASH],	[AC_DEFUN([_LT_AC_PROG_ECHO_BACKSLASH])])
m4_ifndef([_LT_AC_SHELL_INIT],		[AC_DEFUN([_LT_AC_SHELL_INIT])])
m4_ifndef([_LT_AC_SYS_LIBPATH_AIX],	[AC_DEFUN([_LT_AC_SYS_LIBPATH_AIX])])
m4_ifndef([_LT_PROG_LTMAIN],		[AC_DEFUN([_LT_PROG_LTMAIN])])
m4_ifndef([_LT_AC_TAGVAR],		[AC_DEFUN([_LT_AC_TAGVAR])])
m4_ifndef([AC_LTDL_ENABLE_INSTALL],	[AC_DEFUN([AC_LTDL_ENABLE_INSTALL])])
m4_ifndef([AC_LTDL_PREOPEN],		[AC_DEFUN([AC_LTDL_PREOPEN])])
m4_ifndef([_LT_AC_SYS_COMPILER],	[AC_DEFUN([_LT_AC_SYS_COMPILER])])
m4_ifndef([_LT_AC_LOCK],		[AC_DEFUN([_LT_AC_LOCK])])
m4_ifndef([AC_LIBTOOL_SYS_OLD_ARCHIVE],	[AC_DEFUN([AC_LIBTOOL_SYS_OLD_ARCHIVE])])
m4_ifndef([_LT_AC_TRY_DLOPEN_SELF],	[AC_DEFUN([_LT_AC_TRY_DLOPEN_SELF])])
m4_ifndef([AC_LIBTOOL_PROG_CC_C_O],	[AC_DEFUN([AC_LIBTOOL_PROG_CC_C_O])])
m4_ifndef([AC_LIBTOOL_SYS_HARD_LINK_LOCKS], [AC_DEFUN([AC_LIBTOOL_SYS_HARD_LINK_LOCKS])])
m4_ifndef([AC_LIBTOOL_OBJDIR],		[AC_DEFUN([AC_LIBTOOL_OBJDIR])])
m4_ifndef([AC_LTDL_OBJDIR],		[AC_DEFUN([AC_LTDL_OBJDIR])])
m4_ifndef([AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH], [AC_DEFUN([AC_LIBTOOL_PROG_LD_HARDCODE_LIBPATH])])
m4_ifndef([AC_LIBTOOL_SYS_LIB_STRIP],	[AC_DEFUN([AC_LIBTOOL_SYS_LIB_STRIP])])
m4_ifndef([AC_PATH_MAGIC],		[AC_DEFUN([AC_PATH_MAGIC])])
m4_ifndef([AC_PROG_LD_GNU],		[AC_DEFUN([AC_PROG_LD_GNU])])
m4_ifndef([AC_PROG_LD_RELOAD_FLAG],	[AC_DEFUN([AC_PROG_LD_RELOAD_FLAG])])
m4_ifndef([AC_DEPLIBS_CHECK_METHOD],	[AC_DEFUN([AC_DEPLIBS_CHECK_METHOD])])
m4_ifndef([AC_LIBTOOL_PROG_COMPILER_NO_RTTI], [AC_DEFUN([AC_LIBTOOL_PROG_COMPILER_NO_RTTI])])
m4_ifndef([AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE], [AC_DEFUN([AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE])])
m4_ifndef([AC_LIBTOOL_PROG_COMPILER_PIC], [AC_DEFUN([AC_LIBTOOL_PROG_COMPILER_PIC])])
m4_ifndef([AC_LIBTOOL_PROG_LD_SHLIBS],	[AC_DEFUN([AC_LIBTOOL_PROG_LD_SHLIBS])])
m4_ifndef([AC_LIBTOOL_POSTDEP_PREDEP],	[AC_DEFUN([AC_LIBTOOL_POSTDEP_PREDEP])])
m4_ifndef([LT_AC_PROG_EGREP],		[AC_DEFUN([LT_AC_PROG_EGREP])])
m4_ifndef([LT_AC_PROG_SED],		[AC_DEFUN([LT_AC_PROG_SED])])
m4_ifndef([_LT_CC_BASENAME],		[AC_DEFUN([_LT_CC_BASENAME])])
m4_ifndef([_LT_COMPILER_BOILERPLATE],	[AC_DEFUN([_LT_COMPILER_BOILERPLATE])])
m4_ifndef([_LT_LINKER_BOILERPLATE],	[AC_DEFUN([_LT_LINKER_BOILERPLATE])])
m4_ifndef([_AC_PROG_LIBTOOL],		[AC_DEFUN([_AC_PROG_LIBTOOL])])
m4_ifndef([AC_LIBTOOL_SETUP],		[AC_DEFUN([AC_LIBTOOL_SETUP])])
m4_ifndef([_LT_AC_CHECK_DLFCN],		[AC_DEFUN([_LT_AC_CHECK_DLFCN])])
m4_ifndef([AC_LIBTOOL_SYS_DYNAMIC_LINKER],	[AC_DEFUN([AC_LIBTOOL_SYS_DYNAMIC_LINKER])])
m4_ifndef([_LT_AC_TAGCONFIG],		[AC_DEFUN([_LT_AC_TAGCONFIG])])
m4_ifndef([AC_DISABLE_FAST_INSTALL],	[AC_DEFUN([AC_DISABLE_FAST_INSTALL])])
m4_ifndef([_LT_AC_LANG_CXX],		[AC_DEFUN([_LT_AC_LANG_CXX])])
m4_ifndef([_LT_AC_LANG_F77],		[AC_DEFUN([_LT_AC_LANG_F77])])
m4_ifndef([_LT_AC_LANG_GCJ],		[AC_DEFUN([_LT_AC_LANG_GCJ])])
m4_ifndef([AC_LIBTOOL_LANG_C_CONFIG],	[AC_DEFUN([AC_LIBTOOL_LANG_C_CONFIG])])
m4_ifndef([_LT_AC_LANG_C_CONFIG],	[AC_DEFUN([_LT_AC_LANG_C_CONFIG])])
m4_ifndef([AC_LIBTOOL_LANG_CXX_CONFIG],	[AC_DEFUN([AC_LIBTOOL_LANG_CXX_CONFIG])])
m4_ifndef([_LT_AC_LANG_CXX_CONFIG],	[AC_DEFUN([_LT_AC_LANG_CXX_CONFIG])])
m4_ifndef([AC_LIBTOOL_LANG_F77_CONFIG],	[AC_DEFUN([AC_LIBTOOL_LANG_F77_CONFIG])])
m4_ifndef([_LT_AC_LANG_F77_CONFIG],	[AC_DEFUN([_LT_AC_LANG_F77_CONFIG])])
m4_ifndef([AC_LIBTOOL_LANG_GCJ_CONFIG],	[AC_DEFUN([AC_LIBTOOL_LANG_GCJ_CONFIG])])
m4_ifndef([_LT_AC_LANG_GCJ_CONFIG],	[AC_DEFUN([_LT_AC_LANG_GCJ_CONFIG])])
m4_ifndef([AC_LIBTOOL_LANG_RC_CONFIG],	[AC_DEFUN([AC_LIBTOOL_LANG_RC_CONFIG])])
m4_ifndef([_LT_AC_LANG_RC_CONFIG],	[AC_DEFUN([_LT_AC_LANG_RC_CONFIG])])
m4_ifndef([AC_LIBTOOL_CONFIG],		[AC_DEFUN([AC_LIBTOOL_CONFIG])])
m4_ifndef([_LT_AC_FILE_LTDLL_C],	[AC_DEFUN([_LT_AC_FILE_LTDLL_C])])
m4_ifndef([_LT_REQUIRED_DARWIN_CHECKS],	[AC_DEFUN([_LT_REQUIRED_DARWIN_CHECKS])])
m4_ifndef([_LT_AC_PROG_CXXCPP],		[AC_DEFUN([_LT_AC_PROG_CXXCPP])])
m4_ifndef([_LT_PREPARE_SED_QUOTE_VARS],	[AC_DEFUN([_LT_PREPARE_SED_QUOTE_VARS])])
m4_ifndef([_LT_PROG_ECHO_BACKSLASH],	[AC_DEFUN([_LT_PROG_ECHO_BACKSLASH])])
m4_ifndef([_LT_PROG_F77],		[AC_DEFUN([_LT_PROG_F77])])
m4_ifndef([_LT_PROG_FC],		[AC_DEFUN([_LT_PROG_FC])])
m4_ifndef([_LT_PROG_CXX],		[AC_DEFUN([_LT_PROG_CXX])])

