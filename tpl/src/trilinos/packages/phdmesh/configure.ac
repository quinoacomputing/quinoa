# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
#       phdMesh : Parallel Heterogneous Dynamic unstructured Mesh
#                 Copyright (2007) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------



# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(phdmesh, 1.1d, hcedwar@sandia.gov)

echo "----------------------------------------"
echo "Running phdMesh Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([src/util/ParallelComm.cpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)
# Configure should create src/phdmesh_config.h from src/phdmesh_config.h.in

AM_CONFIG_HEADER(src/phdmesh_config.h:src/phdmesh_config.h.in)

# Allow users to specify their own "install" command.  If none is specified,
# the default is install-sh found in the config subdirectory.
                                                                                
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)
                                                                                
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake

#  - Required version of automake.
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI
# phdmesh has a custom macro for MPI
if test "X${HAVE_PKG_MPI}" = "Xyes"; then
   AC_DEFINE(PHDMESH_HAS_MPI,,[define if we want to use MPI])
fi

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

# ------------------------------------------------------------------------
# Determine how to build C programs using POSIX
# ------------------------------------------------------------------------

ACX_PTHREAD(,AC_DEFINE(NO_PTHREADS,,[define if pthreads is not available]))

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------
# This needs to be added for phdMesh

AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)

if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" != "Xyes"; then
  AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  New_Package requires at least one or the other.]])
fi

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

# Use ExodusII libraries conditionally for mesh IO
TAC_ARG_ENABLE_FEATURE_SUB2( phdmesh, exodus, [Enable ExodusII support for phdmesh mesh IO], PHDMESH_HAS_SNL_EXODUSII, no)
#AM_CONDITIONAL(HAVE_THREADS, [test "X$ac_cv_use_pdmesh_has_threads" != "Xno"])

# ------------------------------------------------------------------------
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( phdmesh, tests, [Make phdMesh tests buildable with 'make tests'], PHDMESH_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_phdmesh_tests" != "Xno")

#TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
#TAC_ARG_ENABLE_FEATURE_SUB_CHECK( new_package, examples, [Make New_Package examples buildable with 'make examples'], NEW_PACKAGE_EXAMPLES)
#AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_new_package_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
#AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

TAC_ARG_ENABLE_FEATURE(libcheck, [Check for some third-party libraries.  (Cannot be disabled unless tests and examples are also disabled.)], LIBCHECK, yes)

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

# ------------------------------------------------------------------------
# Checks for libraries
# ------------------------------------------------------------------------

# If tests, examples and libcheck are disabled, we don't have to check
# for these libraries.

# #np# -
# If a package does not have tests or examples, the corresponding check(s)
# should be pulled out of the "if" statement below.
#if test "X$ac_cv_use_phdmesh_tests" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then
#if test "X$ac_cv_use_new_package_tests" != "Xno" || test "X$ac_cv_use_new_package_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

#fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
AC_CONFIG_FILES([
		Makefile
		Makefile.export.phdmesh
		src/Makefile
		tests/Makefile
		tests/gears/Makefile
		tests/txblas/Makefile
                tests/util/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "---------------------------------------------"
echo "Finished Running phdMesh Configure Script"
echo "---------------------------------------------"
