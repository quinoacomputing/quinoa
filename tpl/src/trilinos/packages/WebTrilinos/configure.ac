# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
#                WebTrilinos: A Web Interface to Trilinos
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

#  #np# - package name, version number, and e-mail address below
AC_INIT(WebTrilinos, 1.3d, trilinos-help@software.sandia.gov)

echo "----------------------------------------"
echo "Running WebTrilinos Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR([common/index.html])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

AC_CONFIG_AUX_DIR(config)

AM_CONFIG_HEADER(common/WebTrilinos_config.h:common/WebTrilinos_config.h.in)

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

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
AC_PROG_F77(f77 g77 f90 xlf90)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
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
# Checks if tests and examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( webtrilinos, tests, [Make WebTrilinos tests buildable with 'make tests'], WEBTRILINOS_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_webtrilinos" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( webtrilinos, examples, [Make WebTrilinos examples buildable with 'make examples'], WEBTRILINOS_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_webtrilinos_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

# ------------------------------------------------------------------------
# Enable optional python wrappers
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_FEATURE_SUB( new_package, python, [Build New_Package Python wrappers (requires --enable-python, which is NOT set by default)], NEW_PACKAGE_PYTHON, yes)
AM_CONDITIONAL(BUILD_NEW_PACKAGE_PYTHON, test "X$ac_cv_use_new_package_python" != "Xno")
TAC_ENABLE_PYTHON

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

# find or choose to supply valid blas and lapack libraries
TAC_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
TAC_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex.  --with-lapack=/usr/path/lib/libcxml.a)]])])

# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.
# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS

# ------------------------------------------------------------------------
# Checks for system services
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------
AC_CONFIG_FILES([
		Makefile
                Makefile.export
                c++/Makefile
                c++/epetra/Makefile
                c++/amesos/Makefile
                c++/aztecoo/Makefile
                c++/ifpack/Makefile
                c++/ml/Makefile
                python/Makefile
                test/Makefile
		])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running WebTrilinos Configure Script"
echo "----------------------------------------"
