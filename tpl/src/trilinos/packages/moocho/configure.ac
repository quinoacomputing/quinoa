# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# @HEADER
# ************************************************************************
# 
# Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
#                  Copyright (2003) Sandia Corporation
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
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ------------------------------------------------------------------------
# Initialization 
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

#  #np# - package name, version number, and e-mail address below
AC_INIT(moocho, Development Version, rabartl@sandia.gov)

echo "----------------------------------------"
echo "Running moocho Configure Script"
echo "----------------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

# #np# - file name below
AC_CONFIG_SRCDIR([src/MoochoUtilities/src/SpecialASCII.hpp])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.

#
#   I wimped out here.  I still think that ../../config can work 
#   for any package with a frim dependency on Trilinos.
#   But, using ../../config puts some ../../config files into a list
#   of files that it looked like it might delete under certain circumstances,
#   or perhaps put into some tar ball.  Anyway, it scared me off.
#   In order to get rid of the local config directory, you have to 
#   reference the config directory in Makefile.am and in bootstrap
#   (which should use aclocal -I ../../config)
#
AC_CONFIG_AUX_DIR(config)

AM_CONFIG_HEADER(src/Moocho_Config.h:src/Moocho_Config.h.in)

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

#  #np# - list required automake version below
AM_INIT_AUTOMAKE(1.10 no-define tar-ustar)

# Specify required version of autoconf.

AC_PREREQ(2.61)

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

#  #np# - can eliminate compiler checks below if your package does not use the
#         language corresponding to the check
# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ c++ cxx)
AC_PROG_F77(f77 g77 gfortran f90 xlf90 f95)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package 
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
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
# Checks for header files
# ------------------------------------------------------------------------

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------

# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for special package flags
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  moocho, ma28, MOOCHO_MA28, ENABLE_MOOCHO_MA28, yes,
  [Build in support for HSL solver MA28 into Moocho.],
  [Build in support for HSL solver MA28 into Moocho.  Can be overridden with --disable-moocho-ma28.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  moocho, thyra, MOOCHO_THYRA, ENABLE_MOOCHO_THYRA, yes,
  [Build the Thyra adapater support in Moocho.],
  [Build the Thyra adapater support in Moocho.  Can be overridden with --disable-moocho-thyra.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  moocho, stratimikos, MOOCHO_STRATIMIKOS, ENABLE_MOOCHO_STRATIMIKOS, yes,
  [Build the Stratimikos adapater support in Moocho.],
  [Build the Stratimikos adapater support in Moocho.  Can be overridden with --disable-moocho-stratimikos.])

TAC_ARG_ENABLE_CAN_USE_PACKAGE(
  moocho, rythmos, MOOCHO_RYTHMOS, ENABLE_MOOCHO_RYTHMOS, yes,
  [Build the Rythmos adapater support in Moocho.],
  [Build the Rythmos adapater support in Moocho.  Can be overridden with --disable-moocho-rythmos.])

TAC_ARG_INFORM_SUBPACKAGE_USED(
  thyra, epetraext, THYRA_EPETRAEXT, ENABLE_THYRA_EPETRAEXT, yes,
  [See thyra/configure --help],
  [See thyra/configure --help])

# ------------------------------------------------------------------------
# Checks if examples should be built
# ------------------------------------------------------------------------

TAC_ARG_ENABLE_FEATURE(tests, [Make tests for all Trilinos packages buildable with 'make tests'], TESTS, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( moocho, tests, [Make Moocho tests buildable with 'make tests'], MOOCHO_TESTS)
AM_CONDITIONAL(BUILD_TESTS, test "X$ac_cv_use_moocho_tests" != "Xno")

TAC_ARG_ENABLE_FEATURE(examples, [Make examples for all Trilinos packages buildable with 'make examples'], EXAMPLES, yes)
TAC_ARG_ENABLE_FEATURE_SUB_CHECK( moocho, examples, [Make Moocho examples buildable with 'make examples'], MOOCHO_EXAMPLES)
AM_CONDITIONAL(BUILD_EXAMPLES, test "X$ac_cv_use_moocho_examples" != "Xno")

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though SUB_TEST and SUB_EXAMPLE will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
AM_CONDITIONAL(SUB_TEST, test "X$ac_cv_use_sub_test" = "Xyes")
AM_CONDITIONAL(SUB_EXAMPLE, test "X$ac_cv_use_sub_example" = "Xyes")

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
if test "X$ac_cv_use_moocho_tests" != "Xno" || test "X$ac_cv_use_moocho_examples" != "Xno" || test "X$ac_cv_use_libcheck" != "Xno"; then

# #np# - 
# To add any libraries that you want to check for, replace pow below 
# with the name of an entry point in that library.
# checks for libraries now
dnl Replace `main' with a function in -lm:
#  #np# These lines can be removed if the math lib is not needed.
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

fi
# end of the list of libraries that don't need to be checked for if
# tests and examples are disabled.

# ------------------------------------------------------------------------
# Checks for compiler characteristics
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Checks for Makefile.export related systems
# ------------------------------------------------------------------------
TAC_ARG_ENABLE_EXPORT_MAKEFILES(yes)

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS

#  ----------------------------------------------------
#  Needs this because Dec handles templates differently
#  ----------------------------------------------------
#  #np# - Only uncomment the following lines if the package
#         creates a subdirectory called cxx_repository for
#         template object files on DEC.  (HINT: try to build on 
#         DEC/CPLANT and see if template-related errors occur)

# -- DEC CXX --
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
  case $target in
    alpha*)
      if test $CXX = cxx; then
        TRILINOS_TEMPLATE_OBJS="cxx_repository/*.o"
      fi
    ;;
  esac
fi

AC_MSG_CHECKING([for template-specific object files to be linked])
if test "X$TRILINOS_TEMPLATE_OBJS" = "X"; then
 AC_MSG_RESULT(none)
else
 AC_MSG_RESULT($TRILINOS_TEMPLATE_OBJS)
fi

AC_SUBST(TRILINOS_TEMPLATE_OBJS)

# ------------------------------------------------------------------------
# Checks for library functions
# ------------------------------------------------------------------------

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
    Makefile.export.moocho
    example/ExampleNLPBanded/Makefile
    example/ExampleNLPBanded/external/Makefile
    example/ExampleNLPBanded/external/Makefile.external.stub
    example/ExampleNLPDirect/Makefile
    example/ExampleNLPFirstOrder/Makefile
    example/Makefile
    example/NLPWBCounterExample/Makefile
    example/NLPWBCounterExample/external/Makefile
    example/NLPWBCounterExample/external/Makefile.external.stub
    example/generate-sample-output.pl.stub
    src/Makefile
    test/DenseLinAlgPack/Makefile
    test/IterationPack/Makefile
    test/Makefile
    test/MoochoUtilities/Makefile
    thyra/Makefile
    thyra/example/Makefile
    thyra/example/NLPThyraEpetraAdvDiffReactOpt/Makefile
    thyra/example/NLPThyraEpetraAdvDiffReactOpt/external/Makefile
    thyra/example/NLPThyraEpetraAdvDiffReactOpt/external/Makefile.external.stub
    thyra/example/NLPThyraEpetraModelEval4DOpt/Makefile
    thyra/example/NLPThyraEpetraModelEval4DOpt/_test_lumped_param_model.py.stub
    thyra/example/NLPThyraEpetraModelEval4DOpt/external/Makefile
    thyra/example/NLPThyraEpetraModelEval4DOpt/external/Makefile.external.stub
    thyra/example/DiagonalTransientInversion/Makefile
    thyra/example/DiagonalQuadraticResponseOnlyOpt/Makefile
    thyra/example/generate-sample-output.pl.stub
    thyra/src/Makefile
    ])

AC_OUTPUT()

# Bye World!
echo "----------------------------------------"
echo "Finished Running Moocho Configure Script"
echo "----------------------------------------"
