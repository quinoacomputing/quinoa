# Every configure script must call AC_INIT before doing anything else.
# AC_INIT (package, version, [bug-report], [tarname])
AC_INIT([H5Part], [1.6.6], [h5part@lists.psi.ch], H5Part)

# Ensure that a recent enough version of Autoconf is being used. 
# If the version of Autoconf being used to create configure is earlier than version, 
# print an error message to the standard error output and do not create configure.
AC_PREREQ(2.59)

# should be called right after AC_INIT.
# configure scripts can create a C header file containing `#define' directives. 
# The AC_CONFIG_HEADERS macro selects this kind of output.
AC_CONFIG_HEADERS(config.h)

AC_CONFIG_MACRO_DIR([m4])

# AM_INIT_AUTOMAKE is required to use autoconf with automake
AM_INIT_AUTOMAKE

BUILD_LIBS='libH5Part'
BUILD_TESTS='test'

###############################################################################
################# --enable-xxx and --with-xxx Argument ########################
###############################################################################

AC_ARG_ENABLE(
	[debug],
	[AC_HELP_STRING([--enable-debug],
	[Compile with debug flags [default=no]])],
	[USE_DEBUG=$enableval])

AC_ARG_ENABLE(
	[fortran],
	[AC_HELP_STRING([--enable-fortran],
	[Compile the Fortran interface [default=no]])],
	[USE_FORTRAN=$enableval])

AC_ARG_ENABLE(
	[parallel],
	[AC_HELP_STRING([--enable-parallel],
	[Compile the MPI/IO interface [default=no]])],
	[USE_PARALLEL=$enableval])

AC_ARG_ENABLE(
	[tools],
        [AC_HELP_STRING([--enable-tools],
	[Compile h5part tools [default=no]])],
        [USE_TOOLS=$enableval])

AC_ARG_WITH(
	[hdf5],
	[AC_HELP_STRING([--with-hdf5],
	[path to HDF5 installation [default=""]])],
	[HDF5PATH=$withval], [HDF5PATH=""])

###############################################################################
############### PATH SERACH FUNCTION - to be used later... ####################
###############################################################################
# /*@@
#   @routine    CCTK_Search
#   @date       Wed Jul 21 11:16:35 1999
#   @author     Tom Goodale
#   @desc
#   Used to search for something in various directories
#   @enddesc
#@@*/

PATH_Search() {
	eval  $1=""
	if test $# -lt 4 ; then
		h5part_basedir=""
	else
		h5part_basedir="$4/"
	fi
	for h5part_place in $2; do
		AC_MSG_CHECKING([looking in $h5part_place ... ])
		if test -r "$h5part_basedir$h5part_place/$3" ; then
			AC_MSG_RESULT([found])
			eval $1="$h5part_place"
			break
		fi
		AC_MSG_RESULT([no])
	done
	return
}


###############################################################################
############# MISC SETTINGS INCLUDING C & C++ COMPILER SETTING ################
###############################################################################
# Compute the canonical host-system type variable, host, and its three
# individual parts host_cpu, host_vendor, and host_os.
AC_CANONICAL_HOST

AC_PROG_MAKE_SET

# Determine a C/C++ compiler to use. 
# If CC is not already set in the environment, check for gcc and cc, then
# for other C compilers. 
# Set output variable CC to the name of the compiler found.
if test "x$USE_PARALLEL" = "xyes"; then
   	CCOMPILERS="mpicc cc"
   	CXXCOMPILERS="mpic++ mpicxx CC"
else
	CCOMPILERS="pgcc pathcc icc gcc cc_r cc"
   	CXXCOMPILERS="pgCC pathCC icpc g++"
fi

AC_PROG_CC($CCOMPILERS)

# Only need C++ for some of the tools
#if test "x$USE_TOOLS" = "xyes"; then
	AC_PROG_CXX($CXXCOMPILERS)
#fi

# Set output variable INSTALL to the path of a BSD-compatible install program,
# if one is found in the current PATH. 
# Otherwise, set INSTALL to `dir/install-sh -c`
AC_PROG_INSTALL

AC_PROG_AWK

# Default prefix for bindir, etc... (eg >> ./build/bin)
AC_PREFIX_DEFAULT(`pwd`/build)

# AC_DEFINE_UNQUOTED (variable, value, [description])
# Define the C preprocessor variable variable to value
# Use this macro instead of AC_DEFINE when variable or value is a shell variable. 
AC_DEFINE_UNQUOTED(MY_BUILD_VENDOR, "$host_vendor", "")
AC_DEFINE_UNQUOTED(MY_BUILD_CPU, "$host_cpu", "")
AC_DEFINE_UNQUOTED(MY_BUILD_OS, "$host_os", "")
AC_DEFINE_UNQUOTED(MY_GNUNAME, "${host_cpu}-${host_vendor}-${host_os}", "")
AC_DEFINE_UNQUOTED(MY_UNAME, "$uname", "")


###############################################################################
######################## CONFIGURE LINE OPTIONS ###############################
###############################################################################

AC_MSG_CHECKING([if debug is enabled])

if test "X$USE_DEBUG" = "Xyes"; then
    AC_MSG_RESULT([yes])
    CFLAGS="$CFLAGS -g"
    CXXFLAGS="$CXXFLAGS -g"
    FFLAGS="$FFLAGS -g"
else
    AC_MSG_RESULT([no])
fi

############################ fortran enabled ##################################
AC_MSG_CHECKING([if fortran interface enabled])
if test "X$USE_FORTRAN" = "Xyes"; then
	AC_MSG_RESULT([yes])
	LIB_FORTRAN="libH5PartF.la"
        BUILD_LIBS="$BUILD_LIBS libH5PartF"
	BUILD_TESTS="$BUILD_TESTS testf"

	if test "X$USE_PARALLEL" = "Xyes"; then
		AC_PROG_FC(mpif90 mpif77 ftn)
	else
		AC_PROG_FC(pgf90 pathf90 ifort ftn xlf_r g95 g90 gfortran)
	fi

	if test -z "$FC" ; then
		AC_MSG_ERROR([Cannot find a Fortran compiler!])
		exit 1
	fi

	if test $FC = "g90"; then
		FFLAGS="${FFLAGS} -fno-second-underscore"
	fi

	if test $FC = "g95"; then
		FFLAGS="${FFLAGS} -fno-second-underscore"
	fi

	AC_MSG_CHECKING([symbol convention in object files])
	`cd src && rm -f TestUnderscore.o TestUnderscoreC.o TestUnderscore`
	`cd src && ${FC} ${FFLAGS} -c TestUnderscore.f`
	`cd src && ${CC} ${CFLAGS} -c TestUnderscoreC.c`
	`cd src && ${FC} ${FFLAGS} -o TestUnderscore TestUnderscore.o TestUnderscoreC.o`

	if test -f src/TestUnderscore ; then
		UNDERSCORE_H=Underscore.h
		`cd src && ./TestUnderscore > Underscore.h`
		AC_MSG_RESULT([ok])
	else
		AC_MSG_RESULT([nok])
		AC_MSG_ERROR([Cannot determine the symbon convention for Fortran object files!])
		exit 1
	fi

else
	AC_MSG_RESULT([no])
fi


# Disable shared libraries by default: can be enabled with --enable-shared
LT_INIT([disable-shared])
AC_PROG_LIBTOOL


######################## parallel interface enabled ###########################
AC_MSG_CHECKING([if parallel interface enabled])
if test "X$USE_PARALLEL" = "Xyes"; then
	AC_MSG_RESULT([yes])

	CFLAGS="${CFLAGS} -DPARALLEL_IO -DMPICH_IGNORE_CXX_SEEK"
	CXXFLAGS="${CXXFLAGS} -DPARALLEL_IO -DMPICH_IGNORE_CXX_SEEK"
	FFLAGS="${FFLAGS} -DPARALLEL_IO"

	AC_MSG_CHECKING([if we can compile MPI code without setting flags])
	AC_TRY_LINK([#include "mpi.h"], [
		MPI_Comm comm;
		int n;
		MPI_Comm_size( comm, &n ); ], 
		[AC_MSG_RESULT([yes]); r='yes'], [AC_MSG_RESULT([no]); r='no'] )

	if test "X$r" = "Xno"; then
		AC_MSG_ERROR([MPI wrapper can't compile or link MPI program! Please set the INCLUDE and LIBS variables manually.])
		exit 1
	fi

else  # --enable-parallel=no
	AC_MSG_RESULT([no])
fi


######################## tools enabled ###########################
AC_MSG_CHECKING([whether tools are enabled])
if test "X$USE_TOOLS" = "Xyes"; then
    AC_MSG_RESULT([yes])
    BUILD_TOOLS="h5pAttrib h5pToGNUplot homdynToH5p"
    # tools + parallel
    if test "X$USE_PARALLEL" = "Xyes"; then
        BUILD_TOOLS="${BUILD_TOOLS} H5PartBench H5BlockBench"
    fi
else
    AC_MSG_RESULT([no])
fi


###############################################################################
######################### PATH CHECKING & SETTING #############################
###############################################################################

AC_MSG_CHECKING([for HDF5 root ])
AC_MSG_RESULT([])
if test -n "${HDF5PATH}" ; then
	P=${HDF5PATH}
elif test -n "${HDF5ROOT}"; then
	P=${HDF5ROOT}
elif test -n "${HDF5HOME}" ; then
	P=${HDF5HOME}
elif test -n "${HDF5_DIR}" ; then
	P=${HDF5_DIR}
else
	P=''
	P="$P /usr"
	P="$P /usr/local"
	P="$P /usr/local/hdf5"
	P="$P /usr/local/packages/hdf5"
	P="$P /apps/hdf5"
	P="$P /opt/hdf5"
fi
PATH_Search HDF5ROOT "$P" include/hdf5.h
if test -z "$HDF5ROOT"; then
	AC_MSG_ERROR([Cannot find an HDF5 library!])
	exit 1
fi

INCLUDES="$INCLUDES -I$HDF5ROOT/include"
LDFLAGS="$LDFLAGS -L$HDF5ROOT/lib"
LIBS="$LIBS -lhdf5"

AC_MSG_CHECKING([if we need to link to libsz ])
if test -n "$HDF5ROOT"; then
	if test -f $HDF5ROOT/lib/libsz.a; then
		AC_MSG_RESULT([yes])
		LDFLAGS="$LDFLAGS -L$HDF5ROOT/lib"
                LIBS="$LIBS -lsz"
	else
		AC_MSG_RESULT([no])
	fi
fi

LIBS="$LIBS -lz -lm"


# Checks for header files.
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow strchr strdup])

###############################################################################
############## EXPORTING VARIABLES & CREATING OUTPUT FILES ####################
###############################################################################
# AC_SUBST (variable, [value])
# Create an output variable from a shell variable.  Make AC_OUTPUT substitute
# the variable variable into output files (typically one or more `Makefile's).
# This means that AC_OUTPUT will replace instances of `@variable@' in input
# files with the value that the shell variable variable has when AC_OUTPUT is
# called.  This value of variable should not contain literal newlines.  If
# value is given, in addition assign it to variable. 
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(FFLAGS)
AC_SUBST(INCLUDES)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)
AC_SUBST(UNDERSCORE_H)
AC_SUBST(LIB_FORTRAN)
AC_SUBST(BUILD_TESTS)
AC_SUBST(BUILD_TOOLS)

# Make AC_OUTPUT create each `file' by copying an input file (by default `file.in'),
# substituting the output variable values.
AC_CONFIG_FILES([
	Makefile
	src/Makefile
	doc/Makefile
	test/Makefile
	tools/Makefile
])

AC_OUTPUT

###############################################################################
########################## PRINTING SUMMARY ###################################
###############################################################################
AC_MSG_RESULT([ ])
AC_MSG_RESULT([Summary:])
AC_MSG_RESULT([ ])
AC_MSG_RESULT([Host OS:             $host_os])
AC_MSG_RESULT([Host CPU:            $host_cpu])
AC_MSG_RESULT([Host vendor:         $host_vendor])
AC_MSG_RESULT([Build libraries:     $BUILD_LIBS])
AC_MSG_RESULT([Build test programs: $BUILD_TESTS])
AC_MSG_RESULT([Build tools:         $BUILD_TOOLS])
AC_MSG_RESULT([CC =                 $CC])
AC_MSG_RESULT([CXX =                $CXX])
AC_MSG_RESULT([FC =                 $FC])
AC_MSG_RESULT([CFLAGS =             $CFLAGS])
AC_MSG_RESULT([CXXFLAGS =           $CXXFLAGS])
AC_MSG_RESULT([FFLAGS =             $FFLAGS])
AC_MSG_RESULT([INCLUDES =           $INCLUDES])
AC_MSG_RESULT([LDFLAGS =            $LDFLAGS])
AC_MSG_RESULT([LIBS =               $LIBS])
AC_MSG_RESULT([HDF5ROOT =           $HDF5ROOT])
AC_MSG_RESULT([ ])
