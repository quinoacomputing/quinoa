#################################################################################
# Check if the compiler defines a specific preprocessor macro
# Arguments:
#  - preprocessor define to check for
#  - action upon success
#  - action upon failure
#################################################################################
AC_DEFUN([SNL_TRY_COMPILER_DEFINE], [
AC_COMPILE_IFELSE([
AC_LANG_PROGRAM( [[#ifndef $1
  choke me
#endif]], []) ],
[$2],[$3])
])

#######################################################################################
# Check for compiler-specific flags.
# Sets the following environmental variables:
#   SNL_CXX_SPECIAL : Any compiler-specific flags which must be specified
#   SNL_CXX_32BIT   : Flag to force compilation of 32-bit code
#   SNL_CXX_64BIT   : Flag to force compilation of 64-bit code
#   SNL_CXX_DEBUG   : Debug flag
#   SNL_CXX_OPTIMIZE: Flag for default optimization level
#   SNL_CXX_FAST    : Fastest optimization, produces machine-specific code
#######################################################################################
AC_DEFUN([SNL_DETECT_CXX], [
AC_REQUIRE([AC_PROG_CXX])

# Detect compiler 
AC_MSG_CHECKING([for known c++ compilers])
# Autoconf does G++ for us
if test x$GXX = xyes; then
  cxx_compiler=GNU
# Search for other compiler types
# For efficiency, limit checks to relevant OSs
else
  cxx_compiler=unknown
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  case "$target_os" in
    aix*)
      SNL_TRY_COMPILER_DEFINE([__IBMCPP__],[cxx_compiler=VisualAge])
      ;;
    solaris*|sunos*)
      SNL_TRY_COMPILER_DEFINE([__SUNPRO_CC],[cxx_compiler=SunWorkshop])
      ;;
    irix*)
      SNL_TRY_COMPILER_DEFINE([__sgi],[cxx_compiler=MIPSpro])
      ;;
    linux*)
      SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
      SNL_TRY_COMPILER_DEFINE([__IBMCPP__],[cxx_compiler=VisualAge])
      SNL_TRY_COMPILER_DEFINE([__DECCXX_VER],[cxx_compiler=Compaq])
      ;;
    hpux*)
      SNL_TRY_COMPILER_DEFINE([__HP_aCC],[cxx_compiler=HP])
      ;;
    windows*)
      SNL_TRY_COMPILER_DEFINE([__MSC_VER],[cxx_compiler=VisualStudio])
      SNL_TRY_COMPILER_DEFINE([__INTEL_COMPILER],[cxx_compiler=Intel])
      SNL_TRY_COMPILER_DEFINE([__DECCXX_VER],[cxx_compiler=Compaq])
      SNL_TRY_COMPILER_DEFINE([__BORLANDC__],[cxx_compiler=Borland])
      SNL_TRY_COMPILER_DEFINE([__CYGWIN__],[cxx_compiler=Cygwin])
      SNL_TRY_COMPILER_DEFINE([__MINGW32__],[cxx_compiler=MinGW])
      ;;
  esac
  AC_LANG_RESTORE
fi
AC_MSG_RESULT([$cxx_compiler])
if test "x$cxx_compiler" = "xunknown"; then
  AC_MSG_WARN([Unrecognized C++ compiler: $CXX])
fi

# Now decide special compiler flags using compiler/cpu combination
AC_MSG_CHECKING([for known compiler/OS combinations])
case "$cxx_compiler:$host_cpu" in
  GNU:sparc*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="-Wall -pipe"
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE=-O3
    SNL_CXX_FAST='-O3 -ffast-math'
    ;;
  GNU:powerpc*)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="-Wall -pipe"
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE=-O3
    SNL_CXX_FAST='-O3 -ffast-math'
    ;;
  GNU:i?86)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="-Wall -pipe"
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE="-O3 -mtune=generic -march=pentium3"
    SNL_CXX_FAST='-O3 -ffast-math -march=native -mfpmath=sse'
    ;;
  GNU:x86_64)
    SNL_CXX_32BIT=-m32
    SNL_CXX_64BIT=-m64
    SNL_CXX_SPECIAL="-Wall -pipe"
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE="-O3 -mtune=generic"
    SNL_CXX_FAST='-O3 -ffast-math -march=native -mfpmath=sse'
    ;;
  GNU:mips*)
    SNL_CXX_32BIT="-mips32 -mabi=32"
    SNL_CXX_64BIT="-mips64 -mabi=64"
    SNL_CXX_SPECIAL="-Wall -pipe"
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE=-O3
    SNL_CXX_FAST='-O3 -ffast-math'
    ;;
  VisualAge:*)
    SNL_CXX_32BIT=-q32
    SNL_CXX_64BIT=-q64
    # Do V5.0 namemangling for compatibility with ACIS, and enable RTTI
    SNL_CXX_SPECIAL="-qrtti=all -qnamemangling=v5"
    AR="ar -X 32_64"
    NM="nm -B -X 32_64"
    SNL_CXX_FAST="${SNL_CXX_OPTIMIZE}"
    ;;
  MIPSpro:mips)
    SNL_CXX_32BIT=-n32
    SNL_CXX_64BIT=-64
    SNL_CXX_SPECIAL=-LANG:std
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE=-O2
    SNL_CXX_FAST="${SNL_CXX_OPTIMIZE}"
    ;;
  MIPSpro:*)
    SNL_CXX_SPECIAL=-LANG:std
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE=-O2
    SNL_CXX_FAST="${SNL_CXX_OPTIMIZE}"
    ;;
  SunWorkshop:sparc*)
    SNL_CXX_32BIT=-xarch=generic
    SNL_CXX_64BIT=-xarch=generic64
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE=-O2
    SNL_CXX_FAST=-fast
    ;;
  *)
    SNL_CXX_DEBUG=-g
    SNL_CXX_OPTIMIZE=-O2
    SNL_CXX_FAST=-O2
    ;;
esac
AC_MSG_RESULT([$cxx_compiler:$host_cpu])
]) # end SNL_CXX_FLAGS

#######################################################################################
# Test if compiler accepts flags.
# Arguments : flag, success-action, fail-action
# Example:
#  SNL_ADD_CXX_FLAG( [-g], [AM_CXXFLAGS="$AM_CXXFLAGS -g"], [] )
#
#######################################################################################
AC_DEFUN([SNL_ADD_CXX_FLAG], [
  AC_REQUIRE([AC_PROG_CXX])
  AC_MSG_CHECKING([if compiler accepts flag: $1])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS $1"
  AC_TRY_COMPILE( [], [], [snl_tmp_success=yes], [snl_tmp_success=no] )
  CXXFLAGS="$OLD_CXXFLAGS"
  AC_LANG_RESTORE
  if test "x$snl_tmp_success" = "xyes"; then
    AC_MSG_RESULT([yes])
    $2
  else
    AC_MSG_RESULT([no])
    $3
  fi
] )

#######################################################################################
# Test a set of CXX flags.  
# Arguments:
#  1) name of shell variable containing list of CXX flags to test
#  2) variable to append valid CXX flags to
#  3) option name to use in warning/error output
#  4) "one" : input is one flag, possibly containng spaces, error if flag doesn't work
#     "yes" : error if no valid flags in list
#     "no"  : silent
#      *    : treat input as list of space-separated flags, warn if no valid flags
#######################################################################################
AC_DEFUN([SNL_ADD_CXX_FLAGS], [
  SNL_CXX_FLAGS_ERROR="$4"
  SNL_TMP_CXX_FLAGS=
  if test "xone" = "x$4"; then
    SNL_CXX_FLAGS_ERROR=yes
    SNL_ADD_CXX_FLAG( [${$1}], [SNL_TMP_CXX_FLAGS="$SNL_TMP_CXX_FLAGS ${$1}"] )
  else
    for flag in ${$1}; do
      SNL_ADD_CXX_FLAG( [$flag], [SNL_TMP_CXX_FLAGS="$SNL_TMP_CXX_FLAGS $flag"] )
    done
  fi
  if test "x" = "x$SNL_TMP_CXX_FLAGS"; then
    if test "xyes" = "x$4"; then
      AC_MSG_ERROR( [Don't know how to $3 on this platform.  Try specifying CXXFLAGS manually.] )
    elif test "xno" != "x$4"; then
      AC_MSG_WARN( [Don't know how to $3 on this platform.  Try specifying CXXFLAGS manually.] )
    fi
  else
    $2="${$2} $SNL_TMP_CXX_FLAGS"
  fi
])

AC_DEFUN([SNL_REQUIRED_CXX_FLAGS],[
dnl  AC_REQUIRE(SNL_DETECT_CXX)
  SNL_ADD_CXX_FLAGS([SNL_CXX_SPECIAL],[AM_CXXFLAGS],[],[no])
])
AC_DEFUN([SNL_CXX_DEBUG_SYMBOLS],[
dnl  AC_REQUIRE(SNL_DETECT_CXX)
  SNL_ADD_CXX_FLAGS([SNL_CXX_DEBUG],[CXXFLAGS],[enable debug symbols],[yes])
])
AC_DEFUN([SNL_CXX_COMPILE_OPTIMIZED],[
dnl  AC_REQUIRE(SNL_DETECT_CXX)
  SNL_ADD_CXX_FLAGS([SNL_CXX_OPTIMIZE],[CXXFLAGS],[enable compiler optimization])
])
AC_DEFUN([SNL_CXX_COMPILE_FAST],[
dnl  AC_REQUIRE(SNL_DETECT_CXX)
  SNL_ADD_CXX_FLAGS([SNL_CXX_FAST],[CXXFLAGS],[enable extra compiler optimization])
])
AC_DEFUN([SNL_CXX_COMPILE_32BIT], [
dnl  AC_REQUIRE(SNL_DETECT_CXX)
  SNL_ADD_CXX_FLAGS([SNL_CXX_32BIT],[CXXFLAGS],[force 32-bit build],[one])
])
AC_DEFUN([SNL_CXX_COMPILE_64BIT], [
dnl  AC_REQUIRE(SNL_DETECT_CXX)
  SNL_ADD_CXX_FLAGS([SNL_CXX_64BIT],[CXXFLAGS],[force 64-bit build],[one])
])

