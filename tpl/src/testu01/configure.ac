#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# AC_PREREQ(2.57)
AC_INIT(TestU01, 1.2.3, lecuyer@iro.umontreal.ca)
AC_CONFIG_SRCDIR([mylib/num.c])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([include/config.h])
AC_CONFIG_HEADER([include/gdefconf.h])

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_CHECK_PROG(have_mathematica,math, yes, no) 
if test "x$have_mathematica" = "xyes" ; then
   AC_DEFINE([HAVE_MATHEMATICA], 1, [Define when Mathematica is installed])
fi

# Checks for libraries.
AC_CHECK_LIBM

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/utsname.h unistd.h windows.h gmp.h stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_BIGENDIAN
AC_CHECK_TYPES([long long, uint32_t, uint8_t])

# Checks for library functions.
# AC_FUNC_MALLOC
# AC_FUNC_REALLOC
LIBS="$LIBM $LIBS"
AC_CHECK_FUNCS([random erf lgamma log1p])

AC_CONFIG_FILES([Makefile
                 include/Makefile mylib/Makefile probdist/Makefile
                 testu01/Makefile examples/Makefile doc/Makefile param/Makefile])
AC_SUBST([LIBM])
AC_OUTPUT
