dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl edited by Jim Willenbring <jmwille@sandia.gov> to check for sgecon
dnl rather than cheev because by default (as of 8-13-2002) Trilinos
dnl does not build the complex portions of the lapack library.  Edited
dnl again on 5-13-2004 to check for dgecon instead of sgecon.
dnl Edited by Jim Willenbring on 4-17-2006 to stop looking for LAPACK if 
dnl a specific LAPACK library specified by a user cannot be used.

AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(dgecon)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $dgecon in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($dgecon, [acx_lapack_ok=yes], [user_spec_lapack_failed=yes])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# If the user specified a LAPACK library that could not be used we will
# halt the search process rather than risk finding a LAPACK library that
# the user did not specify.

if test "x$user_spec_lapack_failed" != xyes; then

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($dgecon, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($lapack, $dgecon,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

fi # If the user specified library wasn't found, we skipped the remaining
   # checks.

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK
