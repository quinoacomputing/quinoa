######################################################################################
# Extract the value of a variable from a makefile fragment.
# Arguments:
#   - The path of the makefile fragment
#   - The name of the makefile variable
#   - Action on success (the shell variable "make_val" will contain the value
#         of the variable)
#   - Action on failure
#######################################################################################
AC_DEFUN([SNL_MAKE_INC_VAR], [
make_val=
snl_makefile="snl_check.mak"
snl_resultfile="snl_check.out"
rm -f $snl_makefile

if test ! -f $1 ; then
  AC_MSG_WARN([File not found: $1])
  $4
else
cat >$snl_makefile <<SNL_END_OF_MAKEFILE
default:
	@echo "\$($2)" > $snl_resultfile

include $1
SNL_END_OF_MAKEFILE
if make -f $snl_makefile > /dev/null 2>&1; then
  make_val=`cat $snl_resultfile`
  rm -f $snl_makefile
  rm -f $snl_resultfile
  $3
else
  rm -f $snl_makefile
  rm -f $snl_resultfile
  $4
fi
fi
])

##############################################################################
# Optinally configure support for an ITAPS API
#
# Arguments:
#   1) The interface name in lowercase (e.g. imesh)
#   2) The interface name in uppercase (e.g. IMESH)
#   3) The interface name in preferred case (e.g. iMesh)
#   4) Interface required for this interface (all caps)
#   5) Second required interface (all caps)
#
# Sets and declares substitutions for:
#    ENABLE_IFACE=yes/no
#    IFACE_DEFS=empty/path to defs file
#    IFACE_LIBS=empty/ld args for linking implementation
#    IFACE_INCL=cpp flag for include path for API header
#    ENABLE_ITAPS=yes
#    IBASE_INCL=cpp flag for locating iBase.h
# Sets AUTOMAKE conditionals for ENABLE_IFACE 
##############################################################################

AC_DEFUN([ITAPS_API], [

AC_ARG_VAR([$2_DEFS],[Path to $3-Defs.inc])

# If no value is set, default to 'disabled'
if test "x${$2_DEFS}" = "x"; then
  $2_DEFS=no
fi

AC_ARG_WITH( [$1], 
[AC_HELP_STRING([--with-$1=<path>],
  [Path to $3-Defs.inc (enable $3)])],
  [test "x${$2_DEFS}" = "xno" && $2_DEFS=$withval])

# If option is disabled
if test "x${$2_DEFS}" = "xno"; then
  $2_DEFS=
  ENABLE_$2=no
  $2_LIBS=
  $2_INCL=
elif test "x${$2_DEFS}" = "xyes"; then
  AC_MSG_ERROR([Path to $3-Defs.inc not specified.  Use with-val or $2_DEFS])
else
  # Verify that we have any required interfaces
  if test "x" != "x$4"; then
    if test "x$ENABLE_$4" != "xyes"; then
      AC_MSG_ERROR([Cannot enable $3 w/out $4])
    fi
    if test "x" != "x$5"; then
      if test "x$ENABLE_$5" != "xyes"; then
        AC_MSG_ERROR([Cannot enable $3 w/out $5])
      fi
    fi
  fi

  # Verify that Defs file exists
  AC_CHECK_FILE( [${$2_DEFS}], [], [AC_MSG_ERROR([Invalid --with-$1 option or $2_DEFS value])] )

  # Extract values from Makefile stub
  SNL_MAKE_INC_VAR( [${$2_DEFS}], [$2_LIBS], [$2_LIBS="$make_val"] )
  SNL_MAKE_INC_VAR( [${$2_DEFS}], [$2_INCLUDES], [$2_INCL="$make_val"] )

  # Set up testing environment
  AC_LANG_PUSH([C++])
  old_LIBS="$LIBS"
  LIBS="$LIBS ${$2_LIBS}"
  ALL_INCL="${$2_INCL}"
  if test "x" != "x$4"; then
    ALL_INCL="$ALL_INCL ${$4_INCL}"
    LIBS="$LIBS ${$4_LIBS}"
    if test "x" != "x$5"; then
      ALL_INCL="$ALL_INCL ${$5_INCL}"
      LIBS="$LIBS ${$5_LIBS}"
    fi
  fi
  old_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS ${ALL_INCL}"
  
  # If common header hasn't been located yet, check for it
  if test "x$IBASE_INCL" == "x"; then
    ENABLE_ITAPS=yes
    AC_CHECK_HEADER([iBase.h],[IBASE_INCL="$ALL_INCL"])
  fi
  
  # Check for existence of interface definition header
  AC_CHECK_HEADER([$3.h],[],[AC_MSG_ERROR([Broken $3: $3.h not found with ${$2_INCL} from ${$2_DEFS}])])
  # Check library.  Note: we cannot use AC_CHECK_LIB because
  # of Fortran name-mangling crap in ITAPS headers
  AC_MSG_CHECKING([for $3_dtor in $3 library])
  AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([#include <$3.h>],[$3_dtor(($3_Instance)0,(int*)0)])],
   [AC_MSG_RESULT([yes])],
   [AC_MSG_RESULT([no])
    AC_MSG_CHECKING([for $3_destroy in $3 library])
    AC_LINK_IFELSE(
     [AC_LANG_PROGRAM([#include <$3.h>],[$3_destroy(($3_Instance)0,(int*)0)])],
     [AC_MSG_RESULT([yes])],
     [AC_MSG_RESULT([no])
      AC_MSG_CHECKING([for $3_destroyPartitionAll in $3 library])
      AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([#include <$3.h>],[$3_destroyPartitionAll((iMesh_Instance)0,($3_PartitionHandle)0,(int*)0)])],
       [AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT([no])
        AC_MSG_ERROR([$3 library unsuable: $2_LIBS=\"${$2_LIBS}\" from ${$2_DEFS}])]
       ]
      )
    )
   ]
  )
   
   # Resore environment
   LIBS="$old_LIBS"
   CPPFLAGS="$old_CPPFLAGS"
   AC_LANG_POP([C++])
   
   ENABLE_$2=yes
   ENABLE_IBASE=yes
fi

AC_SUBST(ENABLE_$2)
AC_SUBST($2_DEFS)
AC_SUBST($2_LIBS)
AC_SUBST($2_INCL)
AM_CONDITIONAL([ENABLE_$2],[test "xyes" = "x$ENABLE_$2"])

])

