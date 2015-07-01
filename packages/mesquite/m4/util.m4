#
# Allow an required program to be specified on the command line.
# If it is not specified or no value is specified for it, attempt
# to detech it.  
# Usage:
#   MSQ_CHK_PROG_WITH( varname, progname, none-val, description, [extra_path] )
#     varname     - The name of the variable in which to store the result
#                   AC_SUBST is called for the specified variable name.
#     progname    - The name of the executable to be specified/found.
#     none-val    - Value to use for program if it cannot be found.
#                   If not specified, configure will abort if the program wasn't found.
#     description - Help string describing option
#     extra_path  - Additional directories to append to search path
#
AC_DEFUN([MSQ_CHK_PROG_WITH], [
  msq_tmp_path=$PATH
  if test -n "$5"; then
    msq_tmp_path="${PATH}:$5"
  fi
  
  AC_ARG_WITH([$2], [$4], [
      # If some relavent CL option was given
    if test "$withval" = "yes"; then
      AC_MSG_ERROR("--with-$2 requires an argument")
    elif test "$withval" = "no"; then
      if test -z "$3"; then
        AC_MSG_ERROR("Cannot build without $2")
      fi
      $1="$3"
      AC_SUBST($1)
    else
      AC_CHECK_FILE( [$withval], [], [AC_MSG_ERROR("$withval does not exist") ])
      $1=$withval
      AC_SUBST($1)
    fi 
    ], [
      # User didn't specify anything, try to detect
    AC_PATH_PROG( $1, $2, $3, [$msq_tmp_path] )
    if test -z "${$1}"; then
      AC_MSG_ERROR("Cannot build without $2")
    fi
  ])
])

