dnl @synopsis TAC_ARG_ENABLE_FEATURE_SUB2(FEATURE_NAME, SUB_FEATURE_NAME, FEATURE_DESCRIPTION, NAME, DEFAULT_VAL)
dnl
dnl This hack gets around the fact that TAC_ARG_ENABLE_FEATURE does not support underscores
dnl in its feature names.  TAC_ARG_ENABLE_FEATURE_SUB2 allows exactly one underscore.  Not great,
dnl but arguably better than supporting no underscores.
dnl
dnl TAC_ARG_ENABLE_FEATURE(feature-sub, [Configure and build feature-sub], FEATURE_SUB2, yes) 
dnl   fails because tac_arg_enable_feature tests for ac_cv_use_feature-sub which gets 
dnl   rejected because the `-' is not allowed in variables.  (AC_ARG_ENABLE sets ac_cv_use_feature_sub
dnl   to avoid this problem.)  Use:
dnl 
dnl TAC_ARG_ENABLE_FEATURE(feature, sub, [Configure and build feature-sub], FEATURE_SUB2, yes) 
dnl   instead.
dnl
dnl Test for --enable-${FEATURE_NAME} and set to DEFAULT_VAL value if feature not specified.
dnl Also calls AC_DEFINE to define ${NAME} if value is not equal to "no"
dnl 
dnl Use this macro to help defining whether or not optional 
dnl features* should compiled.  For example:
dnl
dnl TAC_ARG_ENABLE_FEATURE(epetra, [Configure and build epetra], EPETRA, yes)
dnl 
dnl will test for --enable-epetra when configure is run.  If it is defined 
dnl and not set to "no" or not defined (default is "yes") then EPETRA will
dnl be defined, if --enable-epetra is defined to be "no", EPETRA will not
dnl be defined.
dnl
dnl *NOTE: epetra, aztecoo, komplex, ifpack, and other software found in
dnl subdirectories of Trilinos/packages are "packages" in their own right.
dnl However, these packages are also "features" of the larger package
dnl "Trilinos".  Therefore, when configuring from the Trilinos directory,
dnl it is appropriate to refer to these software packages as "features".
dnl
dnl The only difference between this macro and tac_arg_enable_feature_sub
dnl is that the define is not prefixed with "HAVE"
dnl
dnl This file was based on tac_arg_enable_feature_sub.m4 by Ken Stanley
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_ENABLE_FEATURE_SUB2],
[
AC_ARG_ENABLE([$1-$2],
AC_HELP_STRING([--enable-$1-$2],[$3 (default is [$5])]),
ac_cv_use_$1_$2=$enableval, ac_cv_use_$1_$2=$5)

AC_MSG_CHECKING(whether to use [$1-$2])

if test "X$ac_cv_use_$1_$2" != "Xno"; then
  AC_MSG_RESULT(yes)  
  AC_DEFINE([$4],,[Define if want to build $1-$2])
else
  AC_MSG_RESULT(no)
fi
])

