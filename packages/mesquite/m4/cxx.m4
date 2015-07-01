# Determine which macro to use for the function name.  Possiblities are:
# __PRETTY_FUNCTION__ - g++ - Decorated function name (namespace, class, return type, args, etc.)
# __FUNCDNAME__             - MS Visual Studio decorated function name
# __FUNCTION__ - g++, Visual Studio, others?
# __func__     - C99 standard (not C++), newer g++, (Sun Forte with -features=extensions)
# __function__ - ?
# __FUNC__     - SGI? 
# If __func__ does not work, redefine it to be whichever one does. If
# platform doesn't support any, define it to be an empty string.

AC_DEFUN([MSQ_CPLUSPLUS_FUNC], [
AC_MSG_CHECKING([for function-name preprocessor macro])
AC_LANG_PUSH(C++)
msq_cpp_func=

  # Loop through list of possibilities...
for i in __PRETTY_FUNCTION__ __FUNCDNAME__ __func__ __FUNCTION__ __function__ __FUNC__; do
  if test -z "$msq_cpp_func"; then
    AC_TRY_COMPILE(
      [#include <stdio.h>],
      [(void)printf("%s\n", $i);],
      [msq_cpp_func=$i]
    )
  fi
done

  # If none was found...
if test -z "$msq_cpp_func"; then
  AC_MSG_RESULT(none)
  AC_MSG_WARN([MSQ_FUNCTION will be defined as an empty string.]);
  AC_DEFINE(MSQ_FUNCTION, [""], [Define to c++ preprocessor macro for function name])
else 
  AC_MSG_RESULT($msq_cpp_func)
  AC_DEFINE_UNQUOTED(MSQ_FUNCTION, $msq_cpp_func, [Define to c++ preprocessor macro for function name])
fi

AC_LANG_POP(C++)
]) #MSQ_CPLUSPLUS_FUNC






