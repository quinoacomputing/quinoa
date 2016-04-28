AC_PREREQ([2.64])
AC_INIT([libaec], [0.3.2], [rosenhauer@dkrz.de])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/libaec.h])
AC_CONFIG_HEADERS(config/config.h)

LT_INIT

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_C_INLINE
AC_C_RESTRICT

# Checks for library functions.
AC_CHECK_FUNCS([memset strstr])
AC_CHECK_DECLS(__builtin_clzll)

AC_CONFIG_FILES([Makefile         \
                 src/Makefile     \
                 tests/Makefile])

AC_OUTPUT
