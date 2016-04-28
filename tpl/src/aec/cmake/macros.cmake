MACRO(CHECK_CLZLL VARIABLE)
  CHECK_C_SOURCE_COMPILES(
    "int main(int argc, char *argv[])
{return __builtin_clzll(1LL);}"
    ${VARIABLE}
    )
ENDMACRO()

MACRO(CHECK_BSR64 VARIABLE)
  CHECK_C_SOURCE_COMPILES(
    "int main(int argc, char *argv[])
{unsigned long foo; unsigned __int64 bar=1LL;
return _BitScanReverse64(&foo, bar);}"
    ${VARIABLE}
    )
ENDMACRO()

MACRO(FIND_INLINE_KEYWORD)
  #Inspired from http://www.cmake.org/Wiki/CMakeTestInline
  SET(INLINE_TEST_SRC "/* Inspired by autoconf's c.m4 */
static inline int static_foo(){return 0\;}
int main(int argc, char *argv[]){return 0\;}
")
  FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CMakeTestCInline.c
    ${INLINE_TEST_SRC})

  FOREACH(KEYWORD "inline" "__inline__" "__inline")
    IF(NOT DEFINED C_INLINE)
      TRY_COMPILE(C_HAS_${KEYWORD}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeTestCInline.c
        COMPILE_DEFINITIONS "-Dinline=${KEYWORD}"
        )
      IF(C_HAS_${KEYWORD})
        SET(C_INLINE TRUE)
        ADD_DEFINITIONS("-Dinline=${KEYWORD}")
        MESSAGE(STATUS "Inline keyword found - ${KEYWORD}")
      ENDIF(C_HAS_${KEYWORD})
    ENDIF(NOT DEFINED C_INLINE)
  ENDFOREACH(KEYWORD)

  IF(NOT DEFINED C_INLINE)
    ADD_DEFINITIONS("-Dinline=")
    MESSAGE(STATUS "Inline keyword - not found")
  ENDIF(NOT DEFINED C_INLINE)
ENDMACRO(FIND_INLINE_KEYWORD)

MACRO(FIND_RESTRICT_KEYWORD)
  SET(RESTRICT_TEST_SRC "/* Inspired by autoconf's c.m4 */
int foo (int * restrict ip){return ip[0]\;}
int main(int argc, char *argv[]){int s[1]\;
int * restrict t = s\; t[0] = 0\; return foo(t)\;}
")

  FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CMakeTestCRestrict.c
    ${RESTRICT_TEST_SRC})

  FOREACH(KEYWORD "restrict" "__restrict" "__restrict__" "_Restrict")
    IF(NOT DEFINED C_RESTRICT)
      TRY_COMPILE(C_HAS_${KEYWORD}
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeTestCRestrict.c
        COMPILE_DEFINITIONS "-Drestrict=${KEYWORD}"
        )
      IF(C_HAS_${KEYWORD})
        SET(C_RESTRICT TRUE)
        ADD_DEFINITIONS("-Drestrict=${KEYWORD}")
        MESSAGE(STATUS "Restrict keyword found - ${KEYWORD}")
      ENDIF(C_HAS_${KEYWORD})
    ENDIF(NOT DEFINED C_RESTRICT)
  ENDFOREACH(KEYWORD)

  IF(NOT DEFINED C_RESTRICT)
    ADD_DEFINITIONS("-Drestrict=")
    MESSAGE(STATUS "Restrict keyword - not found")
  ENDIF(NOT DEFINED C_RESTRICT)
ENDMACRO(FIND_RESTRICT_KEYWORD)
