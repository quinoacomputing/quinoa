MACRO(MAKEFILE_VAL file name varname)

SET(MAKE_STUB ${CMAKE_CURRENT_BINARY_DIR}/make.stub)
SET(MAKE_OUT ${CMAKE_CURRENT_BINARY_DIR}/make.out)

FILE(WRITE ${MAKE_STUB}
"
default:
	echo \${${name}} >make.out

include ${file}
")

# Should use ${CMAKE_MAKE_PROGRAM}, or is that set to
# something other than nmake when using visual studio?
IF(WIN32)
  EXECUTE_PROCESS(WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMAND nmake /f ${MAKE_STUB})
ELSE(WIN32)
  EXECUTE_PROCESS(WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMAND make -f ${MAKE_STUB})
ENDIF(WIN32) 

FILE(READ ${MAKE_OUT} ${varname})
FILE(REMOVE ${MAKE_STUB} ${MAKE_OUT})

ENDMACRO(MAKEFILE_VAL)
