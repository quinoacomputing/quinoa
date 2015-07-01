
#
# On Mac OS X Gnu compilers, add dynamic lookup for undefined symbols
# to the pytrilinos library and PyTrilinos extension modules
SET(EXTRA_LINK_ARGS "")
IF(APPLE)
  IF((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    APPEND_SET(EXTRA_LINK_ARGS "-undefined dynamic_lookup")
  ENDIF()
ENDIF(APPLE)


#
# Core PyTrilinos headers and sources
APPEND_SET(HEADERS
  PySundanceCallback.hpp
  PySundanceCellPredicate.hpp
  PySundanceLinearSolver.hpp
  PyTeuchos_Utils.hpp
  PlayaHack.hpp
  )

APPEND_SET(SOURCES
  PySundanceCallback.cpp
  PySundanceCellPredicate.cpp
  PySundanceLinearSolver.cpp
  PyTeuchos_Utils.cpp
  )

# Get the list of all packages, assumed to be in order of dependency
# (most-dependent packages last). The Trilinos cmake system specifies 
# that the master list of packages be in this order, so this assumption
# should be valid unless there's an error in the Trilinos package file.
SET(PACK_LIST ${Trilinos_PACKAGES})

# Remove from the package list all packages after ${TRIBITS_PACKAGE} in the dependency
# chain. This way each package can create its own minimalist export makefile
# with no upstream libraries. 
SET(TMP_PACK_LIST)
SET(SKIP FALSE)
FOREACH(PACK_NAME ${PACK_LIST})
  IF(NOT SKIP)
    LIST(APPEND TMP_PACK_LIST ${PACK_NAME})
    IF(TRIBITS_PACKAGE STREQUAL ${PACK_NAME})
      SET(SKIP TRUE)
    ENDIF()
  ENDIF()
ENDFOREACH()

SET(PACK_LIST ${TMP_PACK_LIST})


# Reverse the order of the package list, letting us loop 
# from most-dependent to least-dependent. 
LIST(REVERSE PACK_LIST)

# Create empty variables for the lists of libraries, TPL libs, and include
# paths. 
SET(LIB_LIST)
SET(TPL_LIB_LIST)
SET(TPL_INCLUDE_LIST)

# ------------------------------------
# --------- Main loop over packages. 
# ------------------------------------
FOREACH(PACK_NAME ${PACK_LIST})

  # Skip packages that haven't been enabled. 
  IF(Trilinos_ENABLE_${PACK_NAME})
    
    # ------- library handling -------------

    # Append this package's libraries to the list.
    LIST(APPEND LIB_LIST ${${PACK_NAME}_LIBRARIES})

    # ------ required TPL handling ---------

    # Get the required TPLs for this package.
    SET(REQ_TPLS ${${PACK_NAME}_LIB_REQUIRED_DEP_TPLS})
    # The TPL list is in reverse order for some reason. 
    IF(REQ_TPLS)
      LIST(REVERSE REQ_TPLS)
    ENDIF()
    # Append each required TPL lib to the main TPL lib list. Likewise the
    # include paths. 
    FOREACH(TPL_NAME ${REQ_TPLS})
      LIST(APPEND TPL_LIB_LIST ${TPL_${TPL_NAME}_LIBRARIES})
      LIST(APPEND TPL_INCLUDE_LIST ${TPL_${TPL_NAME}_INCLUDE_DIRS})
    ENDFOREACH()

    # ------ optional TPL handling ---------

    # Get the optional TPLs for this package.
    SET(OPT_TPLS ${${PACK_NAME}_LIB_OPTIONAL_DEP_TPLS})
    # The TPL list is in reverse order for some reason. 
    IF(OPT_TPLS)
      LIST(REVERSE OPT_TPLS)
    ENDIF()

    # Append each enabled optional TPL lib to the main TPL lib list. 
    # Likewise the include paths. 
    FOREACH(TPL_NAME ${OPT_TPLS})
      # Skip TPLs that haven't been anabled.
      IF(TPL_ENABLE_${TPL_NAME}) 
        LIST(APPEND TPL_LIB_LIST ${TPL_${TPL_NAME}_LIBRARIES})
        LIST(APPEND TPL_INCLUDE_LIST ${TPL_${TPL_NAME}_INCLUDE_DIRS})
      ENDIF()
    ENDFOREACH()
  ENDIF()
ENDFOREACH()

PRINT_VAR(LIB_LIST)
PRINT_VAR(TPL_LIB_LIST)
LIST(APPEND LIB_LIST ${TPL_LIB_LIST})
LIST(APPEND LIB_LIST libgfortran.so)
PRINT_VAR(LIB_LIST)

#
# Define the target for the swig runtime header
SET(SWIG_RUNTIME swigpyrun.h)
ADD_CUSTOM_COMMAND(
  OUTPUT ${SWIG_RUNTIME}
  COMMAND ${SWIG_EXECUTABLE} -python -external-runtime
  )

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(AFTER ${Sundance_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${Sundance_LIBRARY_DIRS})


#
# Loop over the PyTrilinos-enabled Trilinos modules and define the
# swig-generated PyTrilinos extension modules
SET(MODULES PySundance)

PRINT_VAR(Sundance_LIBRARIES)




# Create the __init__.py file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/__init__.src
  ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)



SET_SOURCE_FILES_PROPERTIES(PySundance.i PROPERTIES CPLUSPLUS ON)

SWIG_ADD_MODULE(PySundance python PySundance.i ${SOURCES})
SWIG_LINK_LIBRARIES(PySundance ${LIB_LIST})

#SET_TARGET_PROPERTIES(PySundance PROPERTIES LINK_FLAGS "${EXTRA_LINK_ARGS}")



