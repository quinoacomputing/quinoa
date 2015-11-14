

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIREPETRA ${CMAKE_CURRENT_SOURCE_DIR}/epetra)
SET_AND_INC_DIRS(DIRTEUCHOS ${CMAKE_CURRENT_SOURCE_DIR}/teuchos)
SET_AND_INC_DIRS(DIRAMESOS ${CMAKE_CURRENT_SOURCE_DIR}/amesos)
SET_AND_INC_DIRS(DIRAZTECOO ${CMAKE_CURRENT_SOURCE_DIR}/aztecoo)
SET_AND_INC_DIRS(DIRGALERI ${CMAKE_CURRENT_SOURCE_DIR}/galeri)
SET_AND_INC_DIRS(DIRIFPACK ${CMAKE_CURRENT_SOURCE_DIR}/ifpack)
SET_AND_INC_DIRS(DIRPLIRIS ${CMAKE_CURRENT_SOURCE_DIR}/pliris)

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

# Set sources and headers

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_Version.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_enums.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_flex_enums.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_table_man.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_external_utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_enum_wrappers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_exceptions.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_utils.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_utils_templ.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_test_utils.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_Table.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_TableRepos.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_difficult.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CNOX_Interface.hpp
  )

APPEND_SET(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_TableRepos.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_utils.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_test_utils.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_table_man.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_external_utils.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CTrilinos_difficult.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CNOX_Interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CNOX_Driver.cpp
  )

APPEND_GLOB(HEADERS ${DIREPETRA}/*.hpp ${DIREPETRA}/*.h
                    ${DIRTEUCHOS}/*.hpp ${DIRTEUCHOS}/*.h)
APPEND_GLOB(SOURCES ${DIREPETRA}/*.cpp ${DIRTEUCHOS}/*.cpp)

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  APPEND_GLOB(HEADERS ${DIRAMESOS}/*.hpp ${DIRAMESOS}/*.h)
  APPEND_GLOB(SOURCES ${DIRAMESOS}/*.cpp)

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_AztecOO)

  APPEND_GLOB(HEADERS ${DIRAZTECOO}/*.hpp ${DIRAZTECOO}/*.h)
  APPEND_GLOB(SOURCES ${DIRAZTECOO}/*.cpp)

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  APPEND_GLOB(HEADERS ${DIRGALERI}/*.hpp ${DIRGALERI}/*.h)
  APPEND_GLOB(SOURCES ${DIRGALERI}/*.cpp)

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  APPEND_GLOB(HEADERS ${DIRIFPACK}/*.hpp ${DIRIFPACK}/*.h)
  APPEND_GLOB(SOURCES ${DIRIFPACK}/*.cpp)

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Pliris)
  IF (TPL_ENABLE_MPI)

    APPEND_GLOB(HEADERS ${DIRPLIRIS}/*.hpp ${DIRPLIRIS}/*.h)
    APPEND_GLOB(SOURCES ${DIRPLIRIS}/*.cpp)

  ENDIF()
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  ctrilinos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
