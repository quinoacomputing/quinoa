INCLUDE(TribitsLibraryMacros)

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

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/../../../Trilinos_version.h
  )

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  )

# Store Fortran modules in a dedicated directory.
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
INCLUDE_DIRECTORIES("${CMAKE_Fortran_MODULE_DIRECTORY}")

# Install *.mod from the module directory to the include directory.
INSTALL(
  DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/"
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
  COMPONENT ${PACKAGE_NAME}
  FILES_MATCHING PATTERN "*.mod"
  )

# Set sources for link-layer build

APPEND_SET(SOURCES
  ForTrilinos_enum_wrappers.F90
  ForTrilinos_enums.F90
  ForTrilinos_table_man.F90
  ForTrilinos_external_utils.F90
  ForTrilinos_utils.F90
  ${DIREPETRA}/forepetra.F90
  ${DIRTEUCHOS}/forteuchos.F90
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  APPEND_SET(SOURCES
    ${DIRAMESOS}/foramesos.F90
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_AztecOO)

  APPEND_SET(SOURCES
    ${DIRAZTECOO}/foraztecoo.F90
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  APPEND_SET(SOURCES
    ${DIRGALERI}/forgaleri.F90
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  APPEND_SET(SOURCES
    ${DIRIFPACK}/forifpack.F90
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Pliris)
  IF (TPL_ENABLE_MPI)

    APPEND_SET(SOURCES
      ${DIRPLIRIS}/forpliris.F90
      )

  ENDIF()
ENDIF()

# Set additional sources for object-oriented build
IF(${PACKAGE_NAME}_ENABLE_OBJECT_ORIENTED)
  APPEND_SET(SOURCES
    ForTrilinos_hermetic.F90
    ForTrilinos_ref_counter.F90
    ForTrilinos_universal.F90
    ForTrilinos_error.F90
    ForTrilinos_assertion_utility.F90
    )

  APPEND_GLOB(SOURCES ${DIREPETRA}/FEpetra*.F90
                      ${DIRTEUCHOS}/FTeuchos*.F90)

  IF (${PACKAGE_NAME}_ENABLE_Amesos)

    APPEND_GLOB(SOURCES ${DIRAMESOS}/FAmesos*.F90)

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_AztecOO)

    APPEND_GLOB(SOURCES ${DIRAZTECOO}/FAztecOO*.F90)

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Galeri)

    APPEND_GLOB(SOURCES ${DIRGALERI}/FGaleri*.F90)

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    APPEND_GLOB(SOURCES ${DIRIFPACK}/FIfpack*.F90)

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Pliris)
    IF (TPL_ENABLE_MPI)

      APPEND_GLOB(SOURCES ${DIRPLIRIS}/FPliris*.F90)

    ENDIF()
  ENDIF()

ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  fortrilinos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#SET_TARGET_PROPERTIES(Trilinos_fortrilinos PROPERTIES LINKER_LANGUAGE CXX)

#EXPORT (TARGETS Trilinos_fortrilinos APPEND FILE ForTrilinos_utils.mod)

#INSTALL (EXPORT DESTINATION ${Trilinos_INSTALL_INCLUDE_DIR})

SET_TARGET_PROPERTIES(fortrilinos PROPERTIES LINKER_LANGUAGE Fortran)

