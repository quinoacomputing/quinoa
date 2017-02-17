# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


INCLUDE(TribitsAddExecutableTestHelpers)
INCLUDE(TribitsAddTestHelpers)
INCLUDE(TribitsGeneralMacros)

INCLUDE(PrintVar)
INCLUDE(AppendSet)
INCLUDE(ParseVariableArguments)


#
# TRIBITS_ADD_EXECUTABLE(...): Function that adds a test/example executable.
#
# TRIBITS_ADD_EXECUTABLE(
#   <execName>
#   SOURCES <src1> <src2> ...
#   [CATEGORIES <category1>  <category2> ...]
#   [HOST <host1> <host2> ...]
#   [XHOST <host1> <host2> ...]
#   [HOSTTYPE <hosttype1> <hosttype2> ...]
#   [XHOSTTYPE <hosttype1> <hosttype2> ...]
#   [NOEXEPREFIX ]
#   [NOEXESUFFIX ]
#   [DIRECTORY <dir> ]
#   [DEPLIBS <lib1> <lib2> ... ]
#   [COMM [serial] [mpi] ]
#   [LINKER_LANGUAGE [C|CXX|Fortran] ]
#   [ADD_DIR_TO_NAME ]
#   [DEFINES <-DSOMEDEFINE>]
#   [INSTALLABLE]
#   )
# 

FUNCTION(TRIBITS_ADD_EXECUTABLE EXE_NAME)

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRIBITS_ADD_EXECUTABLE: ${EXE_NAME} ${ARGN}")
  ENDIF()
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    "SOURCES;CATEGORIES;HOST;XHOST;HOSTTYPE;XHOSTTYPE;DIRECTORY;DEPLIBS;COMM;LINKER_LANGUAGE;DEFINES"
    #options
    "NOEXEPREFIX;NOEXESUFFIX;ADD_DIR_TO_NAME;INSTALLABLE"
    ${ARGN}
    )

  #
  # B) Exclude building the test executable based on some several criteria
  #

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_CATEGORIES(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  SET(ADD_THE_TEST FALSE)
  TRIBITS_ADD_TEST_PROCESS_HOST_HOSTTYPE(ADD_THE_TEST)
  IF (NOT ADD_THE_TEST)
    RETURN()
  ENDIF()

  TRIBITS_PROCESS_COMM_ARGS(ADD_SERIAL_EXE  ADD_MPI_EXE  ${PARSE_COMM})
  IF (NOT ADD_SERIAL_EXE AND NOT ADD_MPI_EXE)
    RETURN()
  ENDIF()

  #
  # C) Add the executable
  #
  
  SET (EXE_SOURCES)
  SET(EXE_BINARY_NAME ${EXE_NAME})
  
  # If requested create a modifier for the name that will be inserted between
  # the package name and the given name or exe_name for the test
  IF(PARSE_ADD_DIR_TO_NAME)
    SET(DIRECTORY_NAME "")
    TRIBITS_CREATE_NAME_FROM_CURRENT_SOURCE_DIRECTORY(DIRECTORY_NAME)
    SET(EXE_BINARY_NAME ${DIRECTORY_NAME}_${EXE_BINARY_NAME})
  ENDIF()
  
  IF(DEFINED PACKAGE_NAME AND NOT PARSE_NOEXEPREFIX)
    SET(EXE_BINARY_NAME ${PACKAGE_NAME}_${EXE_BINARY_NAME})
  ENDIF()

  # If exe is in subdirectory prepend that dir name to the source files
  IF(PARSE_DIRECTORY ) 
    FOREACH( SOURCE_FILE ${PARSE_SOURCES} )
      IF(IS_ABSOLUTE ${SOURCE_FILE})
        SET (EXE_SOURCES ${EXE_SOURCES} ${SOURCE_FILE})
      ELSE()
        SET (EXE_SOURCES ${EXE_SOURCES} ${PARSE_DIRECTORY}/${SOURCE_FILE})
      ENDIF()
    ENDFOREACH( )
  ELSE()
    FOREACH( SOURCE_FILE ${PARSE_SOURCES} )
      SET (EXE_SOURCES ${EXE_SOURCES} ${SOURCE_FILE})
    ENDFOREACH( )
  ENDIF()

  FOREACH(DEPLIB ${PARSE_DEPLIBS})
    IF (${DEPLIB}_INCLUDE_DIRS)
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Adding include directories ${DEPLIB}_INCLUDE_DIRS ...")
        #PRINT_VAR(${DEPLIB}_INCLUDE_DIRS)
      ENDIF()
      INCLUDE_DIRECTORIES(${${DEPLIB}_INCLUDE_DIRS})
    ENDIF()
  ENDFOREACH()

  IF (PARSE_DEFINES)
    ADD_DEFINITIONS(${PARSE_DEFINES})
  ENDIF()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("TRIBITS_ADD_EXECUTABLE: ADD_EXECTUABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})")
  ENDIF()
  ADD_EXECUTABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})
  APPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${EXE_BINARY_NAME})

  IF(PARSE_NOEXESUFFIX AND NOT WIN32)
    SET_TARGET_PROPERTIES(${EXE_BINARY_NAME} PROPERTIES SUFFIX "")
  ENDIF()

  IF(PARSE_LINKER_LANGUAGE)
    SET(LINKER_LANGUAGE ${PARSE_LINKER_LANGUAGE})
  ELSEIF (${PROJECT_NAME}_ENABLE_CXX)
    SET(LINKER_LANGUAGE CXX)
  ELSEIF(${PROJECT_NAME}_ENABLE_C)
    SET(LINKER_LANGUAGE C)
  ELSE()
    SET(LINKER_LANGUAGE)
  ENDIF()

  IF (LINKER_LANGUAGE)
    SET_PROPERTY(TARGET ${EXE_BINARY_NAME} APPEND PROPERTY
      LINKER_LANGUAGE ${LINKER_LANGUAGE})
  ENDIF()

  SET(LINK_LIBS)

  # First, add in the passed in dependent libraries
  IF (PARSE_DEPLIBS)
    APPEND_SET(LINK_LIBS ${PARSE_DEPLIBS})
  ENDIF()
  # 2009/01/09: rabartl: Above, I moved the list of dependent
  # libraries first to get around a problem with test-only libraries
  # creating multiple duplicate libraries on the link line with
  # CMake.

  # Second, add the package's own regular libraries
  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    APPEND_SET(LINK_LIBS ${${PACKAGE_NAME}_LIBRARIES})
  ELSE()
    APPEND_SET(LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})
  ENDIF()

  # Call INCLUDE_DIRECTORIES() and LINK_DIRECTORIES(...) for upstream
  # dependent Packages and TPLs and accumulate the list of libraries that will
  # need to be linked to.

  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    AND NOT ${PACKAGE_NAME}_INCLUDE_DIRS
    )
    # No libraries have been added for this package so
    # add the upstream package and TPL includes and libraries
    TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS) 
    TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  LIB  LINK_LIBS)
  ENDIF()

  TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
    ${PACKAGE_NAME}  TEST  LINK_LIBS) 
  
  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
      ${PACKAGE_NAME}  TEST  LINK_LIBS)
  ELSE()
    APPEND_SET(LINK_LIBS ${${PACKAGE_NAME}_INSTALLATION_TPL_LIBRARIES})
  ENDIF()

  # Last, add last_lib to get extra link options on the link line
  IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    APPEND_SET(LINK_LIBS last_lib)
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(LINK_LIBS)
  ENDIF()

  TARGET_LINK_LIBRARIES(${EXE_BINARY_NAME} ${LINK_LIBS})

  IF ("${CMAKE_VERSION}" VERSION_GREATER "2.8.4")
    ASSERT_DEFINED(${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
    IF (${PROJECT_NAME}_LINK_SEARCH_START_STATIC)
      #MESSAGE("${EXE_BINARY_NAME}: Adding property LINK_SEARCH_START_STATIC")
      SET_PROPERTY(TARGET ${EXE_BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
    ENDIF()
  ENDIF()
  

  IF(PARSE_DIRECTORY)
    SET_TARGET_PROPERTIES( ${EXE_BINARY_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${PARSE_DIRECTORY} )
  ENDIF()

  SET_PROPERTY(TARGET ${EXE_BINARY_NAME} APPEND PROPERTY
    LABELS ${PARENT_PACKAGE_NAME})

  IF(${PROJECT_NAME}_INSTALL_EXECUTABLES AND PARSE_INSTALLABLE)
    INSTALL(
      TARGETS ${EXE_BINARY_NAME}
      EXPORT ${PROJECT_NAME}
        DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}
      COMPONENT ${PACKAGE_NAME}
    )
  ENDIF()
ENDFUNCTION()


#
# Setup include directories and library dependencies
#

IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  MESSAGE("TribitsAddExecutable.cmake")
  PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
  PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
ENDIF()

IF (NOT PACKAGE_ADD_EXECUTABLE_UNIT_TESTING)
  INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
    ${${PACKAGE_NAME}_INCLUDE_DIRS})
  SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS
    ${${PACKAGE_NAME}_LIBRARY_DIRS})
ENDIF()
