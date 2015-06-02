# Attempt to find tools required for code coverage analysis

find_program( GCOV_PATH gcov )
find_program( LCOV_PATH lcov )
find_program( GENHTML_PATH genhtml )
find_program( GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/tests )

# Code coverage analysis only supported if
#  * all prerequisites found,
#  * with GNU and Clang,
#  * on debug builds
if ( GCOV_PATH AND
     LCOV_PATH AND
     GENHTML_PATH AND
     GCOVR_PATH AND
     ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
       "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") AND
     (CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR
      CMAKE_BUILD_TYPE STREQUAL "COVERAGE") )

  message(STATUS "Code coverage analysis enabled: compiler:${CMAKE_CXX_COMPILER_ID}, build:${CMAKE_BUILD_TYPE}, gcov:${GCOV_PATH}, lcov:${LCOV_PATH}, genhtml:${GENHTML_PATH}, gcovr:${GCOVR_PATH}")

  # Enable code coverage analysis
  SET(ENABLE_CODE_COVERAGE ON)

  # Set compiler and linker flags enabling code coverage analysis
  SET(CMAKE_CXX_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
  SET(CMAKE_C_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage"
    CACHE STRING "Flags used by the C compiler during coverage builds."
    FORCE )
  SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used for linking binaries during coverage builds."
    FORCE )
  SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
    ""
    CACHE STRING "Flags used by the shared libraries linker during coverage builds."
    FORCE )

  # Make compiler and linker flags enabling code coverage analysis available
  MARK_AS_ADVANCED(
    ENABLE_CODE_COVERAGE
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )

  # Only include code coverage cmake functions if all prerequsites are met
  include(CodeCoverage)

else()

  message(STATUS "Code coverage analysis disabled as not all prerequisites found: compiler:${CMAKE_CXX_COMPILER_ID} (only supported with GNU and Clang), build:${CMAKE_BUILD_TYPE} (only supported with DEBUG), gcov:${GCOV_PATH}, lcov:${LCOV_PATH}, genhtml:${GENHTML_PATH}, gcovr:${GCOVR_PATH}")

endif()
