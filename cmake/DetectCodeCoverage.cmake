################################################################################
#
# \file      cmake/DetectCodeCoverage.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Detect prerequesites for code coverage analysis
# \date      Thu 02 Mar 2017 10:30:26 AM MST
#
################################################################################

# Attempt to find tools required for code coverage analysis

if( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  find_program( GCOV ${CMAKE_CURRENT_SOURCE_DIR}/../script/llvm-gcov.sh )
elseif( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
  find_program( GCOV gcov )
endif()

find_program( LCOV lcov )
find_program( GENHTML genhtml )
find_program( SED sed )

# Code coverage analysis only supported if all prerequisites (gcov, lcov,
# genhtml) found using the GNU compiler and only on debug builds.
if ( GCOV AND
     LCOV AND
     GENHTML AND
     SED AND
     ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
       CMAKE_CXX_COMPILER_ID STREQUAL "GNU" ) )

  message(STATUS "Code coverage analysis enabled: compiler:${CMAKE_CXX_COMPILER_ID}, build:${CMAKE_BUILD_TYPE}, gcov:${GCOV}, lcov:${LCOV}, genhtml:${GENHTML}, sed:${SED}")

  # Enable code coverage analysis
  SET(CODE_COVERAGE ON)

  # Make flag enabling code coverage analysis available in parent cmake scope
  MARK_AS_ADVANCED(CODE_COVERAGE)

  # Only include code coverage cmake functions if all prerequsites are met
  include(CodeCoverage)

else()

  message(STATUS "Code coverage analysis disabled. Not all prerequisites found: compiler:${CMAKE_CXX_COMPILER_ID} (only supported with GNU and Clang), build:${CMAKE_BUILD_TYPE} (only supported with DEBUG), gcov:${GCOV}, lcov:${LCOV}, genhtml:${GENHTML}, sed:${SED}")

endif()
