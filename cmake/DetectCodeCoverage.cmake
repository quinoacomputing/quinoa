################################################################################
#
# \file      cmake/DetectCodeCoverage.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Detect prerequesites for code coverage analysis
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

# Code coverage analysis only supported if all prerequisites found and the user
# has requested it via the cmake variable COVERAGE=on..
if ( COVERAGE AND
     GCOV AND
     LCOV AND
     GENHTML AND
     SED AND
     ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
       CMAKE_CXX_COMPILER_ID STREQUAL "GNU" ) )

  message(STATUS "Code coverage analysis enabled: All prerequisites found: COVERAGE:${COVERAGE}, compiler:${CMAKE_CXX_COMPILER_ID}, gcov:${GCOV}, lcov:${LCOV}, genhtml:${GENHTML}, sed:${SED}")

  # Enable code coverage analysis.
  SET(CODE_COVERAGE ON)

  # Make flag enabling code coverage analysis available in parent cmake scope
  MARK_AS_ADVANCED(CODE_COVERAGE)

  # Only include code coverage cmake functions if all prerequsites are met
  include(CodeCoverage)

else()

  message(STATUS "Code coverage analysis disabled. Not all prerequisites found: COVERAGE:${COVERAGE}, compiler:${CMAKE_CXX_COMPILER_ID}, gcov:${GCOV}, lcov:${LCOV}, genhtml:${GENHTML}, sed:${SED}")

endif()
