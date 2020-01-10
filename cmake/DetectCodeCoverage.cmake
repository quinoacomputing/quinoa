################################################################################
#
# \file      DetectCodeCoverage.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Detect prerequesites for code coverage analysis
#
################################################################################

# Attempt to find tools required for code coverage analysis

find_program( GCOV gcov )
find_program( FASTCOV fastcov.py )
find_program( GENHTML genhtml )
find_program( SED sed )

# Code coverage analysis only supported if all prerequisites found and the user
# has requested it via the cmake variable COVERAGE=on.
if ( COVERAGE AND
     GCOV AND
     FASTCOV AND
     GENHTML AND
     SED AND
     CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )

  message(STATUS "Code coverage analysis enabled with compiler:${CMAKE_CXX_COMPILER_ID}, gcov:${GCOV}, fastcov:${FASTCOV}, genhtml:${GENHTML}, sed:${SED}")

  # Enable code coverage analysis.
  SET(CODE_COVERAGE ON)

  # Make flag enabling code coverage analysis available in parent cmake scope
  MARK_AS_ADVANCED(CODE_COVERAGE)

  # Only include code coverage cmake functions if all prerequsites are met
  include(CodeCoverage)

else()

  message(STATUS "Code coverage analysis disabled. Not all prerequisites found: COVERAGE (must be true):${COVERAGE}, compiler (must be GNU):${CMAKE_CXX_COMPILER_ID}, and the following must be found: gcov:${GCOV}, fastcov:${FASTCOV}, genhtml:${GENHTML}, sed:${SED}")

endif()
