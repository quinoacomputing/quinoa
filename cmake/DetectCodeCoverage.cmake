# Attempt to find tools required for code coverage analysis
find_program( GCOV gcov )
find_program( LCOV lcov )
find_program( GENHTML genhtml )

# Code coverage analysis only supported if all prerequisites (gcov, lcov,
# genhtml) found using the GNU compiler and only on debug builds.
if ( GCOV AND
     LCOV AND
     GENHTML AND
     CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
     CMAKE_BUILD_TYPE STREQUAL "DEBUG" )

  message(STATUS "Code coverage analysis enabled: compiler:${CMAKE_CXX_COMPILER_ID}, build:${CMAKE_BUILD_TYPE}, gcov:${GCOV}, lcov:${LCOV}, genhtml:${GENHTML}")

  # Enable code coverage analysis
  SET(CODE_COVERAGE ON)

  # Make flag enablin code coverage analysis available
  MARK_AS_ADVANCED(CODE_COVERAGE)

  # Only include code coverage cmake functions if all prerequsites are met
  include(CodeCoverage)

else()

  message(STATUS "Code coverage analysis disabled. Not all prerequisites found: compiler:${CMAKE_CXX_COMPILER_ID} (only supported with GNU), build:${CMAKE_BUILD_TYPE} (only supported with DEBUG), gcov:${GCOV}, lcov:${LCOV}, genhtml:${GENHTML}")

endif()
