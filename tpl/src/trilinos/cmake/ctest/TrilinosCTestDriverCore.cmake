#
# Backward compatibility script for Trilinos CTest/CDash drivers
#

get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Set project defaults (only the defaults for variables set in
# TribitsCTestDriverCore.cmake)

SET(TRIBITS_PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
#MESSAGE("TRIBITS_PROJECT_ROOT = '${TRIBITS_PROJECT_ROOT}'")

SET(Trilinos_REPOSITORY_LOCATION_DEFAULT
  "software.sandia.gov:/space/git/Trilinos")
SET(Trilinos_REPOSITORY_LOCATION_NIGHTLY_DEFAULT
  "software.sandia.gov:/space/git/nightly/Trilinos")
SET(CTEST_DROP_SITE_COVERAGE_DEFAULT
  "testing.sandia.gov")
SET(CTEST_DROP_LOCATION_COVERAGE_DEFAULT
  "/extended/cdash/submit.php?project=Trilinos")

# Many of the existing scripts use the variable TRILINOS_CMAKE_DIR, so we set
# it here.
SET(TRILINOS_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR}/../)

INCLUDE("${TRIBITS_PROJECT_ROOT}/cmake/tribits/ctest/TribitsCTestDriverCore.cmake")

macro(TRILINOS_CTEST_DRIVER)
  TRIBITS_CTEST_DRIVER()
endmacro()
