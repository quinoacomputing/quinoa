################################################################################
#
# \file      CodeCoverage.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Setup target for code coverage analysis
#
################################################################################

# ##############################################################################
# Function to add code coverage target
#
# setup_target_for_coverage( <suite> <path> <targetname> <testrunner>
#                            [TESTRUNNER_ARGS ...]
#                            [DEPENDS dep1 dep2 ... ] )
#
# Mandatory arguments:
# --------------------
#
# suite - Test suite name to be displayed in HTML report title.
#
# path - Path to prepend to where the report is generated:
# <path>${targetname}/index.html.
#
# targetname - The name of the code coverage target. The HTML report on code
# coverage is generated at the path <path>/${targetname}/index.html.
#
# testrunner - Command line of the test runner.
#
# Optional arguments:
# -------------------
#
# TESTRUNNER_ARGS arg1 arg2 ... - Optional arguments to test runner. Pass them
# in list form, e.g.: "-v;-g;group" for passing '-v -g group'. Default: "".
#
# DEPENDS dep1 dep2 ... - Optional dependencies added to test coverage target.
# Default: "". Here all dependencies should be given that should be covered by
# the test suite the coverage is being setup for, as well as those that are
# required for successfully building the tests and the test runner.
#
# ##############################################################################
FUNCTION(SETUP_TARGET_FOR_COVERAGE suite path targetname testrunner)

  set(multiValueArgs TESTRUNNER_ARGS DEPENDS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  IF(NOT LCOV)
    MESSAGE(FATAL_ERROR "lcov not found! Aborting...")
  ENDIF()

  IF(NOT GENHTML)
    MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
  ENDIF()

  IF(NOT SED)
    MESSAGE(FATAL_ERROR "sed not found! Aborting...")
  ENDIF()

  # Set shortcut for output: path/target
  set(OUTPUT ${path}/${targetname})
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${path})

  # Convert test runner args to (semi-colon-separated cmake) list of arguments
  string(REPLACE " " ";" ARG_TESTRUNNER_ARGS "${ARG_TESTRUNNER_ARGS}")

  # Setup code coverage target
  ADD_CUSTOM_TARGET(${targetname}
    # Cleanup lcov
    COMMAND ${LCOV} --gcov-tool ${GCOV} --directory . --zerocounters
    # Capture initial state yielding zero coverage baseline
    COMMAND ${LCOV} --gcov-tool ${GCOV} --capture --initial --directory . --output-file ${OUTPUT}.base.info
    # Run test suite
    COMMAND ${testrunner} ${ARG_TESTRUNNER_ARGS}
    # Capture lcov counters
    COMMAND ${LCOV} --gcov-tool ${GCOV}  --capture --rc lcov_branch_coverage=1 --directory . --output-file ${OUTPUT}.test.info
    # Combine trace files
    COMMAND ${LCOV} --gcov-tool ${GCOV} --rc lcov_branch_coverage=1 --add-tracefile ${OUTPUT}.base.info --add-tracefile ${OUTPUT}.test.info --output-file ${OUTPUT}.total.info
    # Filter out unwanted files
    COMMAND ${LCOV} --gcov-tool ${GCOV} --rc lcov_branch_coverage=1 --remove ${OUTPUT}.total.info "*/tests/*" "*/c++/*" "*/include/*" "*/boost/*" "*/charm/*" "*.decl.h" "*.def.h" "*/STDIN" "*/openmpi/*" "*/pstreams/*" "*/Random123/*" "*/pegtl/*" "*/tut/*" "*/highwayhash/*" "*/moduleinit*" --output-file ${OUTPUT}.filtered.info
    # Copy over report customization files for genhtml
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/quinoa.gcov.css
            ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/quinoa.lcov.prolog
            ${CMAKE_BINARY_DIR}
    # Generate HTML report
    COMMAND ${GENHTML} --legend --branch-coverage --demangle-cpp --css-file quinoa.gcov.css --ignore-errors source --html-prolog quinoa.lcov.prolog --title "${GIT_SHA1}" -o ${OUTPUT} ${OUTPUT}.filtered.info
    # Customize page headers in generated html to own
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/LCOV - code coverage report/Quinoa ${suite} test code coverage report/g" {} +
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/<td class=\"headerItem\">Test:<\\/td>/<td class=\"headerItem\">Commit:<\\/td>/g" {} +
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/Quinoa_v\\\(.*\\\)-\\\(.*\\\)-g\\\(.*\\\)<\\/td>/<a href=\"https:\\/\\/github.com\\/quinoacomputing\\/quinoa\\/commit\\/\\3\">Quinoa_v\\1-\\2-g\\3<\\/a><\\/td>/g" {} +
    # Cleanup intermediate data
    COMMAND ${CMAKE_COMMAND} -E remove ${OUTPUT}.base.info ${OUTPUT}.test.info ${OUTPUT}.total.info ${OUTPUT}.filtered.info
    # Set work directory for target
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # Echo what is being done
    COMMENT "Quinoa ${suite} test code coverage report"
    VERBATIM USES_TERMINAL
  )

  # Make test coverage target dependent on optional dependencies passed in using
  # keyword DEPENDS
  add_dependencies(${targetname} ${ARG_DEPENDS})

  # Output code coverage target enabled
  string(REPLACE ";" " " ARGUMENTS "${ARG_TESTRUNNER_ARGS}")
  message(STATUS "Enabling code coverage target '${targetname}' tested by '${testrunner} ${ARGUMENTS}', dependencies {${ARG_DEPENDS}}, report at ${OUTPUT}/index.html")

ENDFUNCTION()

# ##############################################################################
# Function to add code coverage target for all individual coverage targets
# included. This is very similar to setup_target_for_coverage(), defined above,
# but compiles a coverage report for all the individual test suites. The
# TESTRUNNER_ARGS are hard-coded.
#
# setup_target_for_all_coverage( <suite> <path> <targetname>
#                                [DEPENDS dep1 dep2 ... ] )
#
# Mandatory arguments:
# --------------------
#
# suite - Test suite name to be displayed in HTML report title.
#
# path - Path to prepend to where the report is generated:
# <path>${targetname}/index.html.
#
# targetname - The name of the code coverage target. The HTML report on code
# coverage is generated at the path <path>/${targetname}/index.html.
#
# unittestrunner - Command line of the unittest harness runner.
#
# unittestrunner_ncpus_arg - Command line argument specifying the number of cpus
# argument of the unittest harness runner.
#
# Optional arguments:
# -------------------
#
# DEPENDS dep1 dep2 ... - Optional dependencies added to test coverage target.
# Default: "". Here all dependencies should be given that should be covered by
# the test suite the coverage is being setup for, as well as those that are
# required for successfully building the tests and the test runner.
#
# ##############################################################################
FUNCTION(SETUP_TARGET_FOR_ALL_COVERAGE suite path targetname unittestrunner
                                       unittestrunner_ncpus_arg)

  set(multiValueArgs DEPENDS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  IF(NOT LCOV)
    MESSAGE(FATAL_ERROR "lcov not found! Aborting...")
  ENDIF()

  IF(NOT GENHTML)
    MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
  ENDIF()

  IF(NOT SED)
    MESSAGE(FATAL_ERROR "sed not found! Aborting...")
  ENDIF()

  # Set shorcut for output: path/target
  set(OUTPUT ${path}/${targetname})
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${path})

  # Setup code coverage target
  ADD_CUSTOM_TARGET(${targetname}
    # Cleanup lcov
    COMMAND ${LCOV} --gcov-tool ${GCOV} --directory . --zerocounters
    # Capture initial state yielding zero coverage baseline
    COMMAND ${LCOV} --gcov-tool ${GCOV} --capture --initial --directory . --output-file ${OUTPUT}.base.info
    # Run all test suites
    COMMAND ${unittestrunner} ${unittestrunner_ncpus_arg} ${PROCESSOR_COUNT} Main/${UNITTEST_EXECUTABLE} -v
    COMMAND ${CMAKE_CTEST_COMMAND} -j${PROCESSOR_COUNT}
    # Capture lcov counters
    COMMAND ${LCOV} --gcov-tool ${GCOV} --capture --rc lcov_branch_coverage=1 --directory . --output-file ${OUTPUT}.test.info
    # Combine trace files
    COMMAND ${LCOV} --gcov-tool ${GCOV} --rc lcov_branch_coverage=1 --add-tracefile ${OUTPUT}.base.info --add-tracefile ${OUTPUT}.test.info --output-file ${OUTPUT}.total.info
    # Filter out unwanted files
    COMMAND ${LCOV} --gcov-tool ${GCOV} --rc lcov_branch_coverage=1 --remove ${OUTPUT}.total.info "*/tests/*" "*/c++/*" "*/include/*" "*/boost/*" "*/charm/*" "*.decl.h" "*.def.h" "*/STDIN" "*/openmpi/*" "*/pstreams/*" "*/Random123/*" "*/pegtl/*" "*/tut/*" "*/highwayhash/*" "*/moduleinit*" --output-file ${OUTPUT}.filtered.info
    # Copy over report customization files for genhtml
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/quinoa.gcov.css
            ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/quinoa.lcov.prolog
            ${CMAKE_BINARY_DIR}
    # Generate HTML report
    COMMAND ${GENHTML} --legend --branch-coverage --demangle-cpp --css-file quinoa.gcov.css --ignore-errors source --html-prolog quinoa.lcov.prolog --title "${GIT_SHA1}" -o ${OUTPUT} ${OUTPUT}.filtered.info
    # Customize page headers in generated html to own
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/LCOV - code coverage report/Quinoa ${suite} test code coverage report/g" {} +
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/<td class=\"headerItem\">Test:<\\/td>/<td class=\"headerItem\">Commit:<\\/td>/g" {} +
    COMMAND find ${OUTPUT} -type f -exec ${SED} -i "s/Quinoa_v\\\(.*\\\)-\\\(.*\\\)-g\\\(.*\\\)<\\/td>/<a href=\"https:\\/\\/github.com\\/quinoacomputing\\/quinoa\\/commit\\/\\3\">Quinoa_v\\1-\\2-g\\3<\\/a><\\/td>/g" {} +
    # Cleanup intermediate data
    COMMAND ${CMAKE_COMMAND} -E remove ${OUTPUT}.base.info ${OUTPUT}.test.info ${OUTPUT}.total.info ${OUTPUT}.filtered.info
    # Set work directory for target
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # Echo what is being done
    COMMENT "Quinoa ${suite} test code coverage report"
    VERBATIM USES_TERMINAL
  )

  # Make test coverage target dependent on optional dependencies passed in using
  # keyword DEPENDS
  add_dependencies(${targetname} ${ARG_DEPENDS})

  # Output code coverage target enabled
  message(STATUS "Enabling code coverage target '${targetname}' generating full coverage, dependencies {${ARG_DEPENDS}}, report at ${OUTPUT}/index.html")

ENDFUNCTION()
