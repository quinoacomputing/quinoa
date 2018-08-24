################################################################################
#
# \file      cmake/add_regression_test.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Function used to add a regression test to the ctest test suite
#
################################################################################

# ##############################################################################
# Function used to add a regression test to the ctest test suite
#
# add_regression_test( <test_name> <executable>
#                      [NUMPES n]
#                      [INPUTFILES file1 file2 ...]
#                      [ARGS arg1 arg2 ...]
#                      [LABELS label1 label2 ...]
#                      [TEXT_DIFF_PROG txtdiff]
#                      [TEXT_BASELINE stat1.std stat2.std ...]
#                      [TEXT_RESULT stat1.txt stat2.txt ...]
#                      [TEXT_DIFF_PROG_CONF ndiff.cfg]
#                      [BIN_DIFF_PROG bindiff]
#                      [BIN_BASELINE stat1.std stat2.std ...]
#                      [BIN_RESULT stat1.bin stat2.bin ...]
#                      [BIN_DIFF_PROG_CONF exodiff.cfg]
#                      [FILECONV_PROG fileconv]
#                      [FILECONV_INPUT arg1 arg2 ...]
#                      [FILECONV_RESULT out.0.exo out.1.exo ...]
#                      [POSTPROCESS_PROG exec]
#                      [POSTPROCESS_PROG_ARGS arg1 arg2 ...]
#                      [POSTPROCESS_PROG_OUTPUT file]
#
# Mandatory arguments:
# --------------------
#
# <test_name> - Name of the test. Note that "<executable>" will be prefixed and
#               "_PEs{NUMPES}" will be postfixed to the test name. See also
#               argument NUMPES below.
#
# <executable> - Name of the executable to test.
#
# Optional arguments:
# -------------------
#
# NUMPES n - The number PEs to pass to charmrun, i.e., +pn. Default: 1.
#
# INPUTFILES file1 file2 ... - Input files required for the test. This list of
# files includes files needed for running the test, e.g., control file
# (containing user input), mesh file, etc., as well as file required for
# evaluating the test, e.g., diff program configuration file. All input files
# are soft-linked from the source dir to the build dir. Default: "".
#
# ARGS arg1 arg2 ... - Arguments to pass to executable tested. Default: "".
#
# LABELS label1 label2 ... - Optional labels associated with the test.
# Default: "${executable}".
#
# TEXT_DIFF_PROG txtdiff - Diff program used for textual diffs. Default:
# numdiff.
#
# TEXT_BASELINE stat1.std stat2.std ... - Textual file(s) containing the known
# good solutions. If empty, no textual diff is performed. Default: "". Note
# that the number of baseline filenames must equal the number of result files.
#
# TEXT_RESULT stat1.txt stat2.txt ... - Textual file(s) produced by the test to
# be tested. If empty, no textual diff is performed. Default: "". Note that the
# number of baseline filenames must equal the number of result files.
#
# TEXT_DIFF_PROG_CONF ndiff.cfg - Textual diff program configuration file.
# Default: "".
#
# TEXT_DIFF_PROG_ARGS arg1 arg2 ... - Textual diff program arguments.
# Default: "".
#
# BIN_DIFF_PROG bindiff - Diff program used for binary diffs. Default: exodiff.
#
# BIN_DIFF_PROG_ARGS arg1 arg2 ... - Binary diff program arguments. Default: "".
#
# BIN_BASELINE stat1.std stat2.std ... - Binary file(s) containing the known
# good solutions. If empty, no binary diff is performed. Default: "". Note
# that the number of baseline filenames must equal the number of result files.
#
# BIN_RESULT stat1.bin stat2.bin ... - Binary file(s) produced by the test to
# be tested. If empty, no binary diff is performed. Default: "". Note that the
# number of baseline filenames must equal the number of result files.
#
# BIN_DIFF_PROG_CONF exodiff1.cfg exodiff2.cfg ... - Binary diff program
# configuration file(s). Default: "".
#
# FILECONV_PROG fileconv - File conversion program to convert between field
# output files. Default: fileconv.
#
# FILECONV_INPUT arg1 arg2 ... - File converter program input files.
# Default: "".
#
# FILECONV_RESULT out.0.exo out.1.exo ... - Output files produced by the
# optional file converter step. Default: "".
#
# POSTPROCESS_PROG exec - Optional postprocess executable to run after test
# run. Default: "".
#
# POSTPROCESS_PROG_ARGS arg1 arg2 ... - Arguments to pass to POSTPROCESS_PROG.
# Default: "".
#
# POSTPROCESS_PROG_OUTPUT file - Filename to save the results of the
# postprocessor program. Default: "".
#
# Author: J. Bakosi
#
# ##############################################################################
function(ADD_REGRESSION_TEST test_name executable)

  set(oneValueArgs NUMPES TEXT_DIFF_PROG BIN_DIFF_PROG TEXT_DIFF_PROG_CONF
                   FILECONV_PROG POSTPROCESS_PROG POSTPROCESS_PROG_OUTPUT)
  set(multiValueArgs INPUTFILES ARGS TEXT_BASELINE TEXT_RESULT BIN_BASELINE
                     BIN_RESULT LABELS POSTPROCESS_PROG_ARGS BIN_DIFF_PROG_ARGS
                     TEXT_DIFF_PROG_ARGS BIN_DIFF_PROG_CONF FILECONV_RESULT
		     FILECONV_INPUT)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # Will collect test properties
  set(test_properties)

  # Set number of processing elements
  set(NUMPES 1)
  if (ARG_NUMPES)
    set(NUMPES ${ARG_NUMPES})
    list(APPEND test_properties PROCESSORS ${ARG_NUMPES})
  endif()

  # Add labels to test
  set(TEST_LABELS ${executable})        # ${executable} is always a label
  if (ARG_LABELS)
    list(APPEND TEST_LABELS ${ARG_LABELS})
  endif()
  # prepare test labels to pass as cmake script arguments
  set(ARG_LABELS LABELS ${TEST_LABELS})
  string(REPLACE ";" " " ARG_LABELS "${ARG_LABELS}")

  # Set textual diff tool
  set(TEXT_DIFF_PROG ${NUMDIFF_EXECUTABLE})
  if (ARG_TEXT_DIFF_PROG)
    set(TEXT_DIFF_PROG ${ARG_TEXT_DIFF_PROG})
  endif()
  # Set binary diff tool
  set(BIN_DIFF_PROG ${EXODIFF_EXECUTABLE})
  if (ARG_BIN_DIFF_PROG)
    set(BIN_DIFF_PROG ${ARG_BIN_DIFF_PROG})
  endif()
  # Set file converter tool
  # FILECONV_EXECUTABLE points to fileconv, refer src/Main/FileConv.C
  set(FILECONV_PROG ${FILECONV_EXECUTABLE})

  # Prefix executable and append NUMPES to test name
  set(test_name "${executable}:${test_name}_pe${NUMPES}")

  # Do sainity check on and prepare to pass as cmake script arguments the
  # filenames of text baseline(s) and text result(s)
  if(ARG_TEXT_BASELINE OR ARG_TEXT_RESULT)
    # Make sure the number of result and baseline files are equal
    list(LENGTH ARG_TEXT_BASELINE nbaseline)
    list(LENGTH ARG_TEXT_RESULT nresult)
    if (NOT nbaseline EQUAL nresult)
      message(FATAL_ERROR
              "Number of text baselines and number of results must be equal.")
    endif()
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_TEXT_BASELINE "${ARG_TEXT_BASELINE}")
    string(REPLACE ";" " " ARG_TEXT_RESULT "${ARG_TEXT_RESULT}")
  endif()

  # Do sainity check on and prepare to pass as cmake script arguments the
  # filenames of binary baseline(s) and binary result(s)
  if(ARG_BIN_BASELINE OR ARG_BIN_RESULT)
    # Make sure the number of result and baseline files are equal
    list(LENGTH ARG_BIN_BASELINE nbaseline)
    list(LENGTH ARG_BIN_RESULT nresult)
    if (NOT nbaseline EQUAL nresult)
      message(FATAL_ERROR
              "Number of binary baselines and number of results must be equal.")
    endif()
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_BIN_BASELINE "${ARG_BIN_BASELINE}")
    string(REPLACE ";" " " ARG_BIN_RESULT "${ARG_BIN_RESULT}")
  endif()

  if(ARG_POSTPROCESS_PROG_ARGS)
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_POSTPROCESS_PROG_ARGS
           "${ARG_POSTPROCESS_PROG_ARGS}")
  endif()

  if(ARG_TEXT_DIFF_PROG_ARGS)
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_TEXT_DIFF_PROG_ARGS "${ARG_TEXT_DIFF_PROG_ARGS}")
  endif()

  # Make exodiff quiet (errors and warnings will still come to output)
  list(APPEND ARG_BIN_DIFF_PROG_ARGS "-q")
  if(ARG_BIN_DIFF_PROG_ARGS)
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_BIN_DIFF_PROG_ARGS "${ARG_BIN_DIFF_PROG_ARGS}")
  endif()

  if(ARG_FILECONV_INPUT)
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_FILECONV_INPUT "${ARG_FILECONV_INPUT}")
  endif()

  # Do sainity check on and prepare to pass as cmake script arguments the
  # filenames of the file converter result(s)
  if(ARG_FILECONV_RESULT)
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_FILECONV_RESULT "${ARG_FILECONV_RESULT}")
  endif()

  if(ARG_BIN_DIFF_PROG_CONF)
    # Convert list to space-separated string for passing as arguments to test
    # runner cmake script below
    string(REPLACE ";" " " ARG_BIN_DIFF_PROG_CONF "${ARG_BIN_DIFF_PROG_CONF}")
  endif()

  # Construct and echo configuration for test being added
  set(msg "Add regression test ${test_name} for ${executable}")
  # Run all regression tests with quiescence detection
  list(APPEND ARG_ARGS "-q")
  if (ARG_ARGS)
    string(REPLACE ";" " " ARGUMENTS "${ARG_ARGS}")
    string(CONCAT msg "${msg}, args: '${ARGUMENTS}'")
  endif()

  # Set and create test run directory
  set(workdir ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  file(MAKE_DIRECTORY ${workdir})

  # Create list of files required to softlink to build directory
  set(reqfiles)
  foreach(file ${ARG_INPUTFILES})
    list(APPEND reqfiles "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
  endforeach()

  # Softlink files required to build directory
  foreach(target ${reqfiles})
    set(LN_COMMAND "ln -sf ${target} ${workdir}")
    exec_program(${LN_COMMAND} OUTPUT_VARIABLE ln_output RETURN_VALUE ln_retval)
    if("${ln_retval}" GREATER 0)
      message(FATAL_ERROR "Problem creating symlink from \"${target}\" to \"${workdir}\":\n${ln_output}")
    endif()
  endforeach()

  # Add the test. See test_runner.cmake for documentation of the arguments.
  add_test(NAME ${test_name}
           COMMAND ${CMAKE_COMMAND}
           -DTEST_NAME=${test_name}
           -DWORKDIR=${workdir}
           -DRUNNER_REQUIRED=${RUNNER_REQUIRED}
           -DRUNNER=${RUNNER}
           -DRUNNER_NCPUS_ARG=${RUNNER_NCPUS_ARG}
           -DRUNNER_ARGS=${RUNNER_ARGS}
           -DTEST_EXECUTABLE=${CMAKE_BINARY_DIR}/Main/${executable}
           -DTEST_EXECUTABLE_ARGS=${ARGUMENTS}
           -DTEST_LABELS=${ARG_LABELS}
           -DNUMPES=${NUMPES}
           -DTEXT_DIFF_PROG=${TEXT_DIFF_PROG}
           -DTEXT_DIFF_PROG_ARGS=${ARG_TEXT_DIFF_PROG_ARGS}
           -DTEXT_DIFF_PROG_CONF=${ARG_TEXT_DIFF_PROG_CONF}
           -DTEXT_BASELINE=${ARG_TEXT_BASELINE}
           -DTEXT_RESULT=${ARG_TEXT_RESULT}
           -DBIN_DIFF_PROG=${BIN_DIFF_PROG}
           -DBIN_DIFF_PROG_ARGS=${ARG_BIN_DIFF_PROG_ARGS}
           -DBIN_DIFF_PROG_CONF=${ARG_BIN_DIFF_PROG_CONF}
           -DBIN_BASELINE=${ARG_BIN_BASELINE}
           -DBIN_RESULT=${ARG_BIN_RESULT}
           -DFILECONV_PROG=${CMAKE_BINARY_DIR}/Main/${FILECONV_PROG}
           -DFILECONV_INPUT=${ARG_FILECONV_INPUT}
           -DFILECONV_RESULT=${ARG_FILECONV_RESULT}
           -DPOSTPROCESS_PROG=${ARG_POSTPROCESS_PROG}
           -DPOSTPROCESS_PROG_ARGS=${ARG_POSTPROCESS_PROG_ARGS}
           -DPOSTPROCESS_PROG_OUTPUT=${ARG_POSTPROCESS_PROG_OUTPUT}
           -P ${TEST_RUNNER}
           WORKING_DIRECTORY ${workdir})

  # build pass regular expression list for test
  set(pass_regexp "")
  # add pass regular expression for numdiff output if needed
  if (ARG_TEXT_BASELINE)
    list(APPEND pass_regexp ".*${test_name}.*PASS")
  endif()
  # add pass regular expression for rngtest output if needed
  if (test_name MATCHES "${RNGTEST_EXECUTABLE}")
    list(APPEND pass_regexp "Failed statistics" "All tests passed")
  endif()
  # add pass regular expression for exodiff output if needed
  if (ARG_BIN_BASELINE)
    list(APPEND pass_regexp "exodiff: Files are the same")
  endif()
  # add pass regular expression for when postprocessor not available, if needed
  if (ENABLE_MESHCONV AND NOT GMSH_FOUND)
    list(APPEND pass_regexp "would be required for this test to be rigorous")
  endif()

  # build fail regular expression list for test
  set(fail_regexp "")
  # add fail regular expression for numdiff output if needed
  if (ARG_TEXT_BASELINE)
    list(APPEND fail_regexp ".*${test_name}.*FAIL")
  endif()
  # add pass regular expression for exodiff output if needed
  if (ARG_BIN_BASELINE)
    list(APPEND fail_regexp "exodiff: Files are different" "exodiff: ERROR")
  endif()
  # add fail regular expression to detect cmake error during test run
  list(APPEND fail_regexp "CMake Error")

  #message("'${test_name}' pass regexp: ${pass_regexp}, fail regexp: ${fail_regexp}")

  # Set test properties and instruct ctest to check textual diff output against
  # the regular expressions specified.
  set_tests_properties(${test_name} PROPERTIES ${test_properties}
                       PASS_REGULAR_EXPRESSION "${pass_regexp}"
                       FAIL_REGULAR_EXPRESSION "${fail_regexp}")

  # Set labels cmake test property. The LABELS built-in cmake property is not
  # passed as part of test_properties above in set_test_properties as
  # TEST_LABELS is a cmake list and passing in lists of lists does not work as
  # expected.
  set_property(TEST ${test_name} PROPERTY LABELS ${TEST_LABELS})

endfunction()
