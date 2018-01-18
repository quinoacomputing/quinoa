################################################################################
#
# \file      cmake/test_runner.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Regression test runner using the cmake scripting language
#
################################################################################

# Regression test runner using the cmake scripting language. See also
# http://www.cmake.org/Wiki/CMake/Language_Syntax.

# Covert string to list of file names of text baseline(s), text result(s), and
# file converter input(s) and output(s)
string(REPLACE " " ";" TEXT_BASELINE "${TEXT_BASELINE}")
string(REPLACE " " ";" TEXT_RESULT "${TEXT_RESULT}")
string(REPLACE " " ";" FILECONV_INPUT "${FILECONV_INPUT}")
string(REPLACE " " ";" FILECONV_RESULT "${FILECONV_RESULT}")
# Covert string to list of file names of binary baseline(s), binary result(s),
# and binary diff program confguration file(s)
string(REPLACE " " ";" BIN_BASELINE "${BIN_BASELINE}")
string(REPLACE " " ";" BIN_RESULT "${BIN_RESULT}")
string(REPLACE " " ";" BIN_DIFF_PROG_CONF "${BIN_DIFF_PROG_CONF}")
# Covert string to list of postprocess program arguments
string(REPLACE " " ";" POSTPROCESS_PROG_ARGS "${POSTPROCESS_PROG_ARGS}")
# Covert string to list of test labels
string(REPLACE " " ";" TEST_LABELS "${TEST_LABELS}")
# Covert string to list of test executable arguments
string(REPLACE " " ";" TEST_EXECUTABLE_ARGS "${TEST_EXECUTABLE_ARGS}")

# Print test runner configuration
message("Test runner configuration:")
message("  TEST_NAME (name of test)                                    : ${TEST_NAME}")
message("  WORKDIR (test run directory)                                : ${WORKDIR}")
message("  RUNNER (used to run parallel and serial jobs inside cmake)  : ${RUNNER}")
message("  RUNNER_NCPUS_ARG (used to specify the number of CPUs)       : ${RUNNER_NCPUS_ARG}")
message("  RUNNER_ARGS (parallel/serial job runner arguments)          : ${RUNNER_ARGS}")
message("  TEST_EXECUTABLE (executable tested)                         : ${TEST_EXECUTABLE}")
message("  TEST_EXECUTABLE_ARGS (executable arguments)                 : ${TEST_EXECUTABLE_ARGS}")
message("  TEST_LABELS (test labels)                                   : ${TEST_LABELS}")
message("  NUMPES (number of processing elements used for test)        : ${NUMPES}")
message("  POSTPROCESS_PROG (executable to run after test)             : ${POSTPROCESS_PROG}")
message("  POSTPROCESS_PROG_ARGS (postprocess program arguments)       : ${POSTPROCESS_PROG_ARGS}")
message("  POSTPROCESS_PROG_OUTPUT (postprocess program output file)   : ${POSTPROCESS_PROG_OUTPUT}")

message("  TEXT_DIFF_PROG (diff tool used for text diffs)              : ${TEXT_DIFF_PROG}")
message("  TEXT_DIFF_PROG_ARGS (text diff tool arguments)              : ${TEXT_DIFF_PROG_ARGS}")
message("  TEXT_DIFF_PROG_CONF (text diff tool configuration file)     : ${TEXT_DIFF_PROG_CONF}")
message("  TEXT_BASELINE (text output known good solution file(s))     : ${TEXT_BASELINE}")
message("  TEXT_RESULT (text output file(s) diffed with good solution) : ${TEXT_RESULT}")

message("  BIN_DIFF_PROG (diff tool used for binary diffs)             : ${BIN_DIFF_PROG}")
message("  BIN_DIFF_PROG_ARGS (binary diff tool arguments)             : ${BIN_DIFF_PROG_ARGS}")
message("  BIN_DIFF_PROG_CONF (binary diff tool configuration file(s)) : ${BIN_DIFF_PROG_CONF}")
message("  BIN_BASELINE (binary output known good solution file(s))    : ${BIN_BASELINE}")
message("  BIN_RESULT (binary output file(s) diffed with good solution): ${BIN_RESULT}")

message("  FILE_CONV_PROG (File conv tool program)                     : ${FILECONV_PROG}")
message("  FILE_CONV_INPUT (File conv tool input file(s))              : ${FILECONV_INPUT}")
message("  FILECONV_RESULT (File conv tool output file(s))             : ${FILECONV_RESULT}")

# Remove previous test output (if any)
if(TEXT_RESULT OR BIN_RESULT OR FILECONV_RESULT OR FILECONV_INPUT)
  message("\nRemoving existing result(s) (if any): ${TEXT_RESULT} ${BIN_RESULT} ${FILECONV_RESULT} ${FILECONV_INPUT}\n")
  file(REMOVE ${TEXT_RESULT} ${BIN_RESULT} ${FILECONV_RESULT} ${FILECONV_INPUT})
endif()

# Configure test run command
set(test_command ${RUNNER} ${RUNNER_NCPUS_ARG} ${NUMPES} ${RUNNER_ARGS}
                 ${TEST_EXECUTABLE} ${TEST_EXECUTABLE_ARGS})
string(REPLACE ";" " " test_command_string "${test_command}")
# Run the test
message("\nRunning test command: '${test_command_string}'\n")
execute_process(COMMAND ${test_command} RESULT_VARIABLE ERROR)

# Check return value from test
if(ERROR)

  message(FATAL_ERROR "Test failed to run: '${test_command_string}' returned error code: ${ERROR}")

else() # Test command ran successfully, attempt to do diffs

  # Run postprocessor if given
  if (POSTPROCESS_PROG)
    set(post_command ${POSTPROCESS_PROG} ${POSTPROCESS_PROG_ARGS})
    string(REPLACE ";" " " post_command_string "${post_command}")
    message("\nRunning postprocessor command: '${post_command_string}'\n")
    execute_process(COMMAND ${post_command} RESULT_VARIABLE ERROR
                    OUTPUT_FILE ${POSTPROCESS_PROG_OUTPUT})
    if(ERROR)
      message(FATAL_ERROR "Postprocessor failed to run: '${post_command_string}' returned error code: ${ERROR}")
    endif()
  elseif(POSTPROCESS_PROG_OUTPUT)
    message(WARNING "Postprocessor not found, but would be required for this test to be rigorous")
  endif()

  # Run fileconv program if args are specified
  if (FILECONV_INPUT)

    # Make sure the number of input and output files are equal
    list(LENGTH FILECONV_INPUT ninput)
    list(LENGTH FILECONV_RESULT nresult)
    if (NOT ninput EQUAL nresult)
      message(FATAL_ERROR
              "Number of input and number of result must be equal.")
    endif()

    # Run fileconv multiple times for all fileconv inputs
    math(EXPR b "0")
    foreach(input IN LISTS FILECONV_INPUT)
      list(GET FILECONV_RESULT ${b} result)
      set(fileconv_command ${FILECONV_PROG} -i ${input} -o ${result})
      string(REPLACE ";" " " fileconv_command_string "${fileconv_command}")
      message("\nRunning file convert command: '${fileconv_command_string}'")
      execute_process(COMMAND ${fileconv_command} RESULT_VARIABLE ERROR )
      if(ERROR)
        message(FATAL_ERROR "File converter failed, returned error code: ${ERROR}")
      endif(ERROR)
      math(EXPR b "${b}+1")
    endforeach(input)

  endif()

  # Do textual diff(s) if
  #  - both TEXT_BASELINE and TEXT_RESULT have been specified and not doing
  #  postprocessing, or
  #  - both TEXT_BASELINE and TEXT_RESULT have been specified and doing
  #  postprocessing (and postprocessing program has been found)
  if( (TEXT_BASELINE AND TEXT_RESULT AND NOT POSTPROCESS_PROG_OUTPUT) OR (TEXT_BASELINE AND TEXT_RESULT AND POSTPROCESS_PROG) )

    # Make sure the number of result and baseline files are equal
    list(LENGTH TEXT_BASELINE nbaseline)
    list(LENGTH TEXT_RESULT nresult)
    if (NOT nbaseline EQUAL nresult)
      message(FATAL_ERROR
              "Number of baselines and number of results must be equal.")
    endif()

    # Do textual diff(s) multiple times diffing matching baseline and result
    math(EXPR b "0")
    foreach(baseline IN LISTS TEXT_BASELINE)
      list(GET TEXT_RESULT ${b} result)
      if (RUNNER_REQUIRED)
        set(runner_prefix ${RUNNER} ${RUNNER_NCPUS_ARG} 1 ${RUNNER_ARGS})
      endif()
      set(text_diff_command ${runner_prefix}
                            ${TEXT_DIFF_PROG} ${TEXT_DIFF_PROG_ARGS}
                            -b -t ${TEST_NAME}
                            ${baseline} ${result} ${TEXT_DIFF_PROG_CONF})
      string(REPLACE ";" " " text_diff_command_string "${text_diff_command}")
      message("\nRunning text diff command: '${text_diff_command_string}'\n")
      execute_process(COMMAND ${text_diff_command} RESULT_VARIABLE ERROR)
      # Check return value from textual diff command
      if(ERROR)
        message(FATAL_ERROR "Textual diff failed to run: '${text_diff_command_string}' returned error code: ${ERROR}")
      endif(ERROR)
      math(EXPR b "${b}+1")
    endforeach(baseline)

  endif()

  # Do binary diff(s) if
  #  - both BIN_BASELINE and BIN_RESULT have been specified and not doing
  #  postprocessing, or
  #  - both BIN_BASELINE and BIN_RESULT have been specified and doing
  #  postprocessing (and postprocessing program has been found)
  if( (BIN_BASELINE AND BIN_RESULT AND NOT POSTPROCESS_PROG_OUTPUT) OR (BIN_BASELINE AND BIN_RESULT AND POSTPROCESS_PROG) )

    # Make sure the number of result and baseline files are equal
    list(LENGTH BIN_BASELINE nbaseline)
    list(LENGTH BIN_RESULT nresult)
    if (NOT nbaseline EQUAL nresult)
      message(FATAL_ERROR
              "Number of baselines and number of results must be equal.")
    endif()

    # If there is only one bin diff program conf, use that for all, if multiple,
    # use one of each
    list(LENGTH BIN_DIFF_PROG_CONF nconf)
    if (NOT nconf EQUAL nresult AND NOT nconf EQUAL 1)
      message(FATAL_ERROR "Number of bin-diff-prog conf files (${nconf}) should either be 1 or it must equal the number of results (${nresult}).")
    endif()

    # Do binary diff(s) multiple times diffing matching baseline and result
    math(EXPR b "0")
    foreach(baseline ${BIN_BASELINE})
      list(GET BIN_RESULT ${b} result)
      if (nconf EQUAL 1)
        list(GET BIN_DIFF_PROG_CONF 0 conf)
      else()
        list(GET BIN_DIFF_PROG_CONF ${b} conf)
      endif()
      if (RUNNER_REQUIRED)
        set(runner_prefix ${RUNNER} ${RUNNER_NCPUS_ARG} 1 ${RUNNER_ARGS})
      endif()
      set(bin_diff_command ${runner_prefix}
                           ${BIN_DIFF_PROG} ${BIN_DIFF_PROG_ARGS}
                           -f ${conf} ${baseline} ${result})
      string(REPLACE ";" " " bin_diff_command_string "${bin_diff_command}")
      message("\nRunning binary diff command: '${bin_diff_command_string}'\n")
      execute_process(COMMAND ${bin_diff_command} RESULT_VARIABLE ERROR
                      ERROR_QUIET)
      # Check return value from binary diff command
      if(ERROR)
        message(FATAL_ERROR "Binary diff returned error code: '${bin_diff_command_string}' returned error code: ${ERROR}")
      endif(ERROR)
      math(EXPR b "${b}+1")
    endforeach(baseline)

  endif()

endif()
