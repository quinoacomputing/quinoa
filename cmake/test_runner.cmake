# Regression test runner using the cmake scripting language. See also
# http://www.cmake.org/Wiki/CMake/Language_Syntax.

# Covert string to list of file names of text baseline(s) and text result(s)
string(REPLACE " " ";" TEXT_BASELINE "${TEXT_BASELINE}")
string(REPLACE " " ";" TEXT_RESULT "${TEXT_RESULT}")

# Print test runner configuration
message("Test runner configuration:")
message("  TEST_NAME (name of test)                                    : ${TEST_NAME}")
message("  WORKDIR (test run directory)                                : ${WORKDIR}")
message("  CHARMRUN (used to run Charm++ executables)                  : ${CHARMRUN}")
message("  MPIRUN_BIND_ARGS (mpirun process binding)                   : ${MPIRUN_BIND_ARGS}")
message("  TEST_EXECUTABLE (executable tested)                         : ${TEST_EXECUTABLE}")
message("  TEST_EXECUTABLE_ARGS (executable arguments)                 : ${TEST_EXECUTABLE_ARGS}")
message("  TEST_LABELS (test labels)                                   : ${TEST_LABLES}")
message("  NUMPES (number of processing elements used for test)        : ${NUMPES}")

message("  TEXT_DIFF_PROG (diff tool used for text diffs)              : ${TEXT_DIFF_PROG}")
message("  TEXT_DIFF_PROG_ARGS (text diff tool arguments)              : ${TEXT_DIFF_PROG_ARGS}")
message("  TEXT_DIFF_PROG_CONF (text diff tool configuration file)     : ${TEXT_DIFF_PROG_CONF}")
message("  TEXT_BASELINE (text output known good solution file(s))     : ${TEXT_BASELINE}")
message("  TEXT_RESULT (text output file(s) diffed with good solution) : ${TEXT_RESULT}")

message("  BIN_DIFF_PROG (diff tool used for binary diffs)             : ${BIN_DIFF_PROG}")
message("  BIN_DIFF_PROG_ARGS (binary diff tool arguments)             : ${BIN_DIFF_PROG_ARGS}")
message("  BIN_DIFF_PROG_CONF (binary diff tool configuration file)    : ${BIN_DIFF_PROG_CONF}")
message("  BIN_BASELINE (binary output known good solution file)       : ${BIN_BASELINE}")
message("  BIN_RESULT (binary output file diffed with good solution)   : ${BIN_RESULT}")

# Remove previous test output (if any)
if(TEXT_RESULT OR BIN_RESULT)
  message("\nRemoving existing result(s) (if any): ${TEXT_RESULT} ${BIN_RESULT}\n")
  file(REMOVE ${TEXT_RESULT} ${BIN_RESULT})
endif()

# Run the test
set(test_command ${CHARMRUN} +p${NUMPES} ${MPIRUN_BIND_ARGS}
                 ${TEST_EXECUTABLE} ${TEST_EXECUTABLE_ARGS})
string(REPLACE ";" " " test_command_string "${test_command}")
message("\nRunning test command: '${test_command_string}'\n")
execute_process(COMMAND ${test_command} RESULT_VARIABLE ERROR)

# Check return value from test
if(ERROR)

  message(FATAL_ERROR "Test failed to run: '${test_command_string}' returned error code: ${ERROR}")

else() # Test command ran successfully, attempt to do diffs

  # Do textual diff(s) if both TEXT_BASELINE and TEXT_RESULT have been specified
  if(TEXT_BASELINE AND TEXT_RESULT)

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
      set(text_diff_command ${TEXT_DIFF_PROG} -b -t ${TEST_NAME}
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

  endif(TEXT_BASELINE AND TEXT_RESULT)

endif()
