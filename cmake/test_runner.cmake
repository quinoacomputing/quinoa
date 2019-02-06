################################################################################
#
# \file      cmake/test_runner.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Regression test runner using the cmake scripting language
#
################################################################################

# Covert string to list of file names of text baseline(s), text result(s), and
# file converter input(s) and output(s)
string(REPLACE " " ";" TEXT_BASELINE "${TEXT_BASELINE}")
string(REPLACE " " ";" TEXT_RESULT "${TEXT_RESULT}")
string(REPLACE " " ";" FILECONV_INPUT "${FILECONV_INPUT}")
string(REPLACE " " ";" FILECONV_RESULT "${FILECONV_RESULT}")
# Covert string to list of file names of binary baseline(s), binary result(s),
# binary diff program argument(s), and binary diff program confguration file(s)
string(REPLACE " " ";" BIN_BASELINE "${BIN_BASELINE}")
string(REPLACE " " ";" BIN_RESULT "${BIN_RESULT}")
string(REPLACE " " ";" BIN_DIFF_PROG_CONF "${BIN_DIFF_PROG_CONF}")
string(REPLACE " " ";" BIN_DIFF_PROG_ARGS "${BIN_DIFF_PROG_ARGS}")
# Covert string to list of postprocess program arguments
string(REPLACE " " ";" POSTPROCESS_PROG_ARGS "${POSTPROCESS_PROG_ARGS}")
# Covert string to list of test labels
string(REPLACE " " ";" TEST_LABELS "${TEST_LABELS}")
# Covert string to list of test executable arguments
string(REPLACE " " ";" TEST_EXECUTABLE_ARGS "${TEST_EXECUTABLE_ARGS}")
# Covert string to list of runner arguments
string(REPLACE " " ";" RUNNER_ARGS "${RUNNER_ARGS}")
# Covert string to list of postfix runner arguments
string(REPLACE " " ";" POSTFIX_RUNNER_ARGS "${POSTFIX_RUNNER_ARGS}")

# Print test runner configuration
message("Test runner configuration:")
message("  TEST_NAME (name of test)                                    : ${TEST_NAME}")
message("  WORKDIR (test run directory)                                : ${WORKDIR}")
message("  RUNNER_REQUIRED (true if an executable runner is required)  : ${RUNNER_REQUIRED}")
message("  RUNNER (used to run parallel and serial jobs inside cmake)  : ${RUNNER}")
message("  RUNNER_NCPUS_ARG (used to specify the number of CPUs)       : ${RUNNER_NCPUS_ARG}")
message("  CHARM_SMP (true/false indicating Charm++ SMP mode)          : ${CHARM_SMP}")
message("  RUNNER_ARGS (parallel/serial job runner arguments)          : ${RUNNER_ARGS}")
message("  POSTFIX_RUNNER_ARGS (postfix job runner arguments)          : ${POSTFIX_RUNNER_ARGS}")
message("  TEST_EXECUTABLE (executable tested)                         : ${TEST_EXECUTABLE}")
message("  TEST_EXECUTABLE_ARGS (executable arguments)                 : ${TEST_EXECUTABLE_ARGS}")
message("  TEST_LABELS (test labels)                                   : ${TEST_LABELS}")
message("  NUMPES (number of processing elements requested for test)   : ${NUMPES}")
message("  NUMNODES (number of logical nodes, in Charm++'s SMP mode)   : ${NUMNODES}")
message("  PPN (number of PEs per logical node, in Charm++'s SMP mode) : ${PPN}")
message("  HARDWARE_NUMPES (number of PEs used in hardware for test)   : ${HARDWARE_NUMPES}")
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

# Set Charm++'s +ppn argument (if configured, used in SMP mode)
if (PPN)
  set(PPN "+ppn;${PPN}")
endif()

# Configure test run command

# In Charm++'s SMP mode, if the runner is mpirun, -n (as RUNNER_NCPUS_ARG)
# specifies the number of compute nodes.
if (CHARM_SMP AND RUNNER MATCHES "mpirun")
  set(test_command ${RUNNER} ${RUNNER_NCPUS_ARG} ${NUMNODES} ${RUNNER_ARGS}
                   ${TEST_EXECUTABLE} ${TEST_EXECUTABLE_ARGS} ${PPN}
                   ${POSTFIX_RUNNER_ARGS})
else()
  set(test_command ${RUNNER} ${RUNNER_NCPUS_ARG} ${NUMPES} ${RUNNER_ARGS}
                   ${TEST_EXECUTABLE} ${TEST_EXECUTABLE_ARGS}
                   ${POSTFIX_RUNNER_ARGS})
endif()

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

    # Set runner for bindiff prog
    if (RUNNER_REQUIRED)
      set(runner_prefix ${RUNNER} ${RUNNER_NCPUS_ARG} 1 ${RUNNER_ARGS})
    endif()

    # Do binary diff(s) multiple times diffing baseline and result.
    #
    # There can be multiple baselines and multiple results if a test is
    # parallel and/or uses overdecomposition and thus the binary diff must be
    # done for a list of baseline + result pairs.
    #
    # In Charm++'s non-SMP mode each baseline must match its corresponding
    # result, i.e., the file numbers must match in lock-step: 0-0, 1-1, etc.
    # However, parallel tests whose baselines have been generated in non-SMP
    # mode (i.e., equivalent to ppn=1 in SMP mode), but run in SMP mode, may
    # yield the same partitions and the same results but with different file
    # numbers (ranks), due to the partitioner labeling the partitions
    # differently. In that case we need to search a matching result for all
    # baselines in the list of results generated by the test. Thus the loop
    # below is a nested double loop, to traverse the Cartesian product of the
    # lists of baselines and results. Note that the Cartesian product is only
    # needed in SMP mode. In non-SMP mode, only the baseline-result pairs are
    # diffed with the same file number, i.e., in lock-step.

    # A baseline file is passed if a matching result is found among the results
    # with any file number. Finding a matching result this way for all
    # baselines, however, is only a necessary condition for passing the entire
    # test. If a baseline or a result is a duplicate (same data with different
    # filenames), due to some error, the test must fail. Thus the sufficient
    # condition for passing the entire test is that if there is exactly one
    # matching baseline for every result, no more no less. This is tested after
    # the nested loop by ensuring that the list of binary output files produced
    # by the test (BINARY_RESULT) and the list of results matched up with the
    # baselines during diffing (matched_results) contain exactly the same
    # filenames, independent of the order.

    set(matched_results)
    set(b "-1")
    foreach(baseline ${BIN_BASELINE})
      math(EXPR b "0${b}+1")

      if (nconf EQUAL 1)
        list(GET BIN_DIFF_PROG_CONF 0 conf)
      else()
        list(GET BIN_DIFF_PROG_CONF ${b} conf)
      endif()

      set(pass FALSE)
      set(matching_result)
      set(baseline_error)
      set(r "-1")
      foreach(result ${BIN_RESULT})
        math(EXPR r "0${r}+1")

        if (NOT CHARM_SMP AND NOT b EQUAL r)
          #message("Charm++ in non-SMP mode: not diffing baseline ${b} with result ${r}")
          continue()
        endif()
        #message("Diffing baseline ${b} (${baseline}) with result ${r} (${result})")

        set(bin_diff_command ${runner_prefix} ${BIN_DIFF_PROG} ${BIN_DIFF_PROG_ARGS} -f ${conf} ${baseline} ${result})
        string(REPLACE ";" " " bin_diff_command_string "${bin_diff_command}")
        #message("Running binary diff command: '${bin_diff_command_string}'")
        execute_process(COMMAND ${bin_diff_command} RESULT_VARIABLE ERROR
                        ERROR_VARIABLE ERROR_OUTPUT OUTPUT_VARIABLE BINDIFF_OUTPUT)

        # remove warnings from exodiff output (this speeds up evaluating the test)
        string(REGEX REPLACE ".*WARNING.*" "" ERROR_OUTPUT "${ERROR_OUTPUT}")

        # A binary diff is passed if there is a test, diffing $baseline with
        # ANY $result, that pass. This allows for a partition of the mesh that
        # yields the same chunks (and the same result) but the ranks numbered
        # differently. This happens if the baseline for a parallel test was
        # saved by running the partitioner (e.g., zoltan in inciter) on, e.g.,
        # 4 ranks (in Charm++'s non-SMP mode), but the test is run on a single
        # rank with 4 threads (in Charm++'s SMP mode).
        if (NOT ERROR AND NOT ERROR_OUTPUT)
          set(pass TRUE)
          set(matching_result ${result})
          list(APPEND matched_results ${result})
        endif()

        # Save return value from binary diff command if any
        if(ERROR)
          set(baseline_error "${bin_diff_command_string}, error: ${ERROR}")
          if (ERROR_OUTPUT)
            set(baseline_error "${baseline_error}, output: ${BINDIFF_OUTPUT}, error output: ${ERROR_OUTPUT}")
          endif()
        endif()

      endforeach(result)

      # Echo pass message if test passed, echo output, error, and error output
      # if failed.
      if (pass)
        message(STATUS "Binary diff found match for '${baseline}' in '${matching_result}'")
      else()
        string(REPLACE "\n" " " baseline_error_out "${baseline_error}")
        message(FATAL_ERROR "Binary diff command failed: ${baseline_error_out}")
      endif()

    endforeach(baseline)

    # Cross-reference the lists of baselines and results: ensure all results
    # are matched up with a baseline and all baselines are matched up with a
    # result. Since the above double loop on baselines and results only ensures
    # searching for ANY match of a baseline to a result, we still have to check
    # if there was a match for all baselines and a match for all results. This
    # is done by comparing the list of baselines and matched results. If the
    # two lists are equal (independent of the order of their elements), the
    # test, containing multiple baselines and results, passed.
    foreach(m ${matched_results})
      list(FIND BIN_RESULT ${m} i)
      if (${i} EQUAL -1)
        message("${m} has not been matched to any of the baselines (${BIN_RESULT})")
      endif()
    endforeach()
    foreach(b ${BIN_RESULT})
      list(FIND matched_results ${b} i)
      if (${i} EQUAL -1)
        message("${b} has not been matched to any of the results (${matched_results})")
      endif()
    endforeach()

  endif()

endif()
