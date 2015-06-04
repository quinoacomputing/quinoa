# Add code coverage target
# Param _targetname The name of new the custom make target
# Param _testrunner The name of the target which runs the tests.
# Param _outputname HTML report is generated in _outputname/index.html
# Param _ncpus Nunmber of CPU cores testrunner should use
# Optional fourth parameter is passed as arguments to _testrunner. Pass them in
# list form, e.g.: "-v;-g;group" for -v -g group.
FUNCTION(SETUP_TARGET_FOR_COVERAGE _targetname _testrunner _outputname _ncpus)

  IF(NOT LCOV)
    MESSAGE(FATAL_ERROR "lcov not found! Aborting...")
  ENDIF()

  IF(NOT GENHTML)
    MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
  ENDIF()

  # Setup code coverage target
  ADD_CUSTOM_TARGET(${_targetname}
    # Cleanup lcov
    ${LCOV} --directory . --zerocounters
    # Capture initial state yielding zero coverage baseline
    COMMAND ${LCOV} --capture --initial --directory . --output-file ${_outputname}.base.info
    # Run test suite
    COMMAND ./charmrun +p${_ncpus} ${_testrunner} ${ARGV4}
    # Capture lcov counters
    COMMAND ${LCOV} --capture --rc lcov_branch_coverage=1 --directory . --output-file ${_outputname}.test.info
    # Combine trace files
    COMMAND ${LCOV} --rc lcov_branch_coverage=1 --add-tracefile ${_outputname}.base.info --add-tracefile ${_outputname}.test.info --output-file ${_outputname}.total.info
    # Filter out unwanted files
    COMMAND ${LCOV} --rc lcov_branch_coverage=1 --remove ${_outputname}.total.info 'UnitTest/tests/*' 'c++/*' 'boost/*' 'charm/*' '*.decl.h' '*.def.h' 'openmpi/*' 'pstreams/*' 'pegtl/*' 'tut/*' 'moduleinit*' --output-file ${_outputname}.filtered.info
    # Copy over report customization files for genhtml
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/quinoa.gcov.css
            ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/../doc/quinoa.lcov.prolog
            ${CMAKE_BINARY_DIR}
    # Generate HTML report
    COMMAND ${GENHTML} --legend --branch-coverage --demangle-cpp --css-file quinoa.gcov.css --html-prolog quinoa.lcov.prolog --title "${GIT_REFSPEC}:${GIT_SHA1}" -o ${_outputname} ${_outputname}.filtered.info
    # Customize page headers in generated html to own
    COMMAND find ${_outputname} -type f -print | xargs file | grep text | cut -f1 -d: | xargs sed -i 's/LCOV - code coverage report/Quinoa unit test code coverage report/g'
    COMMAND find ${_outputname} -type f -print | xargs file | grep text | cut -f1 -d: | xargs sed -i 's/<td class="headerItem">Test:<\\/td>/<td class="headerItem">Commit:<\\/td>/g'
    # Cleanup intermediate data
    COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.base.info ${_outputname}.test.info ${_outputname}.total.info ${_outputname}.filtered.info
    # Set work directory for target
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # Echo what is being done
    COMMENT "Code coverage report"
  )

  # Show info where to find the report
  ADD_CUSTOM_COMMAND(TARGET ${_targetname} POST_BUILD COMMAND ;
    COMMENT "Code coverage report at ./${_outputname}/index.html"
  )

ENDFUNCTION()
