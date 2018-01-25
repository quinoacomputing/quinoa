################################################################################
#
# \file      cmake/CppCheck.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
# \brief     Setup target for code coverage analysis
#
################################################################################

find_program( SED sed )
find_program( CPPCHECK cppcheck )
find_program( CPPCHECK_HTMLREPORT cppcheck-htmlreport )

if(SED AND CPPCHECK AND CPPCHECK_HTMLREPORT)

  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/cppcheck)

  string(REGEX REPLACE "Quinoa_v(.*)-(.*)-g(.*)" "\\3" SHA1 ${GIT_SHA1})

  # Setup cppcheck static analysis target
  ADD_CUSTOM_TARGET(cppcheck
    # Run cppcheck static analysis
    COMMAND ${CPPCHECK} --xml --xml-version=2 --enable=all
            ${CMAKE_CURRENT_SOURCE_DIR} 2> doc/cppcheck/cppcheck-report.xml
    # Generate html output
    COMMAND ${CPPCHECK_HTMLREPORT} --file=doc/cppcheck/cppcheck-report.xml
            --report-dir=doc/cppcheck --source-dir=.
    # Customize page headers in generated html
    COMMAND find . -type f -print | xargs file | cut -f1 -d: | xargs ${SED} -i 's/project name/<a href="https:\\/\\/github.com\\/quinoacomputing\\/quinoa\\/commit\\/${SHA1}">${GIT_SHA1}<\\/a>/g'
    # Set work directory for target
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # Echo what is being done
    COMMENT "Quinoa cppcheck static analysis report"
  )

  # Output code coverage target enabled
  string(REPLACE ";" " " ARGUMENTS "${ARG_TESTRUNNER_ARGS}")
  message(STATUS "Enabling cppcheck static analysis target 'cppcheck', report at ${CMAKE_BINARY_DIR}/doc/cppcheck/index.html")

endif()
