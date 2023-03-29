################################################################################
#
# \file      CppCheck.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Setup target for code coverage analysis
#
################################################################################

find_program( CPPCHECK cppcheck )
find_program( CPPCHECK_HTMLREPORT cppcheck-htmlreport )

if(CPPCHECK AND CPPCHECK_HTMLREPORT)

  ADD_CUSTOM_TARGET(cppcheck
    # Run cppcheck static analysis
    COMMAND ${CPPCHECK} --inline-suppr --enable=all --force
                        --error-exitcode=1 -j${PROCESSOR_COUNT}
                        -I${PROJECT_SOURCE_DIR}/Base
                        -I${PROJECT_SOURCE_DIR}/Control
                        -I${PROJECT_SOURCE_DIR}/NoWarning
                        -I${PROJECT_BINARY_DIR}/Main
                        ${CMAKE_CURRENT_SOURCE_DIR}
    # Set work directory for target
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # Echo what is being done
    COMMENT "Quinoa cppcheck static analysis report"
    VERBATIM USES_TERMINAL
  )
  # Output code coverage target enabled
  message(STATUS "Enabling cppcheck static analysis target 'cppcheck'")

  find_program( FILEFIND find )
  find_program( SED sed )

  if(FILEFIND AND SED)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/cppcheck)
    string(REGEX REPLACE "Quinoa_v(.*)-(.*)-g(.*)" "\\3" SHA1 ${GIT_SHA1})
    ADD_CUSTOM_TARGET(cppcheck-xml
      # Run cppcheck static analysis
      COMMAND ${CPPCHECK} --inline-suppr --enable=all --force
                          --xml --xml-version=2 -j${PROCESSOR_COUNT}
                          -I${PROJECT_SOURCE_DIR}/Base
                          -I${PROJECT_SOURCE_DIR}/Control
                          -I${PROJECT_SOURCE_DIR}/NoWarning
                          -I${PROJECT_BINARY_DIR}/Main
                          ${CMAKE_CURRENT_SOURCE_DIR}
                          2> doc/cppcheck/cppcheck-report.xml
      # Generate html output
      COMMAND ${CPPCHECK_HTMLREPORT} --file=doc/cppcheck/cppcheck-report.xml
              --report-dir=doc/cppcheck --source-dir=.
      # Customize page headers in generated html
      COMMAND ${FILEFIND} doc/cppcheck -type f -exec ${SED} -i "s/project name/<a href=https:\\/\\/github.com\\/quinoacomputing\\/quinoa\\/commit\\/${SHA1}>${GIT_SHA1}<\\/a>/g" {} +
      # Set work directory for target
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      # Echo what is being done
      COMMENT "Quinoa cppcheck-xml static analysis report"
      VERBATIM USES_TERMINAL
    )
    # Output code coverage target enabled
    message(STATUS "Enabling cppcheck static analysis target 'cppcheck-xml', report at ${CMAKE_BINARY_DIR}/doc/cppcheck/index.html")
  endif()

endif()
